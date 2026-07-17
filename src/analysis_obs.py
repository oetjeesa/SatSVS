import os
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from numpy.linalg import norm
from math import sin, cos, asin, degrees, radians
import xarray as xr
from astropy import time
from astropy.coordinates import get_sun, ITRS
import pandas as pd

# Project modules
from constants import R_EARTH
from analysis import AnalysisBase, AnalysisObs, AnalysisPlot3D, make_map_cyl, make_map_polar
import misc_fn
import logging_svs as ls


# from multiprocessing import Process, Value, Array, RawArray

class AnalysisObsSwathConical(AnalysisBase, AnalysisObs, AnalysisPlot3D):

    def __init__(self):
        super().__init__()
        self.polar_view = None
        self.revisit = None
        self.statistic = None
        self.user_pos_ecf = None
        self.user_metric = None
        self.save_output = None
        self.earth_angle_swath = None
        self.init_3d()
        self.sat_pos_hist = None  # Satellite ECI positions for the 3D plot
        self.swath_edges = None  # Left/right swath edge points (2D ribbon map + 3D plot)

    def read_config(self, node):
        if node.find('PolarView') is not None:
            self.polar_view = float(node.find('PolarView').text)
        if node.find('Revisit') is not None:
            self.revisit = misc_fn.str2bool(node.find('Revisit').text)
        if node.find('Statistic') is not None:
            self.statistic = node.find('Statistic').text.lower()
        if node.find('SaveOutput') is not None:
            self.save_output = node.find('SaveOutput').text.lower()
        self.read_config_3d(node)

    def before_loop(self, sm):
        self.det_angles_from_swath_before_loop(sm)
        self.user_pos_ecf = np.zeros((len(sm.users),4))
        self.user_metric = np.zeros((len(sm.users), sm.num_epoch), dtype=np.uint8)
        for idx_user, user in enumerate(sm.users):
            self.user_pos_ecf[idx_user,0:3] = user.pos_ecf
            self.user_pos_ecf[idx_user,3] = norm(self.user_pos_ecf[idx_user,0:3])
        # Swath edge/orbit histories: the edges feed the smooth 2D ribbon map
        # and, together with the ECI positions, the optional 3D globe render
        self.sat_pos_hist = np.zeros((len(sm.satellites), sm.num_epoch, 3))
        self.swath_edges = np.zeros((len(sm.satellites), sm.num_epoch, 2, 3))

    def det_angles_from_swath_before_loop(self, sm):
        for satellite in sm.satellites:
            idx_found = 0
            for idx, constellation in enumerate(sm.constellations):
                if satellite.constellation_id == constellation.constellation_id:
                    idx_found = idx
            const = sm.constellations[idx_found]
            sat_altitude = satellite.kepler.semi_major_axis - R_EARTH
            if const.obs_swath_stop is not None:  # if swath defined by swath length rather than incidence
                satellite.obs_swath_stop = const.obs_swath_stop  #
                satellite.obs_inci_angle_stop = misc_fn.incl_from_swath(
                    const.obs_swath_stop, R_EARTH, sat_altitude)
            else:
                satellite.obs_inci_angle_stop = const.obs_inci_angle_stop
            alfa_critical = asin(R_EARTH / (R_EARTH + sat_altitude))  # If incidence angle shooting off Earth -> error
            if satellite.obs_inci_angle_stop > alfa_critical:
                ls.logger.error(f'Incidence angle stop: {degrees(satellite.obs_inci_angle_stop)} ' +
                                f'larger than critical angle {round(degrees(alfa_critical),1)}')
                exit()

    def in_loop(self, sm):
        # Computed by angle distance point and satellite ground point
        # Just 10% faster if done by checking normal euclidean distance
        for idx_sat, satellite in enumerate(sm.satellites):
            self.det_angles_from_swath_in_loop(satellite)
            self.user_metric[:,sm.cnt_epoch] = \
                misc_fn.check_users_from_nadir(self.user_metric, self.user_pos_ecf, satellite.pos_ecf,
                                               self.earth_angle_swath, sm.cnt_epoch)
            # Cross-track swath extremes: the subsatellite ground point
            # rotated about the along-track horizontal axis by the swath
            # earth angle, in both directions. The orbit history is kept in
            # ECI (plot_3d draws the inertial path at the final epoch).
            r_hat = satellite.pos_ecf / norm(satellite.pos_ecf)
            ground = r_hat * misc_fn.earth_radius_lat(satellite.lla[0])
            axis = np.array(satellite.vel_ecf) - np.dot(satellite.vel_ecf, r_hat) * r_hat
            self.swath_edges[idx_sat, sm.cnt_epoch, 0] = \
                misc_fn.rot_vec_vec(ground, axis, self.earth_angle_swath)
            self.swath_edges[idx_sat, sm.cnt_epoch, 1] = \
                misc_fn.rot_vec_vec(ground, axis, -self.earth_angle_swath)
            self.sat_pos_hist[idx_sat, sm.cnt_epoch] = satellite.pos_eci

    def det_angles_from_swath_in_loop(self, satellite):
        satellite.det_lla()
        r_earth = misc_fn.earth_radius_lat(satellite.lla[0])
        sat_altitude = norm(satellite.pos_ecf) - r_earth
        if satellite.obs_swath_stop is not None:  # if swath defined by swath length rather than incidence
            satellite.obs_inci_angle_stop = misc_fn.incl_from_swath(satellite.obs_swath_stop, r_earth, sat_altitude)
        radius = misc_fn.det_swath_radius(sat_altitude, satellite.obs_inci_angle_stop, r_earth)
        self.earth_angle_swath = misc_fn.earth_angle_beta(radius, r_earth)

    def export2nc(self, sm, file_name):
        user3d_data = np.zeros((len(sm.user_latitudes),len(sm.user_longitudes),sm.num_epoch), dtype=np.uint8)
        for idx_usr, user in enumerate(sm.users):
            idx_lat = np.searchsorted(sm.user_latitudes,degrees(user.lla[0])).flatten()
            idx_lon = np.searchsorted(sm.user_longitudes,degrees(user.lla[1])).flatten()
            user3d_data[int(idx_lat),int(idx_lon),:] = self.user_metric[idx_usr,:]
        da = xr.DataArray(user3d_data,
                          dims=('lat', 'lon', 'time_mjd'),
                          coords={'lat': sm.user_latitudes,
                                  'lon': sm.user_longitudes,
                                  'time_mjd': self.times_mjd},
                          name='swath_coverage')
        da.to_netcdf(file_name)


    def after_loop(self, sm):

        if self.save_output=='numpy':
            np.save(sm.output_path('user_cov_swath'), self.user_metric)  # Save to numpy array
        if self.save_output=='netcdf':
            self.export2nc(sm, sm.output_path('user_cov_swath.nc'))  # Save to netcdf file

        write_swath_coverage_csv(self, sm)
        self.plot_swath_coverage(sm, self.swath_edges, self.polar_view)

        if self.revisit:
            self.plot_swath_revisit(sm, self.user_metric, self.statistic, self.polar_view)
            self.plot_swath_revisit_latitude(sm, self.user_metric)

        if self.plot_3d:
            plot_swath_3d_from_analysis(self, sm)

        if self.mp4:
            import plot_movie
            plot_movie.movie_ribbons_2d(sm, self.swath_edges,
                                        sm.output_path(self.type + '_2d.mp4'))
            self.render_movie_3d(sm, sm.satellites, self.sat_pos_hist,
                                 swath_edges=self.swath_edges)


def write_swath_coverage_csv(analysis, sm):
    """Shared data dump of the swath analyses: per user grid point the number
    of epochs inside the swath (the flag history behind revisit/coverage)."""
    lons = [degrees(user.lla[1]) for user in sm.users]
    lats = [degrees(user.lla[0]) for user in sm.users]
    analysis.write_csv(sm, ['lon_deg', 'lat_deg', 'epochs_in_swath'],
                       np.column_stack([lons, lats, analysis.user_metric.sum(axis=1)]))


def plot_swath_3d_from_analysis(analysis, sm):
    """Shared <Plot3D> hook of the swath analyses: render the recorded swath
    edge histories as 3D ribbons on the textured globe (needs pyvista)."""
    p3d = analysis._plot_3d_module()
    if p3d is None:
        return
    p3d.plot_swath_3d(sm, sm.satellites, analysis.sat_pos_hist,
                      analysis.swath_edges,
                      sm.output_path(analysis.type + '_3d.png'),
                      **analysis._kwargs_3d())


class AnalysisObsSwathPushBroom(AnalysisBase, AnalysisObs, AnalysisPlot3D):

    def __init__(self):
        super().__init__()
        self.polar_view = None
        self.revisit = None
        self.statistic = None
        self.planes = np.zeros((4,3))
        self.user_pos_ecf = None
        self.user_metric = None
        self.save_output = None
        self.init_3d()
        self.sat_pos_hist = None  # Satellite ECI positions for the 3D plot
        self.swath_edges = None  # Left/right swath edge points (2D ribbon map + 3D plot)

    def read_config(self, node):
        if node.find('PolarView') is not None:
            self.polar_view = float(node.find('PolarView').text)
        if node.find('Revisit') is not None:
            self.revisit = misc_fn.str2bool(node.find('Revisit').text)
        if node.find('Statistic') is not None:
            self.statistic = node.find('Statistic').text.lower()
        if node.find('SaveOutput') is not None:
            self.save_output = node.find('SaveOutput').text.lower()
        self.read_config_3d(node)

    def before_loop(self, sm):
        # Get the incidence angles for each of the satelllites
        self.det_angles_from_swath_before_loop(sm)
        self.user_pos_ecf = np.zeros((len(sm.users),3))  # User position in ECF
        self.user_metric = np.zeros((len(sm.users), sm.num_epoch), dtype=np.uint8)  # Range
        # self.shared_array = RawArray('i', len(sm.users))
        for idx_user, user in enumerate(sm.users):
            self.user_pos_ecf[idx_user,:] = user.pos_ecf
        # Swath edge/orbit histories: the edges feed the smooth 2D ribbon map
        # and, together with the ECI positions, the optional 3D globe render
        self.sat_pos_hist = np.zeros((len(sm.satellites), sm.num_epoch, 3))
        self.swath_edges = np.zeros((len(sm.satellites), sm.num_epoch, 2, 3))

    def det_angles_from_swath_before_loop(self, sm):
        # Previous-epoch swath edge points per satellite (see in_loop)
        self._prev_edges = np.zeros((len(sm.satellites), 2, 3))
        for satellite in sm.satellites:
            idx_found = 0
            for idx, constellation in enumerate(sm.constellations):
                if satellite.constellation_id == constellation.constellation_id:
                    idx_found = idx
            const = sm.constellations[idx_found]
            sat_altitude = satellite.kepler.semi_major_axis - R_EARTH
            if const.obs_swath_start is not None:  # if swath defined by swath length rather than incidence
                satellite.obs_swath_start = const.obs_swath_start  # Copy over from constellation
                satellite.obs_inci_angle_start = misc_fn.incl_from_swath(
                    const.obs_swath_start, R_EARTH, sat_altitude)
            else:
                satellite.obs_inci_angle_start = const.obs_inci_angle_start
            if const.obs_swath_stop is not None:  # if swath defined by swath length rather than incidence
                satellite.obs_swath_stop = const.obs_swath_stop  # Copy over from constellation
                satellite.obs_inci_angle_stop = misc_fn.incl_from_swath(
                    const.obs_swath_stop, R_EARTH, sat_altitude)
            else:
                satellite.obs_inci_angle_stop = const.obs_inci_angle_stop
            alfa_critical = asin(R_EARTH / (R_EARTH + sat_altitude))  # If incidence angle shooting off Earth -> error
            if np.abs(satellite.obs_inci_angle_start) > alfa_critical:
                ls.logger.error(f'Incidence angle start: {degrees(satellite.obs_inci_angle_start)} ' +
                                f'larger than critical angle {round(degrees(alfa_critical),1)}')
                exit()
            if np.abs(satellite.obs_inci_angle_stop) > alfa_critical:
                ls.logger.error(f'Incidence angle stop: {degrees(satellite.obs_inci_angle_stop)} ' +
                                f'larger than critical angle {round(degrees(alfa_critical),1)}')
                exit()

    def in_loop(self, sm):

        for idx_sat, satellite in enumerate(sm.satellites):
            r_earth = self.det_angles_from_swath_in_loop(satellite)
            point_vec1 = misc_fn.rot_vec_vec(-satellite.pos_ecf, np.array(satellite.vel_ecf),
                                             -satellite.obs_inci_angle_start)  # minus for right looking, plus for left
            point_vec2 = misc_fn.rot_vec_vec(-satellite.pos_ecf, np.array(satellite.vel_ecf),
                                             -satellite.obs_inci_angle_stop)  # minus for right looking, plus for left
            intersect, p1b, satellite.p1 = misc_fn.line_sphere_intersect(
                satellite.pos_ecf, satellite.pos_ecf + point_vec1, r_earth, np.zeros(3))
            intersect, p2b, satellite.p2 = misc_fn.line_sphere_intersect(
                satellite.pos_ecf, satellite.pos_ecf + point_vec2, r_earth, np.zeros(3))
            # 4 Planes of pyramid need to be carefully chosen with normal outwards of
            # pyramid. The previous-epoch edge points are kept per analysis
            # (self._prev_edges), NOT on the satellite: with several push-broom
            # analyses in one run a shared satellite attribute would already be
            # overwritten within the epoch, degenerating the swath pyramid
            prev1 = self._prev_edges[idx_sat, 0]
            prev2 = self._prev_edges[idx_sat, 1]
            self.planes[0,:] = misc_fn.plane_normal(satellite.p1, satellite.p2)
            self.planes[1,:] = misc_fn.plane_normal(prev2, prev1)
            self.planes[2,:] = misc_fn.plane_normal(satellite.p2, prev2)
            self.planes[3,:] = misc_fn.plane_normal(prev1, satellite.p1)
            self._prev_edges[idx_sat, 0] = satellite.p1
            self._prev_edges[idx_sat, 1] = satellite.p2
            if sm.cnt_epoch > 0:  # Now valid point 3 and 4
                # misc_fn.check_users_in_plane(
                #      self.user_pos_ecf, self.planes, self.shared_array)
                self.user_metric[:,sm.cnt_epoch] = misc_fn.check_users_in_plane(self.user_metric, self.user_pos_ecf,
                                                                                self.planes, sm.cnt_epoch)
            # The two line-of-sight ground intersections are the swath edges;
            # the orbit history is kept in ECI (plot_3d draws the inertial
            # path at the final epoch)
            self.swath_edges[idx_sat, sm.cnt_epoch, 0] = satellite.p1
            self.swath_edges[idx_sat, sm.cnt_epoch, 1] = satellite.p2
            self.sat_pos_hist[idx_sat, sm.cnt_epoch] = satellite.pos_eci

    def det_angles_from_swath_in_loop(self, satellite):

        satellite.det_lla()
        r_earth = misc_fn.earth_radius_lat(satellite.lla[0])
        sat_altitude = norm(satellite.pos_ecf) - r_earth
        if satellite.obs_swath_start is not None:  # if swath defined by swath length rather than incidence
            satellite.obs_inci_angle_start = misc_fn.incl_from_swath(satellite.obs_swath_start, r_earth, sat_altitude)
        if satellite.obs_swath_stop is not None:  # if swath defined by swath length rather than incidence
            satellite.obs_inci_angle_stop = misc_fn.incl_from_swath(satellite.obs_swath_stop, r_earth, sat_altitude)
        return r_earth

    def export2nc(self, sm, file_name):
        user3d_data = np.zeros((len(sm.user_latitudes),len(sm.user_longitudes),sm.num_epoch), dtype=np.uint8)
        for idx_usr, user in enumerate(sm.users):
            idx_lat = np.searchsorted(sm.user_latitudes,degrees(user.lla[0])).flatten()
            idx_lon = np.searchsorted(sm.user_longitudes,degrees(user.lla[1])).flatten()
            user3d_data[int(idx_lat),int(idx_lon),:] = self.user_metric[idx_usr,:]
        da = xr.DataArray(user3d_data,
                          dims=('lat', 'lon', 'time_mjd'),
                          coords={'lat': sm.user_latitudes,
                                  'lon': sm.user_longitudes,
                                  'time_mjd': self.times_mjd},
                          name='swath_coverage')
        da.to_netcdf(file_name)

    def after_loop(self, sm):

        if self.save_output=='numpy':
            np.save(sm.output_path('user_cov_swath'), self.user_metric)  # Save to numpy array
        if self.save_output=='netcdf':
            self.export2nc(sm, sm.output_path('user_cov_swath.nc'))  # Save to netcdf file

        write_swath_coverage_csv(self, sm)
        self.plot_swath_coverage(sm, self.swath_edges, self.polar_view)

        if self.revisit:
            self.plot_swath_revisit(sm, self.user_metric, self.statistic, self.polar_view)
            self.plot_swath_revisit_latitude(sm, self.user_metric)

        if self.plot_3d:
            plot_swath_3d_from_analysis(self, sm)

        if self.mp4:
            import plot_movie
            plot_movie.movie_ribbons_2d(sm, self.swath_edges,
                                        sm.output_path(self.type + '_2d.mp4'))
            self.render_movie_3d(sm, sm.satellites, self.sat_pos_hist,
                                 swath_edges=self.swath_edges)


class AnalysisObsSzaPushBroom(AnalysisBase, AnalysisPlot3D): # In very early stages, runs but very slow

    # Tried it but solar angle computation makes this way too slow...
    # Just kept it not to loose the effort...

    def __init__(self):
        super().__init__()
        self.polar_view = None
        self.statistic = None
        self.planes = np.zeros((4,3))
        self.user_pos_ecf = None
        self.user_pos_lla = None
        self.user_metric = None
        self.save_output = None
        self.init_3d()

    def read_config(self, node):
        if node.find('PolarView') is not None:
            self.polar_view = float(node.find('PolarView').text)
        if node.find('Statistic') is not None:
            self.statistic = node.find('Statistic').text.lower()
        if node.find('SaveOutput') is not None:
            self.save_output = node.find('SaveOutput').text.lower()
        self.read_config_3d(node)

    def before_loop(self, sm):
        # Get the incidence angles for each of the satelllites
        self.det_angles_from_swath_before_loop(sm)
        self.user_pos_ecf = np.zeros((len(sm.users),3))  # User position in ECF
        self.user_pos_lla = np.zeros((len(sm.users),3))  # User lat,lon,alt in radians,m
        self.user_metric = np.zeros((len(sm.users), sm.num_epoch), dtype=float)  # Range
        # self.shared_array = RawArray('i', len(sm.users))
        for idx_user, user in enumerate(sm.users):
            self.user_pos_ecf[idx_user,:] = user.pos_ecf
            self.user_pos_lla[idx_user, :] = user.lla
        self.before_loop_3d(sm)

    def det_angles_from_swath_before_loop(self, sm):
        # Previous-epoch swath edge points per satellite (see in_loop)
        self._prev_edges = np.zeros((len(sm.satellites), 2, 3))
        for satellite in sm.satellites:
            idx_found = 0
            for idx, constellation in enumerate(sm.constellations):
                if satellite.constellation_id == constellation.constellation_id:
                    idx_found = idx
            const = sm.constellations[idx_found]
            sat_altitude = satellite.kepler.semi_major_axis - R_EARTH
            if const.obs_swath_start is not None:  # if swath defined by swath length rather than incidence
                satellite.obs_swath_start = const.obs_swath_start  # Copy over from constellation
                satellite.obs_inci_angle_start = misc_fn.incl_from_swath(
                    const.obs_swath_start, R_EARTH, sat_altitude)
            else:
                satellite.obs_inci_angle_start = const.obs_inci_angle_start
            if const.obs_swath_stop is not None:  # if swath defined by swath length rather than incidence
                satellite.obs_swath_stop = const.obs_swath_stop  # Copy over from constellation
                satellite.obs_inci_angle_stop = misc_fn.incl_from_swath(
                    const.obs_swath_stop, R_EARTH, sat_altitude)
            else:
                satellite.obs_inci_angle_stop = const.obs_inci_angle_stop
            alfa_critical = asin(R_EARTH / (R_EARTH + sat_altitude))  # If incidence angle shooting off Earth -> error
            if np.abs(satellite.obs_inci_angle_start) > alfa_critical:
                ls.logger.error(f'Incidence angle start: {degrees(satellite.obs_inci_angle_start)} ' +
                                f'larger than critical angle {round(degrees(alfa_critical),1)}')
                exit()
            if np.abs(satellite.obs_inci_angle_stop) > alfa_critical:
                ls.logger.error(f'Incidence angle stop: {degrees(satellite.obs_inci_angle_stop)} ' +
                                f'larger than critical angle {round(degrees(alfa_critical),1)}')
                exit()

    def in_loop(self, sm):

        epoch = time.Time(sm.time_mjd, format='mjd')
        epoch.delta_ut1_utc = 0.0  # avoid getting IERS outside range error
        for idx_sat, satellite in enumerate(sm.satellites):
            r_earth = self.det_angles_from_swath_in_loop(satellite)
            point_vec1 = misc_fn.rot_vec_vec(-satellite.pos_ecf, np.array(satellite.vel_ecf),
                                             -satellite.obs_inci_angle_start)  # minus for right looking, plus for left
            point_vec2 = misc_fn.rot_vec_vec(-satellite.pos_ecf, np.array(satellite.vel_ecf),
                                             -satellite.obs_inci_angle_stop)  # minus for right looking, plus for left
            intersect, p1b, satellite.p1 = misc_fn.line_sphere_intersect(
                satellite.pos_ecf, satellite.pos_ecf + point_vec1, r_earth, np.zeros(3))
            intersect, p2b, satellite.p2 = misc_fn.line_sphere_intersect(
                satellite.pos_ecf, satellite.pos_ecf + point_vec2, r_earth, np.zeros(3))
            # 4 Planes of pyramid need to be carefully chosen with normal outwards of
            # pyramid. The previous-epoch edge points are kept per analysis
            # (self._prev_edges), NOT on the satellite: with several push-broom
            # analyses in one run a shared satellite attribute would already be
            # overwritten within the epoch, degenerating the swath pyramid
            prev1 = self._prev_edges[idx_sat, 0]
            prev2 = self._prev_edges[idx_sat, 1]
            self.planes[0,:] = misc_fn.plane_normal(satellite.p1, satellite.p2)
            self.planes[1,:] = misc_fn.plane_normal(prev2, prev1)
            self.planes[2,:] = misc_fn.plane_normal(satellite.p2, prev2)
            self.planes[3,:] = misc_fn.plane_normal(prev1, satellite.p1)
            self._prev_edges[idx_sat, 0] = satellite.p1
            self._prev_edges[idx_sat, 1] = satellite.p2
            if sm.cnt_epoch > 0:  # Now valid point 3 and 4
                # misc_fn.check_users_in_plane(
                #      self.user_pos_ecf, self.planes, self.shared_array)
                self.user_metric[:, sm.cnt_epoch] = misc_fn.check_users_in_plane(self.user_metric, self.user_pos_ecf,
                                                                                 self.planes, sm.cnt_epoch)
                self.user_metric[:, sm.cnt_epoch] = misc_fn.det_sza_fast(self.user_metric, self.user_pos_lla,
                                                                        epoch, sm.cnt_epoch)
        self.in_loop_3d(sm)

    def det_angles_from_swath_in_loop(self, satellite):

        satellite.det_lla()
        r_earth = misc_fn.earth_radius_lat(satellite.lla[0])
        sat_altitude = norm(satellite.pos_ecf) - r_earth
        if satellite.obs_swath_start is not None:  # if swath defined by swath length rather than incidence
            satellite.obs_inci_angle_start = misc_fn.incl_from_swath(satellite.obs_swath_start, r_earth, sat_altitude)
        if satellite.obs_swath_stop is not None:  # if swath defined by swath length rather than incidence
            satellite.obs_inci_angle_stop = misc_fn.incl_from_swath(satellite.obs_swath_stop, r_earth, sat_altitude)
        return r_earth

    def export2nc(self, sm, file_name):
        user3d_data = np.zeros((len(sm.user_latitudes),len(sm.user_longitudes),sm.num_epoch), dtype=np.uint8)
        for idx_usr, user in enumerate(sm.users):
            idx_lat = np.searchsorted(sm.user_latitudes,degrees(user.lla[0])).flatten()
            idx_lon = np.searchsorted(sm.user_longitudes,degrees(user.lla[1])).flatten()
            user3d_data[int(idx_lat),int(idx_lon),:] = self.user_metric[idx_usr,:]
        da = xr.DataArray(user3d_data,
                          dims=('lat', 'lon', 'time_mjd'),
                          coords={'lat': sm.user_latitudes,
                                  'lon': sm.user_longitudes,
                                  'time_mjd': self.times_mjd},
                          name='swath_coverage')
        da.to_netcdf(file_name)

    def after_loop(self, sm):

        if self.save_output=='numpy':
            np.save(sm.output_path('user_cov_swath'), self.user_metric)  # Save to numpy array
        if self.save_output=='netcdf':
            self.export2nc(sm, sm.output_path('user_cov_swath.nc'))  # Save to netcdf file

        self.plot_sza_coverage(sm, self.user_metric, self.polar_view)

        if self.plot_3d:
            points = []
            for idx_user, user in enumerate(sm.users):
                in_swath = self.user_metric[idx_user, np.nonzero(self.user_metric[idx_user, :])]
                if in_swath.size and in_swath.mean() >= 1:
                    points.append([degrees(user.lla[1]), degrees(user.lla[0]), in_swath.mean()])
            self.render_3d_points(sm, np.array(points), 'Solar Zenith Angle Mean [deg]',
                                  point_size=10)

    def plot_sza_coverage(self, sm, user_metric, polar_view):
        plot_points = np.zeros((len(sm.users), 3))
        for idx_user, user in enumerate(sm.users):
            if idx_user % 1000 == 0:
                ls.logger.info(f'User sza coverage {user.user_id} of {len(sm.users)}')
            if user_metric[idx_user, :].any():  # Any value bigger than 0
                sza_stat =  user_metric[idx_user,np.nonzero(user_metric[idx_user,:])].mean()
                if sza_stat >= 1:
                    plot_points[idx_user, :] = [degrees(user.lla[1]), degrees(user.lla[0]), sza_stat]
        plot_points = plot_points[~np.all(plot_points == 0, axis=1)]  # Clean up empty rows
        self.write_csv(sm, ['lon_deg', 'lat_deg', 'mean_sza_deg'], plot_points)
        if polar_view is not None:
            fig, ax = make_map_polar(polar_view)
            self.decorate_map2d(sm, ax)
        else:
            fig, ax = make_map_cyl()
            self.decorate_map2d(sm, ax)
        sc = ax.scatter(plot_points[:,0], plot_points[:,1], s=12, marker='o', cmap=plt.cm.jet,
                        c=plot_points[:,2], alpha=.3, transform=ccrs.PlateCarree())
        cb = plt.colorbar(sc, ax=ax, shrink=0.85)
        cb.set_label('Solar Zenith Angle Mean [deg]', fontsize=10)
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()



class AnalysisObsSzaSubSat(AnalysisBase, AnalysisPlot3D):

    def __init__(self):
        super().__init__()
        self.polar_view = None
        self.user_metric = None
        self.epoch = None
        self.save_output = None  # default if SaveOutput not in config
        self.range_lat = [-90, 91, 10]  # default if RangeLatitude not in config
        self.init_3d()

    def read_config(self, node):
        if node.find('PolarView') is not None:
            self.polar_view = float(node.find('PolarView').text)
        if node.find('SaveOutput') is not None:
            self.save_output = node.find('SaveOutput').text.lower()
        if node.find('RangeLatitude') is not None:
            self.range_lat = [int(i) for i in node.find('RangeLatitude').text.split(',')]
        self.read_config_3d(node)

    def before_loop(self, sm):
        # Get the incidence angles for each of the satelllites
        self.user_metric = np.zeros((sm.num_epoch,4))
        self.epoch = time.Time(sm.time_mjd, format='mjd')
        self.epoch.delta_ut1_utc = 0.0  # avoid getting IERS outside range error
        self.before_loop_3d(sm)

    def in_loop(self, sm):

        self.epoch = time.Time(sm.time_mjd, format='mjd')

        for satellite in sm.satellites:
            satellite.det_lla()
            self.user_metric[sm.cnt_epoch, 2] = misc_fn.det_sza([degrees(satellite.lla[0]),degrees(satellite.lla[1])], self.epoch)
            if (self.user_metric[sm.cnt_epoch, 2] != 0) :
                self.user_metric[sm.cnt_epoch, 0] = degrees(satellite.lla[0])
                self.user_metric[sm.cnt_epoch, 1] = degrees(satellite.lla[1])
                self.user_metric[sm.cnt_epoch, 3] = self.times_f_doy[sm.cnt_epoch]
        self.in_loop_3d(sm)

    def after_loop(self, sm):

        if self.plot_3d:
            points = self.user_metric[~np.all(self.user_metric == 0, axis=1)]
            self.render_3d_points(sm, points[:, [1, 0, 2]], 'Solar Zenith Angle [deg]')

        self.plot_sza_subsat(sm, self.user_metric, self.polar_view)
        self.plot_sza_latitude(sm, self.user_metric, self.polar_view, range(self.range_lat[0],self.range_lat[1],self.range_lat[2]))
        self.plot_sza_latitude_year(sm, self.user_metric, self.polar_view, range(self.range_lat[0],self.range_lat[1],self.range_lat[2]))

    def plot_sza_subsat(self, sm, user_metric, polar_view):

        self.user_metric = self.user_metric[~np.all(self.user_metric == 0, axis=1)]
        self.write_csv(sm, ['lat_deg', 'lon_deg', 'sza_deg', 'doy'], self.user_metric)
        if polar_view is not None:
            fig, ax = make_map_polar(polar_view)
            self.decorate_map2d(sm, ax)
        else:
            fig, ax = make_map_cyl()
            self.decorate_map2d(sm, ax)
        sc = ax.scatter(self.user_metric[:, 1], self.user_metric[:, 0], s=12, marker='o',
                        cmap=plt.cm.jet, c=self.user_metric[:, 2], alpha=.5,
                        transform=ccrs.PlateCarree())
        cb = plt.colorbar(sc, ax=ax, shrink=0.85)
        cb.set_label('Solar Zenith Angle [deg]', fontsize=10)
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

    def plot_sza_latitude(self, sm, user_metric, polar_view, range_lat):

        self.user_metric = self.user_metric[~np.all(self.user_metric == 0, axis=1)] # Clean up night values where Sun was not visible

        df = pd.DataFrame({'Latitude': self.user_metric[:, 0], 'Longitude': self.user_metric[:, 1],
                                'SZA' : self.user_metric[:,2], 'DOY' : self.user_metric[:, 3]})

        step_size = range_lat[1]-range_lat[0]
        results = np.zeros((len(range_lat),2))
        for i, lat in enumerate(range_lat):
            results[i,0] = lat
            results[i,1] = df[(df.Latitude>lat-step_size/2) & (df.Latitude<lat+step_size/2)].SZA.mean()

        self.write_csv(sm, ['lat_deg', 'mean_sza_deg'], results, suffix='lat')
        fig = plt.figure(figsize=(10, 5))
        plt.plot(results[:,0],results[:,1])
        plt.xlabel('Latitude [deg]')
        plt.ylabel('Solar Zenith Angle [deg]')
        plt.grid()
        plt.savefig(sm.output_path(self.type + '_lat.png'))
        plt.show()

        if self.save_output=='numpy':
            np.save(sm.output_path('user_sza_latitude'), results)  # Save to numpy array


    def plot_sza_latitude_year(self, sm, user_metric, polar_view, range_lat):

        self.user_metric = self.user_metric[~np.all(self.user_metric == 0, axis=1)]  # Clean up night values where Sun was not visible

        df = pd.DataFrame({'Latitude': self.user_metric[:, 0], 'Longitude': self.user_metric[:, 1],
                           'SZA': self.user_metric[:, 2], 'DOY': self.user_metric[:, 3]})

        step_size = range_lat[1]-range_lat[0]
        fig = plt.figure(figsize=(10, 5))

        for i, lat in enumerate(range_lat):
            df2 = df[(df.Latitude > lat - step_size / 2) & (df.Latitude < lat + step_size / 2)]
            # Daily-mean SZA plotted at the day-of-year bin midpoints (the
            # positional index would collapse short runs onto x=0, giving an
            # empty-looking plot); markers keep single-day runs visible
            df3 = df2.groupby(pd.cut(df2["DOY"], np.arange(0, 367, 1)),
                              observed=True).SZA.mean().reset_index().dropna()
            plt.plot([interval.mid for interval in df3.DOY], df3.SZA, 'o-',
                     markersize=3, linewidth=1.0, label=str(lat))

        plt.legend()
        plt.xlabel('DOY [-]')
        plt.ylabel('Solar Zenith Angle [deg]')
        plt.grid()
        plt.savefig(sm.output_path(self.type + '_lat_year.png'))
        plt.show()


class AnalysisObsAoiRevisit(AnalysisObsSwathPushBroom):
    """Revisit and coverage build-up statistics over an area of interest.
    The AOI is the configured user segment - typically Type Polygon (a grid
    clipped to an inline polygon or a shapefile), but a regional Grid works
    too. The instrument is the push-broom swath defined in the
    <Constellation> block, evaluated with the same machinery as
    obs_swath_push_broom. Produces a map of the revisit statistic per AOI
    grid point (zoomed to the AOI), the fraction of the AOI covered versus
    time (fill-up curve), and aggregate revisit numbers in the log."""

    def __init__(self):
        super().__init__()
        self.statistic = 'max'

    def read_config(self, node):
        if node.find('Statistic') is not None:
            self.statistic = node.find('Statistic').text.lower()

    def after_loop(self, sm):
        if not sm.users:
            ls.logger.error(f'{self.type} needs a user segment (Type Polygon or '
                            f'Grid) as the area of interest. No plot produced.')
            return
        lons = np.array([degrees(user.lla[1]) for user in sm.users])
        lats = np.array([degrees(user.lla[0]) for user in sm.users])
        covered = self.user_metric > 0  # (num_user, num_epoch)

        # Per-point pass count and revisit gaps (hours between successive passes)
        num_passes = np.array([int(c[0]) + int((c[1:] & ~c[:-1]).sum()) for c in covered])
        gaps_list = self.revisit_gaps_hours(self.user_metric, sm.time_step)
        stat_fn = {'min': np.min, 'mean': np.mean, 'max': np.max,
                   'std': np.std, 'median': np.median}.get(self.statistic, np.max)
        stat_gap = np.array([stat_fn(g) if len(g) else np.nan for g in gaps_list])
        mean_gap = np.array([np.mean(g) if len(g) else np.nan for g in gaps_list])
        max_gap = np.array([np.max(g) if len(g) else np.nan for g in gaps_list])

        # AOI fill-up: fraction of the AOI points seen at least once up to t
        fraction = (np.cumsum(covered, axis=1) > 0).mean(axis=0)
        times = np.asarray(self.times_f_doy)
        name = getattr(sm.users[0], 'name', '') or 'AOI'
        ls.logger.info(f'{self.type}: {name} {len(sm.users)} grid points, '
                       f'{fraction[-1] * 100:.1f}% covered at the end of the run')
        for level in (0.5, 0.9, 0.99):
            reached = np.flatnonzero(fraction >= level)
            if reached.size:
                hours = reached[0] * sm.time_step / 3600.0
                ls.logger.info(f'{self.type}: {level * 100:.0f}% of the AOI covered '
                               f'after {hours:.1f} h')
        if np.isfinite(mean_gap).any():
            ls.logger.info(f'{self.type}: revisit over the AOI: mean of mean gaps '
                           f'{np.nanmean(mean_gap):.1f} h, worst max gap '
                           f'{np.nanmax(max_gap):.1f} h')

        # Map of the revisit statistic, zoomed to the AOI with a margin
        # (regional map with automatic gridline spacing instead of the fixed
        # 60/30 degree global-map locators)
        fig = plt.figure(figsize=(8, 6), layout='constrained')
        ax = plt.axes(projection=ccrs.PlateCarree())
        margin = 5.0
        ax.set_extent([lons.min() - margin, lons.max() + margin,
                       max(lats.min() - margin, -90.0), min(lats.max() + margin, 90.0)],
                      ccrs.PlateCarree())
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray')
        gl.top_labels = False
        gl.right_labels = False
        ax.coastlines()
        self.decorate_map2d(sm, ax)
        if np.isfinite(stat_gap).any():
            sc = ax.scatter(lons, lats, s=14, marker='s', cmap=plt.cm.jet, c=stat_gap,
                            transform=ccrs.PlateCarree())
            cb = plt.colorbar(sc, ax=ax, shrink=0.85)
            cb.set_label(f'{self.statistic.capitalize()} Revisit Interval [hours]',
                         fontsize=10)
        else:  # No point with two passes yet: show the pass count instead
            sc = ax.scatter(lons, lats, s=14, marker='s', cmap=plt.cm.jet, c=num_passes,
                            transform=ccrs.PlateCarree())
            cb = plt.colorbar(sc, ax=ax, shrink=0.85)
            cb.set_label('Number of Passes [-]', fontsize=10)
        ax.set_title(f'{name}: revisit over the area of interest')
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        # AOI coverage fill-up curve
        fig = plt.figure(figsize=(10, 5))
        plt.plot((times - times[0]) * 24.0, fraction * 100.0, 'b-')
        plt.xlabel('Elapsed time [hours]')
        plt.ylabel('AOI covered at least once [%]')
        plt.ylim(0, 105)
        plt.grid()
        plt.savefig(sm.output_path(self.type + '_coverage.png'))
        plt.show()

        self.write_csv(sm, ['lon_deg', 'lat_deg', 'num_passes', 'mean_gap_hours',
                            'max_gap_hours'],
                       np.column_stack([lons, lats, num_passes, mean_gap, max_gap]))
        self.write_csv(sm, ['doy', 'covered_fraction'],
                       np.column_stack([times, fraction]), suffix='coverage')


class AnalysisObsTargetImaging(AnalysisBase):
    """Imaging opportunities over a list of point targets for an agile
    satellite: a target can be imaged when it lies within MaxOffNadir degrees
    of the satellite nadir direction (the pointing agility cone) and,
    optionally, while the Sun is at least MinSunElevation degrees above the
    target horizon (optical imaging daylight constraint). Reports the
    opportunity windows per target (start, duration, best off-nadir angle),
    the per-target opportunity counts and revisit gaps, and a map of the
    targets coloured by opportunity count."""

    def __init__(self):
        super().__init__()
        self.constellation_id = 0  # Optional selection (0 = all satellites)
        self.max_off_nadir = 0.0  # Mandatory, pointing agility cone [rad]
        self.min_sun_elevation = None  # Optional daylight constraint [rad]
        self.target_names = []
        self.target_lla = []  # Per target [lat, lon] in radians
        self.target_pos_ecf = None
        self.metric = None  # (num_target, num_epoch) best off-nadir [deg], -1 = not imaged

    def read_config(self, node):
        if node.find('ConstellationID') is not None:
            self.constellation_id = int(node.find('ConstellationID').text)
        self.max_off_nadir = radians(float(node.find('MaxOffNadir').text))
        if node.find('MinSunElevation') is not None:
            self.min_sun_elevation = radians(float(node.find('MinSunElevation').text))
        for target in node.findall('Target'):  # "Name, lat_deg, lon_deg"
            name, lat, lon = [v.strip() for v in target.text.split(',')]
            self.target_names.append(name)
            self.target_lla.append([radians(float(lat)), radians(float(lon))])
        if node.find('TargetFile') is not None:  # CSV lines "Name, lat_deg, lon_deg"
            with open(misc_fn.resolve_path(node.find('TargetFile').text)) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    name, lat, lon = [v.strip() for v in line.split(',')]
                    self.target_names.append(name)
                    self.target_lla.append([radians(float(lat)), radians(float(lon))])

    def before_loop(self, sm):
        num_target = len(self.target_lla)
        self.target_pos_ecf = np.zeros((num_target, 3))
        for i, (lat, lon) in enumerate(self.target_lla):
            self.target_pos_ecf[i] = misc_fn.lla2xyz(np.array([lat, lon, 0.0]))
        self.metric = np.full((num_target, sm.num_epoch), -1.0)
        ls.logger.info(f'{self.type}: {num_target} targets, max off-nadir '
                       f'{degrees(self.max_off_nadir):.1f} deg')

    def in_loop(self, sm):
        # Optional daylight constraint: one Sun direction (ECF) per epoch,
        # elevation over the target horizon from a dot product per target
        sun_ok = np.ones(len(self.target_pos_ecf), dtype=bool)
        if self.min_sun_elevation is not None:
            epoch = time.Time(sm.time_mjd, format='mjd')
            epoch.delta_ut1_utc = 0.0
            sun_itrs = get_sun(epoch).transform_to(ITRS(obstime=epoch))
            sun_dir = np.array([sun_itrs.x.value, sun_itrs.y.value, sun_itrs.z.value])
            sun_dir /= norm(sun_dir)
            for i, (lat, lon) in enumerate(self.target_lla):
                up = np.array([cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat)])
                sun_ok[i] = np.dot(up, sun_dir) >= sin(self.min_sun_elevation)

        for satellite in sm.satellites:
            if self.constellation_id > 0 and \
                    satellite.constellation_id != self.constellation_id:
                continue
            sat_pos = np.asarray(satellite.pos_ecf, dtype=float)
            norm_sat = norm(sat_pos)
            for i, target in enumerate(self.target_pos_ecf):
                if not sun_ok[i]:
                    continue
                to_target = target - sat_pos
                off_nadir = misc_fn.angle_two_vectors(-sat_pos, to_target,
                                                      norm_sat, norm(to_target))
                # Within the agility cone and above the target horizon
                if off_nadir <= self.max_off_nadir and np.dot(target, -to_target) > 0:
                    off_nadir_deg = degrees(off_nadir)
                    if self.metric[i, sm.cnt_epoch] < 0 or \
                            off_nadir_deg < self.metric[i, sm.cnt_epoch]:
                        self.metric[i, sm.cnt_epoch] = off_nadir_deg

    def after_loop(self, sm):
        times = np.asarray(self.times_f_doy)
        opportunities = []  # [target_id, doy_start, duration_min, min_off_nadir_deg]
        target_rows = []  # [target_id, lat, lon, num_opps, total_min, mean/max gap_h]
        for i, name in enumerate(self.target_names):
            imaged = self.metric[i] >= 0
            edges = np.diff(np.concatenate(([False], imaged, [False])).astype(np.int8))
            starts = np.flatnonzero(edges == 1)
            ends = np.flatnonzero(edges == -1) - 1
            for s, e in zip(starts, ends):
                opportunities.append([i + 1, times[s], (e - s + 1) * sm.time_step / 60.0,
                                      self.metric[i, s:e + 1].min()])
            gaps_h = ((starts[1:] - ends[:-1] - 1) * sm.time_step / 3600.0
                      if len(starts) > 1 else np.array([]))
            target_rows.append([i + 1, degrees(self.target_lla[i][0]),
                                degrees(self.target_lla[i][1]), len(starts),
                                imaged.sum() * sm.time_step / 60.0,
                                np.mean(gaps_h) if gaps_h.size else np.nan,
                                np.max(gaps_h) if gaps_h.size else np.nan])
            gap_str = (f'mean/max gap {gaps_h.mean():.1f}/{gaps_h.max():.1f} h'
                       if gaps_h.size else 'single or no opportunity')
            ls.logger.info(f'{self.type}: target {name}: {len(starts)} opportunities, '
                           f'{imaged.sum() * sm.time_step / 60.0:.1f} min total, {gap_str}')
        target_rows = np.array(target_rows)

        fig, ax = make_map_cyl()
        self.decorate_map2d(sm, ax)
        sc = ax.scatter(target_rows[:, 2], target_rows[:, 1], s=60, marker='^',
                        cmap=plt.cm.jet, c=target_rows[:, 3], edgecolors='black',
                        linewidths=0.5, transform=ccrs.PlateCarree(), zorder=5)
        for i, name in enumerate(self.target_names):
            ax.annotate(name, (target_rows[i, 2], target_rows[i, 1]),
                        xytext=(4, 4), textcoords='offset points', fontsize=8)
        cb = plt.colorbar(sc, ax=ax, shrink=0.85)
        cb.set_label('Number of Imaging Opportunities [-]', fontsize=10)
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        self.write_csv(sm, ['target_id', 'doy_start', 'duration_min', 'min_off_nadir_deg'],
                       opportunities)
        self.write_csv(sm, ['target_id', 'lat_deg', 'lon_deg', 'num_opportunities',
                            'total_duration_min', 'mean_gap_hours', 'max_gap_hours'],
                       target_rows, suffix='targets')