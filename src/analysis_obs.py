import os
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from numpy.linalg import norm
from math import sin, cos, asin, degrees, radians
import xarray as xr
from astropy import time
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
        self.sat_pos_hist = None  # Satellite ECF positions for the 3D plot
        self.swath_edges = None  # Left/right swath edge points for the 3D plot

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
        if self.plot_3d:
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
            if self.plot_3d:
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
                                  'time_mjd': sm.analysis.times_mjd},
                          name='swath_coverage')
        da.to_netcdf(file_name)


    def after_loop(self, sm):

        if self.save_output=='numpy':
            np.save('../output/user_cov_swath', self.user_metric)  # Save to numpy array
        if self.save_output=='netcdf':
            self.export2nc(sm, '../output/user_cov_swath.nc')  # Save to netcdf file

        self.plot_swath_coverage(sm, self.user_metric, self.polar_view)

        if self.revisit:
            self.plot_swath_revisit(sm, self.user_metric, self.statistic, self.polar_view)

        if self.plot_3d:
            plot_swath_3d_from_analysis(self, sm)


def plot_swath_3d_from_analysis(analysis, sm):
    """Shared <Plot3D> hook of the swath analyses: render the recorded swath
    edge histories as 3D ribbons on the textured globe (needs pyvista)."""
    p3d = analysis._plot_3d_module()
    if p3d is None:
        return
    p3d.plot_swath_3d(sm, sm.satellites, analysis.sat_pos_hist,
                      analysis.swath_edges,
                      '../output/' + analysis.type + '_3d.png',
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
        self.sat_pos_hist = None  # Satellite ECF positions for the 3D plot
        self.swath_edges = None  # Left/right swath edge points for the 3D plot

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
        if self.plot_3d:
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
            # 4 Planes of pyramid need to be carefully chosen with normal outwards of pyramid
            self.planes[0,:] = misc_fn.plane_normal(satellite.p1, satellite.p2)
            self.planes[1,:] = misc_fn.plane_normal(satellite.p4, satellite.p3)
            self.planes[2,:] = misc_fn.plane_normal(satellite.p2, satellite.p4)
            self.planes[3,:] = misc_fn.plane_normal(satellite.p3, satellite.p1)
            satellite.p3 = satellite.p1.copy()  # Copy for next run, without .copy() python just refers to same list
            satellite.p4 = satellite.p2.copy()
            if sm.cnt_epoch > 0:  # Now valid point 3 and 4
                # misc_fn.check_users_in_plane(
                #      self.user_pos_ecf, self.planes, self.shared_array)
                self.user_metric[:,sm.cnt_epoch] = misc_fn.check_users_in_plane(self.user_metric, self.user_pos_ecf,
                                                                                self.planes, sm.cnt_epoch)
            if self.plot_3d:
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
                                  'time_mjd': sm.analysis.times_mjd},
                          name='swath_coverage')
        da.to_netcdf(file_name)

    def after_loop(self, sm):

        if self.save_output=='numpy':
            np.save('../output/user_cov_swath', self.user_metric)  # Save to numpy array
        if self.save_output=='netcdf':
            self.export2nc(sm, '../output/user_cov_swath.nc')  # Save to netcdf file

        self.plot_swath_coverage(sm, self.user_metric, self.polar_view)

        if self.revisit:
            self.plot_swath_revisit(sm, self.user_metric, self.statistic, self.polar_view)

        if self.plot_3d:
            plot_swath_3d_from_analysis(self, sm)


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
        for satellite in sm.satellites:
            r_earth = self.det_angles_from_swath_in_loop(satellite)
            point_vec1 = misc_fn.rot_vec_vec(-satellite.pos_ecf, np.array(satellite.vel_ecf),
                                             -satellite.obs_inci_angle_start)  # minus for right looking, plus for left
            point_vec2 = misc_fn.rot_vec_vec(-satellite.pos_ecf, np.array(satellite.vel_ecf),
                                             -satellite.obs_inci_angle_stop)  # minus for right looking, plus for left
            intersect, p1b, satellite.p1 = misc_fn.line_sphere_intersect(
                satellite.pos_ecf, satellite.pos_ecf + point_vec1, r_earth, np.zeros(3))
            intersect, p2b, satellite.p2 = misc_fn.line_sphere_intersect(
                satellite.pos_ecf, satellite.pos_ecf + point_vec2, r_earth, np.zeros(3))
            # 4 Planes of pyramid need to be carefully chosen with normal outwards of pyramid
            self.planes[0,:] = misc_fn.plane_normal(satellite.p1, satellite.p2)
            self.planes[1,:] = misc_fn.plane_normal(satellite.p4, satellite.p3)
            self.planes[2,:] = misc_fn.plane_normal(satellite.p2, satellite.p4)
            self.planes[3,:] = misc_fn.plane_normal(satellite.p3, satellite.p1)
            satellite.p3 = satellite.p1.copy()  # Copy for next run, without .copy() python just refers to same list
            satellite.p4 = satellite.p2.copy()
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
                                  'time_mjd': sm.analysis.times_mjd},
                          name='swath_coverage')
        da.to_netcdf(file_name)

    def after_loop(self, sm):

        if self.save_output=='numpy':
            np.save('../output/user_cov_swath', self.user_metric)  # Save to numpy array
        if self.save_output=='netcdf':
            self.export2nc(sm, '../output/user_cov_swath.nc')  # Save to netcdf file

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
        if polar_view is not None:
            fig, ax = make_map_polar(polar_view)
        else:
            fig, ax = make_map_cyl()
        sc = ax.scatter(plot_points[:,0], plot_points[:,1], s=12, marker='o', cmap=plt.cm.jet,
                        c=plot_points[:,2], alpha=.3, transform=ccrs.PlateCarree())
        cb = plt.colorbar(sc, ax=ax, shrink=0.85)
        cb.set_label('Solar Zenith Angle Mean [deg]', fontsize=10)
        plt.savefig('../output/'+self.type+'.png')
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
        if polar_view is not None:
            fig, ax = make_map_polar(polar_view)
        else:
            fig, ax = make_map_cyl()
        sc = ax.scatter(self.user_metric[:, 1], self.user_metric[:, 0], s=12, marker='o',
                        cmap=plt.cm.jet, c=self.user_metric[:, 2], alpha=.5,
                        transform=ccrs.PlateCarree())
        cb = plt.colorbar(sc, ax=ax, shrink=0.85)
        cb.set_label('Solar Zenith Angle [deg]', fontsize=10)
        plt.savefig('../output/' + self.type + '.png')
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

        fig = plt.figure(figsize=(10, 5))
        plt.plot(results[:,0],results[:,1])
        plt.xlabel('Latitude [deg]')
        plt.ylabel('Solar Zenith Angle [deg]')
        plt.grid()
        plt.savefig('../output/' + self.type + '_lat.png')
        plt.show()

        if self.save_output=='numpy':
            np.save('../output/user_sza_latitude', results)  # Save to numpy array


    def plot_sza_latitude_year(self, sm, user_metric, polar_view, range_lat):

        self.user_metric = self.user_metric[~np.all(self.user_metric == 0, axis=1)]  # Clean up night values where Sun was not visible

        df = pd.DataFrame({'Latitude': self.user_metric[:, 0], 'Longitude': self.user_metric[:, 1],
                           'SZA': self.user_metric[:, 2], 'DOY': self.user_metric[:, 3]})

        step_size = range_lat[1]-range_lat[0]
        fig = plt.figure(figsize=(10, 5))

        for i, lat in enumerate(range_lat):
            print('*** Analysing latitude: '+str(lat))
            df2 = df[(df.Latitude > lat - step_size / 2) & (df.Latitude < lat + step_size / 2)]
            df3 = df2.groupby(pd.cut(df2["DOY"], np.arange(0, 365, 1))).SZA.mean().reset_index().dropna()
            plt.plot(df3.index, df3.SZA, label=str(lat))

        plt.legend()
        plt.xlabel('DOY [-]')
        plt.ylabel('Solar Zenith Angle [deg]')
        plt.grid()
        plt.savefig('../output/' + self.type + '_lat_year.png')
        plt.show()