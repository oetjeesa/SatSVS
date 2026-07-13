import os
import numpy as np
import pandas as pd
from math import degrees, radians
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# Import project modules
import misc_fn
import logging_svs as ls
from constants import PI
from analysis import AnalysisBase, AnalysisPlot3D, make_map_cyl, map_pcolormesh, get_user_grid_shape


class AnalysisCovDepthOfCoverage(AnalysisBase, AnalysisPlot3D):

    def __init__(self):
        super().__init__()
        self.sat_metric = None  # Per-satellite metric memory (num_sat, num_epoch, 3)
        self.init_3d()

    def read_config(self, node):
        self.read_config_3d(node)

    def before_loop(self, sm):
        # Columns: lat [deg], lon [deg], number of stations in view
        self.sat_metric = np.zeros((sm.num_sat, sm.num_epoch, 3))
        self.before_loop_3d(sm)

    def in_loop(self, sm):
        for idx_sat, satellite in enumerate(sm.satellites):
            satellite.det_lla()
            self.sat_metric[idx_sat, sm.cnt_epoch, 0] = degrees(satellite.lla[0])
            self.sat_metric[idx_sat, sm.cnt_epoch, 1] = degrees(satellite.lla[1])
            self.sat_metric[idx_sat, sm.cnt_epoch, 2] = len(satellite.idx_stat_in_view)
        self.in_loop_3d(sm)

    def after_loop(self, sm):
        fig, ax = make_map_cyl(figsize=(10, 4))
        sc = None
        for idx_sat in range(sm.num_sat):
            metric = self.sat_metric[idx_sat]
            sc = ax.scatter(metric[:, 1], metric[:, 0], cmap='RdYlBu',
                            c=metric[:, 2], vmin=0, vmax=len(sm.stations),
                            transform=ccrs.PlateCarree())
        plt.colorbar(sc, ax=ax, shrink=0.85)
        for station in sm.stations:
            ax.plot(degrees(station.lla[1]), degrees(station.lla[0]), 'r^',
                    transform=ccrs.PlateCarree())
        ax.text(50, 80, 'Red triangles: station locations', transform=ccrs.PlateCarree())
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        times = np.asarray(self.times_f_doy)
        self.write_csv(sm, ['doy', 'sat_id', 'lat_deg', 'lon_deg', 'num_stations_in_view'],
                       np.vstack([np.column_stack([times, np.full(sm.num_epoch, satellite.sat_id),
                                                   self.sat_metric[idx_sat]])
                                  for idx_sat, satellite in enumerate(sm.satellites)]))

        if self.plot_3d:
            points = np.vstack([m[:, [1, 0, 2]] for m in self.sat_metric])
            self.render_3d_points(sm, points, 'Number of stations in view [-]',
                                  cmap='RdYlBu', clim=(0, len(sm.stations)))

        if self.mp4:
            import plot_movie
            plot_movie.movie_track_2d(sm, self.sat_metric,
                                      sm.output_path(self.type + '_2d.mp4'),
                                      color_index=2, clim=(0, len(sm.stations)),
                                      label='Number of stations in view [-]')
            self.render_movie_3d(sm, sm.satellites, self.sat_pos_hist_3d,
                                 track_latlon=self.sat_metric[:, :, 0:2])


class AnalysisCovGroundTrack(AnalysisBase, AnalysisPlot3D):

    def __init__(self):
        super().__init__()
        self.constellation_id = 0  # Mandatory
        self.satellite_id = 0  # Optional
        self.sat_metric = None  # Per-satellite metric memory (num_sat, num_epoch, 5)
        self.init_3d()

    def read_config(self, node):
        if node.find('ConstellationID') is not None:
            self.constellation_id = int(node.find('ConstellationID').text)
        if node.find('SatelliteID') is not None:
            self.satellite_id = int(node.find('SatelliteID').text)
        self.read_config_3d(node)

    def _selected(self, satellite):
        if satellite.constellation_id != self.constellation_id:
            return False
        return self.satellite_id == 0 or satellite.sat_id == self.satellite_id

    def before_loop(self, sm):
        # Columns: lat [deg], lon [deg], ECI position x, y, z [m] (the 3D
        # plot draws the inertial orbit path oriented at the final epoch)
        self.sat_metric = np.zeros((sm.num_sat, sm.num_epoch, 5))

    def in_loop(self, sm):
        for idx_sat, satellite in enumerate(sm.satellites):
            if self._selected(satellite):
                satellite.det_lla()
                self.sat_metric[idx_sat, sm.cnt_epoch, 0] = degrees(satellite.lla[0])
                self.sat_metric[idx_sat, sm.cnt_epoch, 1] = degrees(satellite.lla[1])
                self.sat_metric[idx_sat, sm.cnt_epoch, 2:5] = satellite.pos_eci

    def after_loop(self, sm):

        fig, ax = make_map_cyl()
        if self.satellite_id > 0:  # Only for one satellite
            for idx_sat, satellite in enumerate(sm.satellites):
                if self._selected(satellite):
                    y, x = self.sat_metric[idx_sat, :, 0], self.sat_metric[idx_sat, :, 1]
                    ax.plot(x, y, 'r.', transform=ccrs.PlateCarree())
        else:
            for idx_sat, satellite in enumerate(sm.satellites):
                y, x = self.sat_metric[idx_sat, :, 0], self.sat_metric[idx_sat, :, 1]
                ax.plot(x, y, '+', label=str(satellite.sat_id), transform=ccrs.PlateCarree())
            ax.legend(fontsize=8)
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        times = np.asarray(self.times_f_doy)
        self.write_csv(sm, ['doy', 'sat_id', 'lat_deg', 'lon_deg'],
                       np.vstack([np.column_stack([times, np.full(sm.num_epoch, satellite.sat_id),
                                                   self.sat_metric[idx_sat, :, 0:2]])
                                  for idx_sat, satellite in enumerate(sm.satellites)
                                  if self._selected(satellite)]))

        if self.mp4:
            import plot_movie
            idx_selected = [i for i, s in enumerate(sm.satellites) if self._selected(s)]
            plot_movie.movie_track_2d(sm, self.sat_metric[idx_selected],
                                      sm.output_path(self.type + '_2d.mp4'))
            self.render_movie_3d(sm, [sm.satellites[i] for i in idx_selected],
                                 self.sat_metric[idx_selected][:, :, 2:5],
                                 track_latlon=self.sat_metric[idx_selected][:, :, 0:2])

        if self.plot_3d:
            p3d = self._plot_3d_module()
            if p3d is None:
                return
            idx_selected = [i for i, s in enumerate(sm.satellites) if self._selected(s)]
            p3d.plot_ground_track_3d(sm, [sm.satellites[i] for i in idx_selected],
                                     [self.sat_metric[i] for i in idx_selected],
                                     sm.output_path(self.type + '_3d.png'), **self._kwargs_3d())


class AnalysisCovPassTime(AnalysisBase, AnalysisPlot3D):

    def __init__(self):
        super().__init__()
        self.constellation_id = 0  # Mandatory
        self.statistic = ''  # Mandatory
        self.user_metric = None  # Per-user metric memory (num_user, num_epoch, num_sat)
        self.init_3d()

    def read_config(self, node):
        if node.find('ConstellationID') is not None:
            self.constellation_id = int(node.find('ConstellationID').text)
        if node.find('Statistic') is not None:
            self.statistic = node.find('Statistic').text
        self.read_config_3d(node)

    def before_loop(self, sm):
        self.user_metric = np.full((sm.num_user, sm.num_epoch, sm.num_sat), False, dtype=bool)
        self.before_loop_3d(sm)

    def in_loop(self, sm):
        for idx_user, user in enumerate(sm.users):
            for j in range(len(user.idx_sat_in_view)):
                if sm.satellites[user.idx_sat_in_view[j]].constellation_id == self.constellation_id:
                    self.user_metric[idx_user, sm.cnt_epoch, user.idx_sat_in_view[j]] = True
        self.in_loop_3d(sm)

    def after_loop(self, sm):
        lats, lons = [], []
        time_step = int((self.times_mjd[1] - self.times_mjd[0]) * 86400)
        metric = np.zeros(len(sm.users))
        for idx_usr, user in enumerate(sm.users):
            valid_value_list = []  # Define and clear
            metric_int = self.user_metric[idx_usr].astype(np.int8)
            for idx_sat, satellite in enumerate(sm.satellites):
                # Vectorised run-length pass detection; like the original per-epoch scan,
                # passes still ongoing at the last epoch are not counted
                transitions = np.diff(metric_int[:, idx_sat])
                ends = np.flatnonzero(transitions == -1) + 1  # First epoch after each completed pass
                if ends.size == 0:
                    continue
                starts = np.flatnonzero(transitions == 1) + 1  # First epoch of each pass
                if metric_int[0, idx_sat]:  # Pass already running at the first epoch
                    starts = np.insert(starts, 0, 0)
                # Runs alternate, so starts/ends pair up in order; +1 keeps the original
                # backward-scan convention (in-view epochs + 1)
                lengths = ends - starts[:ends.size] + 1
                valid_value_list.extend((lengths * time_step).tolist())  # Add pass lengths to the list for this user

            if len(valid_value_list) == 0:
                metric[idx_usr] = -1.0
            else:
                if self.statistic == "Min":
                    metric[idx_usr] = np.min(valid_value_list)
                if self.statistic == "Mean":
                    metric[idx_usr] = np.mean(valid_value_list)
                if self.statistic == "Max":
                    metric[idx_usr] = np.max(valid_value_list)
                if self.statistic == "Std":
                    metric[idx_usr] = np.std(valid_value_list)
                if self.statistic == "Median":
                    metric[idx_usr] = np.median(valid_value_list)
            lats.append(degrees(sm.users[idx_usr].lla[0]))
            lons.append(degrees(sm.users[idx_usr].lla[1]))

        self.write_csv(sm, ['lon_deg', 'lat_deg', f'{self.statistic.lower()}_pass_time_s'],
                       np.column_stack([lons, lats, metric]))
        grid_shape = get_user_grid_shape(sm, self.type)
        if grid_shape is None:
            return
        x_new = np.reshape(np.array(lons), grid_shape)
        y_new = np.reshape(np.array(lats), grid_shape)
        z_new = np.reshape(np.array(metric), grid_shape)
        fig, ax = make_map_cyl()
        im1 = map_pcolormesh(ax, x_new, y_new, z_new, cmap=plt.cm.jet)
        cb = plt.colorbar(im1, ax=ax, shrink=0.85, pad=0.02)
        cb.set_label(self.statistic + ' Pass Time Interval [s]', fontsize=10)
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        if self.plot_3d:
            self.render_3d_grid(sm, sm.user_latitudes, sm.user_longitudes, z_new,
                                self.statistic + ' Pass Time Interval [s]')

        if self.mp4:
            import plot_movie
            plot_movie.movie_grid_2d(sm, self.user_metric.sum(axis=2).astype(float),
                                     self.type, 'Satellites of constellation in view [-]',
                                     sm.output_path(self.type + '_2d.mp4'))
            self.render_movie_3d(sm, sm.satellites, self.sat_pos_hist_3d,
                                 grid=(sm.user_latitudes, sm.user_longitudes, z_new,
                                       self.statistic + ' Pass Time Interval [s]',
                                       'jet', None))


class AnalysisCovSatelliteContour(AnalysisBase, AnalysisPlot3D):

    def __init__(self):
        super().__init__()
        self.constellation_id = 0  # Mandatory
        self.satellite_id = 0  # Mandatory
        self.elevation_mask = 0  # Mandatory
        self.idx_found_satellite = 0
        self.idx_found_satellites = []
        self.init_3d()

    def read_config(self, node):
        if node.find('ConstellationID') is not None:
            self.constellation_id = int(node.find('ConstellationID').text)
        if node.find('SatelliteID') is not None:
            self.satellite_id = int(node.find('SatelliteID').text)
        if node.find('ElevationMask') is not None:
            self.elevation_mask = radians(float(node.find('ElevationMask').text))
        self.read_config_3d(node)

    def before_loop(self, sm):
        if self.satellite_id > 0:  # One satellite
            # Find the index of the satellite that is needed
            for idx_sat, satellite in enumerate(sm.satellites):
                if satellite.constellation_id == self.constellation_id and \
                        satellite.sat_id == self.satellite_id:
                    self.idx_found_satellite = idx_sat
                    break
        else:  # Whole constellation
            # Find the index of the satellites that is needed
            for idx_sat, satellite in enumerate(sm.satellites):
                if satellite.constellation_id == self.constellation_id:
                    self.idx_found_satellites.append(idx_sat)
        self.before_loop_3d(sm)

    def in_loop(self, sm):
        self.in_loop_3d(sm)

    def after_loop(self, sm):
        if self.satellite_id > 0:  # One satellite
            idx_selected = [self.idx_found_satellite]
        else:
            idx_selected = self.idx_found_satellites
        fig, ax = make_map_cyl()
        contours = []
        for idx_sat in idx_selected:
            sm.satellites[idx_sat].det_lla()
            contour = misc_fn.sat_contour(sm.satellites[idx_sat].lla, self.elevation_mask)
            contours.append(contour)
            if self.satellite_id > 0:
                ax.plot(contour[:, 1] / PI * 180, contour[:, 0] / PI * 180, 'r.',
                        transform=ccrs.PlateCarree())
            else:
                ax.plot(contour[:, 1] / PI * 180, contour[:, 0] / PI * 180, '.',
                        label='Satellite ID: '+str(sm.satellites[idx_sat].sat_id),
                        transform=ccrs.PlateCarree())
                ax.legend()
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        self.write_csv(sm, ['sat_id', 'lat_deg', 'lon_deg'],
                       np.vstack([np.column_stack([np.full(len(contour), sm.satellites[idx_sat].sat_id),
                                                   np.degrees(contour)])
                                  for idx_sat, contour in zip(idx_selected, contours)]))

        if self.plot_3d:
            selected = [sm.satellites[i] for i in idx_selected]
            self.render_3d_contours(sm, contours, satellites=selected)


class AnalysisCovSatelliteHighest(AnalysisBase, AnalysisPlot3D):

    def __init__(self):
        super().__init__()
        self.statistic = ''  # Mandatory
        self.constellation_id = 0  # Mandatory
        self.user_metric = None  # Per-user metric memory (num_user, num_epoch)
        self.init_3d()

    def read_config(self, node):
        if node.find('Statistic') is not None:
            self.statistic = node.find('Statistic').text
        if node.find('ConstellationID') is not None:
            self.constellation_id = int(node.find('ConstellationID').text)
        self.read_config_3d(node)

    def before_loop(self, sm):
        self.user_metric = np.zeros((sm.num_user, sm.num_epoch))
        self.before_loop_3d(sm)

    def in_loop(self, sm):
        for idx_user, user in enumerate(sm.users):
            best_satellite_value = -1
            for idx_sat in range(len(user.idx_sat_in_view)):
                if sm.satellites[user.idx_sat_in_view[idx_sat]].constellation_id == self.constellation_id:
                    elevation = degrees(sm.usr2sp[idx_user][user.idx_sat_in_view[idx_sat]].elevation)
                    if elevation > best_satellite_value:
                        best_satellite_value = elevation
            self.user_metric[idx_user, sm.cnt_epoch] = best_satellite_value
        self.in_loop_3d(sm)

    def after_loop(self, sm):
        metric, lats, lons = [], [], []
        for idx_usr, user in enumerate(sm.users):
            if self.statistic == 'Min':
                metric.append(np.min(self.user_metric[idx_usr]))
            if self.statistic == 'Mean':
                metric.append(np.mean(self.user_metric[idx_usr]))
            if self.statistic == 'Max':
                metric.append(np.max(self.user_metric[idx_usr]))
            if self.statistic == 'Std':
                metric.append(np.std(self.user_metric[idx_usr]))
            if self.statistic == 'Median':
                metric.append(np.median(self.user_metric[idx_usr]))
            lats.append(degrees(user.lla[0]))
            lons.append(degrees(user.lla[1]))
        self.write_csv(sm, ['lon_deg', 'lat_deg', f'{self.statistic.lower()}_max_elevation_deg'],
                       np.column_stack([lons, lats, metric]))
        grid_shape = get_user_grid_shape(sm, self.type)
        if grid_shape is None:
            return
        x_new = np.reshape(np.array(lons), grid_shape)
        y_new = np.reshape(np.array(lats), grid_shape)
        z_new = np.reshape(np.array(metric), grid_shape)
        fig, ax = make_map_cyl()
        im1 = map_pcolormesh(ax, x_new, y_new, z_new, cmap=plt.cm.jet)
        cb = plt.colorbar(im1, ax=ax, shrink=0.85, pad=0.02)
        cb.set_label(self.statistic + ' of Max Elevation satellites in view [deg]', fontsize=10)
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        if self.plot_3d:
            self.render_3d_grid(sm, sm.user_latitudes, sm.user_longitudes, z_new,
                                self.statistic + ' Max Elevation in view [deg]')

        if self.mp4:
            import plot_movie
            plot_movie.movie_grid_2d(sm, self.user_metric, self.type,
                                     'Max elevation in view [deg]',
                                     sm.output_path(self.type + '_2d.mp4'))
            self.render_movie_3d(sm, sm.satellites, self.sat_pos_hist_3d,
                                 grid=(sm.user_latitudes, sm.user_longitudes, z_new,
                                       self.statistic + ' Max Elevation in view [deg]',
                                       'jet', None))


class AnalysisCovSatellitePvt(AnalysisBase):

    def __init__(self):
        super().__init__()
        self.constellation_id = 0  # Mandatory
        self.satellite_id = 0  # Mandatory
        self.file_orbits = None

    def read_config(self, node):
        if node.find('ConstellationID') is not None:
            self.constellation_id = int(node.find('ConstellationID').text)
        if node.find('SatelliteID') is not None:
            self.satellite_id = int(node.find('SatelliteID').text)

    def before_loop(self, sm):
        # One file handle for the whole run; opening per satellite per epoch dominated the loop
        self.file_orbits = open(sm.output_path('orbits.txt'), 'w')

    def in_loop(self, sm):
        for satellite in sm.satellites:
            if satellite.constellation_id != self.constellation_id:
                continue
            if self.satellite_id > 0 and satellite.sat_id != self.satellite_id:  # Only one satellite
                continue
            self.file_orbits.write("%13.6f,%d,%13.6f,%13.6f,%13.6f,%13.6f,%13.6f,%13.6f\n"
                                   % (sm.time_mjd, satellite.sat_id,
                                      satellite.pos_eci[0], satellite.pos_eci[1], satellite.pos_eci[2],
                                      satellite.vel_eci[0], satellite.vel_eci[1], satellite.vel_eci[2]))

    def after_loop(self, sm):
        self.file_orbits.close()
        self.file_orbits = None
        if os.path.getsize(sm.output_path('orbits.txt')) == 0:
            ls.logger.error(f'No satellite matched ConstellationID {self.constellation_id} / '
                            f'SatelliteID {self.satellite_id}, nothing recorded. Available: ' +
                            ', '.join(f'{s.constellation_id}/{s.sat_id}' for s in sm.satellites))
            return
        data = pd.read_csv(sm.output_path('orbits.txt'), sep=',', header=None,
                           names=['RunTime', 'ID', 'x', 'y', 'z', 'x_vel', 'y_vel', 'z_vel'])
        # Plot the configured satellite, or the first recorded one (satellite ids are
        # NORAD catalog numbers for TLE-defined constellations, not 1-based)
        sat_id = self.satellite_id if self.satellite_id > 0 else int(data.ID.iloc[0])
        data2 = data[data.ID == sat_id]
        fig, ax1 = plt.subplots(figsize=(10, 6))
        plt.grid()
        ax1.set_ylabel('Position ECI [m]')
        ax1.plot(self.times_f_doy, data2.x, 'r+-', label='x_pos')
        ax1.plot(self.times_f_doy, data2.y, 'g+-', label='y_pos')
        ax1.plot(self.times_f_doy, data2.z, 'b+-', label='z_pos')
        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
        ax2.set_ylabel('Velocity ECI [m]')  # we already handled the x-label with ax1
        ax2.plot(self.times_f_doy, data2.x_vel, 'm+-', label='x_vel')
        ax2.plot(self.times_f_doy, data2.y_vel, 'y+-', label='y_vel')
        ax2.plot(self.times_f_doy, data2.z_vel, 'k+-', label='z_vel')
        ax1.legend(loc=2); ax2.legend(loc=0)
        plt.xlabel('DOY[-]'); fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        self.write_csv(sm, ['mjd', 'sat_id', 'x_m', 'y_m', 'z_m',
                            'x_vel_ms', 'y_vel_ms', 'z_vel_ms'], data2.values)


class AnalysisCovSatelliteSkyAngles(AnalysisBase):

    def __init__(self):
        super().__init__()
        self.constellation_id = 0  # Mandatory
        self.satellite_id = 0  # Mandatory
        self.idx_found_satellite = 0
        self.user_metric = None  # Metric memory for the first user (num_epoch, 2)

    def read_config(self, node):
        if node.find('ConstellationID') is not None:
            self.constellation_id = int(node.find('ConstellationID').text)
        if node.find('SatelliteID') is not None:
            self.satellite_id = int(node.find('SatelliteID').text)

    def before_loop(self, sm):
        # Find the index of the satellite that is needed
        for i, satellite in enumerate(sm.satellites):
            if satellite.constellation_id == self.constellation_id and \
                    satellite.sat_id == self.satellite_id:
                self.idx_found_satellite = i
                break
        self.user_metric = np.zeros((sm.num_epoch, 2))

    def in_loop(self, sm):
        # Sky angles are plotted for the first user only
        link = sm.usr2sp[0][self.idx_found_satellite]
        if link.elevation > 0:
            self.user_metric[sm.cnt_epoch, 0] = degrees(link.azimuth)
            self.user_metric[sm.cnt_epoch, 1] = degrees(link.elevation)

    def after_loop(self, sm):
        fig, ax1 = plt.subplots(figsize=(10, 6))
        plt.grid()
        plt.subplots_adjust(left=.1, right=.92, top=0.95, bottom=0.07)
        ax1.set_ylabel('Azimuth [deg]')
        ax1.yaxis.label.set_color('red')
        ax1.tick_params(axis='y', colors='red')
        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
        ax2.set_ylabel('Elevation [deg]')
        ax2.yaxis.label.set_color('blue')
        ax2.tick_params(axis='y', colors='blue')
        # Sky angles are plotted for the first user only
        ax1.plot(self.times_f_doy, self.user_metric[:, 0], 'r+', label='Azimuth')
        ax2.plot(self.times_f_doy, self.user_metric[:, 1], 'b+', label='Elevation')
        plt.xlabel('DOY[-]');
        ax1.legend(loc='upper left'); ax2.legend(loc='upper right')
        plt.grid()
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        self.write_csv(sm, ['doy', 'azimuth_deg', 'elevation_deg'],
                       np.column_stack([self.times_f_doy, self.user_metric]))


class AnalysisCovSatelliteVisible(AnalysisBase):

    def __init__(self):
        super().__init__()
        self.user_metric = None  # Per-user metric memory (num_user, num_epoch)

    def read_config(self, node):
        pass

    def before_loop(self, sm):
        self.user_metric = np.zeros((sm.num_user, sm.num_epoch))

    def in_loop(self, sm):
        for idx_user, user in enumerate(sm.users):
            self.user_metric[idx_user, sm.cnt_epoch] = len(user.idx_sat_in_view)

    def after_loop(self, sm):
        fig = plt.figure(figsize=(10, 6))
        plt.subplots_adjust(left=.1, right=.95, top=0.95, bottom=0.07)
        for idx_user, user in enumerate(sm.users):
            plt.plot(self.times_f_doy, self.user_metric[idx_user], '-',
                     label=f'User lat/lon {round(degrees(user.lla[0]),1)} {round(degrees(user.lla[1]),1)}')
        plt.xlabel('DOY[-]'); plt.ylabel('Number of satellites in view [-]')
        plt.grid(); plt.legend()
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        self.write_csv(sm, ['doy'] + [f'user_{user.user_id}' for user in sm.users],
                       np.column_stack([self.times_f_doy, self.user_metric.T]))


class AnalysisCovSatelliteVisibleGrid(AnalysisBase, AnalysisPlot3D):

    def __init__(self):
        super().__init__()
        self.statistic = ''  # Mandatory
        self.user_metric = None  # Per-user metric memory (num_user, num_epoch)
        self.init_3d()

    def read_config(self, node):
        if node.find('Statistic') is not None:
            self.statistic = node.find('Statistic').text
        self.read_config_3d(node)

    def before_loop(self, sm):
        self.user_metric = np.zeros((sm.num_user, sm.num_epoch))
        self.before_loop_3d(sm)

    def in_loop(self, sm):
        for idx_user, user in enumerate(sm.users):
            self.user_metric[idx_user, sm.cnt_epoch] = len(user.idx_sat_in_view)
        self.in_loop_3d(sm)

    def after_loop(self, sm):
        metric, latitudes, longitudes = [], [], []
        for idx_user, user in enumerate(sm.users):
            if self.statistic == 'Min':
                metric.append(np.min(self.user_metric[idx_user]))
            if self.statistic == 'Mean':
                metric.append(np.mean(self.user_metric[idx_user]))
            if self.statistic == 'Max':
                metric.append(np.max(self.user_metric[idx_user]))
            if self.statistic == 'Std':
                metric.append(np.std(self.user_metric[idx_user]))
            if self.statistic == 'Median':
                metric.append(np.median(self.user_metric[idx_user]))
            latitudes.append(degrees(user.lla[0]))
            longitudes.append(degrees(user.lla[1]))
        self.write_csv(sm, ['lon_deg', 'lat_deg', f'{self.statistic.lower()}_satellites_in_view'],
                       np.column_stack([longitudes, latitudes, metric]))
        grid_shape = get_user_grid_shape(sm, self.type)
        if grid_shape is None:
            return
        x_new = np.reshape(np.array(longitudes), grid_shape)
        y_new = np.reshape(np.array(latitudes), grid_shape)
        z_new = np.reshape(np.array(metric), grid_shape)
        fig, ax = make_map_cyl()
        im1 = map_pcolormesh(ax, x_new, y_new, z_new, cmap=plt.cm.jet)
        cb = plt.colorbar(im1, ax=ax, shrink=0.85, pad=0.02)
        cb.set_label(self.statistic + ' Number of satellites in view [-]', fontsize=10)
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        if self.plot_3d:
            self.render_3d_grid(sm, sm.user_latitudes, sm.user_longitudes, z_new,
                                self.statistic + ' Satellites in view [-]')

        if self.mp4:
            import plot_movie
            plot_movie.movie_grid_2d(sm, self.user_metric, self.type,
                                     'Number of satellites in view [-]',
                                     sm.output_path(self.type + '_2d.mp4'))
            self.render_movie_3d(sm, sm.satellites, self.sat_pos_hist_3d,
                                 grid=(sm.user_latitudes, sm.user_longitudes, z_new,
                                       self.statistic + ' Satellites in view [-]',
                                       'jet', None))


class AnalysisCovSatelliteVisibleId(AnalysisBase):

    def __init__(self):
        super().__init__()
        self.constellation_id = 0  # Mandatory
        self.user_metric = None  # Metric memory for the first user (num_epoch, num_sat)

    def read_config(self, node):
        if node.find('ConstellationID') is not None:
            self.constellation_id = int(node.find('ConstellationID').text)

    def before_loop(self, sm):
        # IDs in view are recorded for the first user only
        self.user_metric = np.ones((sm.num_epoch, sm.num_sat)) * np.nan

    def in_loop(self, sm):
        for idx_sat in range(len(sm.users[0].idx_sat_in_view)):
            if sm.satellites[sm.users[0].idx_sat_in_view[idx_sat]].constellation_id == self.constellation_id:
                self.user_metric[sm.cnt_epoch, idx_sat] = sm.satellites[sm.users[0].idx_sat_in_view[idx_sat]].sat_id

    def after_loop(self, sm):
        fig = plt.figure(figsize=(10, 6))
        plt.subplots_adjust(left=.1, right=.95, top=0.95, bottom=0.07)
        plt.plot(self.times_f_doy, self.user_metric, 'r+')
        plt.xlabel('DOY[-]'); plt.ylabel('IDs of satellites in view [-]')
        plt.grid()
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        self.write_csv(sm, ['doy'] + [f'slot_{i}' for i in range(sm.num_sat)],
                       np.column_stack([self.times_f_doy, self.user_metric]))


