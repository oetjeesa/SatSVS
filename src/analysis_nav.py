import os
import numpy as np
import pandas as pd
from math import degrees, radians
from matplotlib import pyplot as plt
# Import project modules
import misc_fn
from constants import PI
from analysis import AnalysisBase, AnalysisPlot3D, make_map_cyl, map_pcolormesh, get_user_grid_shape


class AnalysisNavDOP(AnalysisBase, AnalysisPlot3D):

    def __init__(self):
        super().__init__()
        self.direction = None  # Mandatory
        self.statistic = None  # Mandatory
        self.user_metric = None  # Per-user metric memory (num_user, num_epoch)
        self.init_3d()

    def read_config(self, node):
        if node.find('Direction') is not None:
            self.direction = node.find('Direction').text
        if node.find('Statistic') is not None:
            self.statistic = node.find('Statistic').text
        self.read_config_3d(node)

    def before_loop(self, sm):
        self.user_metric = np.zeros((sm.num_user, sm.num_epoch))
        self.before_loop_3d(sm)

    def in_loop(self, sm):
        for user_idx, user in enumerate(sm.users):
            if len(user.idx_sat_in_view) > 3:
                h_mat = np.ones((len(user.idx_sat_in_view),4))
                for i, idx_sat in enumerate(user.idx_sat_in_view):
                    h_mat[i,0:3] = sm.usr2sp[user_idx][idx_sat].usr2sp_ecf/sm.usr2sp[user_idx][idx_sat].distance
                hth_inv = misc_fn.inverse4by4(np.matmul(np.transpose(h_mat),h_mat))  # Fast implementation
                #hth_inv2 = np.linalg.inv(np.matmul(np.transpose(h_mat),h_mat))
                q = misc_fn.ecef2enu(hth_inv,user.lla[0],user.lla[1])
                if self.direction == "Pos":
                    dop = np.sqrt(q[0, 0] + q[1, 1] + q[2, 2])
                elif self.direction == "Hor":
                    dop = np.sqrt(q[0, 0] + q[1, 1])
                elif self.direction == "Ver":
                    dop = np.sqrt(q[2, 2])
                self.user_metric[user_idx, sm.cnt_epoch] = dop
            else:
                self.user_metric[user_idx, sm.cnt_epoch] = np.nan
        self.in_loop_3d(sm)

    def after_loop(self, sm):
        lats, lons = [], []
        metric = np.zeros(len(sm.users))
        for idx_usr, user in enumerate(sm.users):
            if self.statistic == "Min":
                metric[idx_usr] = np.nanmin(self.user_metric[idx_usr])
            if self.statistic == "Mean":
                metric[idx_usr] = np.nanmean(self.user_metric[idx_usr])
            if self.statistic == "Max":
                metric[idx_usr] = np.nanmax(self.user_metric[idx_usr])
            if self.statistic == "Std":
                metric[idx_usr] = np.nanstd(self.user_metric[idx_usr])
            if self.statistic == "Median":
                metric[idx_usr] = np.nanmedian(self.user_metric[idx_usr])
            lats.append(degrees(sm.users[idx_usr].lla[0]))
            lons.append(degrees(sm.users[idx_usr].lla[1]))

        self.write_csv(sm, ['lon_deg', 'lat_deg',
                            f'{self.statistic.lower()}_{self.direction.lower()}_dop'],
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
        cb.set_label(self.statistic + ' ' + self.direction + ' Dilution of Precision [-]', fontsize=10)
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        if self.plot_3d:
            self.render_3d_grid(sm, sm.user_latitudes, sm.user_longitudes, z_new,
                                self.statistic + ' ' + self.direction + ' DOP [-]')

        if self.mp4:
            import plot_movie
            plot_movie.movie_grid_2d(sm, self.user_metric, self.type,
                                     self.direction + ' Dilution of Precision [-]',
                                     sm.output_path(self.type + '_2d.mp4'))
            self.render_movie_3d(sm, sm.satellites, self.sat_pos_hist_3d,
                                 grid=(sm.user_latitudes, sm.user_longitudes, z_new,
                                       self.statistic + ' ' + self.direction + ' DOP [-]',
                                       'jet', None))


class AnalysisNavAccuracy(AnalysisBase, AnalysisPlot3D):

    def __init__(self):
        super().__init__()
        self.direction = None  # Mandatory
        self.statistic = None  # Mandatory
        self.user_metric = None  # Per-user metric memory (num_user, num_epoch)
        self.init_3d()

    def read_config(self, node):
        if node.find('Direction') is not None:
            self.direction = node.find('Direction').text
        if node.find('Statistic') is not None:
            self.statistic = node.find('Statistic').text
        self.read_config_3d(node)

    def before_loop(self, sm):
        self.user_metric = np.zeros((sm.num_user, sm.num_epoch))
        self.before_loop_3d(sm)

    def in_loop(self, sm):
        for user_idx, user in enumerate(sm.users):
            sum_user_error = 0
            if len(user.idx_sat_in_view) > 3:
                h_mat = np.ones((len(user.idx_sat_in_view),4))
                for i, idx_sat in enumerate(user.idx_sat_in_view):
                    h_mat[i,0:3] = sm.usr2sp[user_idx][idx_sat].usr2sp_ecf/sm.usr2sp[user_idx][idx_sat].distance
                    num_uere = len(sm.satellites[idx_sat].uere_list)
                    el_piece_angle = PI/2/num_uere
                    idx_el = int(np.floor(sm.usr2sp[user_idx][idx_sat].elevation/el_piece_angle))
                    uere = np.power(sm.satellites[idx_sat].uere_list[idx_el],2)
                    sum_user_error += uere
                hth_inv = misc_fn.inverse4by4(np.matmul(np.transpose(h_mat),h_mat))  # Fast implementation
                #hth_inv = np.linalg.inv(np.matmul(np.transpose(h_mat),h_mat))
                q = misc_fn.ecef2enu(hth_inv,user.lla[0],user.lla[1])
                error = 2 * np.sqrt(sum_user_error/len(user.idx_sat_in_view))
                if self.direction == "Pos":
                    acc = error * np.sqrt(q[0, 0] + q[1, 1] + q[2, 2])
                elif self.direction == "Hor":
                    acc = error * np.sqrt(q[0, 0] + q[1, 1])
                elif self.direction == "Ver":
                    acc = error * np.sqrt(q[2, 2])
                self.user_metric[user_idx, sm.cnt_epoch] = acc
            else:
                self.user_metric[user_idx, sm.cnt_epoch] = np.nan
        self.in_loop_3d(sm)

    def after_loop(self, sm):
        lats, lons = [], []
        metric = np.zeros(len(sm.users))
        for idx_usr, user in enumerate(sm.users):
            if self.statistic == "Min":
                metric[idx_usr] = np.nanmin(self.user_metric[idx_usr])
            if self.statistic == "Mean":
                metric[idx_usr] = np.nanmean(self.user_metric[idx_usr])
            if self.statistic == "Max":
                metric[idx_usr] = np.nanmax(self.user_metric[idx_usr])
            if self.statistic == "Std":
                metric[idx_usr] = np.nanstd(self.user_metric[idx_usr])
            if self.statistic == "Median":
                metric[idx_usr] = np.nanmedian(self.user_metric[idx_usr])
            lats.append(degrees(sm.users[idx_usr].lla[0]))
            lons.append(degrees(sm.users[idx_usr].lla[1]))

        self.write_csv(sm, ['lon_deg', 'lat_deg',
                            f'{self.statistic.lower()}_{self.direction.lower()}_accuracy_m'],
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
        cb.set_label(self.statistic + ' ' + self.direction + ' Navigation Accuracy 95% [m]', fontsize=10)
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        if self.plot_3d:
            self.render_3d_grid(sm, sm.user_latitudes, sm.user_longitudes, z_new,
                                self.statistic + ' ' + self.direction + ' Nav Accuracy 95% [m]')

        if self.mp4:
            import plot_movie
            plot_movie.movie_grid_2d(sm, self.user_metric, self.type,
                                     self.direction + ' Navigation Accuracy 95% [m]',
                                     sm.output_path(self.type + '_2d.mp4'))
            self.render_movie_3d(sm, sm.satellites, self.sat_pos_hist_3d,
                                 grid=(sm.user_latitudes, sm.user_longitudes, z_new,
                                       self.statistic + ' ' + self.direction + ' Nav Accuracy 95% [m]',
                                       'jet', None))



