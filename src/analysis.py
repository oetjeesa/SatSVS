import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy.crs as ccrs
from math import sin, cos, asin, degrees, radians
import numpy as np

# Modules from project
import logging_svs as ls


def make_map_cyl(figsize=(10, 5)):
    """
    Global cylindrical (PlateCarree) map with gridlines and coastlines,
    the cartopy equivalent of the previous Basemap(projection='cyl') setup.
    Data plotted on the returned axes must pass transform=ccrs.PlateCarree().
    """
    fig = plt.figure(figsize=figsize)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()
    gl = ax.gridlines(draw_labels=True, xlocs=np.arange(-180., 181., 60.),
                      ylocs=np.arange(-90., 91., 30.), linewidth=0.5, color='gray')
    gl.top_labels = False
    gl.right_labels = False
    ax.coastlines()
    return fig, ax


def make_map_polar(polar_view, figsize=(7, 6)):
    """
    North (polar_view > 0) or south polar stereographic map bounded at latitude
    polar_view with a circular boundary, the cartopy equivalent of the previous
    Basemap npstere/spstere setup.
    """
    fig = plt.figure(figsize=figsize)
    if polar_view > 0:
        ax = plt.axes(projection=ccrs.NorthPolarStereo())
        ax.set_extent([-180, 180, polar_view, 90], ccrs.PlateCarree())
    else:
        ax = plt.axes(projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180, 180, -90, polar_view], ccrs.PlateCarree())
    theta = np.linspace(0, 2 * np.pi, 200)
    circle = mpath.Path(np.vstack([np.sin(theta), np.cos(theta)]).T * 0.5 + [0.5, 0.5])
    ax.set_boundary(circle, transform=ax.transAxes)
    ax.gridlines(xlocs=np.arange(-180., 181., 60.), ylocs=np.arange(-90., 91., 30.),
                 linewidth=0.5, color='gray')
    ax.coastlines()
    return fig, ax


def map_pcolormesh(ax, x, y, z, **kwargs):
    """
    pcolormesh of a lon/lat grid on a cartopy map. vmin/vmax are pinned at draw
    time: cartopy draws cells straddling the +/-180 map edge as a separate
    collection which does not follow a later colorbar re-normalisation, so an
    unpinned norm leaves dark strips at the map edges when all values are equal.
    """
    if 'vmin' not in kwargs and 'vmax' not in kwargs:
        zmin, zmax = np.nanmin(z), np.nanmax(z)
        if np.isfinite(zmin) and np.isfinite(zmax):
            if zmin == zmax:  # degenerate range, expand like matplotlib does
                zmin, zmax = zmin - 0.1, zmax + 0.1
            kwargs['vmin'] = zmin
            kwargs['vmax'] = zmax
    return ax.pcolormesh(x, y, z, shading='auto', transform=ccrs.PlateCarree(), **kwargs)


class AnalysisBase:

    def __init__(self):

        self.times_mjd = []  # Time list of simulation
        self.times_f_doy = []  # Time list of simulation

        self.type = ''

    def read_config(self, node):  # Node in xml element tree
        pass

    def before_loop(self, sm):
        pass

    def in_loop(self, sm):
        pass

    def after_loop(self, sm):
        pass


class AnalysisObs:   # Common methods needed for some OBS analysis

    def plot_swath_coverage(self, sm, user_metric, polar_view):
        plot_points = np.zeros((len(sm.users), 3))
        for idx_user, user in enumerate(sm.users):
            if idx_user % 1000 == 0:
                ls.logger.info(f'User swath coverage {user.user_id} of {len(sm.users)}')
            if user_metric[idx_user, :].any():  # Any value bigger than 0
                num_swaths = len(np.flatnonzero(np.diff(user_metric[idx_user, :]))) / 2
                if num_swaths >= 1:
                    plot_points[idx_user, :] = [degrees(user.lla[1]), degrees(user.lla[0]), num_swaths]
        plot_points = plot_points[~np.all(plot_points == 0, axis=1)]  # Clean up empty rows
        if polar_view is not None:
            fig, ax = make_map_polar(polar_view)
        else:
            fig, ax = make_map_cyl()
        sc = ax.scatter(plot_points[:,0], plot_points[:,1], s=3, marker='o', cmap=plt.cm.jet,
                        c=plot_points[:,2], alpha=.3, transform=ccrs.PlateCarree())
        cb = plt.colorbar(sc, ax=ax, shrink=0.85)
        cb.set_label('Number of passes [-]', fontsize=10)
        plt.subplots_adjust(left=.1, right=.9, top=0.92, bottom=0.1)
        plt.savefig('../output/'+self.type+'.png')
        plt.show()

    def plot_swath_revisit(self, sm, user_metric, statistic, polar_view):
        plot_points = np.zeros((len(sm.users), 3))
        metric = 0
        for idx_user, user in enumerate(sm.users):
            if idx_user % 1000 == 0:
                ls.logger.info(f'User revisit {user.user_id} of {len(sm.users)}')
            gaps = np.diff(np.where(np.diff(user_metric[idx_user,:]) != 0)).flatten()
            if len(gaps) > 1:
                gaps = np.delete(gaps, np.where(gaps == 1), axis=0) * sm.time_step / 3600.0
                if statistic == "min":
                    metric = (np.nanmin(gaps))
                if statistic == "mean":
                    metric = (np.nanmean(gaps))
                if statistic == "max":
                    metric = (np.nanmax(gaps))
                if statistic == "std":
                    metric = (np.nanstd(gaps))
                if statistic == "median":
                    metric = (np.nanmedian(gaps))
                if polar_view is not None:
                    if np.abs(degrees(sm.users[idx_user].lla[0])) > np.abs(polar_view):
                        plot_points[idx_user, :] = [degrees(sm.users[idx_user].lla[1]), degrees(sm.users[idx_user].lla[0]), metric]
                else:
                    plot_points[idx_user, :] = [degrees(sm.users[idx_user].lla[1]), degrees(sm.users[idx_user].lla[0]), metric]
        plot_points = plot_points[~np.all(plot_points == 0, axis=1)]  # Clean up empty rows
        if metric>0:  # Only plot if not empty
            if polar_view is not None:
                fig, ax = make_map_polar(polar_view)
                point_size = 5
            else:
                fig, ax = make_map_cyl()
                point_size = 3
            sc = ax.scatter(plot_points[:,0], plot_points[:,1], s=point_size, cmap=plt.cm.jet,
                            c=plot_points[:,2], vmin=np.nanmin(plot_points[:,2]),
                            vmax=np.nanmax(plot_points[:,2]), transform=ccrs.PlateCarree())
            cb = plt.colorbar(sc, ax=ax, shrink=0.85)
            cb.set_label(statistic.capitalize() + ' Revisit Interval [hours]', fontsize=10)
            plt.subplots_adjust(left=.1, right=.9, top=0.9, bottom=0.1)
            plt.savefig(f'../output/{self.type}_revisit.png')
            plt.show()

