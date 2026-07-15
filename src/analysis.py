import os

import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patheffects as patheffects
import cartopy.crs as ccrs
from math import sin, cos, asin, degrees, radians
import numpy as np

# Modules from project
import logging_svs as ls
import misc_fn


def make_map_cyl(figsize=(10, 5.6)):
    """
    Global cylindrical (PlateCarree) map with gridlines and coastlines,
    the cartopy equivalent of the previous Basemap(projection='cyl') setup.
    Data plotted on the returned axes must pass transform=ccrs.PlateCarree().
    Constrained layout with the figure height a bit above the fixed 2:1 map
    aspect leaves a comfortable white margin above and below the map.
    """
    fig = plt.figure(figsize=figsize, layout='constrained')
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()
    gl = ax.gridlines(draw_labels=True, xlocs=np.arange(-180., 181., 60.),
                      ylocs=np.arange(-90., 91., 30.), linewidth=0.5, color='gray')
    gl.top_labels = False
    gl.right_labels = False
    ax.coastlines()
    return fig, ax


def make_map_polar(polar_view, figsize=(7, 5.8)):
    """
    North (polar_view > 0) or south polar stereographic map bounded at latitude
    polar_view with a circular boundary, the cartopy equivalent of the previous
    Basemap npstere/spstere setup.
    """
    fig = plt.figure(figsize=figsize, layout='constrained')
    # PolarView 90/-90 (bounding latitude at the pole, as documented in readme.md)
    # would give a zero-area extent that crashes cartopy at draw time; show the
    # full hemisphere instead, matching the old Basemap tolerance for this value.
    bounding_lat = polar_view if abs(polar_view) < 90 else 0
    if polar_view > 0:
        ax = plt.axes(projection=ccrs.NorthPolarStereo())
        ax.set_extent([-180, 180, bounding_lat, 90], ccrs.PlateCarree())
    else:
        ax = plt.axes(projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180, 180, -90, bounding_lat], ccrs.PlateCarree())
    theta = np.linspace(0, 2 * np.pi, 200)
    circle = mpath.Path(np.vstack([np.sin(theta), np.cos(theta)]).T * 0.5 + [0.5, 0.5])
    ax.set_boundary(circle, transform=ax.transAxes)
    ax.gridlines(xlocs=np.arange(-180., 181., 60.), ylocs=np.arange(-90., 91., 30.),
                 linewidth=0.5, color='gray')
    ax.coastlines()
    return fig, ax


def get_user_grid_shape(sm, analysis_type):
    """
    (num_lat, num_lon) shape of the user grid for the map-statistics analyses,
    or None (with a clear error log) when the configured user segment is not a
    Grid — those analyses reshape their per-user metric into a lat/lon grid.
    """
    num_lat = sm.users[0].num_lat if sm.users else 0
    num_lon = sm.users[0].num_lon if sm.users else 0
    if num_lat == 0 or num_lat * num_lon != len(sm.users):
        ls.logger.error(f'Analysis {analysis_type} needs a user segment of Type Grid, '
                        f'found {len(sm.users)} user(s) not forming a lat/lon grid. '
                        f'No plot produced.')
        return None
    return num_lat, num_lon


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


class AnalysisMap2D:
    """Mixin of AnalysisBase: shared optional decorations of the 2D world-map
    plots, available on every world-map analysis (all default off):

    <EarthImage>True</EarthImage>            Earth photo as the map background
        (input/earth_texture.jpg, cartopy stock image as fallback)
    <ShowStations>True</ShowStations>        ground stations as red triangles
        (drawn underneath ground tracks) with halo-outlined name labels
    <ShowUsers>True</ShowUsers>              user locations as blue triangles
        (skipped with a warning for large user grids)
    <ShowGroundTrack>True</ShowGroundTrack>  subsatellite ground track(s),
        recorded per epoch, as thin grey lines
    <Coastlines>False</Coastlines>           turn the coastline outlines off
        (default on; e.g. with EarthImage the photo already shows the coasts)

    The config reading and the loop hooks are wired centrally (config.py
    load_simulation, main.py run loops); a map-drawing analysis only calls
    decorate_map2d(sm, ax) once per created map."""

    EARTH_TEXTURE = '../input/earth_texture.jpg'

    def init_map2d(self):
        self.show_stations = False
        self.show_users = False
        self.earth_image = False
        self.show_ground_track = False
        self.coastlines = True
        self._map2d_track = None  # (num_sat, num_epoch, 2) lat/lon [deg]

    def read_config_map2d(self, node):
        if node.find('ShowStations') is not None:
            self.show_stations = misc_fn.str2bool(node.find('ShowStations').text)
        if node.find('ShowUsers') is not None:
            self.show_users = misc_fn.str2bool(node.find('ShowUsers').text)
        if node.find('EarthImage') is not None:
            self.earth_image = misc_fn.str2bool(node.find('EarthImage').text)
        if node.find('ShowGroundTrack') is not None:
            self.show_ground_track = misc_fn.str2bool(node.find('ShowGroundTrack').text)
        if node.find('Coastlines') is not None:
            self.coastlines = misc_fn.str2bool(node.find('Coastlines').text)

    def before_loop_map2d(self, sm):
        if self.show_ground_track:
            self._map2d_track = np.full((sm.num_sat, sm.num_epoch, 2), np.nan)

    def in_loop_map2d(self, sm):
        if self._map2d_track is not None:
            for idx_sat, satellite in enumerate(sm.satellites):
                satellite.det_lla()
                self._map2d_track[idx_sat, sm.cnt_epoch] = [degrees(satellite.lla[0]),
                                                            degrees(satellite.lla[1])]

    def decorate_map2d(self, sm, ax):
        """Draw the enabled decorations on a 2D map; call once right after
        the map is created. The Earth image goes under the analysis data
        (zorder 0), the station/user triangles above the data fields but
        underneath the ground track and analysis lines, and the station name
        labels (halo-outlined for readability on any background) on top."""
        if not self.coastlines:
            import cartopy.mpl.feature_artist as feature_artist
            for artist in [c for c in ax.collections
                           if isinstance(c, feature_artist.FeatureArtist)]:
                artist.remove()
        if self.earth_image:
            texture = misc_fn.resolve_path(self.EARTH_TEXTURE)
            if texture is not None and os.path.isfile(texture):
                ax.imshow(plt.imread(texture), extent=(-180, 180, -90, 90),
                          transform=ccrs.PlateCarree(), origin='upper', zorder=0)
            else:
                ls.logger.warning(f'EarthImage: texture {self.EARTH_TEXTURE} not '
                                  f'found, using the cartopy stock image')
                ax.stock_img()
        if self._map2d_track is not None:
            for idx_sat in range(self._map2d_track.shape[0]):
                lat = self._map2d_track[idx_sat, :, 0]
                lon = self._map2d_track[idx_sat, :, 1]
                used = np.isfinite(lat)
                if not used.any():
                    continue
                x, y = misc_fn.track_with_gaps(lon[used], lat[used])
                # zorder 3: above the data fields and the station triangles,
                # underneath marker-based analysis data (e.g. imaging targets)
                ax.plot(x, y, '-', color='0.35', linewidth=0.7,
                        transform=ccrs.PlateCarree(), zorder=3)
        if self.show_stations:
            for station in sm.stations:
                lon, lat = degrees(station.lla[1]), degrees(station.lla[0])
                # zorder 1.8: above the data fields (pcolormesh/patches at 1)
                # but underneath the ground track and analysis lines (2+)
                ax.plot(lon, lat, '^', color='red', markersize=9,
                        markeredgecolor='black', transform=ccrs.PlateCarree(),
                        zorder=1.8)
                ax.text(lon + 2.0, lat + 2.0, station.station_name, fontsize=8,
                        color='red', transform=ccrs.PlateCarree(), zorder=7,
                        path_effects=[patheffects.withStroke(linewidth=2,
                                                             foreground='white')])
        if self.show_users:
            if len(sm.users) > 200:
                ls.logger.warning(f'ShowUsers: skipped, the user segment has '
                                  f'{len(sm.users)} points (a grid would flood the map)')
            else:
                for user in sm.users:
                    ax.plot(degrees(user.lla[1]), degrees(user.lla[0]), '^',
                            color='blue', markersize=7, markeredgecolor='black',
                            transform=ccrs.PlateCarree(), zorder=1.8)


class AnalysisBase(AnalysisMap2D):

    def __init__(self):

        self.times_mjd = []  # Time list of simulation
        self.times_f_doy = []  # Time list of simulation

        self.type = ''
        self.init_map2d()  # Shared 2D-map decoration flags (AnalysisMap2D)

    def read_config(self, node):  # Node in xml element tree
        pass

    def before_loop(self, sm):
        pass

    def in_loop(self, sm):
        pass

    def after_loop(self, sm):
        pass

    def write_csv(self, sm, columns, data, suffix=None):
        """Numeric dump of the plotted metric next to the PNG, as
        <type>[_suffix].csv with a one-line header. Every analysis dumps the
        data behind its plot so results can be post-processed without a rerun
        (and the regression tests compare numbers instead of pixels)."""
        name = self.type + (f'_{suffix}' if suffix else '') + '.csv'
        data = np.atleast_2d(np.asarray(data, dtype=float))
        if data.size == 0:  # Header-only file: the run produced no data points
            data = data.reshape(0, len(columns))
        np.savetxt(sm.output_path(name), data, delimiter=',', fmt='%.10g',
                   header=','.join(columns), comments='')
        ls.logger.info(f'Written data dump {name} ({data.shape[0]} rows)')


class AnalysisPlot3D:
    """Mixin for analyses that offer a 3D globe version of their world map
    (<Plot3D>True</Plot3D>, rendered by plot_3d.py with pyvista).

    Provides the shared configuration flags (ShowSatellite, ShowOrbit,
    SatelliteModelFile, SatelliteModelScale), the per-epoch recording of the
    satellite positions for the orbit tracks, and guarded wrappers around the
    plot_3d renderers. Call init_3d() in __init__, read_config_3d(node) in
    read_config, before_loop_3d(sm) in before_loop and in_loop_3d(sm) in
    in_loop; then use one of the render_3d_* methods in after_loop."""

    def init_3d(self):
        self.plot_3d = False  # Optional 3D globe plot (needs pyvista)
        self.show_satellite = True  # Draw the 3D satellite model(s)
        self.show_orbit = True  # Draw the orbital track(s)
        self.model_file = None  # Optional satellite mesh (STL/OBJ/PLY/VTK)
        self.model_scale = 200e3  # Satellite model size in m (exaggerated)
        self.mp4 = False  # Optional MP4 movies of the 2D and 3D world maps
        self.earth_clouds = False  # Cloud layer on the 3D globe
        # (num_sat, num_epoch, 3) ECI positions: the orbit is drawn as the
        # inertial path oriented at the final epoch (plot_3d rotates it by the
        # final GMST), which reads as the familiar orbit ellipse instead of the
        # pretzel a ground-repeat orbit traces in the rotating frame
        self.sat_pos_hist_3d = None

    def read_config_3d(self, node):
        if node.find('Plot3D') is not None:
            self.plot_3d = misc_fn.str2bool(node.find('Plot3D').text)
        if node.find('ShowSatellite') is not None:
            self.show_satellite = misc_fn.str2bool(node.find('ShowSatellite').text)
        if node.find('ShowOrbit') is not None:
            self.show_orbit = misc_fn.str2bool(node.find('ShowOrbit').text)
        if node.find('SatelliteModelFile') is not None:
            self.model_file = misc_fn.resolve_path(node.find('SatelliteModelFile').text)
        if node.find('SatelliteModelScale') is not None:
            self.model_scale = float(node.find('SatelliteModelScale').text)
        if node.find('MP4') is not None:
            self.mp4 = misc_fn.str2bool(node.find('MP4').text)
        if node.find('EarthClouds') is not None:
            self.earth_clouds = misc_fn.str2bool(node.find('EarthClouds').text)

    def before_loop_3d(self, sm):
        # The 3D movie needs the ECI position history even without ShowOrbit
        if (self.plot_3d and self.show_orbit) or self.mp4:
            self.sat_pos_hist_3d = np.zeros((len(sm.satellites), sm.num_epoch, 3))

    def in_loop_3d(self, sm):
        if self.sat_pos_hist_3d is not None:
            for idx_sat, satellite in enumerate(sm.satellites):
                self.sat_pos_hist_3d[idx_sat, sm.cnt_epoch] = satellite.pos_eci

    def _plot_3d_module(self):
        try:
            import plot_3d  # Deferred: pyvista only needed for 3D plots
            return plot_3d
        except ImportError:
            ls.logger.error('Plot3D requested but pyvista is not installed '
                            '(pip install pyvista). 3D plot skipped.')
            return None

    def _kwargs_3d(self):
        return dict(model_file=self.model_file, model_scale=self.model_scale,
                    show_satellite=self.show_satellite, show_orbit=self.show_orbit,
                    clouds=self.earth_clouds)

    def render_movie_3d(self, sm, satellites, pos_hist_eci, track_latlon=None,
                        swath_edges=None, grid=None):
        """MP4 of the 3D scene: camera circling the first satellite while the
        satellites fly and the analysis content grows (see plot_3d.movie_3d)."""
        p3d = self._plot_3d_module()
        if p3d is None:
            return
        p3d.movie_3d(sm, satellites, pos_hist_eci,
                     sm.output_path(self.type + '_3d.mp4'),
                     track_latlon=track_latlon, swath_edges=swath_edges, grid=grid,
                     model_file=self.model_file, model_scale=self.model_scale,
                     show_orbit=self.show_orbit, clouds=self.earth_clouds)

    def _sats_and_hist(self, sm, satellites):
        """Satellite subset with the matching rows of the position history."""
        if satellites is None:
            return sm.satellites, self.sat_pos_hist_3d
        if self.sat_pos_hist_3d is None:
            return satellites, None
        idx = [sm.satellites.index(satellite) for satellite in satellites]
        return satellites, self.sat_pos_hist_3d[idx]

    def render_3d_points(self, sm, points_llv, scalar_label, cmap='jet', clim=None,
                         point_size=7.0, satellites=None):
        p3d = self._plot_3d_module()
        if p3d is None:
            return
        sats, hist = self._sats_and_hist(sm, satellites)
        p3d.plot_points_3d(sm, sats, hist, np.asarray(points_llv), scalar_label,
                           sm.output_path(self.type + '_3d.png'), cmap=cmap, clim=clim,
                           point_size=point_size, **self._kwargs_3d())

    def render_3d_grid(self, sm, lats_deg, lons_deg, values, scalar_label,
                       cmap='jet', clim=None, satellites=None):
        p3d = self._plot_3d_module()
        if p3d is None:
            return
        sats, hist = self._sats_and_hist(sm, satellites)
        p3d.plot_grid_3d(sm, sats, hist, lats_deg, lons_deg, values, scalar_label,
                         sm.output_path(self.type + '_3d.png'), cmap=cmap, clim=clim,
                         **self._kwargs_3d())

    def render_3d_contours(self, sm, contours, satellites=None):
        p3d = self._plot_3d_module()
        if p3d is None:
            return
        sats, hist = self._sats_and_hist(sm, satellites)
        p3d.plot_contours_3d(sm, sats, hist, contours,
                             sm.output_path(self.type + '_3d.png'), **self._kwargs_3d())


def swath_ribbon_polygons(edges, max_step_deg=30.0):
    """Filled-polygon outlines (lon, lat arrays in deg inside [-180, 180]) of a
    swath edge history, for drawing smooth swaths on a 2D map.

    edges: (num_epoch, 2, 3) left/right swath edge points on the Earth surface
    in ECF; all-zero rows mark epochs without a swath. The history is cut into
    strips at recording gaps, at along-track jumps larger than max_step_deg
    (coarse time steps) and at the latitude turning points of the orbit (so a
    strip cannot fold back over itself). Longitudes are unwrapped along track
    and each strip is clipped against the map window once per overlapping
    360 deg copy, so date line crossings come out as two clean pieces.
    """
    from shapely.geometry import Polygon, box

    # Clip a hair inside the map window: rings lying exactly on the +/-180
    # projection boundary can crash cartopy's polygon reassembly at draw time
    # (MultiPolygon of a GeometryCollection)
    window = box(-179.999, -89.999, 179.999, 89.999)
    polygons = []
    valid = ~np.all(edges.reshape(len(edges), -1) == 0, axis=1)
    idx_valid = np.flatnonzero(valid)
    if idx_valid.size < 2:
        return polygons
    for run in np.split(idx_valid, np.flatnonzero(np.diff(idx_valid) > 1) + 1):
        if len(run) < 2:
            continue
        e = edges[run].astype(float)
        r = np.linalg.norm(e, axis=2)
        lat = np.degrees(np.arcsin(np.clip(e[:, :, 2] / r, -1.0, 1.0)))
        lon = np.degrees(np.arctan2(e[:, :, 1], e[:, :, 0]))
        # Continuous coordinates across the date line: unwrap the left edge
        # along track and keep the right edge within +/-180 deg of it
        lon_l = np.degrees(np.unwrap(np.radians(lon[:, 0])))
        lon_r = lon_l + (lon[:, 1] - lon_l + 180.0) % 360.0 - 180.0
        lat_l, lat_r = lat[:, 0], lat[:, 1]

        # Along-track step in ground degrees (cos scales out the meridian
        # convergence, so near-polar epochs are not mistaken for jumps)
        mean_lat = np.radians((lat_l[1:] + lat_l[:-1]) / 2.0)
        step = np.hypot(np.diff(lon_l) * np.cos(mean_lat), np.diff(lat_l))
        dlat = np.diff(lat_l)
        turning = set((np.flatnonzero(np.sign(dlat[1:]) * np.sign(dlat[:-1]) < 0) + 1).tolist())
        strips, current = [], [0]
        for k in range(1, len(run)):
            if step[k - 1] > max_step_deg:  # Jump: start a new strip
                strips.append(current)
                current = [k]
            else:
                current.append(k)
                if k in turning:  # Fold point: new strip sharing the vertex
                    strips.append(current)
                    current = [k]
        strips.append(current)

        for strip in strips:
            if len(strip) < 2:
                continue
            s = np.asarray(strip)
            px = np.concatenate([lon_l[s], lon_r[s][::-1]])  # Left edge out, right edge back
            py = np.concatenate([lat_l[s], lat_r[s][::-1]])
            for shift in range(int(np.ceil((-180.0 - px.max()) / 360.0)),
                               int(np.floor((180.0 - px.min()) / 360.0)) + 1):
                poly = Polygon(np.column_stack([px + 360.0 * shift, py]))
                if not poly.is_valid:
                    poly = poly.buffer(0)
                inter = poly.intersection(window)
                geoms = [inter] if inter.geom_type == 'Polygon' else getattr(inter, 'geoms', [])
                for g in geoms:
                    if g.geom_type == 'Polygon' and not g.is_empty:
                        x, y = g.exterior.xy
                        polygons.append((np.asarray(x), np.asarray(y)))
    return polygons


def fill_ribbon_safe(ax, lon, lat, **style):
    """ax.fill of a lon/lat polygon with the map projection applied up front:
    cartopy's deferred draw-time projection crashes on some degenerate rings
    (polar stereographic ring reassembly can produce a GeometryCollection),
    so project through shapely here, retry with a repaired/simplified
    outline, and let the caller skip the ribbon instead of the whole run
    dying inside savefig. Returns True when the ribbon was drawn."""
    from shapely.geometry import Polygon
    poly = Polygon(np.column_stack([lon, lat]))
    for attempt in (poly, poly.buffer(0).simplify(0.05)):
        if attempt.is_empty:
            continue
        try:
            geom = ax.projection.project_geometry(attempt, ccrs.PlateCarree())
        except Exception:
            continue
        for g in getattr(geom, 'geoms', [geom]):
            if g.geom_type == 'Polygon' and not g.is_empty:
                x, y = g.exterior.xy
                ax.fill(np.asarray(x), np.asarray(y), **style)  # projection coords
        return True
    return False


class AnalysisObs:   # Common methods needed for some OBS analysis

    def plot_swath_coverage(self, sm, swath_edges, polar_view):
        """Swath coverage as smooth semi-transparent ribbons on the map — the
        2D counterpart of the 3D globe render, built from the same left/right
        swath edge histories instead of colouring the discrete user grid
        points. Overlapping passes show darker through the alpha stacking."""
        if polar_view is not None:
            fig, ax = make_map_polar(polar_view)
        else:
            fig, ax = make_map_cyl()
        if isinstance(self, AnalysisMap2D):
            self.decorate_map2d(sm, ax)
        skipped = 0
        for idx_sat in range(swath_edges.shape[0]):
            for lon, lat in swath_ribbon_polygons(swath_edges[idx_sat]):
                if not fill_ribbon_safe(ax, lon, lat, facecolor='orangered',
                                        edgecolor='orangered', linewidth=0.3,
                                        alpha=0.4):
                    skipped += 1
        if skipped:
            ls.logger.warning(f'{self.type}: {skipped} swath ribbon(s) could not '
                              f'be projected onto the map and were skipped')
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

    def revisit_gaps_hours(self, user_metric, time_step):
        """Per-user revisit gap arrays in hours from the (num_user, num_epoch)
        in-swath flag history: intervals between successive coverage
        transitions, dropping the one-epoch intervals (the pass itself)."""
        gaps_list = []
        for idx_user in range(user_metric.shape[0]):
            gaps = np.diff(np.where(np.diff(user_metric[idx_user, :]) != 0)).flatten()
            if len(gaps) > 1:
                gaps = np.delete(gaps, np.where(gaps == 1), axis=0) * time_step / 3600.0
            else:
                gaps = np.array([])
            gaps_list.append(gaps)
        return gaps_list

    def plot_swath_revisit(self, sm, user_metric, statistic, polar_view):
        gaps_list = self.revisit_gaps_hours(user_metric, sm.time_step)
        plot_points = np.zeros((len(sm.users), 3))
        metric = 0
        for idx_user, user in enumerate(sm.users):
            if idx_user % 1000 == 0:
                ls.logger.info(f'User revisit {user.user_id} of {len(sm.users)}')
            gaps = gaps_list[idx_user]
            if len(gaps) > 0:
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
        self.write_csv(sm, ['lon_deg', 'lat_deg', f'{statistic}_revisit_hours'],
                       plot_points, suffix='revisit')
        if metric>0:  # Only plot if not empty
            if polar_view is not None:
                fig, ax = make_map_polar(polar_view)
                point_size = 5
            else:
                fig, ax = make_map_cyl()
                point_size = 3
            if isinstance(self, AnalysisMap2D):
                self.decorate_map2d(sm, ax)
            sc = ax.scatter(plot_points[:,0], plot_points[:,1], s=point_size, cmap=plt.cm.jet,
                            c=plot_points[:,2], vmin=np.nanmin(plot_points[:,2]),
                            vmax=np.nanmax(plot_points[:,2]), transform=ccrs.PlateCarree())
            cb = plt.colorbar(sc, ax=ax, shrink=0.85)
            cb.set_label(statistic.capitalize() + ' Revisit Interval [hours]', fontsize=10)
            plt.savefig(sm.output_path(f'{self.type}_revisit.png'))
            plt.show()

    def plot_swath_revisit_latitude(self, sm, user_metric):
        """Longitude-averaged max and mean revisit time versus latitude in
        days: per user grid point the max/mean revisit gap, averaged over all
        grid points sharing the same latitude."""
        gaps_list = self.revisit_gaps_hours(user_metric, sm.time_step)
        lat_max, lat_mean = {}, {}
        for idx_user, user in enumerate(sm.users):
            gaps = gaps_list[idx_user]
            if len(gaps) == 0:
                continue
            lat = round(degrees(user.lla[0]), 6)
            lat_max.setdefault(lat, []).append(np.max(gaps) / 24.0)
            lat_mean.setdefault(lat, []).append(np.mean(gaps) / 24.0)
        if not lat_max:
            ls.logger.warning('No revisit data found for the latitude profile. No plot produced.')
            return
        lats = np.array(sorted(lat_max))
        max_avg = np.array([np.mean(lat_max[lat]) for lat in lats])
        mean_avg = np.array([np.mean(lat_mean[lat]) for lat in lats])
        self.write_csv(sm, ['lat_deg', 'max_revisit_days', 'mean_revisit_days'],
                       np.column_stack([lats, max_avg, mean_avg]), suffix='revisit_lat')
        fig = plt.figure(figsize=(10, 6))
        plt.plot(lats, max_avg, 'r-', label='Max Revisit Time')
        plt.plot(lats, mean_avg, 'b-', label='Mean Revisit Time')
        plt.title('Longitude-averaged Max and Mean Revisit Time')
        plt.xlabel('Latitude [deg]'); plt.ylabel('Revisit Time [days]')
        plt.xticks(np.arange(-90, 91, 15))
        plt.ylim(bottom=0)
        plt.grid(); plt.legend()
        plt.savefig(sm.output_path(f'{self.type}_revisit_lat.png'))
        plt.show()

