"""
3D visualisation helpers for SatSVS based on PyVista (VTK).

Used by the analyses that accept <Plot3D>True</Plot3D> (all world-map analyses):
they render their result (ground track, swath ribbon, statistic field, point
cloud or visibility contour cap) on a textured 3D Earth (NASA Blue Marble,
input/earth_texture.jpg), together with the satellite orbital track(s) and a 3D
satellite model at the last simulated epoch, nadir-pointing with the solar
panel axis cross-track. A custom satellite mesh (STL/OBJ/PLY/VTK) can be
supplied through <SatelliteModelFile>; otherwise a procedural bus + solar
panel + dish model is built. The satellite is drawn hugely exaggerated
(default 200 km, set by <SatelliteModelScale> in metres) so it is visible at
globe scale.

Surface geometry (tracks, swaths, fields, contours, stations) is in the
Earth-fixed frame so it lines up with the texture. The orbital track is the
inertial (ECI) path oriented at the final epoch — the familiar orbit ellipse
around the Earth of that instant (see _orbit_scene_points). The plot is
rendered off-screen to ../output/<analysis>_3d.png when matplotlib runs
headless (MPLBACKEND=Agg, as the test runner does); otherwise an interactive
window opens first and the screenshot is saved on close.
"""
import os
from math import radians, degrees

import numpy as np

from constants import R_EARTH
import logging_svs as ls

EARTH_TEXTURE = '../input/earth_texture.jpg'


def _make_earth(radius):
    """Textured Earth sphere from a lat/lon structured grid. The seam column is
    duplicated (lon -180 and +180) so the equirectangular texture wraps cleanly,
    with texture coordinates assigned in the grid's own point order."""
    import pyvista as pv
    lon = np.radians(np.linspace(-180.0, 180.0, 361))
    lat = np.radians(np.linspace(-90.0, 90.0, 181))
    lon2, lat2 = np.meshgrid(lon, lat, indexing='ij')
    x = radius * np.cos(lat2) * np.cos(lon2)
    y = radius * np.cos(lat2) * np.sin(lon2)
    z = radius * np.sin(lat2)
    grid = pv.StructuredGrid(x, y, z)
    u = (np.degrees(lon2) + 180.0) / 360.0
    v = (np.degrees(lat2) + 90.0) / 180.0
    # StructuredGrid flattens the input arrays in Fortran order; use the same
    # order so every texture coordinate belongs to its own grid point
    grid.active_texture_coordinates = np.column_stack((u.ravel(order='F'),
                                                       v.ravel(order='F')))
    texture = pv.read_texture(EARTH_TEXTURE)
    return grid, texture


def _make_satellite_meshes(model_file=None):
    """Satellite model in the body frame (+x flight direction, +y solar panel
    axis, +z towards nadir), roughly unit-sized. Returns [(mesh, color), ...]."""
    import pyvista as pv
    if model_file:
        mesh = pv.read(model_file)
        # Normalise user models to unit size around their centre
        mesh = mesh.translate(np.array(mesh.center) * -1.0)
        extent = max(mesh.length, 1e-9)
        mesh = mesh.scale(2.0 / extent)
        return [(mesh, 'lightgray')]
    bus = pv.Cube(x_length=0.9, y_length=0.9, z_length=1.4)
    panel1 = pv.Cube(center=(0.0, 1.6, 0.0), x_length=0.9, y_length=2.2, z_length=0.04)
    panel2 = pv.Cube(center=(0.0, -1.6, 0.0), x_length=0.9, y_length=2.2, z_length=0.04)
    dish = pv.Cone(center=(0.0, 0.0, 0.85), direction=(0.0, 0.0, -1.0),
                   height=0.5, radius=0.45, resolution=24, capping=False)
    return [(bus, 'gold'), (panel1, 'midnightblue'), (panel2, 'midnightblue'),
            (dish, 'lightgray')]


def _body_to_ecf(pos_ecf, vel_ecf):
    """Rotation matrix with columns = body axes in ECF: +z nadir, +x along the
    (horizontal component of the) velocity, +y completing the triad."""
    z_axis = -pos_ecf / np.linalg.norm(pos_ecf)  # Nadir
    vel = vel_ecf if np.linalg.norm(vel_ecf) > 0 else np.array([1.0, 0.0, 0.0])
    x_axis = vel - np.dot(vel, z_axis) * z_axis
    x_axis = x_axis / np.linalg.norm(x_axis)
    y_axis = np.cross(z_axis, x_axis)
    return np.column_stack((x_axis, y_axis, z_axis))


def _track_points(lats_deg, lons_deg, radius):
    lats = np.radians(lats_deg)
    lons = np.radians(lons_deg)
    return np.column_stack((radius * np.cos(lats) * np.cos(lons),
                            radius * np.cos(lats) * np.sin(lons),
                            radius * np.sin(lats)))


def _open_scene(sm):
    """Plotter with black background and the textured Earth already added."""
    import pyvista as pv
    off_screen = os.environ.get('MPLBACKEND', '').lower() == 'agg'
    plotter = pv.Plotter(off_screen=off_screen, window_size=(1600, 1200))
    plotter.set_background('black')
    earth, texture = _make_earth(R_EARTH)
    # Ambient light keeps the dark oceans of the texture visibly distinct from
    # the black space background, so the full Earth disk and its limb stay
    # readable (otherwise geometry on the surface can appear to float in space)
    plotter.add_mesh(earth, texture=texture, smooth_shading=True, ambient=0.3)
    return plotter, off_screen


def _orbit_scene_points(pos_hist_eci, gmst):
    """Rotate an ECI position history into the scene (Earth-fixed) frame as it
    stands at the render epoch: the drawn orbit is the familiar inertial orbit
    path around the Earth of that instant. Drawing the Earth-fixed history
    instead would show the path relative to the rotating Earth, which for e.g.
    a ground-repeat GNSS orbit is a confusing pretzel-shaped loop (that
    information is what the ground track on the surface already shows)."""
    cos_g, sin_g = np.cos(gmst), np.sin(gmst)
    return np.column_stack((cos_g * pos_hist_eci[:, 0] + sin_g * pos_hist_eci[:, 1],
                            -sin_g * pos_hist_eci[:, 0] + cos_g * pos_hist_eci[:, 1],
                            pos_hist_eci[:, 2]))


def _add_orbit_and_model(plotter, satellite, pos_hist, model_file, model_scale,
                         show_satellite=True, show_orbit=True, draw_nadir=True):
    """Orbit polyline (already in scene coordinates), the satellite model at
    the last epoch (nadir pointing) and optionally the nadir line down to the
    subsatellite point (suppressed by the callers for large constellations,
    where the many crossing nadir lines only clutter the scene)."""
    import pyvista as pv
    if show_orbit and pos_hist is not None and len(pos_hist) > 1:
        plotter.add_mesh(pv.lines_from_points(pos_hist), color='cyan', line_width=1)
    if not show_satellite:
        return
    rot = _body_to_ecf(satellite.pos_ecf, satellite.vel_ecf)
    for mesh, color in _make_satellite_meshes(model_file):
        pts = np.asarray(mesh.points, dtype=float) * model_scale
        mesh = mesh.copy()
        mesh.points = (rot @ pts.T).T + satellite.pos_ecf
        plotter.add_mesh(mesh, color=color, smooth_shading=False)
    if draw_nadir:
        nadir_end = satellite.pos_ecf / np.linalg.norm(satellite.pos_ecf) * R_EARTH
        plotter.add_mesh(pv.Line(satellite.pos_ecf, nadir_end),
                         color='yellow', line_width=1)


def _add_stations(plotter, sm):
    import pyvista as pv
    for station in sm.stations:
        marker = pv.Sphere(radius=R_EARTH * 0.008, center=station.pos_ecf)
        plotter.add_mesh(marker, color='magenta')


def _finish_scene(plotter, satellites, off_screen, file_name, content_radius=None):
    """Camera above the (last) satellite looking at the Earth centre, then
    render off-screen or show interactively (screenshot saved on close).

    Near-Earth scenes (LEO) use the default perspective camera. Scenes whose
    content reaches far above the Earth (MEO/GEO constellations) switch to a
    fitted parallel projection: a perspective camera close enough to keep the
    globe large magnifies the orbit arcs near the camera by several times,
    which makes the Earth look far too small compared to the orbits."""
    ref = satellites[-1].pos_ecf if len(satellites) else np.array([R_EARTH * 2, 0, 0])
    r_content = max([R_EARTH, content_radius or 0.0] +
                    [np.linalg.norm(s.pos_ecf) for s in satellites])
    direction = ref / np.linalg.norm(ref)
    if len(satellites):
        # Tilt the view ~30 deg towards the orbit normal: a camera exactly
        # above the satellite looks edge-on onto its orbit plane, so the drawn
        # (inertial) orbit ellipse would degenerate to a straight line
        normal = np.cross(ref, satellites[-1].vel_ecf)
        norm_normal = np.linalg.norm(normal)
        if norm_normal > 0.0:
            direction = direction + 0.58 * normal / norm_normal
            direction = direction / np.linalg.norm(direction)
    if r_content > 1.5 * R_EARTH:
        plotter.camera_position = [(direction * r_content * 4.0).tolist(),
                                   [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
        plotter.enable_parallel_projection()
        plotter.camera.parallel_scale = 1.08 * r_content  # Vertical half-extent
    else:
        cam = direction * max(R_EARTH * 3.2, np.linalg.norm(ref) * 1.7)
        plotter.camera_position = [cam.tolist(), [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
    if off_screen:
        plotter.screenshot(file_name)
    else:
        plotter.show(screenshot=file_name)  # Interactive; saved on close
    plotter.close()
    ls.logger.info(f'Saved 3D plot to {file_name}')


SCALAR_BAR_ARGS = dict(color='white', title_font_size=22, label_font_size=17,
                       vertical=False, position_x=0.25, position_y=0.03, width=0.5)


def plot_ground_track_3d(sm, satellites, file_name,
                         model_file=None, model_scale=200e3,
                         show_satellite=True, show_orbit=True):
    """Render ground track(s), ECF orbit path(s), stations and a 3D satellite
    model on a textured Earth. satellites: list of Satellite objects whose
    metric holds per-epoch [lat_deg, lon_deg, x_ecf, y_ecf, z_ecf]."""
    import pyvista as pv

    plotter, off_screen = _open_scene(sm)
    r_content = 0.0
    for satellite in satellites:
        used = ~np.all(satellite.metric == 0, axis=1)  # Skip never-filled epochs
        if not used.any():
            continue
        lats, lons = satellite.metric[used, 0], satellite.metric[used, 1]

        # Ground track slightly above the surface to avoid z-fighting
        track = pv.lines_from_points(_track_points(lats, lons, R_EARTH * 1.003))
        plotter.add_mesh(track, color='red', line_width=3)

        pos_hist = satellite.metric[used, 2:5]  # ECI history
        if show_orbit:
            r_content = max(r_content, np.linalg.norm(pos_hist, axis=1).max())
        _add_orbit_and_model(plotter, satellite, _orbit_scene_points(pos_hist, sm.time_gmst),
                             model_file, model_scale, show_satellite, show_orbit)
    _add_stations(plotter, sm)
    _finish_scene(plotter, satellites, off_screen, file_name, r_content)


def _slerp_across(edges, n_across):
    """Spherically interpolate between the left and right edge points of every
    epoch so the ribbon follows the Earth curvature instead of cutting through
    it as a straight chord. edges: (n, 2, 3) -> (n, n_across, 3)."""
    left, right = edges[:, 0, :], edges[:, 1, :]
    r_left = np.linalg.norm(left, axis=1)
    r_right = np.linalg.norm(right, axis=1)
    u_left = left / r_left[:, None]
    u_right = right / r_right[:, None]
    omega = np.arccos(np.clip(np.sum(u_left * u_right, axis=1), -1.0, 1.0))
    sin_omega = np.where(np.sin(omega) < 1e-12, 1.0, np.sin(omega))  # Degenerate: width ~ 0
    t = np.linspace(0.0, 1.0, n_across)
    w_left = np.sin((1.0 - t)[None, :] * omega[:, None]) / sin_omega[:, None]
    w_right = np.sin(t[None, :] * omega[:, None]) / sin_omega[:, None]
    narrow = omega < 1e-12
    if narrow.any():  # Fall back to linear weights where the arc vanishes
        w_left[narrow] = (1.0 - t)[None, :]
        w_right[narrow] = t[None, :]
    radius = r_left[:, None] * (1.0 - t)[None, :] + r_right[:, None] * t[None, :]
    direction = w_left[:, :, None] * u_left[:, None, :] + w_right[:, :, None] * u_right[:, None, :]
    return direction * radius[:, :, None]


def plot_swath_3d(sm, satellites, sat_pos_hist, swath_edges, file_name,
                  model_file=None, model_scale=200e3,
                  show_satellite=True, show_orbit=True):
    """Render the swath coverage as a semi-transparent ribbon on the textured
    Earth, plus orbit path(s), stations and the 3D satellite model.

    sat_pos_hist: (num_sat, num_epoch, 3) satellite ECF positions
    swath_edges:  (num_sat, num_epoch, 2, 3) left/right swath edge points on the
                  Earth surface in ECF (for push broom the two line-of-sight
                  ground intersections, for conical the cross-track extremes)
    """
    import pyvista as pv

    N_ACROSS = 9  # Cross-track samples of the ribbon (follows Earth curvature)

    plotter, off_screen = _open_scene(sm)
    # Proper per-fragment ordering of the overlapping translucent ribbons;
    # without it the crossing ascending/descending strips render blotchy
    plotter.enable_depth_peeling(number_of_peels=8)
    for idx_sat, satellite in enumerate(satellites):
        edges = swath_edges[idx_sat]
        used = ~np.all(edges.reshape(len(edges), -1) == 0, axis=1)
        if used.sum() > 1:
            # Curvature-following strip between the edge histories, raised above
            # the globe. The raise ramps slowly over the simulation so that
            # overlapping passes sit at slightly different radii: two translucent
            # layers at the same radius would z-fight into zigzag artifacts.
            ramp = 1.004 + 0.004 * np.linspace(0.0, 1.0, int(used.sum()))
            strip = _slerp_across(edges[used], N_ACROSS) * ramp[:, None, None]
            grid = pv.StructuredGrid(
                np.ascontiguousarray(strip[:, :, 0].T.reshape(N_ACROSS, -1, 1)),
                np.ascontiguousarray(strip[:, :, 1].T.reshape(N_ACROSS, -1, 1)),
                np.ascontiguousarray(strip[:, :, 2].T.reshape(N_ACROSS, -1, 1)))
            plotter.add_mesh(grid, color='orangered', opacity=0.5)
        pos_used = ~np.all(sat_pos_hist[idx_sat] == 0, axis=1)
        _add_orbit_and_model(plotter, satellite,
                             _orbit_scene_points(sat_pos_hist[idx_sat][pos_used], sm.time_gmst),
                             model_file, model_scale, show_satellite, show_orbit)
    _add_stations(plotter, sm)
    r_content = float(np.linalg.norm(sat_pos_hist.reshape(-1, 3), axis=1).max()) \
        if show_orbit else 0.0
    _finish_scene(plotter, satellites, off_screen, file_name, r_content)


def _add_orbits_and_models(plotter, gmst, satellites, sat_pos_hist,
                           model_file, model_scale, show_satellite, show_orbit):
    """Orbit + model for every satellite (sat_pos_hist holds ECI positions and
    may be None). Returns the largest radius of the drawn orbit content."""
    draw_nadir = len(satellites) <= 2
    r_content = 0.0
    for idx_sat, satellite in enumerate(satellites):
        pos_hist = None
        if sat_pos_hist is not None:
            pos_hist = sat_pos_hist[idx_sat]
            pos_hist = pos_hist[~np.all(pos_hist == 0, axis=1)]
            if show_orbit and len(pos_hist):
                r_content = max(r_content, np.linalg.norm(pos_hist, axis=1).max())
            pos_hist = _orbit_scene_points(pos_hist, gmst)
        _add_orbit_and_model(plotter, satellite, pos_hist,
                             model_file, model_scale, show_satellite, show_orbit,
                             draw_nadir)
    return r_content


def plot_points_3d(sm, satellites, sat_pos_hist, points_llv, scalar_label,
                   file_name, cmap='jet', clim=None, point_size=7.0,
                   model_file=None, model_scale=200e3,
                   show_satellite=True, show_orbit=True):
    """Render a value-coloured point cloud on the textured Earth (3D version of
    the scatter world maps). points_llv: (n, 3) array of longitude [deg],
    latitude [deg], value."""
    import pyvista as pv

    plotter, off_screen = _open_scene(sm)
    if len(points_llv):
        cloud = pv.PolyData(_track_points(points_llv[:, 1], points_llv[:, 0],
                                          R_EARTH * 1.003))
        cloud[scalar_label] = points_llv[:, 2]
        plotter.add_mesh(cloud, scalars=scalar_label, cmap=cmap, clim=clim,
                         point_size=point_size, render_points_as_spheres=True,
                         scalar_bar_args=dict(title=scalar_label, **SCALAR_BAR_ARGS))
    r_content = _add_orbits_and_models(plotter, sm.time_gmst, satellites, sat_pos_hist,
                                       model_file, model_scale, show_satellite, show_orbit)
    _add_stations(plotter, sm)
    _finish_scene(plotter, satellites, off_screen, file_name, r_content)


def plot_grid_3d(sm, satellites, sat_pos_hist, lats_deg, lons_deg, values,
                 scalar_label, file_name, cmap='jet', clim=None,
                 model_file=None, model_scale=200e3,
                 show_satellite=True, show_orbit=True):
    """Render a lat/lon grid statistic as a coloured, slightly transparent layer
    draped over the textured Earth (3D version of the pcolormesh world maps).
    values: (num_lat, num_lon) array matching lats_deg/lons_deg."""
    import pyvista as pv

    plotter, off_screen = _open_scene(sm)
    # Subdivide the (typically coarse) user grid to <= 2 deg mesh cells with
    # nearest-cell scalar sampling: flat quads between 10 deg grid points would
    # sag below the globe surface, leaving only patches around the vertices
    lats_deg = np.asarray(lats_deg, dtype=float)
    lons_deg = np.asarray(lons_deg, dtype=float)
    values = np.asarray(values, dtype=float)
    up = max(1, int(np.ceil(max(np.max(np.diff(lats_deg)), np.max(np.diff(lons_deg))) / 2.0)))
    fine_lats = np.linspace(lats_deg[0], lats_deg[-1], (len(lats_deg) - 1) * up + 1)
    fine_lons = np.linspace(lons_deg[0], lons_deg[-1], (len(lons_deg) - 1) * up + 1)
    idx_lat = np.rint(np.interp(fine_lats, lats_deg, np.arange(len(lats_deg)))).astype(int)
    idx_lon = np.rint(np.interp(fine_lons, lons_deg, np.arange(len(lons_deg)))).astype(int)
    fine_values = values[np.ix_(idx_lat, idx_lon)]
    lon2, lat2 = np.meshgrid(np.radians(fine_lons), np.radians(fine_lats), indexing='ij')
    r = R_EARTH * 1.002
    grid = pv.StructuredGrid(r * np.cos(lat2) * np.cos(lon2),
                             r * np.cos(lat2) * np.sin(lon2),
                             r * np.sin(lat2))
    # Same Fortran point order as the meshgrid arrays given to StructuredGrid
    grid[scalar_label] = fine_values.T.ravel(order='F')
    plotter.add_mesh(grid, scalars=scalar_label, cmap=cmap, clim=clim,
                     opacity=0.75, nan_opacity=0.0,
                     scalar_bar_args=dict(title=scalar_label, **SCALAR_BAR_ARGS))
    r_content = _add_orbits_and_models(plotter, sm.time_gmst, satellites, sat_pos_hist,
                                       model_file, model_scale, show_satellite, show_orbit)
    _add_stations(plotter, sm)
    _finish_scene(plotter, satellites, off_screen, file_name, r_content)


def plot_contours_3d(sm, satellites, sat_pos_hist, contours, file_name,
                     model_file=None, model_scale=200e3,
                     show_satellite=True, show_orbit=True):
    """Render satellite visibility contours on the textured Earth as filled,
    semi-transparent spherical caps with a bold outline, so the covered region
    reads unambiguously as lying on the surface. contours: list of (n, 2)
    arrays of [lat_rad, lon_rad] per satellite (same order as satellites, used
    for the colour cycle)."""
    import pyvista as pv
    from matplotlib import cm

    N_RADIAL = 24  # Radial samples of the cap fill (follows Earth curvature)

    plotter, off_screen = _open_scene(sm)
    plotter.enable_depth_peeling(number_of_peels=8)  # Overlapping translucent caps
    colors = cm.tab10(np.linspace(0.0, 0.9, 10))
    for idx, contour in enumerate(contours):
        color = colors[idx % 10][:3].tolist()
        ring = _track_points(np.degrees(contour[:, 0]), np.degrees(contour[:, 1]),
                             R_EARTH * 1.006)
        # Filled cap: spherical fan between the ring's spherical centroid (the
        # subsatellite point for a plain visibility cone) and the ring itself
        centre_dir = ring.mean(axis=0)
        centre_dir = centre_dir / np.linalg.norm(centre_dir)
        edges = np.empty((len(ring), 2, 3))
        edges[:, 0, :] = centre_dir * R_EARTH * 1.004
        edges[:, 1, :] = ring / 1.006 * 1.004
        fan = _slerp_across(edges, N_RADIAL)
        cap = pv.StructuredGrid(
            np.ascontiguousarray(fan[:, :, 0].T.reshape(N_RADIAL, -1, 1)),
            np.ascontiguousarray(fan[:, :, 1].T.reshape(N_RADIAL, -1, 1)),
            np.ascontiguousarray(fan[:, :, 2].T.reshape(N_RADIAL, -1, 1)))
        plotter.add_mesh(cap, color=color, opacity=0.35)
        # Bold closed outline on top of the fill
        outline = np.vstack([ring, ring[:1]])
        plotter.add_mesh(pv.lines_from_points(outline), color=color, line_width=5)
    r_content = _add_orbits_and_models(plotter, sm.time_gmst, satellites, sat_pos_hist,
                                       model_file, model_scale, show_satellite, show_orbit)
    _add_stations(plotter, sm)
    _finish_scene(plotter, satellites, off_screen, file_name, r_content)
