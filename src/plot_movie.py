"""
2D world-map movies (MP4) for the analyses that accept <MP4>True</MP4>.

The world map fills up with the analysis data over the simulation time:
a growing ground track (movie_track_2d), growing semi-transparent swath
ribbons (movie_ribbons_2d) or the instantaneous per-epoch metric field on the
user grid (movie_grid_2d). Movies are written as ../output/<type>_2d.mp4 with
imageio + imageio-ffmpeg (the 3D counterpart movie_3d lives in plot_3d.py and
shares open_writer/frame_epochs below).
"""
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import logging_svs as ls

MAX_FRAMES = 360  # Upper bound on movie frames; epochs are subsampled beyond it


def open_writer(file_name, fps=20):
    """MP4 writer, or None (with a clear error log) when imageio/ffmpeg is
    missing. macro_block_size=8 accepts the map figure's 1000x440 canvas."""
    try:
        import imageio.v2 as imageio
    except ImportError:
        ls.logger.error('MP4 requested but imageio / imageio-ffmpeg is not '
                        'installed (pip install imageio imageio-ffmpeg). '
                        'Movie skipped.')
        return None
    return imageio.get_writer(file_name, fps=fps, codec='libx264', quality=7,
                              macro_block_size=8)


def frame_epochs(num_epoch, max_frames=MAX_FRAMES):
    """Epoch index of every movie frame: all epochs when they fit, otherwise
    evenly subsampled, always ending on the last epoch."""
    step = max(1, int(np.ceil(num_epoch / max_frames)))
    epochs = list(range(step - 1, num_epoch, step))
    if epochs[-1] != num_epoch - 1:
        epochs.append(num_epoch - 1)
    return epochs


def _grab(fig):
    """Current figure canvas as an RGB frame."""
    fig.canvas.draw()
    return np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()


def movie_track_2d(sm, sat_metrics, file_name, color_index=None, cmap='RdYlBu',
                   clim=None, label=None, fps=20):
    """Progressive ground-track movie: the map fills up with the per-epoch
    positions of every satellite. sat_metrics: (num_sat, num_epoch, >=2) with
    [:, :, 0]=lat [deg] and [:, :, 1]=lon [deg]; all-zero rows are skipped.
    With color_index the points are coloured by that metric column."""
    from analysis import make_map_cyl
    writer = open_writer(file_name, fps)
    if writer is None:
        return
    fig, ax = make_map_cyl()
    num_sat, num_epoch = sat_metrics.shape[0], sat_metrics.shape[1]
    filled = [~np.all(sat_metrics[i, :, 0:2] == 0, axis=1) for i in range(num_sat)]
    scatters = []
    for i in range(num_sat):
        if color_index is None:
            sc = ax.scatter(np.empty(0), np.empty(0), s=3, color='red',
                            transform=ccrs.PlateCarree())
        else:
            if clim is None:
                vals = np.concatenate([sat_metrics[j, filled[j], color_index]
                                       for j in range(num_sat)] or [np.zeros(1)])
                clim = (np.nanmin(vals), np.nanmax(vals))
            sc = ax.scatter(np.empty(0), np.empty(0), s=3, c=np.empty(0),
                            cmap=cmap, vmin=clim[0], vmax=clim[1],
                            transform=ccrs.PlateCarree())
        scatters.append(sc)
    if color_index is not None and label:
        cb = fig.colorbar(scatters[0], ax=ax, shrink=0.85)
        cb.set_label(label, fontsize=10)
    for k in frame_epochs(num_epoch):
        for i, sc in enumerate(scatters):
            m = filled[i][:k + 1]
            data = sat_metrics[i, :k + 1][m]
            sc.set_offsets(np.column_stack([data[:, 1], data[:, 0]]))
            if color_index is not None:
                sc.set_array(data[:, color_index])
        ax.set_title(sm.times_str_pre[k], fontsize=10)
        writer.append_data(_grab(fig))
    writer.close()
    plt.close(fig)
    ls.logger.info(f'Saved 2D movie to {file_name}')


def movie_ribbons_2d(sm, swath_edges, file_name, fps=20):
    """Progressive swath movie: the semi-transparent ribbons fill the map up
    to the current epoch. swath_edges: (num_sat, num_epoch, 2, 3) ECF edges."""
    from analysis import make_map_cyl, swath_ribbon_polygons
    writer = open_writer(file_name, fps)
    if writer is None:
        return
    fig, ax = make_map_cyl()
    num_epoch = swath_edges.shape[1]
    artists = []
    for k in frame_epochs(num_epoch):
        for artist in artists:
            artist.remove()
        artists = []
        for idx_sat in range(swath_edges.shape[0]):
            for lon, lat in swath_ribbon_polygons(swath_edges[idx_sat, :k + 1]):
                artists.extend(ax.fill(lon, lat, facecolor='orangered',
                                       edgecolor='orangered', linewidth=0.3,
                                       alpha=0.4, transform=ccrs.PlateCarree()))
        ax.set_title(sm.times_str_pre[k], fontsize=10)
        writer.append_data(_grab(fig))
    writer.close()
    plt.close(fig)
    ls.logger.info(f'Saved 2D movie to {file_name}')


def movie_grid_2d(sm, user_metric, analysis_type, label, file_name,
                  cmap='jet', fps=20):
    """Instantaneous metric-field movie on the user grid: every frame shows
    the per-epoch values of user_metric (num_user, num_epoch) as a world map
    (needs a user segment of Type Grid)."""
    from math import degrees
    from analysis import make_map_cyl, get_user_grid_shape
    grid_shape = get_user_grid_shape(sm, analysis_type)
    if grid_shape is None:
        return
    finite = np.isfinite(user_metric)
    if not finite.any():
        ls.logger.warning(f'No data for the {analysis_type} 2D movie. Skipped.')
        return
    writer = open_writer(file_name, fps)
    if writer is None:
        return
    lats = np.reshape([degrees(u.lla[0]) for u in sm.users], grid_shape)
    lons = np.reshape([degrees(u.lla[1]) for u in sm.users], grid_shape)
    vmin, vmax = np.nanmin(user_metric[finite]), np.nanmax(user_metric[finite])
    dlat = (lats.max() - lats.min()) / max(grid_shape[0] - 1, 1)
    dlon = (lons.max() - lons.min()) / max(grid_shape[1] - 1, 1)
    fig, ax = make_map_cyl()
    num_epoch = user_metric.shape[1]
    im = ax.imshow(user_metric[:, 0].reshape(grid_shape), origin='lower',
                   extent=(lons.min() - dlon / 2, lons.max() + dlon / 2,
                           lats.min() - dlat / 2, lats.max() + dlat / 2),
                   transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax,
                   interpolation='nearest')
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label(label, fontsize=10)
    for k in frame_epochs(num_epoch):
        im.set_data(user_metric[:, k].reshape(grid_shape))
        ax.set_title(sm.times_str_pre[k], fontsize=10)
        writer.append_data(_grab(fig))
    writer.close()
    plt.close(fig)
    ls.logger.info(f'Saved 2D movie to {file_name}')
