import datashader as ds
import pandas as pd
from colorcet import fire
from datashader import transfer_functions as tf
import bokeh.plotting as bp
from datashader.bokeh_ext import InteractiveImage

df = pd.read_csv('/Users/micheltossaint/Documents/86 - Ipython Notebooks/Tutorial/datashader-examples/data/nyc_taxi.csv', usecols=['dropoff_x', 'dropoff_y'])
print(df.head())

canvas = ds.Canvas(plot_width=10000, plot_height=10000,
x_axis_type='linear', y_axis_type='linear')
agg = ds.Canvas().points(df,'dropoff_x', 'dropoff_y')
img = tf.set_background(tf.shade(agg, cmap=fire),None)
ds.utils.export_image(img=img,filename='test', fmt=".png", background=None)

# opts.defaults(opts.WMTS(width=500, height=500))
# tiles = gv.WMTS('https://maps.wikimedia.org/osm-intl/{Z}/{X}/{Y}@2x.png')
# points = gv.operation.project_points(gv.Points(verts, vdims=['z']))
# tiles * datashade(hv.TriMesh((tris, points)), aggregator=ds.mean('z'), precompute=True)

# bp.output_notebook()
# p = bp.figure(tools='pan,wheel_zoom,reset', x_range=(-5,5), y_range=(-5,5), plot_width=500, plot_height=500)
# InteractiveImage(p, image_callback)