
# Define user grid shapefile
import ast
#from shapely.geometry import Point, Polygon
import numpy as np
from geopandas import GeoSeries, GeoDataFrame
import geopandas
import matplotlib.pyplot as plt


# Get the shape-file for NYC
#poly = GeoDataFrame.from_file('../input/example_shp_files/EMSA_draft_AOIs_CleanSeaNet.shp')
#poly = poly[poly.ADMIN=='Denmark']['geometry'].iloc[0]
#poly = poly[poly.ADMIN=='Denmark']['geometry'].iloc[0]
# plt.rcParams["figure.figsize"] = [8,6]
# poly.plot()
# plt.show()
# print(poly.head())
# print(poly['geometry'].iloc[0])
# print(poly.bounds)

# print(poly.iloc[100])

# # Create a Polygon
# coords = list(ast.literal_eval('(-50,-50),(-50,30),(50,40),(50,-50)'))
# poly = Polygon(coords)
# print(poly.bounds)
fig, ax = plt.subplots(figsize=(15, 8))

df = GeoDataFrame.from_file('../input/example_shp_files/EMSA_draft_AOIs_CleanSeaNet.shp')
df.plot(ax=ax, legend=True, column='Monitoring', cmap='Set1')

world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
#world.plot(ax=ax,color='k',alpha=.3)
world.plot(ax=ax,color='white', edgecolor='black')
poly = GeoDataFrame.from_file('../input/example_shp_files/EMSA_draft_AOIs_CleanSeaNet.shp')
poly.plot(ax=ax,alpha=.4)
# plt.show()

# lon_step = .5
# lat_step = .5
#
# for i in range(9):
#     poly1 = GeoDataFrame.from_file('../input/example_shp_files/EMSA_draft_AOIs_CleanSeaNet.shp')
#     # poly1 = GeoDataFrame.from_file('../input/example_shp_files/EMSA_draft_AOIs_CMS.shp')
#     poly = poly1['geometry'].iloc[i]
#     # Create a grid
#     xmin, xmax, ymin, ymax = poly.bounds[0], poly.bounds[2], poly.bounds[1], poly.bounds[3]
#     print(xmin, xmax, ymin, ymax )
#     xx, yy = np.meshgrid(np.arange(xmin,xmax,lon_step), np.arange(ymin,ymax,lat_step))
#     xc = xx.flatten()
#     yc = yy.flatten()
#     # Chck the ones within the polygon
#     pts = GeoSeries([Point(x, y) for x, y in zip(xc, yc)])
#     in_map = np.array([pts.within(poly)]).sum(axis=0)
#     pts = GeoSeries([val for pos,val in enumerate(pts) if in_map[pos]])
#     pts.plot(alpha=0.2,markersize=2,marker='s',label=str(i)+' '+poly1['Monitoring'].iloc[i],ax=ax)
#     plt.text(pts[1].x,pts[1].y,str(i))
# plt.legend()


# colors=['red','blue','yellow','magenta','green','pink','orange','black','gray','purple']
# for i in range(9):
#     df.loc[[i], 'geometry'].plot(ax=ax,color=colors[i],label=str(i)+' ',legend=True, alpha=0.5)
    # poly1 = GeoDataFrame.from_file('../input/example_shp_files/EMSA_draft_AOIs_CMS.shp')
    # poly = df['geometry'].iloc[i]
    # geopandas.boundary.plot(poly, linewidth=2, color='red',ax=ax)
    #poly.plot(alpha=0.2,markersize=2,marker='s',label=str(i)+' '+poly1['Monitoring'].iloc[i],ax=ax)
    #plt.text(poly[1].x,poly[1].y,str(i))
# poly1.plot(ax=ax,alpha=0.1)
# plt.legend()
plt.xlim((-180,180))
plt.ylim((-90,90))
plt.show()

