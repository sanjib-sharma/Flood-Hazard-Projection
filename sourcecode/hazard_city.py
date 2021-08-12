import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import shapefile as shp
import numpy as np
from matplotlib import rcParams
from functions_plot_map import plot_basemap, plot_border
import gdal
from matplotlib.patches import Circle
from matplotlib import cm
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D

# Import PA cities shapefile
sf1 = shp.Reader('../data/shapefile/pa_proj/pa_proj.shp')
recs1 = sf1.records()
shapes1 = sf1.shapes()
Nshp1 = len(shapes1)
# Import PA state shapefile
sf2 = shp.Reader('../data/shapefile/PAState2019_01_Project.shp')
recs2 = sf2.records()
shapes2 = sf2.shapes()
Nshp2 = len(shapes2)
# Import floodmap raster
floodmap = gdal.Open('../flowraster/raster_proj2.tif', gdal.GA_ReadOnly)
data = floodmap.ReadAsArray()
data[data<0] = np.nan
gt = floodmap.GetGeoTransform()
proj = floodmap.GetProjection()
xsize = floodmap.RasterXSize
ysize = floodmap.RasterYSize
xres = gt[1]
yres = gt[5]
## get the edge coordinates and add half the resolution 
## to go to center coordinates
xmin = gt[0] + xres * 0.5
xmax = gt[0] + (xres * xsize) - xres * 0.5
ymin = gt[3] + (yres * ysize) + yres * 0.5
ymax = gt[3] - yres * 0.5
## create a grid of xy coordinates in the original projection
x, y = np.mgrid[xmin:xmax+xres:xres, ymax+yres:ymin:yres]
'''
# Import landcover raster
landcover = gdal.Open('../landcover_raster/lc_wgs84_wmas/w001001.adf', gdal.GA_ReadOnly)
data2 = landcover.ReadAsArray()
data2 = data2.astype(float)
data2[data2 < 583] = np.NaN
data2[data2 > 583] = np.NaN
#data2 = data2.astype(int)
gt = landcover.GetGeoTransform()
proj = landcover.GetProjection()
xsize = landcover.RasterXSize
ysize = landcover.RasterYSize
xres = gt[1]
yres = gt[5]
## get the edge coordinates and add half the resolution 
## to go to center coordinates
xmin = gt[0] + xres * 0.5
xmax = gt[0] + (xres * xsize) - xres * 0.5
ymin = gt[3] + (yres * ysize) + yres * 0.5
ymax = gt[3] - yres * 0.5
## create a grid of xy coordinates in the original projection
x2, y2 = np.mgrid[xmin:xmax+xres:xres, ymax+yres:ymin:yres]
'''
# Import hazard data
hazard = pd.read_csv('../data/hazard/hazardranking.txt', delim_whitespace=True, converters={'ID':lambda x: '00'+x if len(x)==3 else '0'+x if len(x)==4 else x})

# Import exposure data
exposure = pd.read_csv('../data/exposure/exposureranking.txt', delim_whitespace=True, converters={'ID':lambda x: '00'+x if len(x)==3 else '0'+x if len(x)==4 else x})

# Import centroids of city shapefiles
centroids = pd.read_csv('../data/centroids/centroids.csv', usecols=['FIPS_MUN_C', 'Longitude', 'Latitude'], converters={
'FIPS_MUN_C':lambda x: '00'+x if len(x)==3 else '0'+x if len(x)==4 else x}).rename(index=str, columns={'FIPS_MUN_C':'ID'})

# Merge centroids to hazard data
cols = ['ID', 'RATIO2018', 'RATIO2099', 'RANK2018', 'RANK2099']
cols_dict_exp = {'RATIO2018':'exp_2018', 'RATIO2099':'exp_2099', 'RANK2018':'exp_rank_2018', 'RANK2099':'exp_rank_2099'}
cols_dict_haz = {'RATIO2018':'haz_2018', 'RATIO2099':'haz_2099', 'RANK2018':'haz_rank_2018', 'RANK2099':'haz_rank_2099'}
merged = hazard[cols].rename(index=str, columns=cols_dict_haz).merge(exposure[cols].rename(index=str, 
columns=cols_dict_exp), how='left', on='ID').merge(centroids, how= 'left', on='ID')
# Remove repeated cities and boroughs in data
repeated_ids = ['00364', '23304', '25136', '46072', '70352', '77272', '06088']
repeated_indexes = merged[merged['ID'].isin(repeated_ids)].index
repeated_ind_del = [1, 223, 242, 401, 647, 717, 819]
merged = merged.drop(index=repeated_ind_del)
haz99 = merged[['Longitude', 'Latitude', 'haz_2099']].copy()
#haz99.to_csv('haz99_points.csv')

# Subset of merged dataframe that are cities
cities_ID = [i[3] for i in recs1 if i[11]=='CITY']
cities_NM = [i[2] for i in recs1 if i[11]=='CITY']
# remove repeated cities in data
cities_ID.pop(3)
cities_NM.pop(3)

city_df = pd.DataFrame({'ID':cities_ID, 'name':cities_NM})
city_df['name'] = city_df['name'].str.title()
merged_city = merged.merge(city_df, how='left', on='ID').dropna()
haz99_city = merged_city[['ID', 'name', 'Longitude', 'Latitude', 'haz_2099']].copy()
exp99_city = merged_city[['ID', 'name', 'Longitude', 'Latitude', 'exp_2099']].copy()
#haz99_city.to_csv('haz99_city_points.csv')
#exp99_city.to_csv('exp99_city_points.csv')

# Subset of merged dataframe that are boroughs
boro_ID = [i[3] for i in recs1 if i[11]=='BORO']
boro_NM = [i[2] for i in recs1 if i[11]=='BORO']
# remove repeated cities in data
elements_delete = [473, 630, 947, 948, 949, 952]
for i in elements_delete:
	boro_ID.pop(i)
	boro_NM.pop(i)

boro_df = pd.DataFrame({'ID':boro_ID, 'name':boro_NM})
merged_boro = merged.merge(boro_df, how='left', on='ID').dropna()
haz99_boro = merged_boro[['ID', 'name', 'Longitude', 'Latitude', 'haz_2099']].copy()
exp99_boro = merged_boro[['ID', 'name', 'Longitude', 'Latitude', 'exp_2099']].copy()
#haz99_boro.to_csv('haz99_boro_points.csv')
#exp99_boro.to_csv('exp99_boro_points.csv')

'''
# Sort dataframe in the order of the NRID_regions in the shapefile
cities_shp = [recs1[i][3] for i in range(Nshp1)]
merged['ID_shp'] = pd.Categorical(merged['ID_shp'], categories=[int(i) for i in cities_shp], ordered=True)
merged = merged.sort_values('ID_shp')
'''

# Plots maps
# Set up the figure specifications
golden_mean = (np.sqrt( 5 ) - 1.0 ) / 2.0 # Aesthetic ratio
fig_width = 4.75 # width in inches
fig_height = fig_width * golden_mean # height in inches
#fig_height = fig_width
fig_size =  [ fig_width, fig_height ]

params = {'backend': 'ps',
'font.family': 'serif',
'axes.labelsize': 10,
'axes.titlesize': 10,
'legend.fontsize': 5,
'xtick.labelsize': 8,
'ytick.labelsize': 8,
#'text.usetex': True,
'figure.dpi': 320,
'figure.figsize': fig_size}
rcParams.update(params)

fig, ax = plot_basemap(shapes=shapes2, recs=recs2, color='lightgrey', lw=1)
plot_border(ax, shapes2, recs2, zpos=1.5, lw=1, ec='k', fill=False)
blues1 = cm.get_cmap('Blues', 512)
blues2 = ListedColormap(blues1(np.linspace(0.4, 1.0, 256)))
ax.pcolormesh(x, y, data[:,:].T, cmap=blues2, alpha=1)
'''
ax.pcolormesh(x2, y2, data2[:,:].T, cmap=cm.get_cmap('Reds_r'), alpha=0.5)
'''
list_circles = haz99_city[['Latitude', 'Longitude', 'haz_2099']].dropna().values
list_names  = merged_city[merged_city['name'].isin(['Connellsville', 'Johnstown', 'Warren', 'Bradford', 'Lock Haven', 
'Williamsport', 'Sunbury', 'York', 'Wilkes-Barre', 'Bethlehem', 'Titusville'])][['Latitude', 'Longitude', 'name']].values
list_delta_xy = [[-30000, 15000], [0, -15000], [60000, 0], [50000, 0], [0, -15000], [0, 15000], [0, -15000], [0, -15000],
[50000, -15000], [0, 15000], [0, 15000]]
for lat, lon, var in list_circles:
	if var<25:
		c=Circle((lon, lat), radius=5000, alpha=0.8, zorder=2, edgecolor='grey', linewidth=0.3, facecolor='darkgrey')
	elif var<50:
		c=Circle((lon, lat), radius=5000, alpha=0.8, zorder=2, edgecolor='grey', linewidth=0.3, facecolor='grey')
	elif var<75:
		c=Circle((lon, lat), radius=5000, alpha=0.8, zorder=2, edgecolor='grey', linewidth=0.3, facecolor='lightcoral')
	else:
		c=Circle((lon, lat), radius=5000, alpha=0.8, zorder=2, edgecolor='grey', linewidth=0.3, facecolor='red')
	a = ax.add_patch(c)
	for (lt,ln,name) ,(dx,dy) in zip(list_names, list_delta_xy):
		txt = ax.text(ln+dx, lt+dy, name, fontsize=6, horizontalalignment='center', verticalalignment='center')

legend_elements = [Line2D([0], [0], marker='o', color='w', label='0-25', markerfacecolor='darkgrey', markersize=6, 
markeredgecolor='grey'),
Line2D([0], [0], marker='o', color='w', label='25-50', markerfacecolor='grey', markersize=6, 
markeredgecolor='grey'),
Line2D([0], [0], marker='o', color='w', label='50-75', markerfacecolor='lightcoral', markersize=6, 
markeredgecolor='grey'),
Line2D([0], [0], marker='o', color='w', label='75-100', markerfacecolor='red', markersize=6, 
markeredgecolor='grey')]
lgd = ax.legend(handles=legend_elements, loc='lower center', title='Hazard [%]', ncol=4, bbox_to_anchor=(0.4, -0.20, 0.2, 0.05),
fontsize=10, frameon=False, columnspacing=0)
lgd.get_title().set_fontsize('10')
fig.savefig('../maps/final_maps/hazard2099_city_final.png', dpi=600, bbox_inches='tight')


