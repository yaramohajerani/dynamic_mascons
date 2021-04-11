#!/usr/bin/env python
u"""
plot_mascon_placement.py
by Yara Mohajerani

Compare the placement of mascons in different solutions

last Update 04/2021
"""
import os
import pickle
from matplotlib.patches import Polygon
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point,Polygon
from shapely.ops import unary_union
from descartes import PolygonPatch
from scipy.spatial import geometric_slerp
import pygplates

regions = ['karakoram_NW','karakoram_SE','nyainqentangla','alaska']

#-- initialize figure
fig, axs = plt.subplots(2,2,figsize = (11,6), dpi=150)
#-- load in world map for plotting in background
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

#-- loop over regions and plot each
for count,reg in enumerate(regions):
	print(reg)
	#-- set up axis for figure
	axi = int(count/2) # row number
	axj = count%2	# column number
	ax = axs[axi,axj]

	param_file = os.path.join(os.getcwd(),'..','parameters_voronoi_mascons_{0}.txt'.format(reg))
	parameters = {}
	#-- for each paramter file, extract parameters
	fid = open(param_file, 'r')
	for fileline in fid:
		part = fileline.split()
		parameters[part[0]] = part[1]
	fid.close()

	#-- input file for voronoi regions
	input_file = os.path.expanduser(parameters['VORONOI_FILE'])
	with open(input_file, 'rb') as in_file:
		sv = pickle.load(in_file)

	#-- get the plot limits
	lat0,lat1 = np.array(parameters['LAT_LIMITS'].split(','),dtype=float)
	lon0,lon1 = np.array(parameters['LON_LIMITS'].split(','),dtype=float)
	#-- make polygon of the frame to check what point are contained in it
	frame_poly = Polygon([[lon0,lat0],[lon0,lat1],[lon1,lat1],[lon1,lat0]])

	#-- load mascon configuration of interest
	mascon_nums = np.array(parameters['MSCN_NUMS'].split(','),dtype=int)
	mascon_name = parameters['MSCN_NAME']

	#-- load data only on the first iteration
	if count == 0:
		#----------------------------------------------------------------------
		#-- load JPL mascons
		#----------------------------------------------------------------------
		jpl_gdf = gpd.read_file(os.path.expanduser(parameters['JPL_MSCNS']))
		#----------------------------------------------------------------------
		#-- load GSFC mascons
		#----------------------------------------------------------------------
		gsfc_gdf = gpd.read_file(os.path.expanduser(parameters['GSFC_MSCNS']))

		#----------------------------------------------------------------------
		#-- Load RGI polygons
		#----------------------------------------------------------------------
		rgi_dir = os.path.expanduser(parameters['RGI_DIR'])
		rgi = {}
		for key,num,name in zip(['c','w','e','a'],[13,14,15,1],['CentralAsia','SouthAsiaWest','SouthAsiaEast','Alaska']):
			rgi[key] = gpd.read_file(os.path.join(rgi_dir,'{0:02d}_rgi60_{1}'.format(num,name),'{0:02d}_rgi60_{1}.shp'.format(num,name)))

	#----------------------------------------------------------------------
	#-- Make plot of mascon polygons
	#----------------------------------------------------------------------
	#-- plot continents 
	world.plot(ax=ax,fc='none',ec='black',linewidth=1.4,rasterized=True)
	#-- add RGI polygons
	if 'alaska' in reg:
		rgi['a'].plot(ax=ax,alpha=0.3,fc='darkkhaki',ec='darkkhaki',rasterized=True)
	else:
		for r in ['c','w','e']:
			rgi[r].plot(ax=ax,alpha=0.3,fc='darkkhaki',ec='darkkhaki',rasterized=True)
	polys = []
	#-- plot regions
	for i,region in enumerate(sv.regions):
		#-- check if polygon is within plotting frame
		pt_lat,pt_lon = pygplates.PointOnSphere(sv.vertices[region][0]).to_lat_lon()
		pt = Point(pt_lon,pt_lat)
		if pt.within(frame_poly):
			n = len(region)
			t = np.linspace(0,1,10)
			for j in range(n):
				start = sv.vertices[region][j]
				end = sv.vertices[region][(j + 1) % n]
				edge = np.array(geometric_slerp(start, end, t))
				edge_lons,edge_lats = np.ones(len(edge)),np.ones(len(edge))
				for v in range(len(edge)):
					edge_lats[v],edge_lons[v] = pygplates.PointOnSphere(edge[v]).to_lat_lon()
				#-- plot the polygon edges if they don't go around the earth making a long line
				if not (np.count_nonzero(edge_lons<0) > 0 and np.count_nonzero(edge_lons>0) > 0):
					ax.plot(edge_lons,edge_lats,color='limegreen',linewidth=1.)
			if i in mascon_nums:
				lon_list = np.zeros(n)
				lat_list = np.zeros(n)
				for v in range(n):
					lat_list[v],lon_list[v] = pygplates.PointOnSphere(sv.vertices[region][v]).to_lat_lon()
				polys.append(Polygon(list(zip(lon_list,lat_list))))
	#-- combine the polygons for the fixed points
	for i,p in enumerate(polys):
		print("Polygon {0}: area = {1}".format(i,p.area))
	poly_sum = unary_union(polys)
	print("Union area = {0}".format(poly_sum.area))
	pp = PolygonPatch(poly_sum,ec='limegreen',fc='lime',alpha=0.6,zorder=2)
	vpatch = ax.add_patch(pp)
	#-- add patch to legend
	lgd_items = [vpatch]
	lgd_lbls = ['Voronoi Mascons']

	#----------------------------------------------------------------------
	#-- plot JPL mascons
	#----------------------------------------------------------------------
	print(len(jpl_gdf))
	for i in range(len(jpl_gdf)):
		#-- include polygons if intersection is more than a quarter of the mascon area
		if poly_sum.intersection(jpl_gdf['geometry'][i]).area > jpl_gdf['geometry'][i].area/4:
			print("Including JPL Mascon #:",i)
			patch = PolygonPatch(jpl_gdf['geometry'][i],ec='blue',fc='dodgerblue',linewidth=1.5,alpha=0.3,zorder=3)
			jplpatch = ax.add_patch(patch)
	lgd_items.append(jplpatch)
	lgd_lbls.append('JPL Mascons')

	#----------------------------------------------------------------------
	#-- plot GSFC mascons
	#----------------------------------------------------------------------
	for i in range(len(gsfc_gdf)):
		#-- include polygons if intersection is more than half of the mascon area
		# if poly_sum.intersects(gsfc_gdf['geometry'][i]):
		if poly_sum.intersection(gsfc_gdf['geometry'][i]).area > gsfc_gdf['geometry'][i].area/2:
			print("Including GSFC Mascon #:",i)
			patch = PolygonPatch(gsfc_gdf['geometry'][i],ec='darkred',fc='red',linewidth=1.5,alpha=0.2,zorder=4)
			gsfcpatch = ax.add_patch(patch)
	lgd_items.append(gsfcpatch)
	lgd_lbls.append('GSFC Mascons')

	#-- add title and labels
	ax.set_title(reg.replace('_',' ').upper())
	if axi == 1:
		ax.set_xlabel('Longitude (Degrees)')
	if axj == 0:
		ax.set_ylabel('Latitude (Degrees)')
	ax.set_xlim([lon0,lon1])
	ax.set_ylim([lat0,lat1])
	ax.set_aspect('equal')

#-- add legend to bottom
fig.subplots_adjust(bottom=0.12,top=0.95)
fig.legend(lgd_items,lgd_lbls, fancybox=True, framealpha=1.0, loc="lower center", 
			bbox_to_anchor=(0.5, 0.0),ncol=3)
# plt.tight_layout()
#-- save plot
outfile = os.path.join(os.getcwd(),'..','data','paper','mascon_placement_comparison.pdf')
plt.savefig(outfile,format='PDF')
plt.close(fig)