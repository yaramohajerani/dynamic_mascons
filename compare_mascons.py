#!/usr/bin/env python
u"""
compare_mascons.py
by Yara Mohajerani

Compare the placement of mascons in different solutions

last Update 12/2020
"""
import os
import sys
import getopt
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

#------------------------------------------------------------------------------
#-- compare mascon solutions
#------------------------------------------------------------------------------
def compare_mascons(parameters,jpl,gsfc):
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

	#----------------------------------------------------------------------
	#-- load JPL mascons
	#----------------------------------------------------------------------
	if jpl:
		jpl_gdf = gpd.read_file(os.path.expanduser(parameters['JPL_MSCNS']))
	#----------------------------------------------------------------------
	#-- load GSFC mascons
	#----------------------------------------------------------------------
	if gsfc:
		gsfc_gdf = gpd.read_file(os.path.expanduser(parameters['GSFC_MSCNS']))

	#----------------------------------------------------------------------
	#-- Load RGI polygons
	#----------------------------------------------------------------------
	rgi_dir = os.path.expanduser(parameters['RGI_DIR'])
	rgi = {}
	if 'alaska' in mascon_name:
		rgi['a'] = gpd.read_file(os.path.join(rgi_dir,'01_rgi60_Alaska','01_rgi60_Alaska.shp'))
	else:
		for key,num,name in zip(['c','w','e'],[13,14,15],['CentralAsia','SouthAsiaWest','SouthAsiaEast']):
			rgi[key] = gpd.read_file(os.path.join(rgi_dir,'{0:d}_rgi60_{1}'.format(num,name),'{0:d}_rgi60_{1}.shp'.format(num,name)))

	#----------------------------------------------------------------------
	#-- Make plot of mascon polygons
	#----------------------------------------------------------------------
	fig, ax = plt.subplots(1,1,figsize = (9,6),sharey=True, dpi=150)
	#-- plot continents 
	#-- load in world map for plotting in background
	world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
	world.plot(ax=ax,fc='none',ec='black',linewidth=1.4,rasterized=True)
	#-- add RGI polygons
	for r in rgi.keys():
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
				#-- plot the polygon edges
				ax.plot(edge_lons,edge_lats,color='limegreen',linewidth=1.,alpha=0.6)
			if i in mascon_nums:
				lon_list = np.zeros(n)
				lat_list = np.zeros(n)
				for v in range(n):
					lat_list[v],lon_list[v] = pygplates.PointOnSphere(sv.vertices[region][v]).to_lat_lon()
				polys.append(Polygon(list(zip(lon_list,lat_list))))
			#-- if not plotting other mascon solutions, add numbers
			if not (jpl and gsfc):
				cnt_lat,cnt_lon = pygplates.PointOnSphere(sv.points[i]).to_lat_lon()
				shift = int(np.log10(i))/2 if (i > 0) else 0
				cnt_lon -= shift
				cnt_lat -= 0.5
				ax.text(cnt_lon,cnt_lat,i,fontsize=10,color='red',weight='bold')
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
	if jpl:
		for i in range(len(jpl_gdf)):
			#-- include polygons if intersection is more than a quarter of the mascon area
			# if poly_sum.intersects(jpl_gdf['geometry'][i]):
			if poly_sum.intersection(jpl_gdf['geometry'][i]).area > jpl_gdf['geometry'][i].area/4:
				print("Including JPL Mascon #:",i)
				patch = PolygonPatch(jpl_gdf['geometry'][i],ec='blue',fc='dodgerblue',linewidth=1.5,alpha=0.3,zorder=3)
				jplpatch = ax.add_patch(patch)
		lgd_items.append(jplpatch)
		lgd_lbls.append('JPL Mascons')
		# jpl_gdf.plot(ax=ax,fc='none',ec='aqua',linewidth=1.,rasterized=True)

	#----------------------------------------------------------------------
	#-- plot GSFC mascons
	#----------------------------------------------------------------------
	if gsfc:
		for i in range(len(gsfc_gdf)):
			#-- include polygons if intersection is more than half of the mascon area
			# if poly_sum.intersects(gsfc_gdf['geometry'][i]):
			if poly_sum.intersection(gsfc_gdf['geometry'][i]).area > gsfc_gdf['geometry'][i].area/2:
				print("Including GSFC Mascon #:",i)
				patch = PolygonPatch(gsfc_gdf['geometry'][i],ec='darkred',fc='red',linewidth=1.5,alpha=0.2,zorder=4)
				gsfcpatch = ax.add_patch(patch)
		lgd_items.append(gsfcpatch)
		lgd_lbls.append('GSFC Mascons')
		# gsfc_gdf.plot(ax=ax,fc='none',ec='seagreen',linewidth=1.,rasterized=True)

	ax.set_xlabel('Longitude (Degrees)')
	ax.set_ylabel('Latitude (Degrees)')
	ax.set_xlim([lon0,lon1])
	ax.set_ylim([lat0,lat1])
	frame = plt.legend(lgd_items,lgd_lbls,fancybox=True, framealpha=1.0)
	plt.tight_layout()
	#-- save plot
	jpl_str = '_JPL' if jpl else ''
	gsfc_str = '_GSFC' if gsfc else ''
	outfile = os.path.join(os.path.dirname(input_file),'mascon_comparison_{0}{1}{2}.png'.format(mascon_name,jpl_str,gsfc_str))
	plt.savefig(outfile,format='PNG')
	plt.close(fig)

#------------------------------------------------------------------------------
#-- main function
#------------------------------------------------------------------------------
def main():
	#-- read parameter files and commandline arguments
	optlist, input_files=getopt.getopt(sys.argv[1:],'JG',['JPL','GSFC'])

	#-- command line parameters
	jpl = False
	gsfc = False
	for opt, arg in optlist:
		if opt in ('-J','--JPL'):
			jpl = True
		if opt in ('-G','--GSFC'):
			gsfc = True
	
	if not jpl:
		print('Not including JPL mascons')
	if not gsfc:
		print('Not including GSFC mascons')

	if not input_files:
		sys.exit('No parameter file given')
	else:
		parameters = {}
		for infile in input_files:
			#-- for each paramter file, extract parameters
			fid = open(infile, 'r')
			for fileline in fid:
				part = fileline.split()
				parameters[part[0]] = part[1]
			fid.close()
			#-- feed parameters to function to compare mascon solutions
			compare_mascons(parameters,jpl,gsfc)

#------------------------------------------------------------------------------
#-- run main program
#------------------------------------------------------------------------------
if __name__ == '__main__':
	main()
