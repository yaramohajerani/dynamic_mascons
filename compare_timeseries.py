#!/usr/bin/env python
u"""
compare_timeseries.py
by Yara Mohajerani

Compare timeseries of different mascon solutions

last Update 12/2020
"""
import os
import sys
import pickle
import numpy as np
import netCDF4 as nc
import h5py as h5
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point,Polygon
from shapely.ops import unary_union
import pygplates
from gravity_toolkit.convert_calendar_decimal import convert_calendar_decimal

rad_e = 6.371e8  # -- Average Radius of the Earth [cm]

#------------------------------------------------------------------------------
#-- compare mascon solutions
#------------------------------------------------------------------------------
def compare_timeseries(parameters):
	#-- input file for voronoi regions
	input_file = os.path.expanduser(parameters['VORONOI_FILE'])
	with open(input_file, 'rb') as in_file:
		sv = pickle.load(in_file)
	DDEG_RASTER = float(parameters['DDEG_RASTER'])
	#-- read harmonic parameters
	LMAX = int(parameters['LMAX'])
	#-- get output directory
	ddir = os.path.expanduser(parameters['DIRECTORY'])
	#-- smoothing radius
	RAD = int(parameters['RAD'])
	#-- Set up output labels
	DS = '_FL' if (parameters['DESTRIPE'] in ['Y','y']) else ''
	OCN = '_OCN' if parameters['MASCON_OCEAN'] in ['Y','y'] else ''

	#-- load mascon configuration of interest
	mascon_nums = np.array(parameters['MSCN_NUMS'].split(','),dtype=int)
	mascon_name = parameters['MSCN_NAME']
	lbl = '{0}_{1}'.format(mascon_name,parameters['MSCN_NUMS'].replace(',','+'))

	#-- make polygon union of region of interest
	polys = []
	for i in mascon_nums:
		region = sv.regions[i]
		n = len(region)
		lon_list = np.zeros(n)
		lat_list = np.zeros(n)
		for v in range(n):
			lat_list[v],lon_list[v] = pygplates.PointOnSphere(sv.vertices[region][v]).to_lat_lon()
		polys.append(Polygon(list(zip(lon_list,lat_list))))
	
	#-- combine the polygons for the fixed points
	for i,p in enumerate(polys):
		print("Polygon {0}: area = {1}".format(i,p.area))
	poly_sum = unary_union(polys)

	#----------------------------------------------------------------------
	#-- Read the customized mascon timeseries
	#----------------------------------------------------------------------
	mscn_file = os.path.join(ddir,'MASCON_{0}_YLMS_{1:.2f}DEG_AW13_ICE6G_GA{2}_L{3:02d}_r{4:d}km{5}.txt'.format(lbl,DDEG_RASTER,OCN,LMAX,RAD,DS))
	ts = np.loadtxt(mscn_file)
	vr_tdec = ts[:,1]
	vr_mass = ts[:,2]

	#----------------------------------------------------------------------
	#-- Get Corresponding JPL Mascons
	#----------------------------------------------------------------------
	jpl_file = os.path.expanduser(parameters['JPL_INPUT'])
	#-- read file
	fid = nc.Dataset(jpl_file, 'r')
	jpl_mass = fid.variables['lwe_thickness'][:]
	jpl_lons = fid.variables['lon'][:]
	jpl_lats = fid.variables['lat'][:]
	jpl_time_var = fid.variables['time']
	jpl_dtime = nc.num2date(jpl_time_var[:],jpl_time_var.units)
	fid.close()
	
	#-- get grid increments for area calculation later on
	dth = np.radians(np.abs(jpl_lats[1] - jpl_lats[0]))
	dphi = np.radians(np.abs(jpl_lons[1] - jpl_lons[0]))

	#-- also read land mask to mask out ocean mascons (because we are comparing to GSM harmonics)
	landmask_file = os.path.expanduser(parameters['JPL_LANDMASK'])
	#-- read file
	fid = nc.Dataset(landmask_file, 'r')
	jpl_landocean = fid.variables['land_mask'][:]
	fid.close()
	#-- multiply mass change over ocean by 0
	jpl_mass *= jpl_landocean

	# #-- also read the JPL mascon mask
	df_mask =  gpd.read_file(os.path.expanduser(parameters['JPL_MSCNS']))

	#-- initialize jpl arrays
	jpl_ts = np.zeros(len(jpl_dtime))
	jpl_tdec = np.zeros(len(jpl_dtime))
	#-- add up mascons are inside the polygon union
	for m in range(len(df_mask)):
		#-- include polygons if intersection is more than a quarter of the mascon area
		if poly_sum.intersection(df_mask['geometry'][m]).area > df_mask['geometry'][m].area/4:
			print("Including JPL Mascon {0}".format(m))
			#-- sum up grid poitns that are within df_mask['geometry'][m]
			for i in range(len(jpl_lons)):
				for j in range(len(jpl_lats)):
					if df_mask['geometry'][m].contains(Point(jpl_lons[i],jpl_lats[j])):
						#-- calculate grid area for cmWE to Gt convesion
						th = np.radians(90. - jpl_lats[j])
						#-- area in cm^2
						area = (rad_e**2)*dphi*dth*np.sin(th)
						#-- convert cmWE * cm^2 to get grams, then convert to Gt
						jpl_ts += (jpl_mass[:,j,i]*area)/1e15
	#-- convert dates to decimal years
	for t in range(len(jpl_dtime)):
		jpl_tdec[t] = convert_calendar_decimal(jpl_dtime[t].year,jpl_dtime[t].month,DAY=jpl_dtime[t].day,\
			HOUR=jpl_dtime[t].hour,MINUTE=jpl_dtime[t].minute,SECOND=jpl_dtime[t].second)
	
	#----------------------------------------------------------------------
	#-- Read Corresponding Goddard Mascons
	#----------------------------------------------------------------------
	gsfc_file = os.path.expanduser(parameters['GSFC_INPUT'])
	f = h5.File(gsfc_file,'r')
	N = int(f['size']['N_mascons'][()])
	gsfc_mass = np.squeeze(np.array(f['solution']['cmwe']))
	gsfc_lons = np.squeeze(np.array(f['mascon']['lon_center']))
	gsfc_lats = np.squeeze(np.array(f['mascon']['lat_center']))
	gsfc_area = np.squeeze(np.array(f['mascon']['area_km2']))
	gsfc_y,gsfc_doy,gsfc_tdec = np.array(f['time']['yyyy_doy_yrplot_middle'])
	#-- initialize mass timeseries
	gsfc_ts = np.zeros(len(gsfc_tdec))
	for i in range(N):
		#-- count mascon if it's center is in polygon
		if poly_sum.contains(Point(gsfc_lons[i],gsfc_lats[i])):
			#-- convert cmWE * km^2 * 10^10 to get grams, then divide by 10^15 convert to Gt
			#-- which is equivalent to dividing by 10^5
			gsfc_ts += (gsfc_mass[i,:]*gsfc_area[i])/1e5
	f.close()

	#----------------------------------------------------------------------
	#-- Plot Comparison
	#----------------------------------------------------------------------
	fig = plt.figure(1,figsize=(9,6))
	plt.plot(vr_tdec,vr_mass,color='k',label='Voronoi Mascons',zorder=1)
	plt.plot(jpl_tdec,jpl_ts,color='red',label='JPL Mascons',zorder=2)
	plt.plot(gsfc_tdec,gsfc_ts,color='cyan',label='GSFC Mascons',zorder=3)
	plt.axvspan(2017.5, 2018.37, color='lightgray',zorder=4)
	plt.legend()
	plt.xlabel('Time [Decimal Year]')
	plt.ylabel('Mass [Gt]')
	plt.title('Comparing Mass Balance Time-Series for Mascons {0} ({1})'.format(\
		parameters['MSCN_NUMS'].replace(',','+'),mascon_name.replace('_',' ')))
	plt.savefig(mscn_file.replace('.txt','_comparison.pdf'),format='PDF')
	plt.close(fig)
	
#------------------------------------------------------------------------------
#-- main function
#------------------------------------------------------------------------------
def main():
	if len(sys.argv) == 1:
		sys.exit('No paramter file given')
	else:
		#-- read input files
		input_files = sys.argv[1:]
		parameters = {}
		for infile in input_files:
			#-- for each paramter file, extract parameters
			fid = open(infile, 'r')
			for fileline in fid:
				part = fileline.split()
				parameters[part[0]] = part[1]
			fid.close()
			#-- feed parameters to function to compare mascon solutions
			compare_timeseries(parameters)

#------------------------------------------------------------------------------
#-- run main program
#------------------------------------------------------------------------------
if __name__ == '__main__':
	main()
