#!/usr/bin/env python
u"""
plot_timeseries_comparisons.py
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
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point,Polygon
from shapely.ops import unary_union
from PyAstronomy import pyasl
import pygplates
from gravity_toolkit.tsregress import tsregress

rad_e = 6.371e8  # -- Average Radius of the Earth [cm]

regions = ['karakoram_NW','karakoram_SE','nyainqentangla','alaska']

#-- initialize figure
fig, axs = plt.subplots(2,2,figsize = (10,6))

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
	#-- GIA
	if parameters['GIA'] == 'ICE6G-D':
		GIA = 'ICE6G-D_VM5A_O512' 
	elif parameters['GIA'] == 'AW13-ICE6G':
		GIA = 'AW13-ICE6G_GA'
	else:
		GIA = ''

	#-- overwriting time-series txt files
	CLOBBER = True if (parameters['CLOBBER'] in ['Y','y']) else False

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

	#-- create dictionaries for time and mass
	tdec = {}
	mass = {}
	errs = {}
	#----------------------------------------------------------------------
	#-- Read the customized mascon timeseries
	#----------------------------------------------------------------------
	mscn_file = os.path.join(ddir,'MASCON_{0}_YLMS_{1:.2f}DEG_{2}{3}_L{4:02d}_r{5:d}km{6}.txt'.format(lbl,DDEG_RASTER,GIA,OCN,LMAX,RAD,DS))
	ts = np.loadtxt(mscn_file)
	tdec['vor'] = ts[:,1]
	mass['vor'] = ts[:,2]
	errs['vor'] = ts[:,3]
	#-- also read leakage error
	leak_file = os.path.join(ddir,'MASCON_{0}_LEAKAGE_ERROR_{1:.2f}DEG_SKERNEL{2}_L{3:d}_r{4:d}km.csv'.format(lbl,DDEG_RASTER,OCN,LMAX,RAD))
	df_leak = pd.read_csv(leak_file)
	leak_err = float(df_leak['leakage'])

	#----------------------------------------------------------------------
	#-- Get Corresponding JPL Mascons
	#----------------------------------------------------------------------
	#-- if timeseries file already exists, read it. Otherwise create it
	jpl_ts_file = os.path.join(ddir,'MASCON_{0}_JPL_timeseries.txt'.format(lbl))
	if CLOBBER or (not os.path.exists(jpl_ts_file)):
		jpl_file = os.path.expanduser(parameters['JPL_INPUT'])
		#-- read file
		fid = nc.Dataset(jpl_file, 'r')
		jpl_mass = fid.variables['lwe_thickness'][:]
		jpl_err = fid.variables['uncertainty'][:]
		jpl_lons = fid.variables['lon'][:]
		jpl_lats = fid.variables['lat'][:]
		jpl_time_var = fid.variables['time']
		jpl_dtime = nc.num2date(jpl_time_var[:],jpl_time_var.units)
		fid.close()

		#-- shift lons from 0 to 360 to -180 to 180
		ind_lon = np.where(jpl_lons > 180)
		jpl_lons[ind_lon] -= 360
		
		#-- get grid increments for area calculation later on
		dth = np.radians(np.abs(jpl_lats[1] - jpl_lats[0]))
		dphi = np.radians(np.abs(jpl_lons[1] - jpl_lons[0]))

		#-- also read land mask to mask out ocean mascons (because we are comparing to GSM harmonics)
		landmask_file = os.path.expanduser(parameters['JPL_LANDMASK'])
		#-- read file
		fid = nc.Dataset(landmask_file, 'r')
		jpl_landocean = fid.variables['land_mask'][:]
		fid.close()

		#-- scale file
		# scale_file = os.path.expanduser(parameters['JPL_SCALE'])
		# #-- read file
		# fid = nc.Dataset(scale_file, 'r')
		# jpl_scale = fid.variables['scale_factor'][:]
		# fid.close()

		#-- multiply mass change over ocean by 0
		jpl_mass *= jpl_landocean #* jpl_scale

		#-- also read the JPL mascon mask
		df_mask =  gpd.read_file(os.path.expanduser(parameters['JPL_MSCNS']))

		#-- initialize jpl arrays
		mass['jpl'] = np.zeros(len(jpl_dtime))
		tdec['jpl'] = np.zeros(len(jpl_dtime))
		errs['jpl'] = np.zeros(len(jpl_dtime))
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
							mass['jpl'] += (jpl_mass[:,j,i]*area)/1e15
							errs['jpl'] += np.abs(jpl_err[:,j,i]*area/1e15)
		# 					errs['jpl'] += (jpl_err[:,j,i]*area/1e15)**2
		# errs['jpl'] = np.sqrt(errs['jpl'])
		#-- convert dates to decimal years
		for t in range(len(jpl_dtime)):
			tdec['jpl'][t] = pyasl.decimalYear(jpl_dtime[t])
		
		#-- write JPL timeseris to file
		fid = open(jpl_ts_file,'w')
		fid.write('Decimal-Years   Mass(Gt)   Errs(Gt)\n')
		for t in range(len(tdec['jpl'])):
			fid.write('{0:.4f}  {1:.4f}  {2:.4f}\n'.format(tdec['jpl'][t],mass['jpl'][t],errs['jpl'][t]))
		fid.close()
	else:
		print("Reading JPL timeseries from file.")
		#-- read file
		jpl_data = np.loadtxt(jpl_ts_file,skiprows=1)
		tdec['jpl'] = jpl_data[:,0]
		mass['jpl'] = jpl_data[:,1]
		errs['jpl'] = jpl_data[:,2]

	#----------------------------------------------------------------------
	#-- Read Corresponding Goddard Mascons
	#----------------------------------------------------------------------
	#-- if timeseries file already exists, read it. Otherwise create it
	gsfc_ts_file = os.path.join(ddir,'MASCON_{0}_GSFC_timeseries.txt'.format(lbl))
	if CLOBBER or (not os.path.exists(gsfc_ts_file)):
		gsfc_file = os.path.expanduser(parameters['GSFC_INPUT'])
		f = h5.File(gsfc_file,'r')
		N = int(f['size']['N_mascons'][()])
		gsfc_mass = np.squeeze(np.array(f['solution']['cmwe']))
		gsfc_lons = np.squeeze(np.array(f['mascon']['lon_center']))
		gsfc_lats = np.squeeze(np.array(f['mascon']['lat_center']))
		gsfc_area = np.squeeze(np.array(f['mascon']['area_km2']))
		gsfc_locs = np.squeeze(np.array(f['mascon']['location']))
		gsfc_leak = np.squeeze(np.array(f['uncertainty']['leakage_2sigma']))
		gsfc_noise = np.squeeze(np.array(f['uncertainty']['noise_2sigma']))
		gsfc_y,gsfc_doy,tdec['gsfc'] = np.array(f['time']['yyyy_doy_yrplot_middle'])
		#-- shift lons from 0 to 360 to -180 to 180
		ind_lon = np.where(gsfc_lons > 180)
		gsfc_lons[ind_lon] -= 360
		#-- initialize mass timeseries
		mass['gsfc'] = np.zeros(len(tdec['gsfc']))
		errs['gsfc'] = np.zeros(len(tdec['gsfc']))
		for i in range(N):
			#-- count mascon if it's center is in polygon and if it's on land (comparable to GSM)
			if (gsfc_locs[i] != 90 and poly_sum.contains(Point(gsfc_lons[i],gsfc_lats[i]))):
				#-- convert cmWE * km^2 * 10^10 to get grams, then divide by 10^15 convert to Gt
				#-- which is equivalent to dividing by 10^5
				mass['gsfc'] += (gsfc_mass[i,:]*gsfc_area[i])/1e5
				errs['gsfc'] += np.sqrt(gsfc_leak[i]**2 + gsfc_noise[i,:]**2)*gsfc_area[i]/1e5
		f.close()

		#-- write GSFC timeseris to file
		fid = open(gsfc_ts_file,'w')
		fid.write('Decimal-Years   Mass(Gt)\n')
		for t in range(len(tdec['gsfc'])):
			fid.write('{0:.4f}  {1:.4f}  {2:.4f}\n'.format(tdec['gsfc'][t],mass['gsfc'][t],errs['gsfc'][t]))
		fid.close()
	else:
		print("Reading GSFC time-series from file.")
		#-- read file
		gsfc_data = np.loadtxt(gsfc_ts_file,skiprows=1)
		tdec['gsfc'] = gsfc_data[:,0]
		mass['gsfc'] = gsfc_data[:,1]
		errs['gsfc'] = gsfc_data[:,2]

	#----------------------------------------------------------------------
	#-- Perform regression for common period
	#----------------------------------------------------------------------
	#-- list of solution keys
	sols = tdec.keys()
	#-- get min and max time
	tmin = np.max([np.min(tdec[s]) for s in sols])
	tmax = np.min([np.max(tdec[s]) for s in sols])
	com_lbl = '{0:.1f}-{1:.1f}'.format(tmin,tmax)
	#-- get indices of common period
	ind = {}
	for s in sols:
		ind[s] = np.squeeze(np.nonzero(
					(tdec[s] > tmin) &
					(tdec[s] < tmax)
					))
		#-- remove common mean
		mass[s] -= np.mean(mass[s][ind[s]])
	"""
	df = {}
	for s in sols:
		df[s] = {}
		fit = tsregress(tdec[s][ind[s]],mass[s][ind[s]],ORDER=1,CYCLES=[0.5,1])
		df[s]['Trend ({0}) Gt/yr'.format(com_lbl)] = fit['beta'][1]
		df[s]['Trend Err ({0}) Gt/yr'.format(com_lbl)] = fit['error'][1]
		df[s]['Seasonal ({0}) Gt'.format(com_lbl)] = np.sqrt(fit['beta'][4]**2 + fit['beta'][5]**2)
		df[s]['Seasonal Err ({0}) Gt'.format(com_lbl)] = np.sqrt(fit['error'][4]**2 + fit['error'][5]**2)
		#-- also get trend for 2015 onwards
		ind15 = np.squeeze(np.nonzero(tdec[s]>2015))
		fit = tsregress(tdec[s][ind15],mass[s][ind15],ORDER=1,CYCLES=[0.5,1])
		df[s]['Trend (2015+) Gt/yr'.format(com_lbl)] = fit['beta'][1]
		df[s]['Trend Err (2015+) Gt/yr'.format(com_lbl)] = fit['error'][1]

	#-- write regression results to file
	df = pd.DataFrame(df)
	df.to_csv(mscn_file.replace('.txt','_comparison_regession.csv'))
	"""
	#-- calculate error timeseries for the voronoi timeseries
	err_line = np.sqrt(errs['vor']**2 + (0.01*leak_err*mass['vor'])**2)
	#----------------------------------------------------------------------
	#-- Plot Comparison
	#----------------------------------------------------------------------
	#-- plot voronoi mascons
	ax.plot(tdec['vor'],mass['vor'],color='limegreen',zorder=3)
	ax.fill_between(tdec['vor'],mass['vor']-err_line,y2=mass['vor']+err_line,
					alpha=0.2,color='limegreen',zorder=3)
	#-- plot JPL mascons
	ax.plot(tdec['jpl'],mass['jpl'],color='blue',zorder=1)
	ax.fill_between(tdec['jpl'],mass['jpl']-errs['jpl'],y2=mass['jpl']+errs['jpl'],
					alpha=0.2,color='blue',zorder=1)
	#-- plot GSFC mascons
	ax.plot(tdec['gsfc'],mass['gsfc'],color='red',zorder=2)
	ax.fill_between(tdec['gsfc'],mass['gsfc']-errs['gsfc'],y2=mass['gsfc']+errs['gsfc'],
					alpha=0.2,color='red',zorder=2)
	ax.axvspan(2017.5, 2018.37, color='palegoldenrod',zorder=4)

	if axi == 1:
		ax.set_xlabel('Time [Decimal Year]')
	if axj == 0:
		ax.set_ylabel('Mass [Gt]')
	#-- add title and labels
	ax.set_title(reg.replace('_',' ').upper())

#-- add legend to bottom
plt.plot([],[],color='limegreen',label='Voronoi Mascons')
plt.plot([],[],color='blue',label='JPL Mascons')
plt.plot([],[],color='red',label='GSFC Mascons')
fig.subplots_adjust(bottom=0.14,top=0.95)
fig.legend(fancybox=True, framealpha=1.0,
			loc="lower center", bbox_to_anchor=(0.5, 0.0),ncol=3)
#-- save file
outfile = os.path.join(os.getcwd(),'..','data','paper','mascon_timeseries_comparison.pdf')
plt.savefig(outfile,format='PDF')
plt.close(fig)
