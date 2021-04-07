#!/usr/bin/env python
u"""
calc_leakage_ratio.py

Calculate first-order estimate of leakage by comparing
kernel to the harmonic representation of the region of interest
"""
import os
import sys
import numpy as np
from gravity_toolkit.ncdf_read import ncdf_read
from gravity_toolkit.ncdf_write import ncdf_write
from gravity_toolkit.ncdf_read_stokes import ncdf_read_stokes
from gravity_toolkit.read_love_numbers import read_love_numbers
from gravity_toolkit.harmonic_summation import harmonic_summation
from gravity_toolkit.plm_mohlenkamp import plm_mohlenkamp

#------------------------------------------------------------------------------
#-- Calculate Leakage ratio
#------------------------------------------------------------------------------
def calc_leakage(parameters):
	#-- read grid parameters
	DDEG_RASTER = float(parameters['DDEG_RASTER'])
	#-- read harmonic parameters
	LMAX = int(parameters['LMAX'])
	MMAX = int(parameters['MMAX'])
	LMIN = int(parameters['LMIN'])
	#-- smoothing radius
	RAD = int(parameters['RAD'])
	#-- Set up output labels
	OCN = '_OCN' if parameters['MASCON_OCEAN'] in ['Y','y'] else ''
	#-- read love numbers from file in parameter file
	love_file = os.path.expanduser(parameters['LOVE_FILE'])
	hl,kl,ll = read_love_numbers(love_file,REFERENCE='CF')

	#-- load mascon configuration of interest
	mascon_nums = np.array(parameters['MSCN_NUMS'].split(','),dtype=int)
	mascon_name = parameters['MSCN_NAME']
	lbl = '{0}_{1}'.format(mascon_name,parameters['MSCN_NUMS'].replace(',','+'))
	
	#-- check if difference file exists
	diff_file = os.path.join(os.path.expanduser(parameters['DIRECTORY']),\
		'MASCON_{0}_YLMS_{1:.2f}DEG_SKERNEL{2}_L{3:d}_r{4:d}km_DIFFERENCE.nc'\
		.format(lbl,DDEG_RASTER,OCN,LMAX,RAD))

	# #-- If file exists, read it. If not, calculate it and save it to file for next time and/or reference.
	# if os.path.exists(diff_file):
	# 	print('Reading difference file...')
	# 	diff_data = ncdf_read(diff_file,DATE=False,ATTRIBUTES=False,TITLE=False)
	# 	diff = diff_data['data']
	# 	lons = diff_data['lon']
	# 	lats = diff_data['lat']
	# else:
		# print('Calculating difference file...')
	#-- read kernel
	kern_file = os.path.join(os.path.expanduser(parameters['DIRECTORY']),\
		'MASCON_{0}_YLMS_{1:.2f}DEG_SKERNEL{2}_L{3:d}_r{4:d}km.nc'\
		.format(lbl,DDEG_RASTER,OCN,LMAX,RAD))
	kern_data = ncdf_read(kern_file,DATE=False,ATTRIBUTES=False,TITLE=False)
	
	#-- extra kernel fields and make grid
	lons = kern_data['lon']
	lats = kern_data['lat']
	kern = kern_data['data']
	glon,glat = np.meshgrid(lons,lats)

	#-- create Legendre polynomials
	th = (90.0 - lats)*np.pi/180.0
	plm = plm_mohlenkamp(LMAX, np.cos(th))

	#-- convert harmonics to mass unit before spatial conversion
	#-- Average Density of the Earth [g/cm^3]
	rho_e = 5.517
	#-- Average Radius of the Earth [cm]
	rad_e = 6.371e8
	l = np.arange(0,LMAX+1)
	dfactor = rho_e*rad_e*(2.0*l+1.0)/(1.0 +kl[l])/3.0
	for m in range(0,MMAX+1):#-- MMAX+1 to include MMAX
		for l in range(m,LMAX+1):#-- LMAX+1 to include LMAX
			plm[l,m,:] = plm[l,m,:]*dfactor[l]

	#-- read mascon harmonics
	Ylms = {}
	harmonic_dir = os.path.join(os.path.dirname(os.path.expanduser(parameters['DIRECTORY'])),'harmonics')
	for i in mascon_nums:
		har_file = os.path.join(harmonic_dir,'mascon_{0:d}_Ylms_L{1:d}_{2:.2f}deg.nc'.format(i,LMAX,DDEG_RASTER))
		Ylms[i] = ncdf_read_stokes(har_file,DATE=False,ATTRIBUTES=False)

	#-- convert the harmonics to spatial domain and sum up mascons
	harm_sum = np.zeros(glon.shape)
	for i in mascon_nums:
		harm_sum += harmonic_summation(Ylms[i]['clm'],Ylms[i]['slm'],lons,lats,LMAX=LMAX,MMAX=MMAX,PLM=plm).transpose()

	#-- take difference between kernel and harmonic representation
	diff = kern - harm_sum

	#-- save difference to file
	ncdf_write(diff,lons,lats,0,FILENAME=diff_file,DATE=False,UNITS='unitless',LONGNAME='Sensitivity_Kernel_DIFFERENCE')

	#-- calculate Mean Percent Error (https://en.wikipedia.org/wiki/Mean_percentage_error)
	#- first calculate a map of percent errors
	percent_error = np.abs(diff/harm_sum) * 100
	#-- calculate area-weighted mean
	th_grid = (90.0 - glat)*np.pi/180.0
	area_scale = np.sin(th_grid)
	mpe = np.abs(np.sum(area_scale*percent_error)/np.sum(area_scale))
	print("Mean Percentage Error: {0:.1f}".format(mpe))

	#-- calculate mean ratio
	ratio = kern/harm_sum * 100
	mean_ratio = np.sum(area_scale*ratio)/np.sum(area_scale)
	print("Mean Ratio: {0:.1f}".format(mean_ratio))
	
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
			#-- feed parameters to function to calculate leakage
			calc_leakage(parameters)

#------------------------------------------------------------------------------
#-- run main program
#------------------------------------------------------------------------------
if __name__ == '__main__':
	main()