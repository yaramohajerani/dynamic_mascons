#!/usr/bin/env python
u"""
calc_gridded_mascon_harmonics.py
by Yara Mohajerani

Calculate harmonic equivalent of GSFC and JPL mascons

Last Update 07/2021
"""
#-- load required modules
import os
import sys
import getopt
import numpy as np
import geopandas as gpd
#-- also import gravity toolkit modules
from gravity_toolkit.plm_holmes import plm_holmes
from gravity_toolkit.ncdf_stokes import ncdf_stokes
from gravity_toolkit.read_love_numbers import read_love_numbers
from gravity_toolkit.gen_spherical_cap import gen_spherical_cap


#------------------------------------------------------------------------------
#-- create harmonics for given voronoi regions
#------------------------------------------------------------------------------
def calc_harmonics(parameters,source=''):
	#-- input mascon grid
	if source == 'GSFC':
		grid_dir = os.path.expanduser(parameters['GSFC_MSCNS'])
	elif source == 'JPL':
		grid_dir = os.path.expanduser(parameters['JPL_MSCNS'])
	else:
		sys.exit('Data source not recognized (JPL or GSFC).')
	df_mask =  gpd.read_file(grid_dir)

	#-- read harmonic parameters
	LMAX = int(parameters['LMAX'])
	MMAX = int(parameters['MMAX'])
	LMIN = int(parameters['LMIN'])
	#-- read love numbers from file in parameter file
	love_file = os.path.expanduser(parameters['LOVE_FILE'])
	hl,kl,ll = read_love_numbers(love_file,REFERENCE='CF')
	#-- get output directory and create it if it doesn't exist
	out_dir = os.path.join(os.path.dirname(grid_dir),'harmonics')
	print('Output directory: ', out_dir)
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)

	# get all coordinates and calculate legendre polynomials
	if source == 'JPL':
		lons = df_mask['Longitude'].values
		lats = df_mask['Latitude'].values
	else:
		lons = df_mask['lon_center'].values
		lats = df_mask['lat_center'].values

	#-- get legendre polynomials
	th = (90.0 - lats)*np.pi/180.0#-- Colatitude in radians
	plm,dplm = plm_holmes(LMAX,np.cos(th))

	# make list of all files
	index_fid = open(os.path.join(out_dir,'mascon_Ylms_index_L{0:02d}_{1}.txt'.format(LMAX,source)),'w')
	# loop through mascons and save harmonics
	for i,m in enumerate(df_mask['id']):
		#-- create harmonics for a spherical cap
		if source == 'JPL':
			Ylms = gen_spherical_cap(1, lons[i], lats[i], LMAX=LMAX, MMAX=MMAX, \
				RAD_KM=df_mask['Radius'][i], UNITS=1, PLM=plm[:,:,i], LOVE=(hl,kl,ll))
		else:
			Ylms = gen_spherical_cap(1, lons[i], lats[i], LMAX=LMAX, MMAX=MMAX, \
				RAD_CAP=np.sqrt(df_mask['area_deg'][i])/2, UNITS=1, PLM=plm[:,:,i], LOVE=(hl,kl,ll))
		
		#-- save harmonics to file
		outfile = os.path.join(out_dir,'mascon_{0:d}_Ylms_L{1:02d}_{2}.nc'.format(int(m),LMAX,source))
		ncdf_stokes(np.array(Ylms.clm),np.array(Ylms.slm),np.arange(LMAX+1),np.arange(LMAX+1),0,0,FILENAME=outfile,DATE=False)
		#-- add file to index list
		index_fid.write('{0}\n'.format(outfile))
	index_fid.close()

	
#------------------------------------------------------------------------------
#-- main function
#------------------------------------------------------------------------------
def main():
	#-- read parameter files and commandline arguments
	optlist, input_files = getopt.getopt(sys.argv[1:],'S:',['source='])

	#-- command line parameters
	source = ''
	for opt, arg in optlist:
		if opt in ('-S','--source'):
			source = arg.upper()
		if opt in ('-G','--GSFC'):
			gsfc = True

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
			calc_harmonics(parameters,source=source)


#------------------------------------------------------------------------------
#-- run main program
#------------------------------------------------------------------------------
if __name__ == '__main__':
	main()
