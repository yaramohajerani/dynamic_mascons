#!/usr/bin/env python
u"""
calc_voronoi_harmonics.py
by Yara Mohajerani

Read the voronoi regions produced by "create_voronoi_regions.py" 
and create harmonic representations of mascons

Last Update 12/2020
"""
#-- load required modules
import os
import sys
import pickle
import numpy as np
import pandas as pd
import importlib.util
from scipy.spatial import SphericalVoronoi,geometric_slerp
#-- import pygplates (https://www.gplates.org/docs/pygplates/pygplates_getting_started.html#installation)
import pygplates
#-- also import gravity toolkit modules
from gravity_toolkit.gen_stokes import gen_stokes
from gravity_toolkit.ncdf_stokes import ncdf_stokes
from gravity_toolkit.ncdf_write import ncdf_write
from gravity_toolkit.ncdf_read_stokes import ncdf_read_stokes
from gravity_toolkit.read_love_numbers import read_love_numbers
from gravity_toolkit.harmonic_summation import harmonic_summation
from gravity_toolkit.plm_mohlenkamp import plm_mohlenkamp

#------------------------------------------------------------------------------
#-- create harmonics for given voronoi regions
#------------------------------------------------------------------------------
def calc_harmonics(parameters):
	#-- input file for voronoi regions
	input_file = os.path.expanduser(parameters['VORONOI_FILE'])
	with open(input_file, 'rb') as in_file:
		sv = pickle.load(in_file)
	#-- read grid parameters
	DDEG_RASTER = float(parameters['DDEG_RASTER'])
	#-- read harmonic parameters
	LMAX = int(parameters['LMAX'])
	MMAX = int(parameters['MMAX'])
	LMIN = int(parameters['LMIN'])
	#-- read love numbers from file in parameter file
	love_file = os.path.expanduser(parameters['LOVE_FILE'])
	hl,kl,ll = read_love_numbers(love_file,REFERENCE='CF')
	#-- get output directory
	out_dir = os.path.expanduser(parameters['OUT_DIRECTORY'])

	#-----------------------------------------------------------------------------
	#-- Now we make a grid so we can rasterize the polygons onto a grid
	#-----------------------------------------------------------------------------
	lons = np.arange(-180,180+DDEG_RASTER,DDEG_RASTER)
	lats = np.arange(-90,90+DDEG_RASTER,DDEG_RASTER)
	glat,glon = np.meshgrid(lats,lons)
	#-- flatten the grid for easier processing
	glon = glon.flatten()
	glat = glat.flatten()

	#-- create Legendre polynomials
	th = (90.0 - lats)*np.pi/180.0
	plm = plm_mohlenkamp(LMAX, np.cos(th))

	#-----------------------------------------------------------------------------
	#-- Convert voronoi regions into harmonics and savbe to file
	#-----------------------------------------------------------------------------
	if parameters['MAKE_HARMONICS'] in ['Y','y']:
		index_fid = open(os.path.join(out_dir,'mascon_Ylms_index_L{0:02d}_{1:.2f}deg.txt'.format(LMAX,DDEG_RASTER)),'w')
		Ylms = {}
		for i,region in enumerate(sv.regions):
			#-- intitialize mascon data as all zeros (no mass)
			mass = np.zeros(len(glon))
			#-- get vertices of region
			reg_vert = sv.vertices[region]
			#-- create polygon from vertices on surface of sphere
			poly = pygplates.PolygonOnSphere(reg_vert)
			#-- loop through grid points to identify which lie inside polygon
			for j in range(len(mass)):
				if poly.is_point_in_polygon((glat[j],glon[j])):
					mass[j] = 1
			#-- now convert the mass field to its harmonic representation
			Ylms[i] = gen_stokes(mass.reshape(len(lons),len(lats)),lons,lats,LMIN=LMIN,LMAX=LMAX,MMAX=MMAX,UNITS=1,PLM=plm,LOVE=(hl,kl,ll))
			#-- save harmonics to file
			outfile = os.path.join(out_dir,'mascon_{0:d}_Ylms_L{1:02d}_{2:.2f}deg.nc'.format(i,LMAX,DDEG_RASTER))
			ncdf_stokes(np.array(Ylms[i]['clm']),np.array(Ylms[i]['slm']),np.arange(LMAX+1),np.arange(LMAX+1),0,0,FILENAME=outfile,DATE=False)
			#-- add file to index list
			index_fid.write('{0}\n'.format(outfile))
		index_fid.close()
	else:
		Ylms = {}
		for i,region in enumerate(sv.regions):
			Ylms[i] = ncdf_read_stokes(os.path.join(out_dir,'mascon_{0:d}_Ylms_L{1:02d}_{2:.2f}deg.nc'.format(i,LMAX,DDEG_RASTER)),DATE=False,ATTRIBUTES=False)

	#-----------------------------------------------------------------------------
	#-- Convert harmonics back to spatial and save to file
	#-----------------------------------------------------------------------------
	if parameters['MAKE_SPATIAL'] in ['Y','y']:
		DDEG_OUT = float(parameters['DDEG_OUT'])
		lons = np.arange(0,360+DDEG_OUT,DDEG_OUT)
		lats = np.arange(-90,90+DDEG_OUT,DDEG_OUT)

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

		#-- convert each mascon back to spatial and save as netCDF
		for i in range(len(Ylms)):
			sp = harmonic_summation(Ylms[i]['clm'],Ylms[i]['slm'],lons,lats,LMAX=LMAX,MMAX=MMAX,PLM=plm)
			#-- save to file
			outfile = os.path.join(out_dir,'mascon_{0:d}_spatial_L{1:02d}_{2:.2f}rasterDeg_{3:.2f}outDeg.nc'.format(i,LMAX,DDEG_RASTER,DDEG_OUT))
			ncdf_write(np.transpose(sp),lons,lats,i,FILENAME=outfile,DATE=False,UNITS='cmWE',LONGNAME='Mascon {0:d}'.format(i))

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
			#-- feed parameters to function to create harmonics
			calc_harmonics(parameters)

#------------------------------------------------------------------------------
#-- run main program
#------------------------------------------------------------------------------
if __name__ == '__main__':
	main()