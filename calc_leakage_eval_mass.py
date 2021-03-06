#!/usr/bin/env python
u"""
calc_leakage_eval_mass.py

Calculate first-order estimate of leakage by comparing
kernel to the harmonic representation of the region of interest
"""
import os
import re
import sys
import numpy as np
from gravity_toolkit.ncdf_read_stokes import ncdf_read_stokes
from gravity_toolkit.ncdf_read_stokes import ncdf_read_stokes
from gravity_toolkit.read_love_numbers import read_love_numbers

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

	#-- read mascon harmonics
	harmonic_dir = os.path.expanduser(parameters['HARMONIC_DIRECTORY'])
	harm_index = os.path.join(harmonic_dir,'mascon_Ylms_index_L{0:d}_{1:.2f}deg.txt'.format(LMAX,DDEG_RASTER))
	fid = open(harm_index,'r')
	file_list = fid.readlines()
	file_list = [f.replace('\n','') for f in file_list if f.endswith('\n')]
	fid.close()
	
	#-- get regex pattern for extracting mascon numbers
	reg = re.compile(r'mascon_(\d+)_Ylms_L{0:d}_{1:.2f}deg.nc'.format(LMAX,DDEG_RASTER))
	Ylms = {}
	nums = []
	for f in file_list:
		fname = os.path.basename(f)
		n = int(reg.search(fname).group(1))
		#-- read the harmonics
		Ylms[n] = ncdf_read_stokes(f,DATE=False)
		nums.append(n)

	#-- read kernels 
	kern = {'clm':np.zeros((LMAX+1,MMAX+1)),'slm':np.zeros((LMAX+1,MMAX+1))}
	for m in mascon_nums:
		kern_file = os.path.join(os.path.expanduser(parameters['DIRECTORY']),\
			'MASCON_{0}_YLMS_{1:.2f}DEG_SKERNEL_CLM{2}_L{3:d}_r{4:d}km.nc'\
			.format(n,DDEG_RASTER,OCN,LMAX,RAD))
		kern_m = ncdf_read_stokes(kern_file,DATE=False)
		kern['clm'] += kern_m['clm'][:,:]
		kern['slm'] += kern_m['slm'][:,:]
	
	#-- convert harmonics to mass unit before spatial conversion
	#-- Average Density of the Earth [g/cm^3]
	rho_e = 5.517
	#-- Average Radius of the Earth [cm]
	rad_e = 6.371e8
	l = np.arange(0,LMAX+1)
	dfactor = rho_e*rad_e*(2.0*l+1.0)/(1.0 +kl[l])/3.0
	for m in range(MMAX+1):
		kern['clm'][:,m] *= dfactor
		kern['slm'][:,m] *= dfactor

	#-- outside 
	inside = 0
	outside = 0
	for n in nums:
		if n in mascon_nums:
			inside += np.sum(kern['clm'][:,:]*Ylms[n]['clm'][:,:] + kern['slm'][:,:]*Ylms[n]['slm'][:,:])
		else:
			outside += np.sum(kern['clm'][:,:]*Ylms[n]['clm'][:,:] + kern['slm'][:,:]*Ylms[n]['slm'][:,:])
	#-- get ratio
	leakage = np.abs(outside/inside) * 100

	#-- save to file
	outfile = os.path.join(os.path.expanduser(parameters['DIRECTORY']),\
			'MASCON_{0}_LEAKAGE_ERROR_{1:.2f}DEG_SKERNEL{2}_L{3:d}_r{4:d}km.csv'\
			.format(lbl,DDEG_RASTER,OCN,LMAX,RAD))
	fid = open(outfile,'w')
	fid.write('region,inside,outside,leakage\n')
	fid.write('{0},{1:.2E},{2:.2E},{3:.1f}'.format(lbl,inside,outside,leakage))
	fid.close()

	#-- also display results
	print('inside: {0:.2E}, outside: {1:.2E}'.format(inside,outside))
	print('{0} Leakage ratio (%): {1:.1f}'.format(lbl,leakage))
	
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