#!/usr/bin/env python
u"""
combine_mascons.py
by Yara Mohajerani

Combine the mascon timeseries for fixed points and plot

Last Update 04/2021
"""
#-- load required modules
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
#-- create sensitivity kernels for given voronoi harmonics
#------------------------------------------------------------------------------
def combine_mascons(parameters):
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

	#-- load mascon configuration of interest
	mascon_nums = np.array(parameters['MSCN_NUMS'].split(','),dtype=int)
	mascon_name = parameters['MSCN_NAME']
	out_lbl = '{0}_{1}'.format(mascon_name,parameters['MSCN_NUMS'].replace(',','+'))

	#----------------------------------------------------------------------
	#-- Read and sum up mascon timeseries corresponding to fixed points
	#----------------------------------------------------------------------
	mascon = {}
	for i in mascon_nums:
		#- read the netcdf files
		mscn_file = os.path.join(ddir,'MASCON_{0:d}_YLMS_{1:.2f}DEG_{2}{3}_L{4:02d}_r{5:d}km{6}.txt'.format(i,DDEG_RASTER,GIA,OCN,LMAX,RAD,DS))
		mascon[i] = np.loadtxt(mscn_file)
	#-- sum up the kernels
	mons = np.array(mascon[i][:,0],dtype=int)
	tdec = mascon[i][:,1]
	mass = np.zeros(len(tdec))
	err = np.zeros(len(tdec))
	area = 0
	for i in mascon_nums:
		mass += mascon[i][:,2]
		err += mascon[i][:,3]**2
		area += mascon[i][0,4]
	#-- RMS of errors
	err = np.sqrt(err)

	#----------------------------------------------------------------------
	#-- write sum to file
	#----------------------------------------------------------------------
	outfile = os.path.join(ddir,'MASCON_{0}_YLMS_{1:.2f}DEG_{2}{3}_L{4:02d}_r{5:d}km{6}.txt'.format(out_lbl,DDEG_RASTER,GIA,OCN,LMAX,RAD,DS))
	fid = open(outfile,'w')
	for i in range(len(tdec)):
		fid.write('{0:03d} {1:12.4f} {2:16.10f} {3:16.10f} {4:16.5f}\n'.format(mons[i],tdec[i],mass[i],err[i],area))
	fid.close()

	#----------------------------------------------------------------------
	#-- plot total mass
	#----------------------------------------------------------------------
	#-- get index of FO gap for plotting
	fo = np.squeeze(np.nonzero(mons==198))

	#-- make plot
	fig = plt.figure(1,figsize=(15,9))
	plt.plot(tdec[:fo],mass[:fo]-mass.mean(),color='blue')
	plt.plot(tdec[fo:],mass[fo:]-mass.mean(),color='blue')
	plt.fill_between(tdec[:fo],mass[:fo]-mass.mean()-err[:fo],y2=mass[:fo]-mass.mean()+err[:fo],color='deepskyblue',alpha=0.5)
	plt.fill_between(tdec[fo:],mass[fo:]-mass.mean()-err[fo:],y2=mass[fo:]-mass.mean()+err[fo:],color='deepskyblue',alpha=0.5)
	#-- Fill FO gap with vertical bar
	#-- index of FO gap
	#ind_gap = np.where((data[p][:,0]>186) & (data[p][:,0]<=198))
	plt.axvspan(2017.5, 2018.37, color='lightgray',zorder=1,alpha=0.8)
	plt.xlabel('Time [Decimal Years]',fontsize=14)
	plt.ylabel('Mass [Gt]',fontsize=14)
	plt.title('Total Mass Balance Time-Series of Mascons {0}'.format(out_lbl),fontsize=14)

	plt.savefig(outfile.replace('.txt','.png'),format='PNG')
	plt.close(fig=fig)

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
			#-- feed parameters to function to combine and plot kernels
			combine_mascons(parameters)

#------------------------------------------------------------------------------
#-- run main program
#------------------------------------------------------------------------------
if __name__ == '__main__':
	main()