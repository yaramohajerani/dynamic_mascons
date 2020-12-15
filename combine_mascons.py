#!/usr/bin/env python
u"""
combine_mascons.py
by Yara Mohajerani

Combine the mascon timeseries for fixed points and plot

Last Update 12/2020
"""
#-- load required modules
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
#-- create sensitivity kernels for given voronoi harmonics
#------------------------------------------------------------------------------
def combine_mascons(parameters):
	#-- read fixed-point coordinates
	coord_file = os.path.expanduser(parameters['COORD_FILE'])
	df = pd.read_csv(coord_file)
	lons=np.array(df['LONS'])
	lats=np.array(df['LATS'])
	DDEG_RASTER = float(parameters['DDEG_RASTER'])
	#-- read harmonic parameters
	LMAX = int(parameters['LMAX'])
	#-- get output directory
	ddir = os.path.expanduser(parameters['MSCN_DIRECTORY'])
	#-- get the spacing for the original generator points (to get indices of fixed points)
	eps = float(parameters['EPSILON'])
	#-- smoothing radius
	RAD = int(parameters['RAD'])
	#-- Set up output labels
	DS = '_FL' if (parameters['DESTRIPE'] in ['Y','y']) else ''
	OCN = '_OCN' if parameters['MASCON_OCEAN'] in ['Y','y'] else ''

	#----------------------------------------------------------------------
	#-- recalculate initial generator grid to get indices of fixed points
	#----------------------------------------------------------------------
	#-- colatitude and longtiude lists in radians
	phis = np.radians(90-lats)
	thetas = np.radians(lons)
	#-- fill the rest of the coordinates (initial generators)
	phi_list = np.concatenate([phis,np.arange(eps,np.min(phis),eps),np.arange(np.max(phis)+eps,np.pi,eps)])
	theta_list = np.concatenate([thetas,np.arange(eps,np.min(thetas),eps),np.arange(np.max(thetas)+eps,2*np.pi,eps)])

	#-- keep track of the index of the fixed points
	ind = np.zeros(len(lons),dtype=int)
	for i in range(1,len(lons)):
		ind[i] = i*(len(phi_list)+1)

	#----------------------------------------------------------------------
	#-- Read and sum up mascon timeseries corresponding to fixed points
	#----------------------------------------------------------------------
	mascon = {}
	for i in ind:
		#- read the netcdf files
		mscn_file = os.path.join(ddir,'MASCON_{0:d}_YLMS_{1:.2f}DEG_AW13_ICE6G_GA{2}_L{3:02d}_r{4:d}km{5}.txt'.format(i,DDEG_RASTER,OCN,LMAX,RAD,DS))
		mascon[i] = np.loadtxt(mscn_file)
	#-- sum up the kernels
	mons = np.array(mascon[i][:,0],dtype=int)
	tdec = mascon[i][:,1]
	mass = np.zeros(len(tdec))
	err = np.zeros(len(tdec))
	area = 0
	out_lbl = ''
	for i in ind:
		mass += mascon[i][:,2]
		err += mascon[i][:,3]**2
		area += mascon[i][0,4]
		out_lbl += '{0:d}+'.format(i)
	#-- RMS of errors
	err = np.sqrt(err)
	#-- remove last '+' from output label
	out_lbl = out_lbl[:-1]

	#----------------------------------------------------------------------
	#-- write sum to file
	#----------------------------------------------------------------------
	outfile = os.path.join(ddir,'MASCON_{0}_YLMS_{1:.2f}DEG_AW13_ICE6G_GA{2}_L{3:02d}_r{4:d}km{5}.txt'.format(out_lbl,DDEG_RASTER,OCN,LMAX,RAD,DS))
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
	plt.title('Total Mass Balance Time-Series of Mascons Corresponding to Fixed Points',fontsize=14)

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