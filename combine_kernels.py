#!/usr/bin/env python
u"""
combine_kernels.py
by Yara Mohajerani

Combine the sensitivity kernels of the sum of the 'fixed points'
and produce netcdf and png outputs

Last Update 12/2020
"""
#-- load required modules
import os
import sys
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
#-- also import gravity toolkit modules
from gravity_toolkit.ncdf_write import ncdf_write
from gravity_toolkit.ncdf_read import ncdf_read

#------------------------------------------------------------------------------
#-- create sensitivity kernels for given voronoi harmonics
#------------------------------------------------------------------------------
def combine_kernels(parameters):
	DDEG_RASTER = float(parameters['DDEG_RASTER'])
	#-- read harmonic parameters
	LMAX = int(parameters['LMAX'])
	#-- get output directory
	ddir = os.path.expanduser(parameters['DIRECTORY'])
	#-- smoothing radius
	RAD = int(parameters['RAD'])
	#-- ocn redistribution label
	OCN = '_OCN' if parameters['MASCON_OCEAN'] in ['Y','y'] else ''

	#-- load mascon configuration of interest
	mascon_nums = np.array(parameters['MSCN_NUMS'].split(','),dtype=int)
	mascon_name = parameters['MSCN_NAME']
	out_lbl = '{0}_{1}'.format(mascon_name,parameters['MSCN_NUMS'].replace(',','+'))

	#----------------------------------------------------------------------
	#-- Read and sum up kernels corresponding to fixed points
	#----------------------------------------------------------------------
	kerns = {}
	for i in mascon_nums:
		#- read the netcdf files
		kern_file = os.path.join(ddir,'MASCON_{0:d}_YLMS_{1:.2f}DEG_SKERNEL{2}_L{3:02d}_r{4:d}km.nc'.format(i,DDEG_RASTER,OCN,LMAX,RAD))
		kerns[i] = ncdf_read(kern_file,DATE=False,TITLE=False)
	#-- sum up the kernels
	kern_sum = kerns[mascon_nums[0]]['data']
	for i in mascon_nums[1:]:
		kern_sum += kerns[i]['data']
	#-- get grid for saving combined sensitivity kernel
	glat = kerns[mascon_nums[0]]['lat']
	glon = kerns[mascon_nums[0]]['lon']

	#----------------------------------------------------------------------
	#-- write kernel sum to file
	#----------------------------------------------------------------------
	outfile = os.path.join(ddir,'MASCON_{0}_YLMS_{1:.2f}DEG_SKERNEL_OCN_L{2:02d}_r{3:d}km.nc'.format(out_lbl,DDEG_RASTER,LMAX,RAD))
	ncdf_write(kern_sum,glon,glat,0,FILENAME=outfile,DATE=False,UNITS='unitless',LONGNAME='Sensitivity_Kernel')

	#----------------------------------------------------------------------
	#-- plot summed kernel
	#----------------------------------------------------------------------
	#-- load in world map for plotting in background
	world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
	#-- plot summed kernel
	fig, ax = plt.subplots(1,figsize = (10,6),dpi=100)
	klim = np.max(np.abs(kern_sum))*0.95
	c = ax.contourf(glon,glat,kern_sum,cmap='bwr',levels=np.linspace(-klim,klim,16))
	#-- use an axis divider for the colorbar
	drx = make_axes_locatable(ax)
	cax = drx.append_axes("right", size="5%", pad=0.1)
	cbar = fig.colorbar(c,cax=cax)
	cbar.set_label('Sensitivity Kernel (min:{0:.1f}, max:{1:.1f})'.format(np.min(kern_sum),np.max(kern_sum)))
	world.plot(ax=ax,alpha=0.3,fc='none',ec='k',linewidth=1.2,rasterized=True)
	plt.tight_layout()
	plt.savefig(outfile.replace('.nc','.png'),format='PNG')
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
			combine_kernels(parameters)

#------------------------------------------------------------------------------
#-- run main program
#------------------------------------------------------------------------------
if __name__ == '__main__':
	main()