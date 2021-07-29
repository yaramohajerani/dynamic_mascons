#!/usr/bin/env python
u"""
plot_voronoi_regions.py
by Yara Mohajerani

Plot the output of create_voronoi_regions

Last Update: 07/2021
"""
import os
import sys
import pickle
import random
import numpy as np
import pandas as pd
import astropy.coordinates as ac
from scipy.spatial.transform import Rotation as R
from plot_configuration_html import plot_html
#-- import pygplates (https://www.gplates.org/docs/pygplates/pygplates_getting_started.html#installation)
import pygplates

#-- configurations for unit sphere on which voronoi regions are constucted
r = 1 
origin = [0,0,0]


#------------------------------------------------------------------------------
#-- calculate voronoi regions based on given fixed points
#------------------------------------------------------------------------------
def calc_regions(parameters):
	#---------------------------------------------------------------
	# Read configuration to plot
	#---------------------------------------------------------------
	input_file = os.path.expanduser(parameters['VORONOI_FILE'])
	with open(input_file, 'rb') as in_file:
		sv = pickle.load(in_file)

	#---------------------------------------------------------------
	# Set up initial generator grid to get fixed point indices
	#---------------------------------------------------------------
	rotate = True if parameters['ROTATE'].upper() in ['TRUE','Y'] else False
	print('Rotate:', rotate)
	#-- read fixed-point coordinates
	coord_file = os.path.expanduser(parameters['COORD_FILE'])
	ddir = os.path.dirname(coord_file)
	df = pd.read_csv(coord_file)
	if rotate:
		lons_orig = np.array(df['LONS'])
		lats_orig = np.array(df['LATS'])
		#-- get reference point coordinates for calcuting distances to reference point
		lat0_orig = np.mean(lats_orig)
		lon0_orig = np.mean(lons_orig)

		# make rotation matrices to rotate fixed point to North Pole
		ry = R.from_euler('y',  -(90-lat0_orig), degrees=True)
		rz = R.from_euler('z', -lon0_orig, degrees=True)

		# make a Cartesian vector fo fixed points and rotate to new frame
		n_fixed = len(lons_orig)
		lats = [None]*n_fixed
		lons = [None]*n_fixed
		for i in range(n_fixed):
			xyz = ac.spherical_to_cartesian(1,np.radians(lats_orig[i]),np.radians(lons_orig[i]))
			v = [k.value for k in xyz]
			#-- rotate coordinates and get new fixed point
			rot_xyz = ry.apply(rz.apply(v))
			rot_latlon = ac.cartesian_to_spherical(rot_xyz[0], rot_xyz[1], rot_xyz[2])
			lats[i] = np.degrees(rot_latlon[1].value)
			lons[i] = np.degrees(rot_latlon[2].value)
		lats = np.array(lats)
		lons = np.array(lons)

		#-- refernce points after rotation
		lat0 = np.mean(lats)
		lon0 = np.mean(lons)
		print(f'Original lon {lon0_orig} lat {lat0_orig}. Transfomed lon {lon0:.2f} lat {lat0:.2f}.')
	else:
		# read fixed points
		lons = np.array(df['LONS'])
		lats = np.array(df['LATS'])
		# reference point
		lat0 = np.mean(lats)
		lon0 = np.mean(lons)

	#-- get grid interval
	eps = float(parameters['EPSILON'])

	#-- colatitude and longtiude lists in radians
	phis = np.radians(90-lats)
	thetas = np.radians(lons)
	# change to -180 to 180 range.
	if (thetas > np.pi).any():
		thetas -= 2*np.pi
	if rotate:
		#-- fill the rest of the coordinates (initial generators)
		theta_list = np.concatenate( [thetas, \
			np.arange(-np.pi, np.min(thetas),eps),\
			np.arange(np.max(thetas)+eps, np.pi, eps)])
		if len(lons)==1:
			print('Only 1 given fixed point')
			phi_list = np.arange(np.max(phis)+eps, np.pi, eps)
			a1,a2 = np.meshgrid(phi_list,theta_list)
			a1 = np.concatenate(( phis, a1.flatten() ))
			a2 = np.concatenate(( thetas, a2.flatten() ))
		else:
			print('Multiple fixed points.')
			phi_list = np.concatenate([phis, np.arange(np.max(phis)+eps, np.pi, eps) ] )
			a1,a2 = np.meshgrid(phi_list,theta_list)
			a1 = a1.flatten()
			a2 = a2.flatten()
	else:
		phi_list = np.concatenate([phis,np.arange(eps,np.min(phis),eps),np.arange(np.max(phis)+eps,np.pi,eps)])
		theta_list = np.concatenate([thetas,np.arange(eps,np.min(thetas),eps),np.arange(np.max(thetas)+eps,2*np.pi,eps)])
		a1,a2 = np.meshgrid(phi_list,theta_list)
		a1 = a1.flatten()
		a2 = a2.flatten()

	# points = np.array([[k.value for k in ac.spherical_to_cartesian(1,phi,theta)] for phi,theta in zip(a1,a2)])
	points = np.array([[r*np.sin(phi)*np.cos(theta),r*np.sin(phi)*np.sin(theta),r*np.cos(phi)] for phi,theta in zip(a1,a2)])
	print('Epsilon: {0:.2f}, Number of Points: {1:d}'.format(eps,len(points)))

	#-- keep track of the index of the fixed points
	ind = np.zeros(len(lons),dtype=int)
	for i in range(1,len(lons)):
		ind[i] = i*(len(phi_list)+1)

	#---------------------------------------------------------------
	#-- Calculate centroids
	#---------------------------------------------------------------
	centroids = np.zeros((len(sv.regions),3))
	for i,region in enumerate(sv.regions):
		reg_vert = sv.vertices[region]
		#-- create polygon on surface of sphere
		poly = pygplates.PolygonOnSphere(reg_vert)
		#-- get centroid
		# centroids[i] = np.array(poly.get_interior_centroid().to_xyz())
		centroids[i] = np.array(poly.get_boundary_centroid().to_xyz())
	# #-- set the centroid of the fixed points back to their original values
	# for i in ind:
	# 	centroids[i] = points[i]

	#---------------------------------------------------------------
	#-- plot the final configuration and save to html
	#---------------------------------------------------------------
	# make color map
	random.seed(1)
	rcol = lambda: random.randint(0,255)
	colors = ['#%02X%02X%02X' % (rcol(), rcol(), rcol()) for i in range(len(sv.regions))]
	plot_html(centroids, sv, ind, colors=colors, outfile=os.path.join(ddir,\
		f'{os.path.basename(input_file)}_spherical_voronoi_regions.html'))

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
			calc_regions(parameters)

#------------------------------------------------------------------------------
#-- run main program
#------------------------------------------------------------------------------
if __name__ == '__main__':
	main()