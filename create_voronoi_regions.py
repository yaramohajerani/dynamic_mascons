#!/usr/bin/env python
u"""
create_voronoi_regions.py
by Yara Mohajerani

Get coordinates of fixed points from user and create
self-adjusting non-uniform voronoi regions around them.

Save voronoi regions as scipy SphericalVoronoi objects
and save an interactive html figure of the regions

Last Update: 12/2020
"""
import os
import sys
import copy
import pickle
import numpy as np
import pandas as pd
import geopandas as gpd
import plotly.graph_objs as go
from scipy.spatial import SphericalVoronoi,geometric_slerp
import astropy.coordinates as ac
from scipy.spatial.transform import Rotation as R
#-- import pygplates (https://www.gplates.org/docs/pygplates/pygplates_getting_started.html#installation)
import pygplates

#-- configurations for unit sphere on which voronoi regions are constucted
r = 1 
origin = [0,0,0]


#------------------------------------------------------------------------------
#-- create plot
#------------------------------------------------------------------------------
def plot_html(new_centroids,new_sv,ind,outfile=''):
	#-- initialize list to be plotted	
	data = []
	
	#-- ploy polygons
	for region in new_sv.regions:
		n = len(region)
		t = np.linspace(0,1,50)
		poly_xy = np.zeros((n,3))
		for i in range(n):
			poly_xy[i] = new_sv.vertices[region][i]
		polys = go.Mesh3d(x=poly_xy[:,0], y=poly_xy[:,1], z=poly_xy[:,2], opacity=1.0)	
		data.append(polys)

	#-- plot generators
	gens = go.Scatter3d(x=new_centroids[:, 0],y=new_centroids[:, 1],z=new_centroids[:, 2],mode='markers',marker={'size': 3,'opacity': 1.0,'color': 'black'},name='Generator Points')
	fixed = [None]*len(ind)
	for i in range(len(ind)):
		fixed[i] = go.Scatter3d(x=[new_centroids[ind[i], 0]],y=[new_centroids[ind[i], 1]],z=[new_centroids[ind[i], 2]],mode='markers',\
			marker={'size': 5,'opacity': 1.0,'color': 'red'},name='Fixed Point {0:d}'.format(i))

	data += [gens] + fixed

	#-- also plot world map
	world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
	for i in range(len(world)):
		if world['geometry'][i].geom_type == 'MultiPolygon':
			for j in range(len(world['geometry'][i])):
				lons,lats = world['geometry'][i][j].exterior.coords.xy
				#-- make line on sphere
				pl = pygplates.PolylineOnSphere(list(zip(lats,lons)))
				xyz = pl.to_xyz_array()
				cline = go.Scatter3d(x=xyz[:,0],y=xyz[:,1],z=xyz[:,2],mode='lines',line={'width': 6,'color': 'darkgreen'},name='Land',showlegend=False)
				data.append(cline)
		else:
			lons,lats = world['geometry'][i].exterior.coords.xy
			#-- make line on sphere
			pl = pygplates.PolylineOnSphere(list(zip(lats,lons)))
			xyz = pl.to_xyz_array()
			cline = go.Scatter3d(x=xyz[:,0],y=xyz[:,1],z=xyz[:,2],mode='lines',line={'width': 6,'color': 'darkgreen'},name='Land',showlegend=False)
			data.append(cline)

	#-- configure layout
	layout = go.Layout(margin={'l': 0, 'r': 0, 'b': 0, 't': 0})

	fig = go.Figure(data=data, layout=layout)
	
	#-- remove axes
	fig.update_scenes(xaxis_visible=False, yaxis_visible=False, zaxis_visible=False)

	#-- save to file
	fig.write_html(outfile)


#------------------------------------------------------------------------------
#-- calculate voronoi regions based on given fixed points
#------------------------------------------------------------------------------
def calc_regions(parameters):
	#---------------------------------------------------------------
	# Set up initial generator grid
	#---------------------------------------------------------------
	n_iter = int(parameters['N_ITERATIONS'])
	r_const = float(parameters['R_CONSTANT'])
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
	if rotate:
		#-- fill the rest of the coordinates (initial generators)
		theta_list = np.concatenate( [thetas, \
			np.arange(-np.pi, np.min(thetas),eps),\
			np.arange(np.max(thetas)+eps, np.pi, eps)])
		# theta_list = np.concatenate( [thetas, np.arange(np.max(thetas)+eps, 2*np.pi, eps)])
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

	#-- test indices
	sph_coords = list(zip(a1,a2))
	print('CHECK:')
	print("Co-Latitude of fixed points (radians): ", a1[:len(lats)])
	print("  Longitude of fixed points (radians): ", a2[:len(lons)])
	for i in range(len(lons)):
		print("Recovered Point {0:d} (ind {1:02d}): {2}".format(i, ind[i], sph_coords[ind[i]]))

	#---------------------------------------------------------------
	#-- create initial spherical voronoi tessellations
	#---------------------------------------------------------------
	sv = SphericalVoronoi(points, r, origin)
	sv.sort_vertices_of_regions()

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
	#-- set the centroid of the fixed points back to their original values
	for i in ind:
		centroids[i] = points[i]

	#---------------------------------------------------------------
	#-- plot the initial configuration and save to html
	#---------------------------------------------------------------
	if rotate:
		# first shift back to true lat/lon before plotting
		ry_recover = R.from_euler('y', 90-lat0_orig, degrees=True)
		rz_recover = R.from_euler('z', lon0_orig, degrees=True)
		true_centroids = np.zeros((len(sv.regions),3))
		for i in range(len(centroids)):
			# Make sure to perform rotation in reverse order of first rotation
			true_centroids[i] = rz_recover.apply(ry_recover.apply(centroids[i]))
			if i in ind:
				rec_latlon = ac.cartesian_to_spherical(true_centroids[i,0], true_centroids[i,1], true_centroids[i,2])
				rec_lat = np.degrees(rec_latlon[1].value)
				rec_lon = np.degrees(rec_latlon[2].value)
				print('Recovered fixed point: ', rec_lat, rec_lon)

		# make regions in true lat/lon
		sv_latlon = SphericalVoronoi(true_centroids, r, origin)
		sv_latlon.sort_vertices_of_regions()
		plot_html(true_centroids, sv_latlon, ind, \
			outfile=os.path.join(ddir,'spherical_voronoi_regions_0.html'))
	else:
		plot_html(centroids, sv, ind, \
			outfile=os.path.join(ddir,'spherical_voronoi_regions_0.html'))

	#---------------------------------------------------------------
	# use the great circle distance between a given centroid and a
	#  reference point to shift the centroids
	#---------------------------------------------------------------
	#-- convert to point on unit sphere
	ref_pt = pygplates.PointOnSphere([lat0,lon0])

	#-- Now iteratively use centroids as new generators to see updated configuration
	new_centroids = copy.copy(centroids)
	for _ in range(n_iter):
		new_sv = SphericalVoronoi(new_centroids, r, origin)
		new_sv.sort_vertices_of_regions()
		#-- calculate the centroids
		new_centroids = np.zeros((len(new_sv.regions),3))
		for i,region in enumerate(new_sv.regions):
			reg_vert = new_sv.vertices[region]
			#-- create polygon on surface of sphere
			poly = pygplates.PolygonOnSphere(reg_vert)
			#-- get centroid
			# cc = poly.get_interior_centroid()
			cc = poly.get_boundary_centroid()
			#-- get great-circle distance to reference point
			dist = cc.distance(cc,ref_pt)
			#-- convert to lat/lon
			clat,clon = np.array(cc.to_lat_lon())
			#-- shift with respect to reference point
			#-- the shift ratio is inversely proportional to the distance
			ratio = r_const*(np.pi-dist)
			if not rotate:
				clon_shift = clon-ratio*(clon-lon0)
			else:
				clon_shift = copy.deepcopy(clon)
			clat_shift = clat-ratio*(clat-lat0)
			#-- convert shifted centroid from lat/lon to xyz on unit sphere
			tmp_pt = pygplates.PointOnSphere([clat_shift, clon_shift])
			new_centroids[i] = np.array(tmp_pt.to_xyz())
		#-- set the centroid of the fixed points back to their original values
		for i in ind:
			new_centroids[i] = points[i]
	
	#-- shift back to true lat/lon after final iteration
	if rotate:
		for i in range(len(new_centroids)):
			# Make sure to perform rotation in reverse order of first rotation
			new_centroids[i] = rz_recover.apply(ry_recover.apply(new_centroids[i]))
			if i in ind:
				rec_latlon = ac.cartesian_to_spherical(new_centroids[i,0], new_centroids[i,1], new_centroids[i,2])
				rec_lat = np.degrees(rec_latlon[1].value)
				rec_lon = np.degrees(rec_latlon[2].value)
				print('Recovered fixed point (Final): ', rec_lat, rec_lon)

	#-- perform final tesselation based on new centroids
	new_sv = SphericalVoronoi(new_centroids, r, origin)
	new_sv.sort_vertices_of_regions()

	#---------------------------------------------------------------
	#-- plot the final configuration and save to html
	#---------------------------------------------------------------
	plot_html(new_centroids,new_sv,ind,outfile=os.path.join(ddir,'spherical_voronoi_regions.html'))

	#---------------------------------------------------------------
	#-- Finally, save spherical voronoi object to file
	#---------------------------------------------------------------
	with open(os.path.join(ddir,'sv_obj'), 'wb') as out_file:
		pickle.dump(new_sv, out_file)

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