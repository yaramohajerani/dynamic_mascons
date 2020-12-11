#!/usr/bin/env python
"""u
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
import getopt
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from scipy.spatial import SphericalVoronoi,geometric_slerp
#-- import pygplates from local installation (change path accordingly or add permanently to paths)
sys.path.append(os.path.expanduser('~/pygplates_rev28_python37_MacOS64'))
import pygplates

#-- configurations for unit sphere on which voronoi regions are constucted
r = 1 
origin = [0,0,0]

#------------------------------------------------------------------------------
#-- create plot
#------------------------------------------------------------------------------
def plot_html(new_centroids,new_sv,ind):
	#-- plot generators
	gens = go.Scatter3d(x=new_centroids[:, 0],y=new_centroids[:, 1],z=new_centroids[:, 2],mode='markers',marker={'size': 3,'opacity': 0.8,'color': 'blue'},name='Generator Points')
	fixed = [None]*len(ind)
	for i in range(len(ind)):
		fixed[i] = go.Scatter3d(x=[new_centroids[ind[i], 0]],y=[new_centroids[ind[i], 1]],z=[new_centroids[ind[i], 2]],mode='markers',\
			marker={'size': 5,'opacity': 0.8,'color': 'red'},name='Fixed Point {0:d}'.format(i))

	data = [gens] + fixed

	for region in new_sv.regions:
		n = len(region)
		t = np.linspace(0,1,50)
		for i in range(n):
			start = new_sv.vertices[region][i]
			end = new_sv.vertices[region][(i + 1) % n]
			result = np.array(geometric_slerp(start, end, t))
			edge = go.Scatter3d(x=result[..., 0],y=result[..., 1],z=result[..., 2],mode='lines',line={'width': 1,'color': 'black'},name='region edge',showlegend=False)
			data.append(edge)

	#-- configure layout
	layout = go.Layout(margin={'l': 0, 'r': 0, 'b': 0, 't': 0},title={
			'text': "Final Generator Setup",
			'y':0.9,
			'x':0.5,
			'xanchor': 'center',
			'yanchor': 'top'})

	fig = go.Figure(data=data, layout=layout)
	
	#-- save to file
	fig.write_html(os.path.join(os.getcwd(),'imgs','spherical_voronoi_regions.html'))


#------------------------------------------------------------------------------
#-- calculate voronoi regions based on given fixed points
#------------------------------------------------------------------------------
def calc_regions(lons=[],lats=[],eps=0,n_iter=0,r_const=0):
	#---------------------------------------------------------------
	# Set up initial generator grid
	#---------------------------------------------------------------
	#-- colatitude and longtiude lists in radians
	phis = np.radians(90-lats)
	thetas = np.radians(lons)
	#-- fill the rest of the coordinates (initial generators)
	phi_list = np.concatenate([phis,np.arange(eps,np.min(phis),eps),np.arange(np.max(phis)+eps,np.pi,eps)])
	theta_list = np.concatenate([thetas,np.arange(eps,np.min(thetas),eps),np.arange(np.max(thetas)+eps,2*np.pi,eps)])
	a1,a2 = np.meshgrid(phi_list,theta_list)
	points = np.array([[r*np.sin(phi)*np.cos(theta),r*np.sin(phi)*np.sin(theta),r*np.cos(phi)] for phi,theta in zip(a1.flatten(),a2.flatten())])

	#-- keep track of the index of the fixed points
	ind = np.zeros(len(lons),dtype=int)
	for i in range(1,len(lons)):
		ind[i] = i*(len(phi_list)+1)

	#-- test indices
	tmp_coords = list(zip(a1.flatten(),a2.flatten()))
	print('CHECK:')
	print("Co-Latitude of fixed points (radians): ",phi_list[:len(lats)])
	print("  Longitude of fixed points (radians): ",theta_list[:len(lons)])
	for i in range(len(lons)):
		print("Recovered Point {0:d} (ind {1:02d}): {2}".format(i, ind[i] ,tmp_coords[ind[i]]))

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
	# use the great circle distance between a given centroid and a
	#  reference point to shift the centroids
	#---------------------------------------------------------------
	#-- get reference point coordinates for calcuting distances to reference point
	lat0 = np.mean(lats)
	lon0 = np.mean(lons)
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
			#-- shift with res[ect to reference point
			#-- the shift ratio is inversely proportional to the distance
			ratio = r_const*(np.pi-dist)
			clon_shift = clon-ratio*(clon-lon0)
			clat_shift = clat-ratio*(clat-lat0)
			#-- convert shifted centroid from lat/lon to xyz on unit sphere
			tmp_pt = pygplates.PointOnSphere([clat_shift,clon_shift])
			new_centroids[i] = np.array(tmp_pt.to_xyz())
		#-- set the centroid of the fixed points back to their original values
		for i in ind:
			new_centroids[i] = points[i]
	#-- perform final tesselation based on new centroids
	new_sv = SphericalVoronoi(new_centroids, r, origin)
	new_sv.sort_vertices_of_regions()

	#---------------------------------------------------------------
	#-- plot the final configuration and save to html
	#---------------------------------------------------------------
	plot_html(new_centroids,new_sv,ind)

	#---------------------------------------------------------------
	#-- Finally, save spherical voronoi object to file
	#---------------------------------------------------------------
	with open(os.path.join(os.getcwd(),'data','sv_obj_test'), 'wb') as out_file:
		pickle.dump(new_sv, out_file)

#------------------------------------------------------------------------------
#-- main function
#------------------------------------------------------------------------------
def main():
	#-- read user-defiend fixed points
	optlist, input_files=getopt.getopt(sys.argv[1:],'E:I:R',['epsilon=','iterations=','ratio='])
	eps = 0.5
	n_iter = 50
	r_const = 0.025
	for opt, arg in optlist:
		if opt in ('-E','--epsilon'):
			eps = float(arg)
		elif opt in ('I','--iterations'):
			n_iter = int(arg)
		elif opt in ('R','--ratio'):
			r_const = float(arg)
	if not input_files:
		sys.exit("No fixed points specified.")
	else:
		for infile in input_files:
			df = pd.read_csv(infile)
			calc_regions(lons=np.array(df['LONS']),lats=np.array(df['LATS']),\
				eps=eps,n_iter=n_iter,r_const=r_const)

#------------------------------------------------------------------------------
#-- run main program
#------------------------------------------------------------------------------
if __name__ == '__main__':
	main()