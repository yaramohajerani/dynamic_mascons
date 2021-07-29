#!/usr/bin/env python
u"""
module for plotting global mascon configuration as an interactive HTML
"""
import random
import numpy as np
import geopandas as gpd
import plotly.graph_objs as go
#-- import pygplates (https://www.gplates.org/docs/pygplates/pygplates_getting_started.html#installation)
import pygplates

#-- configurations for unit sphere on which voronoi regions are constucted
r = 1 
origin = [0,0,0]

#------------------------------------------------------------------------------
#-- create plot
#------------------------------------------------------------------------------
def plot_html(new_centroids, new_sv, ind, colors=[], outfile=''):
	#-- if  given color range doesn't cover all indices, repeat
	while len(colors) < len(new_sv.regions):
		colors += colors
	
	#-- initialize list to be plotted	
	data = [None]*len(new_sv.regions)
	
	#-- ploy polygons
	for count,region in enumerate(new_sv.regions):
		n = len(region)
		poly_xy = np.zeros((n,3))
		for i in range(n):
			poly_xy[i] = new_sv.vertices[region][i]
		polys = go.Mesh3d(x=poly_xy[:,0], y=poly_xy[:,1], z=poly_xy[:,2], opacity=1.0, color=colors[count])	
		data[count] = polys

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
				cline = go.Scatter3d(x=xyz[:,0],y=xyz[:,1],z=xyz[:,2],mode='lines',line={'width': 2,'color': 'white'},name='Land',showlegend=False)
				data.append(cline)
		else:
			lons,lats = world['geometry'][i].exterior.coords.xy
			#-- make line on sphere
			pl = pygplates.PolylineOnSphere(list(zip(lats,lons)))
			xyz = pl.to_xyz_array()
			cline = go.Scatter3d(x=xyz[:,0],y=xyz[:,1],z=xyz[:,2],mode='lines',line={'width': 2,'color': 'white'},name='Land',showlegend=False)
			data.append(cline)

	#-- configure layout
	layout = go.Layout(margin={'l': 0, 'r': 0, 'b': 0, 't': 0})

	fig = go.Figure(data=data, layout=layout)
	
	#-- remove axes
	fig.update_scenes(xaxis_visible=False, yaxis_visible=False, zaxis_visible=False)

	#-- save to file
	fig.write_html(outfile)
	