#!/usr/bin/env python
u"""
get_scentroids.py
by Yara Mohajerani (11/2020)

Module for finding the weighted centroid of a given Voronoi region
"""
import numpy as np
from calc_weights import calc_weights

def get_weighted_centroid(pts,density=''):
	#-- Now get all points on the discretized boundary
	for i in range(len(pts)-1):
		#-- note as long as the edge is discretized into
		#-- the same number of elements in both directions,
		#-- the points will be on the line
		#-- But the number of points needs to be proportional
		#-- to the length of the line segment so each line
		#-- contributes equally to the average          
		length = np.sqrt(np.sum((pts[i+1]-pts[i])**2))
		n_pts = np.int(length/200)
		edge = {}
		for j in range(3):
			edge[j] = np.linspace(pts[i,j],pts[i+1,j],n_pts)

		#-- add to the total list of points
		try:
			pts = np.concatenate((pts,np.array(list(zip(edge[0],edge[1],edge[2])))), axis=0)
		except:    
			print("Unable to concatenate. 'pts' shape: {0}, array shape: {1}".format(\
				pts.shape,np.array(zip(edge[0],edge[1],edge[2])).shape))

	#-- calculate weights
	weights = calc_weights(pts,density)
	#-- calculate the centroid and return
	C = np.average(pts,axis=0,weights=weights)

	return C