#!/usr/bin/env python
u"""
calc_weights.py
by Yara Mohajerani (11/2020)

Calcalate weights for given points to calculate vornoi centroid 
based on the given density label
"""
import sys
import numpy as np

def calc_weights(pts,density):
	if density in ['uniform','Uniform','uni']:
		weights = np.ones(len(pts))
	else:
		sys.exit('Specified density model does not exist yet.')
	
	return weights

