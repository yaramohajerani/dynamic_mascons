#!/usr/bin/env python
u"""
plot_harmonics.py
by Yara Mohajerani
01/2021

Simple script for plotting the 
harmonic representation of given regions.
"""
import os
import pickle
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.spatial import geometric_slerp
import pygplates

from gravity_toolkit.read_love_numbers import read_love_numbers
from gravity_toolkit.ncdf_read_stokes import ncdf_read_stokes
from gravity_toolkit.plm_mohlenkamp import plm_mohlenkamp
from gravity_toolkit.harmonic_summation import harmonic_summation

#-- directory setup
base_dir = os.path.join(os.getcwd(),'..','data')

#-- load in world map for plotting in background
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

#-- make a list of regions to be plotted
regs = ['karakoram','nyainqentangla','alaska']

#-- set up grid
DDEG = 0.5
lons = np.arange(-180,180+DDEG,DDEG)
lats = np.arange(-90,90+DDEG,DDEG)
glat,glon = np.meshgrid(lats,lons)
LMAX, MMAX = 60,60

#-- create Legendre polynomials
th = (90.0 - lats)*np.pi/180.0
plm = plm_mohlenkamp(LMAX, np.cos(th))

#-- read love numbers required for generating geoid harmonics
hl,kl,ll = read_love_numbers(os.path.expanduser('~/read-GRACE-harmonics/gravity_toolkit/data/love_numbers'),REFERENCE='CF')

#-- convert harmonics to mass unit before spatial conversion
#-- Average Density of the Earth [g/cm^3]
rho_e = 5.517
#-- Average Radius of the Earth [cm]
rad_e = 6.371e8
l = np.arange(0,LMAX+1)
dfactor = rho_e*rad_e*(2.0*l+1.0)/(1.0 +kl[l])/3.0
for m in range(0,MMAX+1):#-- MMAX+1 to include MMAX
    for l in range(m,LMAX+1):#-- LMAX+1 to include LMAX
        plm[l,m,:] = plm[l,m,:]*dfactor[l]

fig, ax = plt.subplots(len(regs),1,figsize = (8,10),sharey=True)#, dpi=100)
for i,(r,k) in enumerate(zip(regs,[70,85,40])):
	#-- First load the saved SphericalVoronoi object
	with open(os.path.join(base_dir,r,'sv_obj'), 'rb') as in_file:
		sv = pickle.load(in_file)
	#-- read harmonics
	m1 = ncdf_read_stokes(os.path.join(base_dir,r,'harmonics','mascon_0_Ylms_L60_0.25deg.nc'),DATE=False,ATTRIBUTES=False)
	m2 = ncdf_read_stokes(os.path.join(base_dir,r,'harmonics','mascon_{0:d}_Ylms_L60_0.25deg.nc'.format(k)),DATE=False,ATTRIBUTES=False)
	#-- calculate the spatial representation of the harmonics
	sp1 = harmonic_summation(m1['clm'],m1['slm'],lons,lats,LMAX=LMAX,MMAX=MMAX,PLM=plm)
	sp2 = harmonic_summation(m2['clm'],m2['slm'],lons,lats,LMAX=LMAX,MMAX=MMAX,PLM=plm)
	levels = np.linspace(-1, 1, 11)
	p = ax[i].contourf(glon,glat,sp1+sp2,cmap='RdBu',levels=levels,extend='both')
	#-- use an axis divider for the colorbar
	drx = make_axes_locatable(ax[i])
	cax = drx.append_axes("right", size="5%", pad=0.1)
	cbar = fig.colorbar(p,cax=cax)#orientation='horizontal',ax=ax[i])
	cbar.set_label('cm water equivalent')
	#-- plot continents 
	world.plot(ax=ax[i],alpha=0.3,fc='none',ec='k',linewidth=1.2,rasterized=True)
	
	#-- also plot the origina polygon
	region = sv.regions[0]
	n = len(region)
	t = np.linspace(0,1,10)
	for j in range(n):
		start = sv.vertices[region][j]
		end = sv.vertices[region][(j + 1) % n]
		edge = np.array(geometric_slerp(start, end, t))
		edge_lons,edge_lats = np.ones(len(edge)),np.ones(len(edge))
		for v in range(len(edge)):
			edge_lats[v],edge_lons[v] = pygplates.PointOnSphere(edge[v]).to_lat_lon()
		#-- plot the polygon edges
		ax[i].plot(edge_lons,edge_lats,color='springgreen',linewidth=1.2)
	
	ax[i].set_xlabel('Longitude (Degrees)')
	ax[i].set_title(r.upper())
	ax[i].set_ylabel('Latitude (Degrees)')
fig.suptitle("Fixed-Point (green) and Sample Far-Field Mascons")
plt.tight_layout()
plt.savefig(os.path.join(base_dir,'paper','harmonic_representation.pdf'),format='PDF')
plt.close(fig)
