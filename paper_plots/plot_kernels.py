#!/usr/bin/env python
u"""
plot_kernels.py
by Yara Mohajerani

multi-panel plot of kernels for paper

Last Update 04/2021
"""
#-- load required modules
import os
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
#-- also import gravity toolkit modules
from gravity_toolkit.ncdf_read import ncdf_read

#-- directory setup
base_dir = os.path.join(os.getcwd(),'..','data')

#-- load in world map for plotting in background
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

# regs = {'karakoram_NW':{'region':'karakoram','label':'NW_18+21','lat_range':[15,45],'lon_range':[60,110]},
# 		'karakoram_NE':{'region':'karakoram','label':'SE_0+30','lat_range':[15,45],'lon_range':[60,110]},
# 		'nyainqentangla':{'region':'nyainqentangla','label':'0+24+27','lat_range':[15,45],'lon_range':[60,110]},
# 		'alaska':{'region':'alaska','label':'50+56+58+170','lat_range':[45,75],'lon_range':[-175,-120]}
# 		#'alaska_SE':{'region':'alaska','label':'east_42+43+163+170','lat_range':[45,75],'lon_range':[-175,-120]}
# 		}

regs = {'karakoram_NW':{'region':'karakoram','label':'NW_1+25','lat_range':[15,45],'lon_range':[60,110],'dir':'rotated_version'},
		'karakoram_NE':{'region':'karakoram','label':'SE_0+57+73','lat_range':[15,45],'lon_range':[60,110],'dir':'rotated_version'},
		'nyainqentangla':{'region':'nyainqentangla','label':'0+24+27','lat_range':[15,45],'lon_range':[60,110],'dir':''},
		'alaska':{'region':'alaska','label':'0+64+78','lat_range':[45,75],'lon_range':[-175,-120],'dir':'rotated_version'}
		#'alaska_SE':{'region':'alaska','label':'east_42+43+163+170','lat_range':[45,75],'lon_range':[-175,-120]}
		}

fig, ax = plt.subplots(2,2,figsize = (10,6))

for c,r in enumerate(regs.keys()):
	i = int(c/2) # row number
	j = c%2	# column number
	#-- read kernel
	kern = ncdf_read(os.path.join(base_dir,regs[r]['region'],regs[r]['dir'], 'MSCNS','MASCON_{0}_{1}_YLMS_0.25DEG_SKERNEL_OCN_L60_r250km.nc'.\
		format(regs[r]['region'],regs[r]['label'])),DATE=False)
	#-- plot kernel
	glon,glat = np.meshgrid(kern['lon'],kern['lat'])
	p = ax[i,j].pcolormesh(glon,glat,kern['data'],cmap='RdYlBu',vmin=-1.6,vmax=1.6,shading='auto',rasterized=True)
	# levels = [-1.5,-1,-0.5,-0.2,0,0.2,0.5,0.8,1,1.5]	#np.linspace(-1.6, 1.6, 17)
	# p = ax[i,j].contourf(glon,glat,kern['data'],levels,cmap='RdBu',extend='both')
	cs = ax[i,j].contour(glon,glat,kern['data'],[1.],colors='lime')#,cmap='RdBu',extend='both')
	#ax[i,j].clabel(cs, inline=True,fontsize=5,inline_spacing=0)
	#-- use an axis divider for the colorbar
	drx = make_axes_locatable(ax[i,j])
	cax = drx.append_axes("right", size="5%", pad=0.1)
	cbar = fig.colorbar(p,cax=cax)#orientation='horizontal',ax=ax[i])
	cbar.set_label('cm water equivalent')
	#-- plot continents 
	world.plot(ax=ax[i,j],alpha=0.3,fc='none',ec='k',linewidth=1.2,rasterized=True)
	#-- set limits
	ax[i,j].set_xlim(regs[r]['lon_range'])
	ax[i,j].set_ylim(regs[r]['lat_range'])
	out_title = r .replace('_',' ').upper()
	ax[i,j].set_title(out_title)
	if j == 0:
		ax[i,j].set_ylabel('Latitude (Degrees)')
	if i == len(regs)-1:
		ax[i,j].set_xlabel('Longitude (Degrees)')
	ax[i,j].text(-0.05, 1.18, chr(c+65), transform=ax[i,j].transAxes,
      fontsize=14, fontweight='bold', va='top', ha='right')
plt.subplots_adjust(wspace=0, hspace=0)
plt.tight_layout()
plt.savefig(os.path.join(base_dir,'paper','kernels_new.pdf'),format='PDF')
plt.close(fig)
