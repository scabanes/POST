#! /usr/bin/env python
from netCDF4 import Dataset as NetCDFFile 
import numpy as np
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt
import matplotlib
from windspharm.standard import VectorWind
from windspharm.tools import prep_data
from windspharm.tools import get_recovery
from mpl_toolkits.basemap import Basemap
import matplotlib.ticker as ticker
import time
import sys
#
import LoadInfos



zefile = 'StatisticalData.nc'
path = 0
#zefile = '/u/scabanes/planeto/TesTsPectra/CI010/StatisticalData.nc'
#zefile = '/u/scabanes/planeto/TesTsPectra/CI007/spectra-iz9/StatisticalData.nc'

nc = NetCDFFile(zefile) # load file

wm = nc.variables['wm'][:] # [time,altitude,lon,lat] Zonal velocity
vort = nc.variables['vort'][:] # [time,altitude,lon,lat] Zonal velocity
lon = nc.variables['lon'][:] # [time,altitude,lon,lat] Zonal velocity
lat = nc.variables['lat'][:] # [time,altitude,lon,lat] Zonal velocity
dsteps = nc.variables['dsteps'][:] #[time]

itmax  = wm.shape[0]
it = itmax-1 # CI013
iz1= 5
iz2 = 11
#it = 7 # CI012

print itmax
print dsteps.shape
print wm.shape

(n,epsilon_dim,epsilon_f,beta_sat,beta_dim,Ck,Cz,E_dim,n_beta,n_R,n_f,n_D) = LoadInfos.LI ( path )

nrows = len(lat)
ncols = len(lon) 

lon_p, lat_p = np.meshgrid(np.linspace(0,360,ncols), np.linspace(-90,90,nrows))
#########################################################################################
#----------------------------------------------------------------------------------------
#				MAP SPHERE ZONAL VELOCITY
#----------------------------------------------------------------------------------------
#-------------------------------------------------------------iz1
fig = plt.subplots()
mmax = 20.0
vmin = -mmax
vmax = mmax-mmax/10
#
vmin2 = -mmax
vmax2 = mmax-mmax/10
# set up map projection
map = Basemap(projection='ortho', lat_0=25, lon_0=115)
# draw lat/lon grid lines every 30 degrees.
map.drawmeridians(np.arange(0, 360, 60))
map.drawparallels(np.arange(-90, 90, 60))
# compute native map projection coordinates of lat/lon grid.
x, y = map(lon_p, lat_p)
# contour data over the map.
vt = np.linspace(vmin, vmax, 110)
cs = map.contourf(x, y,np.transpose(wm[it,iz1,:,:]),vt)
plt.set_cmap('bwr') # a good start: blue to white to red colormap
#plt.set_cmap('jet') # a good start: blue to white to red colormap
cbar = map.colorbar(cs)
cs.set_clim(vmin2, vmax2)
plt.title('zonal velocity at '+ str(dsteps[it])+'at level = '+str(iz1))

#-------------------------------------------------------------iz2
fig = plt.subplots()
mmax = 40.0
vmin = -mmax
vmax = mmax-mmax/10
#
vmin2 = -mmax
vmax2 = mmax-mmax/10
# set up map projection
map = Basemap(projection='ortho', lat_0=25, lon_0=115)
# draw lat/lon grid lines every 30 degrees.
map.drawmeridians(np.arange(0, 360, 60))
map.drawparallels(np.arange(-90, 90, 60))
# compute native map projection coordinates of lat/lon grid.
x, y = map(lon_p, lat_p)
# contour data over the map.
vt = np.linspace(vmin, vmax, 110)
cs = map.contourf(x, y,np.transpose(wm[it,iz2,:,:]),vt)
plt.set_cmap('bwr') # a good start: blue to white to red colormap
#plt.set_cmap('jet') # a good start: blue to white to red colormap
cbar = map.colorbar(cs)
cs.set_clim(vmin2, vmax2)
plt.title('zonal velocity at '+ str(dsteps[it])+'at level = '+str(iz2))

#matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
#csb = map.contour(x, y, sf,10,colors='k')
#########################################################################################
#----------------------------------------------------------------------------------------
#				MAP SPHERE VORTICITY
#----------------------------------------------------------------------------------------
figHn, axFn = plt.subplots()
mmax = 0.00001
vmin = -mmax
vmax = mmax-mmax/10
#
vmin2 = -mmax
vmax2 = mmax-mmax/10
# set up map projection
map = Basemap(projection='ortho', lat_0=25, lon_0=115)
# draw lat/lon grid lines every 30 degrees.
map.drawmeridians(np.arange(0, 360, 60))
map.drawparallels(np.arange(-90, 90, 60))
# compute native map projection coordinates of lat/lon grid.
x, y = map(lon_p, lat_p)
# contour data over the map.
vt = np.linspace(vmin, vmax, 110)
cs = map.contourf(x, y,np.transpose(vort[it,iz1,:,:]),vt)
plt.set_cmap('bwr') # a good start: blue to white to red colormap
#plt.set_cmap('jet') # a good start: blue to white to red colormap
cbar = map.colorbar(cs)
cs.set_clim(vmin2, vmax2)
plt.title('Vorticity at '+ str(dsteps[it]))
#matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
#csb = map.contour(x, y, sf,10,colors='k')

plt.show()

###########################################
vmin = -50.0
vmax = 49.9
vt = np.linspace(vmin, vmax, 110)
cs = plt.contourf(vort[0,0,:,:])
plt.set_cmap('bwr') # a good start: blue to white to red colormap
plt.colorbar(orientation='vertical')
plt.show()





