#! /usr/bin/env python
from netCDF4 import Dataset as NetCDFFile 
import numpy as np
import math
from pylab import *
from decimal import *
# --> Librairie Simon
import windowing
import FourierTransform1D
############################################################################################
# 		    PTS: "Parameters to set"
############################################################################################
zefile = 'StatisticalData-FullTime.nc'
#---------------------------------------------------------------Mode azimutaux
im = 3 		# indice en longitude
#---------------------------------------------------------------Moyenne glissee en temps & truncation
nbT = 3 	# -nbT de donnees enleve a la fin du signal
Tstep = 1 	# time step pour la moyenne glisse
#---------------------------------------------------------------for multiple value of l
nbn = 3 	# va jusqu a l indice nbl
omega_sat = 0.000165121 #rad.s-1
JourSaturn_s = 118.9125*320# jour saturne en seconde est dt=118.9125 fois 320 le nombre de pas en temps pour faire une journee ~(2*pi)/0.000165121 = 38052,00614809
############################################################################################
# 		    Load Data
############################################################################################
#ptimestep= 3805.2 # time step en seconde entre chaque injection
# ..........LOAD data from StatisticalData.nc
nc = NetCDFFile(zefile) # load file
ET = nc.variables['ET'][:] #[time]
#dsteps = nc.variables['dsteps'][:-1000] #[time]
time_counter = nc.variables['time_counter'][:] #[time]
#epsilonf = nc.variables['epsilonf'][:-1000] #[time]
#............SIGNAL TOTAL LENGTH
L_T = (time_counter[-1]-time_counter[0])
N_T = len(time_counter)
print '..........................INITIAL SIGNAL LENGTH'
print 'Signal length : ',N_T 
print 'Signal length in days : ',L_T/JourSaturn_s
print 'Signal length in second : ',L_T
#............SIGNAL USED FOR FFT
LL_T = (time_counter[0:-nbT+0]-time_counter[0])# exclude the value -nbT
L_T = LL_T[-1]
N_T = len(time_counter[0:-nbT+0])
#FFT_YTuk= np.zeros([N_T/2+1,nbT])
print '..........................LENGTH USED FOR FFT'
print 'Signal length : ',N_T 
print 'Signal length in days : ',L_T/JourSaturn_s
print 'Signal length in second : ',L_T
############################################################################################
#sf_r = nc.variables['sf_r'][:-nbT,0,1,0] #[t,z,n,m]
#print sf_r.shape
#(sfr_Tuk) = windowing.WindowingTukey(sf_r,N_T,L_T,1 )#BLM_t[l,m,t]
#( FFT_YTuk,E_YTuk[:,iz,it],modes ) = FourierTransform1D.fft1D ( sfr_Tuk, N_T )
#exit()
############################################################################################
# 		    WRITE NETCEDEF
############################################################################################
# get dimensions from netCDF file
zedim = nc.dimensions.keys()
loncoord = zedim[0]
latcoord = zedim[1]
zcoord = zedim[2]
tcoord = zedim[3]
fcoord = "fq"
mcoord = "m"
ncoord = "n"

presnivs = nc.variables['altitude'][:]
# creer un fichier netcdef --> http://www.ceda.ac.uk/static/media/uploads/ncas-reading-2015/11_create_netcdf_python.pdf
dataset = NetCDFFile('TempModalSpectra'+'-im-'+str(im)+'-upto-n-'+str(nbn)+'.nc', 'w', format='NETCDF3_CLASSIC')
       	#----Dimensions-----------------------------------------------------------------------------
presnivs = dataset.createDimension(zcoord, len(presnivs))
presnivss = dataset.createVariable(zcoord, np.float64, (zcoord,))
#presnivss[:] = presnivs

fqmodes = np.arange(0.,N_T/2+1) # Attention!!! la gamme des modes demarre a 0 et pas 1, lunquist inclu
w = 2.*np.pi*fqmodes/L_T #[2*pi*m/LT] = rad.s-1
frequence = dataset.createDimension(fcoord, N_T/2+1)
frequences = dataset.createVariable('frequences', np.float64, (fcoord,))
frequences[:] = w

#iz_counter = dataset.createDimension(zcoord, None)
#mz_counter = dataset.createDimension(mcoord, None)
nz_counter = dataset.createDimension(ncoord, None)

EFT= dataset.createVariable('EFT', np.float64, (ncoord,fcoord,zcoord))
omega_RHW= dataset.createVariable('omega_RHW', np.float64, (ncoord))



#mmode = dataset.createDimension(mcoord, 4)
#mmodes = dataset.createVariable(mcoord, np.float64, (mcoord,))
#mmodes[:] = 

#nmode = dataset.createDimension(ncoord, 5)
#nmodes = dataset.createVariable(ncoord, np.float64, (ncoord,))
#nmodes[:] = 

niz=len(presnivs)
Trange =  np.arange(0, nbT,Tstep)
#inrange = np.arange(im, nbn+1)
E_YTuk= np.zeros([N_T/2+1,niz,len(Trange)])
#_FFT_= np.zeros([N_T/2+1,niz])
#ToOMEGA = 2.*omega_sat
#omega_RHW[il,im] = (2.*omega_sat*float(im))/(float(il)*(float(il)+1.)) # en rad.s-1
for ni in range(im, nbn+1):
	print ni
	for iz in range(0,niz):
		sf_r = nc.variables['sf_r'][:,iz,ni,im] # [time,altitude,n,m]
		print sf_r.shape
		for it in Trange:#range(0, nbT,Tstep)
			(sfr_Tuk) = windowing.WindowingTukey(sf_r[it:-nbT+it],N_T,L_T,0 )#BLM_t[l,m,t]
			( EFT_YTuk,E_YTuk[:,iz,it],modes ) = FourierTransform1D.fft1D ( sfr_Tuk, N_T )


	EFT[ni,:,:] = np.mean(E_YTuk[:,:,:],2)#[n,fq,iz]
	omega_RHW[ni] = (2.*omega_sat*float(im))/(float(ni)*(float(ni)+1.)) # en rad.s-1
	print 'completed n = '+str(ni)+' on '+str(nbn)
	E_YTuk= np.zeros([N_T/2+1,niz,len(Trange)])

del E_YTuk

dataset.close()

print '///////////////////',dataset.file_format,'named: sfvpData.nc writen///////////////'

print '----------------------------------------------------------------------------------'
print '------ Created file netcdf'
ncData = NetCDFFile('TempModalSpectra'+'-im-'+str(im)+'-upto-n-'+str(nbn)+'.nc')
for varname in ncData.variables.keys():
        var = ncData.variables[varname]
        print varname, var.dtype, var.dimensions, var.shape
print '----------------------------------------------------------------------------------'




