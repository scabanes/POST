from netCDF4 import Dataset as NetCDFFile
import numpy as np
import time
import sys

####################################################################
#                   PTS: "Parameters to set"
####################################################################
#****************************************************************
# This code is preceded of CropData and suceded by 
# statistical_analysis_FullTime. It leads to long time series at
# a single depth layer. All time step are selected and not adjustable.
# -----
# INPUTS 1:
# path =  est le chemin vers les Xhinsts.
# Seront selectionnes les fichiers allant de 
# Xhinst_init : Xhinst_'Xhinst_end'.
# Avec les iterations en temps souhaite:
# time_range_i = premiere iteration souhaitee du fichier Xhinst --> uniquement indicatif ici
# time_range_e = derniere iteration souhaitee du fichier Xhinst --> uniquement indicatif ici
# time_step =    time step souhaitee du fichier Xhinst		--> uniquement indicatif ici
# iz = selectionne le niveau en profondeur depuis le fond de 
#      de l'atmosphere iz=0 a iz=32
# INPUTS 2 :
# Inputs d inormation sur les fichiers Xhinsts.
# fq_Time_sampling = la frequence d enregistrement en jours.
# time_lenght =  nombre d iteration dans un fichier Xhinst
# -----
# OUTPUTS :
# a file uvData-FullTime-istep-XX-nstep-XXX-izX.nc with the velocity components
# u and v, dsteps in days only indicative, time_counter is more 
# precise time in seconds, lat, lon.
#*****************************************************************
path = '/store/cabanes/CI018/TempSeries/'
#-------------------------- PTS : Parameter To Set --
Xhinst_init = 1
Xhinst_end = 6#13#21 # +1 que le numero du dernier fichier Xhinst
iz=7
#--------------------------
le=115
JourSaturn_s = 118.9125*320# jour saturne en seconde est dt=118.9125 fois 320 le nombre de pas en temps pour faire une journee ~(2*pi)/0.000165121 = 38052,00614809
#-------------------------- parametres modulables
#0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,1,1,1,1,1]
#1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,]
time_range_i= np.ones(le,dtype=int)*0 #[0,0,0,0,0,0,0,0,0,0] #
time_range_e=[20,20,20,200,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400]#it18
time_step= np.ones(le,dtype=int)*1
#-------------------------- parmetres fixes
fq_Time_sampling= [0.1,0.1,0.1,0.1,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]#it17
time_lenght= [20,20,20,200,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400]#it18
#####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   Infos processing
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# I am using my knowledge of the sampling frequency and time length
# of each file Xhistins to comple a step days file. It appears however
# that when the sampling is of 0.1 day, a delay of -0.19921875 s exist
# in step_days that differ from time_counter of the Xhistins files. 
# I would rather believe Xhistins. I also have to add a day, 38052 s
# to step_days that does not account for the makestart step to initiate
# the radiative profile at the begening of a simulation.
zefile = path+'Xhistins_1-uvtiz'+str(iz)+'.nc'
#zefile = path+'Xhistins_1-to-3-uvtiz7.nc'

nb_Xhinst = Xhinst_end-Xhinst_init
nb_step=0#len(np.arange(time_range[0],time_range[1],time_step))*(Xhinst_end-Xhinst_init)
step_days = []
step_it = []
if Xhinst_init == 1:
	step_prec= 0
else:
	step_prec = np.array(time_lenght[0:Xhinst_init-2+1])*np.array(fq_Time_sampling[0:Xhinst_init-2+1] )
step_count = np.sum(step_prec)
print step_count
#print len(time_lenght), len(fq_Time_sampling), len(time_range_i)
##------------------------------------------------------------------
iXh = -1 
if Xhinst_end-Xhinst_init > len(fq_Time_sampling):
        print '!!! -----> The readme.time file needs the information for the saving frequency in days of the Xhinst_.nc, size:'
        print len(fq_Time_sampling)
        print 'instead of:'
        print Xhinst_end-Xhinst_init
	exit()
else:
	for iXhinst in range(Xhinst_init-1,Xhinst_end-1):
		iXh = iXh + 1
		nb_step = nb_step + len(np.arange(time_range_i[iXhinst],time_range_e[iXhinst],time_step[iXhinst]))
		if iXh == 0:
                	step_array = (np.arange(time_range_i[iXhinst]+1,time_range_e[iXhinst]+1,time_step[iXhinst])*fq_Time_sampling[iXhinst]) + step_count
                else:
			step_count = step_count + time_lenght[iXhinst-1]*fq_Time_sampling[iXhinst-1]
			step_array = (np.arange(time_range_i[iXhinst]+1,time_range_e[iXhinst]+1,time_step[iXhinst])*fq_Time_sampling[iXhinst]) + step_count

		step_days = np.concatenate((step_days,step_array),axis=0)
		print step_days

# I have to add a day that is the makestart day to initiate the simulation
step_days = step_days+1

del iXh
#exit()
##------------------------------------------------------------------
print 'The infos_AS array ---------- '
print #step_days
print '----------------------------- '
####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                         Data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##-----------------
nc = NetCDFFile(zefile) # load file
# get dimensions from netCDF file
zedim = nc.dimensions.keys()
altdim = zedim[0]
latdim = zedim[1]
londim = zedim[2]
timedim = zedim[3]
axisdim = zedim[4]
#pcoord = "profondeur"
print 'coordonates dimensions'
print zedim
#
print '----------------------------------------------------------------------------------'
print '----- Initial file uploaded'
for varname in nc.variables.keys():
  	var = nc.variables[varname]
        print varname, var.dtype, var.dimensions, var.shape
print '----------------------------------------------------------------------------------'
#
lat_AS = nc.variables[latdim][:] #(128,)
lon_AS = nc.variables[londim][:]
presnivs_AS = nc.variables[altdim][:]
#u_AS = nc.variables['u'][t,Prof(une seule),theta,phi] Zonal velocity
#v_AS = nc.variables['v'][t,Prof,theta,phi] meridional velocity
#-----------------------------------------------------------------------
#                       WRITE NETCDF
#-----------------------------------------------------------------------
# creer un fichier netcdef --> http://www.ceda.ac.uk/static/media/uploads/ncas-reading-2015/11_create_netcdf_python.pdf
dataset = NetCDFFile('uvData-FullTime-istep-'+str(np.sum(step_prec))+'-nstep-'+str(nb_step)+'-iz'+str(iz)+'.nc', 'w', format='NETCDF3_CLASSIC')
       	#----Dimensions-----------------------------------------------------------------------------
presnivs = dataset.createDimension(altdim, len(presnivs_AS))
presnivss = dataset.createVariable(altdim, np.float64, (altdim,))
presnivss[:] = presnivs_AS
	#
#profondeurs = dataset.createDimension(pcoord, len(profondeurs_AS))
#profondeurss = dataset.createVariable(pcoord, np.float64, (pcoord,))
#profondeurss[:] = profondeurs_AS
       #
longitude = dataset.createDimension(londim, len(lon_AS))
longitudes = dataset.createVariable(londim, np.float64, (londim,))
longitudes[:] = lon_AS
       #
latitude = dataset.createDimension(latdim, len(lat_AS))
latitudes = dataset.createVariable(latdim, np.float64, (latdim,))
latitudes[:] = lat_AS
       #
time_counter = dataset.createDimension(timedim, None)
       #
time_counter = dataset.createVariable('time_counter', np.float64, (timedim))
v = dataset.createVariable('v', np.float64, (timedim,altdim,latdim,londim))
u = dataset.createVariable('u', np.float64, (timedim,altdim,latdim,londim))
       #
step_infos = dataset.createDimension('step_infos', len(step_days))
dsteps = dataset.createVariable('dsteps', 'f4', ('step_infos',))
dsteps[:] = step_days
       #
#stepit = dataset.createDimension('stepit', len(step_it))
#StepIt = dataset.createVariable('StepIt', 'i4', ('stepit',))
#StepIt[:] = step_it

nb_t_T = 0
for iXhinst in range(Xhinst_init,Xhinst_end):
	print '--------------  iXhinst  '+str(iXhinst)
	zefile = path+'Xhistins_'+str(iXhinst)+'-uvtiz'+str(iz)+'.nc'
       #
	nc = NetCDFFile(zefile) # load file
	#
	print np.arange(time_range_i[iXhinst-1],time_range_e[iXhinst-1],time_step[iXhinst-1])
	#
	u_AS = nc.variables['u'][:,:,:,:] # [t,Prof,theta,phi] Zonal velocity
	v_AS = nc.variables['v'][:,:,:,:] # [t,Prof,theta,phi] meridional velocity
	time_counter_AS = nc.variables['time_counter'][:]
	#
	#nb_t = len(np.arange(time_range_i[iXhinst-1],time_range_e[iXhinst-1],time_step[iXhinst-1]))
	nb_t = u_AS.shape[0]
	#
	print 'step 1 done.....'
	print v_AS.shape
	time_counter[nb_t_T:nb_t_T+nb_t] = time_counter_AS
	v[nb_t_T:nb_t_T+nb_t,:,:,:] = v_AS
        u[nb_t_T:nb_t_T+nb_t,:,:,:] = u_AS
	print 'step 2 done.....'
	#
	nb_t_T = nb_t_T + nb_t
	
	del v_AS, u_AS

dataset.close()

print '///////////////////',dataset.file_format,'named: sfvpData.nc writen///////////////'

print '----------------------------------------------------------------------------------'
print '------ Created file netcdf'
ncData = NetCDFFile('uvData-FullTime-istep-'+str(np.sum(step_prec))+'-nstep-'+str(nb_step)+'-iz'+str(iz)+'.nc')
for varname in ncData.variables.keys():
        var = ncData.variables[varname]
        print varname, var.dtype, var.dimensions, var.shape
print '----------------------------------------------------------------------------------'


####################################################################
#                       OUTPUTS
####################################################################
# sfvpData-xxx-xxx.nc : 
# A file netcdet that contain velocities, u,v and sf, vp the stream-
# function and the velocity potential. 
