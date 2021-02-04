from netCDF4 import Dataset as NetCDFFile
import numpy as np
import time
import sys

####################################################################
#                   PTS: "Parameters to set"
####################################################################
#****************************************************************
# INPUTS 1:
# path =  est le chemin vers les Xhinsts.
# Seront selectionnes les fichiers allant de 
# Xhinst_init : Xhinst_'Xhinst_end'.
# Avec les iterations en temps souhaite:
# time_range_i = premiere iteration souhaitee du fichier Xhinst
# time_range_e = derniere iteration souhaitee du fichier Xhinst 
# time_step =    time step souhaitee du fichier Xhinst
# iz = selectionne les niveaux en profondeur depuis le fond de 
#      de l'atmosphere iz=0 a iz=32
# INPUTS 2 :
# Inputs d inormation sur les fichiers Xhinsts.
# fq_Time_sampling = la frequence d enregistrement en jours.
# time_lenght =  nombre d iteration dans un fichier Xhinst
#*****************************************************************
path_e = '/store/aspigaplaneto/simu_REF_100days_continued/'
path = '/store/aspigaplaneto/aspigaheat-store/dynamico_giant/saturn/simu_REF_100days/'
#--------------------------
Xhinst_init = 170#200  #150 #313 #292#1  #273#147#1#269
Xhinst_end = 174#370 #202#375   #165#375 #315 #321#21  #360#171#23#266 # +1 que le numero du dernier fichier Xhinst
le=380
#--------------------------
time_range_i=np.ones(le,dtype=int)*0 #[0,0,0,0,0,0,0,0,0,0] #
time_range_e=np.ones(le,dtype=int)*50 #[50,50,50,50,50,50,50] # !! Atention nombre max iteration dans le fichier Xhinst_*
time_step=   np.ones(le,dtype=int)*10#20#[10,1,25,25,25,25,25]
#time_step=   [2,5,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10]
#--------------------------
#iz=[0,1,2,8,9,10,18,19,20,29,30,31] # choisir un niveau en profondeur
#iz=[5,7,9,10,11,13,15,16,17,19,21,23] # choisis pour ICARUS
iz=[4,5,6,7,8,9,10,11,12,13,14,15,17]# troposphere + lower stratosphere
iz=[19,20,21,22,23,24,25,26,27,28,29] # upper stratosphere
iz=[6,7,8,9,10,11,12,13]# troposphere + lower stratosphere --> pour ZENODO
#iz=[15]
iphi=0
fq_Time_sampling= np.ones(le,dtype=int)*20#[20,20,20,20,20,20,20,20,20]
time_lenght= np.ones(le, dtype=int)*50
#*****************************************************************
#####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                   Infos processing
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zefile = path+'Xhistins_170.nc'
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
ccoord = zedim[0]
xcoord = zedim[1]
ycoord = zedim[2]
zcoord = zedim[3]
tcoord = zedim[4]
pcoord = "profondeur"
print 'coordonates length'
print zedim,zedim[3]

nbiz=len(iz) # On est pour le moment limite a un seul iz.
# On forme la martice presniv_AS qui n a qu une seule valeure. 
# C est pour quelle est une longueure.
presnivs_AS = np.zeros([nbiz])
#
print '----------------------------------------------------------------------------------'
print '----- Initial file uploaded'
for varname in nc.variables.keys():
  	var = nc.variables[varname]
        print varname, var.dtype, var.dimensions, var.shape
print '----------------------------------------------------------------------------------'
#
lat_AS = nc.variables[ycoord][:] #(128,)
lon_AS = nc.variables[xcoord][:]
presnivs_AS[0:nbiz] = nc.variables[zcoord][iz]
profondeurs_AS = nc.variables[zcoord][:]
print nc.variables[zcoord][iz]
print profondeurs_AS

#
#u_AS = nc.variables['u'][time_range[0]:time_range[1]:time_step,iz,:,:] # [t,Prof,theta,phi] Zonal velocity
#v_AS = nc.variables['v'][time_range[0]:time_range[1]:time_step,iz,:,:] # [t,Prof,theta,phi] meridional velocity
#-----------------------------------------------------------------------
#                       WRITE NETCDF
#-----------------------------------------------------------------------
# creer un fichier netcdef --> http://www.ceda.ac.uk/static/media/uploads/ncas-reading-2015/11_create_netcdf_python.pdf
dataset = NetCDFFile('uvData-istep-'+str(np.sum(step_prec))+'-nstep-'+str(nb_step)+'-niz-'+str(len(iz))+'.nc', 'w', format='NETCDF3_CLASSIC')
       	#----Dimensions-----------------------------------------------------------------------------
presnivs = dataset.createDimension(zcoord, len(presnivs_AS))
presnivss = dataset.createVariable(zcoord, np.float64, (zcoord,))
presnivss[:] = presnivs_AS
	#
profondeurs = dataset.createDimension(pcoord, len(profondeurs_AS))
profondeurss = dataset.createVariable(pcoord, np.float64, (pcoord,))
profondeurss[:] = profondeurs_AS
       #
longitude = dataset.createDimension(xcoord, len(lon_AS))
longitudes = dataset.createVariable(xcoord, np.float64, (xcoord,))
longitudes[:] = lon_AS
       #
latitude = dataset.createDimension(ycoord, len(lat_AS))
latitudes = dataset.createVariable(ycoord, np.float64, (ycoord,))
latitudes[:] = lat_AS
       #
time_counter = dataset.createDimension(tcoord, None)
       #
time_counter = dataset.createVariable('time_counter', np.float64, (tcoord))
v = dataset.createVariable('v', np.float64, (tcoord,zcoord,ycoord,xcoord))
u = dataset.createVariable('u', np.float64, (tcoord,zcoord,ycoord,xcoord))
	#
v_SMerid = dataset.createVariable('v_SMerid', np.float64, (tcoord,pcoord,ycoord))
u_SMerid = dataset.createVariable('u_SMerid', np.float64, (tcoord,pcoord,ycoord))
      #
#dulat_diss1 = dataset.createVariable('dulat_diss1', np.float64, (tcoord,zcoord,ycoord,xcoord))
#dulon_diss1 = dataset.createVariable('dulon_diss1', np.float64, (tcoord,zcoord,ycoord,xcoord))
#dulat_diss2 = dataset.createVariable('dulat_diss2', np.float64, (tcoord,zcoord,ycoord,xcoord))
#dulon_diss2 = dataset.createVariable('dulon_diss2', np.float64, (tcoord,zcoord,ycoord,xcoord))
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
	zefile = path+'Xhistins_'+str(iXhinst)+'.nc'
	#zefile2 = path+'Xdissip_'+str(iXhinst)+'.nc'
       #
	nc = NetCDFFile(zefile) # load file
	#nc2 = NetCDFFile(zefile2) # load file
	#
	print np.arange(time_range_i[iXhinst-1],time_range_e[iXhinst-1],time_step[iXhinst-1])
	#*********************************** load
	u_AS = nc.variables['u'][time_range_i[iXhinst-1]:time_range_e[iXhinst-1]:time_step[iXhinst-1],iz,:,:] # [t,Prof,theta,phi] Zonal velocity
	v_AS = nc.variables['v'][time_range_i[iXhinst-1]:time_range_e[iXhinst-1]:time_step[iXhinst-1],iz,:,:] # [t,Prof,theta,phi] meridional velocity
	u_SMerid_AS = nc.variables['u'][time_range_i[iXhinst-1]:time_range_e[iXhinst-1]:time_step[iXhinst-1],:,:,iphi] # [t,Prof,theta,phi] Zonal velocity
	v_SMerid_AS = nc.variables['v'][time_range_i[iXhinst-1]:time_range_e[iXhinst-1]:time_step[iXhinst-1],:,:,iphi] # [t,Prof,theta,phi] meridional velocity
	time_counter_AS = nc.variables['time_counter'][time_range_i[iXhinst-1]:time_range_e[iXhinst-1]:time_step[iXhinst-1]]
	#
        #dulat_diss1_AS = 0#nc2.variables['dulat_diss1'][time_range_i[iXh-1]:time_range_e[iXh-1]:time_step[iXh-1],iz,:,:] # [t,Prof,theta,phi] 
        #dulon_diss1_AS = 0#nc2.variables['dulon_diss1'][time_range_i[iXh-1]:time_range_e[iXh-1]:time_step[iXh-1],iz,:,:] # [t,Prof,theta,phi]
        #dulat_diss2_AS = 0#nc2.variables['dulat_diss2'][time_range_i[iXh-1]:time_range_e[iXh-1]:time_step[iXh-1],iz,:,:] # [t,Prof,theta,phi] 
        #dulon_diss2_AS = 0#nc2.variables['dulon_diss2'][time_range_i[iXh-1]:time_range_e[iXh-1]:time_step[iXh-1],iz,:,:] # [t,Prof,theta,phi]
	#
	nb_t = len(np.arange(time_range_i[iXhinst-1],time_range_e[iXhinst-1],time_step[iXhinst-1]))
	#
	print v_AS.shape
	#*********************************** Save
	time_counter[nb_t_T:nb_t_T+nb_t] = time_counter_AS
	v[nb_t_T:nb_t_T+nb_t,0:nbiz,:,:] = v_AS
        u[nb_t_T:nb_t_T+nb_t,0:nbiz,:,:] = u_AS
	#
	v_SMerid[nb_t_T:nb_t_T+nb_t,:,:] = v_SMerid_AS
        u_SMerid[nb_t_T:nb_t_T+nb_t,:,:] = u_SMerid_AS
	#
	#print dulat_diss1_AS.shape
        #dulat_diss1[nb_t_T:nb_t_T+nb_t,0:nbiz,:,:] = dulat_diss1_AS
        #dulon_diss1[nb_t_T:nb_t_T+nb_t,0:nbiz,:,:] = dulon_diss1_AS
        #dulat_diss2[nb_t_T:nb_t_T+nb_t,0:nbiz,:,:] = dulat_diss2_AS
        #dulon_diss2[nb_t_T:nb_t_T+nb_t,0:nbiz,:,:] = dulon_diss2_AS
	#
	nb_t_T = nb_t_T + nb_t
	#nb_t = len(np.arange(time_range[0],time_range[1],time_step))
	#v[nb_t*(iXhinst-1):nb_t*(iXhinst-1)+nb_t,nbiz-1,:,:] = v_AS
	#u[nb_t*(iXhinst-1):nb_t*(iXhinst-1)+nb_t,nbiz-1,:,:] = u_AS
	
	del v_AS, u_AS, u_SMerid_AS, v_SMerid_AS, time_counter_AS

dataset.close()

print '///////////////////',dataset.file_format,'named: StatisticalData.nc writen///////////////'

print '----------------------------------------------------------------------------------'
print '------ Created file netcdf'
ncData = NetCDFFile('uvData-istep-'+str(np.sum(step_prec))+'-nstep-'+str(nb_step)+'-niz-'+str(len(iz))+'.nc')
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
