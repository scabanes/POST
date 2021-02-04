#! /usr/bin/env python
from netCDF4 import Dataset as NetCDFFile 
import numpy as np
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt
import matplotlib
import time
import sys
sys.path.append('/u/scabanes/planeto/SPECTRA/POST/plots')
import LoadInfos
#...........Parametres
fonthesize = 14
matplotlib.rcParams.update({'font.size': 14})
#-------------------------------------------fille.PTS.plots
pp = open('filePTS.plots.infos', "r")
print pp.readline()
print pp.readline()
print pp.readline()
omega_ref = float(pp.readline().split()[1])
pp.readline()
itime = (pp.readline().split()[1:])
itime= [int(i) for i in itime] # from list to array of integer
pp.readline()
iz1 = int(pp.readline().split()[1])
iz2 = int(pp.readline().split()[1])
pp.readline()
nbtime = int(pp.readline().split()[1])
# ..........File to load
zefile = 'StatisticalData.nc'
path=0
nc = NetCDFFile(zefile) # load file
# ..........LOAD data from StatisticalData.nc
hn_zm = nc.variables['hn_zm'][:] # [t,Prof,theta,phi] Zonal velocity
fn_zm = nc.variables['fn_zm'][:] # [t,Prof,theta,phi] Zonal velocity
Emn = nc.variables['Emn'][:,:,:,0:12] # [time,altitude,n,m]
Enmo = nc.variables['Enmo'][:] # [t,Prof,theta,phi] Zonal velocity
Ene = nc.variables['Ene'][:] # [t,Prof,theta,phi] Zonal velocity
hn = nc.variables['hn'][:] # [t,Prof,theta,phi] Zonal velocity
fn = nc.variables['fn'][:] # [t,Prof,theta,phi] Zonal velocity
ET = nc.variables['ET'][:] #[time]
ER = nc.variables['ER'][:] #[time]
EZ = nc.variables['EZ'][:] #[time]
En_rot = nc.variables['En_rot'][:] # [time,altitude,latitude]
En_div = nc.variables['En_div'][:] # [time,altitude,latitude]
dsteps = nc.variables['dsteps'][:] #[time]
epsilonf = nc.variables['epsilonf'][:] #[time]
altitude = nc.variables['altitude'][iz1:iz2] #[time]
n = (np.arange(0.,360,1.))
#
itmax  = fn_zm.shape[0]
izmax  = fn.shape[1]
#............LOAD infos
(n,L_dim,T_dim,SatDay_s,omega_sat,epsilon_dim,epsilon_f,epsilon_i,EnstF_dim,beta_sat,beta_dim,Ck,Cz,E_dim,n_beta,n_R,n_f,n_D) = LoadInfos.LI ( path,ET  )
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#									 ENERGY FLUXES:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (0==0):
	figA = plt.subplots()
	plt.semilogx(n,np.mean(np.mean(fn[nbtime:,iz1:iz2,:],0),0)*(10**(5)), '-', color='black', linewidth=1.5) 
	plt.ylabel(r'$\Pi_n^E  (10^{-5} W kg^{-1})$', fontsize=fonthesize)
	plt.xlabel(r'$n$', fontsize=fonthesize)
	plt.xlim((0,361.))

# -------------------Typpical scales
	plt.axvline(n_beta, color='black', linestyle=':', linewidth=1.,  alpha=1)
	plt.axvline(n_R, color='black', linestyle=':', linewidth=1.,  alpha=1)
	#plt.axvline(2.*n_f, color='red', linestyle='solid', linewidth=1.5)
	plt.axhline(0., color='black', linestyle='-', linewidth=0.5,  alpha=1)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#								      ENSTROPHY FLUXES:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (0==0):
	figA = plt.subplots()
	plt.semilogx(n,np.mean(np.mean(hn[nbtime:,iz1:iz2,:],0),0), '-', color=(0,0.4,0.6), linewidth=1.5) 
	plt.ylabel(r'$\Pi_n^Z  (s^{-3})$', fontsize=fonthesize)
	plt.xlabel(r'$n$', fontsize=fonthesize)
	plt.xlim((0,361.))

# -------------------Typpical scales
	plt.axvline(n_beta, color='black', linestyle=':', linewidth=1.,  alpha=1)
	plt.axvline(n_R, color='black', linestyle=':', linewidth=1.,  alpha=1)
	#plt.axvline(2.*n_f, color='red', linestyle='solid', linewidth=1.5)
	plt.axhline(0., color='black', linestyle='-', linewidth=0.5,  alpha=1)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#								      	        ENERGY:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
en=1
ax,figB = plt.subplots()

#plt.loglog(n[0:-en],np.mean(np.mean(Enmo[nbtime:,iz1:iz2,0:-en],0),0), '-', color='darkred', linewidth=1.5) 

#_______________________________________________________________________________________________________________________
# ------------ Here is how Galperin et al. (2014) https://doi.org/10.1016/j.icarus.2013.08.030 plot the spectra
#plt.loglog(n[0:-en],np.mean(np.mean(Ene[nbtime:,iz1:iz2,0:-en],0),0)*np.pi*L_dim, '-', color='black', linewidth=1.5) 
#plt.loglog(n[0:-en],np.mean(np.mean(En_rot[nbtime:,iz1:iz2,0:-en],0),0)*np.pi*L_dim, '-', color='blue', linewidth=1.) 
#plt.loglog(n[1:],Ck*(epsilon_f**(2./3.))*((n[1:]/(2.*np.pi*L_dim))**(-5./3.)), '--', color='black', linewidth=1.5)

#_______________________________________________________________________________________________________________________
# ------------ Here is how it has to be ploted, by normalizing the energy extracted from flow analysis with [E*2*R] = m^3.s^-2, then we have
# 	       an enrgy density. 2*R likely comes from the integral over m, the zonal modes. It leads to the same energetic amplitude than
#	       the one obtained by Young & Read, NatPhys DOI: 10.1038/NPHYS4227, by a totally different approach, using structures function
#	       in cartesian, see details in supplementary materials - it leads for Jupiter Cassini data to a fit with epsilon = 9.10^-5 
#	       m^3.s^-2. Then to be consistent the fit is E_theoric = Ck eps^2/3 (n/2R)^-5/3, [E_theoric] = m^3.s^-2. Or I can do 
#	       E_theoric/2R if I want everything in m^2 s^-2, which is an energy.
plt.loglog(n[0:-en],np.mean(np.mean(Enmo[nbtime:,iz1:iz2,0:-en],0),0), '-', color='darkred', linewidth=1.5) 
plt.loglog(n[0:-en],np.mean(np.mean(Ene[nbtime:,iz1:iz2,0:-en],0),0), '-', color='black', linewidth=1.5) 
#plt.loglog(n[0:-en],np.mean(np.mean(En_rot[nbtime:,iz1:iz2,0:-en],0),0), '-', color='blue', linewidth=1.) 

plt.loglog(n[1:],(1./(2.*L_dim))*Ck*((epsilon_f)**(2./3.))*((n[1:]/(2.*L_dim))**(-5./3.)), '--', color='black', linewidth=1.5)
plt.loglog(n[1:],(1./(2.*L_dim))*Cz*((beta_sat)**2.)*((n[1:]/(2.*L_dim))**(-5.)), '--', color='darkred', linewidth=1.5) 

#_______________________________________________________________________________________________________________________
# ------------ If I want x-axis in km instead of latitudinal indices
#Lambda = (2*np.pi*L_dim)/n[1:]
#plt.loglog((2.*np.pi*L_dim/1000)/n[0:-en],np.mean(np.mean(Ene[nbtime:,iz1:iz2,0:-en],0),0)*(2.*L_dim), '-', color='black', linewidth=1.5) 
#plt.loglog(Lambda/1000.,Ck*((epsilon_f)**(2./3.))*((n[1:]/(2.*L_dim))**(-5./3.)), '--', color='black', linewidth=1.5)

#_______________________________________________________________________________________________________________________
#													 	   Axis:
#_______________________________________________________________________________________________________________________
#plt.ylabel(r'$E_n^R/(R \epsilon_i)^{2/3}$', fontsize=fonthesize)
plt.ylabel(r'$E_n  (J kg^{-1})$', fontsize=fonthesize)
plt.xlabel(r'$n$', fontsize=fonthesize)
plt.xlim((0,361.))
plt.ylim((5*10**(-5),2*10**3))
#plt.ylim((3*10**(-2),2*10**1))
#plt.gca().invert_xaxis()

#_______________________________________________________________________________________________________________________
#													 	 Scales:
#_______________________________________________________________________________________________________________________
# -------------------Typpical scales
#trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
plt.axvline(n_beta, color='black', linestyle=':', linewidth=1.,  alpha=1)
plt.axvline(n_R, color='black', linestyle=':', linewidth=1.,  alpha=1)
#plt.axvline(2.*n_f, color='red', linestyle='solid', linewidth=1.5)
#plt.fill_between(np.arange(n_D,n_D+100,1), 0, 1, facecolor='grey', alpha=0.5, transform=trans)

print 'Energy total'
print np.sum(np.mean(np.mean(Enmo[nbtime:,iz1:iz2,0:-en],0),0))
print np.sum(np.mean(np.mean(Ene[nbtime:,iz1:iz2,0:-en],0),0))
print np.sum(np.mean(np.mean(Enmo[nbtime:,iz1:iz2,0:-en],0),0))+np.sum(np.mean(np.mean(Ene[nbtime:,iz1:iz2,0:-en],0),0))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#								      Barotropic modes:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ETmn_barotrop = np.mean(Ene[nbtime:,:,0:-en] + Enmo[nbtime:,:,0:-en],1)
ETmn = Ene[nbtime:,:,0:-en] + Enmo[nbtime:,:,0:-en]

nbz = ETmn.shape[1]
ETmn_baroclin = np.zeros((ETmn.shape))
for iz in range(1, nbz+1):
	ETmn_baroclin[:,iz-1,:] = abs(ETmn[:,iz-1,:]-ETmn_barotrop[:,:])
print np.mean(np.mean(ETmn_baroclin[nbtime:,:,0:-en],0),0).shape
print np.mean(ETmn_barotrop[nbtime:,0:-en],0).shape

print (np.mean(np.sum(ETmn_barotrop[nbtime:,1:-en],1),0)/np.mean(np.sum(np.mean(ETmn_baroclin[nbtime:,:,1:-en],1),1),0))

figA = plt.subplots()
plt.loglog(n[0:-en],np.mean(np.mean(ETmn_baroclin[nbtime:,:,:],0),0)/((epsilon_f*L_dim)**(2./3.)), '-', color='blue', linewidth=1.5)
plt.loglog(n[0:-en],np.mean(ETmn_barotrop[nbtime:,:],0)/((epsilon_f*L_dim)**(2./3.)), '-', color='black', linewidth=1.5)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#								      	     ENERGY II:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# same than previously but with different representation
en=1
ax,figB = plt.subplots()


#_______________________________________________________________________________________________________________________
# ------------ Here is how Galperin et al. (2014) https://doi.org/10.1016/j.icarus.2013.08.030 plot the spectra
#plt.loglog(n[0:-en],np.mean(np.mean(Ene[nbtime:,iz1:iz2,0:-en],0),0)*np.pi*L_dim, '-', color='black', linewidth=1.5) 
#plt.loglog(n[0:-en],np.mean(np.mean(En_rot[nbtime:,iz1:iz2,0:-en],0),0)*np.pi*L_dim, '-', color='blue', linewidth=1.) 
#plt.loglog(n[1:],Ck*(epsilon_f**(2./3.))*((n[1:]/(2.*np.pi*L_dim))**(-5./3.)), '--', color='black', linewidth=1.5)

#_______________________________________________________________________________________________________________________
# ------------ Here is how it has to be ploted, by normalizing the energy extracted from flow analysis with [E*2*R] = m^3.s^-2, then we have
# 	       an enrgy density. 2*R likely comes from the integral over m, the zonal modes. It leads to the same energetic amplitude than
#	       the one obtained by Young & Read, NatPhys DOI: 10.1038/NPHYS4227, by a totally different approach, using structures function
#	       in cartesian, see details in supplementary materials - it leads for Jupiter Cassini data to a fit with epsilon = 9.10^-5  
#	       m^3.s^-2. Then to be consistent the fit is E_theoric = Ck eps^2/3 (n/2R)^-5/3, [E_theoric] = m^3.s^-2. Or I can do 
#	       E_theoric/2R if I want everything in m^2 s^-2, which is an energy.
plt.loglog(n[0:-en],np.mean(np.mean(Enmo[nbtime:,iz1:iz2,0:-en]+Ene[nbtime:,iz1:iz2,0:-en],0),0), '-', color='darkred', linewidth=1.5) 

plt.loglog(n[1:],(1./(2.*L_dim))*Ck*((epsilon_f)**(2./3.))*((n[1:]/(2.*L_dim))**(-5./3.)), '--', color='black', linewidth=1.5)
plt.loglog(n[1:],(1./(2.*L_dim))*Cz*((beta_sat)**2.)*((n[1:]/(2.*L_dim))**(-5.)), '--', color='darkred', linewidth=1.5) 

#_______________________________________________________________________________________________________________________
# ------------ If I want x-axis in km instead of latitudinal indices
#Lambda = (2*np.pi*L_dim)/n[1:]
#plt.loglog((2.*np.pi*L_dim/1000)/n[0:-en],np.mean(np.mean(Ene[nbtime:,iz1:iz2,0:-en],0),0)*(2.*L_dim), '-', color='black', linewidth=1.5) 
#plt.loglog(Lambda/1000.,Ck*((epsilon_f)**(2./3.))*((n[1:]/(2.*L_dim))**(-5./3.)), '--', color='black', linewidth=1.5)


#_______________________________________________________________________________________________________________________
#													 	   Axis:
#_______________________________________________________________________________________________________________________
#plt.ylabel(r'$E_n^R/(R \epsilon_i)^{2/3}$', fontsize=fonthesize)
plt.ylabel(r'$E_n  (J kg^{-1})$', fontsize=fonthesize)
plt.xlabel(r'$n$', fontsize=fonthesize)
plt.xlim((0,361.))
plt.ylim((1*10**(-4),2*10**3))
#plt.ylim((3*10**(-2),2*10**1))
#plt.gca().invert_xaxis()

#plt.show()
#exit()






#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#								      	    ENERGY Emn:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figA = plt.subplots()
plt.loglog(n[0:-en],np.mean(np.mean(Ene[nbtime:,iz1:iz2,0:-en],0),0)/E_dim, '-', color='black', linewidth=1.5) 
#plt.loglog(n[0:-en],np.mean(np.mean(Ene[:,iz1:iz2,0:-en],0),0)/((np.mean(epsilonf)*L_dim)**(2./3.)), '-', color='black', linewidth=2.) 
plt.loglog(n[0:-en],np.mean(np.mean(Emn[nbtime:,iz1:iz2,0:-en,0],0),0)/E_dim, '-', color='darkred', linewidth=1.) 
plt.loglog(n[0:-en],np.mean(np.mean(Emn[nbtime:,iz1:iz2,0:-en,1],0),0)/E_dim, ':', color='darkgreen', linewidth=1.) 
plt.loglog(n[0:-en],np.mean(np.mean(Emn[nbtime:,iz1:iz2,0:-en,2],0),0)/E_dim, '-', color='darkgreen', linewidth=1.) 
plt.loglog(n[0:-en],np.mean(np.mean(Emn[nbtime:,iz1:iz2,0:-en,3],0),0)/E_dim, ':', color='darkgreen', linewidth=1.) 
plt.loglog(n[0:-en],np.mean(np.mean(Emn[nbtime:,iz1:iz2,0:-en,4],0),0)/E_dim, '-', color='darkgreen', linewidth=1.) 
plt.loglog(n[0:-en],np.mean(np.mean(Emn[nbtime:,iz1:iz2,0:-en,5],0),0)/E_dim, ':', color='darkgreen', linewidth=1.) 
plt.loglog(n[0:-en],np.mean(np.mean(Emn[nbtime:,iz1:iz2,0:-en,6],0),0)/E_dim, '-', color='darkgreen', linewidth=1.) 
plt.loglog(n[0:-en],np.mean(np.mean(Emn[nbtime:,iz1:iz2,0:-en,7],0),0)/E_dim, ':', color='darkgreen', linewidth=1.) 
plt.loglog(n[0:-en],np.mean(np.mean(Emn[nbtime:,iz1:iz2,0:-en,8],0),0)/E_dim, '-', color='darkgreen', linewidth=1.) 
plt.loglog(n[0:-en],np.mean(np.mean(Emn[nbtime:,iz1:iz2,0:-en,9],0),0)/E_dim, ':', color='darkgreen', linewidth=1.) 
plt.loglog(n[0:-en],np.mean(np.mean(Emn[nbtime:,iz1:iz2,0:-en,10],0),0)/E_dim, '-', color='darkgreen', linewidth=1.) 
# -------------------Theoretical spectra
plt.loglog(n[1:],Ck*((epsilon_f/epsilon_dim)**(2./3.))*(n[1:]**(-5./3.)), '--', color='black', linewidth=1.5)
plt.loglog(n[1:],0.5*Ck*((epsilon_f/epsilon_dim)**(2./3.))*(n[1:]**(-8./3.)), '--', color='darkgreen', linewidth=1.5)
plt.loglog(n[1:],Cz*((beta_sat/beta_dim)**2.)*(n[1:]**(-5.)), '--', color='darkred', linewidth=1.5) 
plt.loglog(n[1:],Ck*(1000.*(epsilon_f/epsilon_dim)**(2./3.))*(n[1:]**(-3.)), '--', color='black', linewidth=1.5)

# -------------------Typpical scales
#trans = mtransforms.blended_transform_factory(figB.transData, figB.transAxes)
#plt.axvline(n_beta, color='grey', linestyle='dashed', linewidth=2.5,  alpha=0.7)
#plt.axvline(n_R, color='grey', linestyle='dashed', linewidth=2.5,  alpha=0.7)
#plt.axvline(n_f, color='red', linestyle='solid', linewidth=1.5)
#plt.fill_between(np.arange(n_D,n_D+100,1), 0, 1, facecolor='grey', alpha=0.5, transform=trans)
plt.axvline(n_beta, color='black', linestyle=':', linewidth=1.,  alpha=1)
plt.axvline(n_R, color='black', linestyle=':', linewidth=1.,  alpha=1)
plt.axvline(2.*n_f, color='red', linestyle='solid', linewidth=1.5)
plt.ylabel(r'$E_{mn}^R/(R \epsilon_i)^{2/3}$', fontsize=fonthesize)
plt.xlabel(r'$n$', fontsize=fonthesize)
plt.xlim((0,361.))
#plt.ylim((1*10**(-9),5*10**1))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#								      	  TOTAL ENERGY:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ct=0.95
if(0==0):
	figC = plt.figure()
	plt.plot( (dsteps*SatDay_s)/T_dim,EZ, '.', color='darkred', linewidth=2.5)# step_days[:-1] 
	plt.plot( (dsteps*SatDay_s)/T_dim,ER, '.', color='black', linewidth=2.5)
	plt.plot( (dsteps*SatDay_s)/T_dim,ET, '.', color='navy', linewidth=2.5)
	plt.plot( (dsteps[nbtime:]*SatDay_s)/T_dim,ER[nbtime:], '*', color=(ct,ct,0.2), markersize=0.9)
	plt.yscale('log')
	plt.ylabel('Ek')
	plt.xlabel('Time in days')

	_epsilonf_ = np.mean(epsilonf[nbtime:])
	print '< epsilonf > = ',_epsilonf_
	figC = plt.figure()
	plt.plot( dsteps,epsilonf, '.', color='darkred', linewidth=2.5)# step_days[:-1] 
	plt.axhline(_epsilonf_, color='red', linestyle='solid', linewidth=1.5)
	#plt.yscale('log')
	plt.ylabel('epsilonf')
	plt.xlabel('Time in days')
	plt.title("< epsilonf > ")


plt.show()


