#! /usr/bin/env python
from netCDF4 import Dataset as NetCDFFile 
import numpy as np
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt
import matplotlib
import time
import sys
import LoadInfos
# ..........File to load
zefile = 'StatisticalData.nc'
nc = NetCDFFile(zefile) # load file
# ..........LOAD data from StatisticalData.nc
hn_zm = nc.variables['hn_zm'][:] # [time,latitude]
fn_zm = nc.variables['fn_zm'][:] # [time,latitude]
Enmo = nc.variables['Enmo'][:] # [time,altitude,latitude]
Ene = nc.variables['Ene'][:] # [time,altitude,latitude]
hn = nc.variables['hn'][:] # [time,altitude,latitude]
fn = nc.variables['fn'][:] #[time,altitude,latitude]
ET = nc.variables['ET'][:] #[time]
ER = nc.variables['ER'][:] #[time]
EZ = nc.variables['EZ'][:] #[time]
dsteps = nc.variables['dsteps'][:] #[time]
epsilonf = nc.variables['epsilonf'][:] #[time]

n = (np.arange(0.,360,1.))
#
itmax  = fn_zm.shape[0]
izmax  = fn.shape[1]
#............LOAD infos
path=0
(n,epsilon_dim,epsilon_f,beta_sat,beta_dim,Ck,Cz,E_dim,n_beta,n_R,n_f,n_D) = LoadInfos.LI ( path  )
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#									 ENERGY FLUXES:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (0==0):
	figA = plt.subplots()
	ax = plt.subplot(321)
	for it in range(0, itmax):
		plt.semilogx(np.mean(fn[it,0:3,:],0)/epsilon_dim, '-', color='black', linewidth=1) 
	plt.title('Energy mean z = 0,1,2 each t')
	plt.grid()

	ax1 = plt.subplot(322)
	for it in range(0,itmax):
		plt.semilogx(np.mean(fn[it,3:6,:],0)/epsilon_dim, '-', color='black', linewidth=1) 
	plt.title('Energy mean  z = 8,9,10 each t')
	plt.grid()

	ax2 = plt.subplot(323)
	for it in range(0, itmax):
		plt.semilogx(np.mean(fn[it,6:9,:],0)/epsilon_dim, '-', color='black', linewidth=1) 
	plt.title('Energy mean  z = 18,19,20 each t')
	plt.grid()

	ax3 = plt.subplot(324)
	for it in range(0, itmax):
		plt.semilogx(np.mean(fn[it,9:11,:],0)/epsilon_dim, '-', color='black', linewidth=1) 
	plt.title('Energy mean z = 29,30,31 each t')
	plt.grid()

	ax4 = plt.subplot(325)
	plt.semilogx(np.mean(np.mean(fn[:,:,:],0),0)/epsilon_dim, '-', color='black', linewidth=1.5) 
	plt.title('Energy mean t & z')
	plt.grid()

# -------------------Typpical scales
ax.axvline(n_R, color='grey', linestyle='dashed', linewidth=2.5,  alpha=0.7)
ax.axvline(n_f, color='red', linestyle='solid', linewidth=1.5)
ax1.axvline(n_R, color='grey', linestyle='dashed', linewidth=2.5,  alpha=0.7)
ax1.axvline(n_f, color='red', linestyle='solid', linewidth=1.5)
ax2.axvline(n_R, color='grey', linestyle='dashed', linewidth=2.5,  alpha=0.7)
ax2.axvline(n_f, color='red', linestyle='solid', linewidth=1.5)
ax3.axvline(n_R, color='grey', linestyle='dashed', linewidth=2.5,  alpha=0.7)
ax3.axvline(n_f, color='red', linestyle='solid', linewidth=1.5)
ax4.axvline(n_R, color='grey', linestyle='dashed', linewidth=2.5,  alpha=0.7)
ax4.axvline(n_f, color='red', linestyle='solid', linewidth=1.5)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#								      	        ENERGY:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
en=1
figB = plt.subplots()
ax = plt.subplot(321)
for it in range(0,itmax):
	#ax.loglog(n[0:-en],np.mean(Enmo[it,0:3,0:-en],1)/E_dim, '-', color='darkred', linewidth=1) 
	ax.loglog(n[0:-en],np.mean(Ene[it,0:3,0:-en],0)/E_dim, '-', color='black', linewidth=1)
# -------------------Theoretical spectra
ax.loglog(n[1:],Ck*((epsilon_f/epsilon_dim)**(2./3.))*(n[1:]**(-5./3.)), '--', color='black', linewidth=2.5)
#ax.loglog(n[1:],Cz*((beta_sat/beta_dim)**2.)*(n[1:]**(-5.)), '--', color='darkred', linewidth=2.5) 
#----------------------------------------
plt.title('Energy mean z = 0,1,2  each t')
ax.grid()

ax1 = plt.subplot(322)
for it in range(0,itmax):
	#ax1.loglog(n[0:-en],np.mean(Enmo[it,0:3,0:-en],1)/E_dim, '-', color='darkred', linewidth=1) 
	ax1.loglog(n[0:-en],np.mean(Ene[it,3:6,0:-en],0)/E_dim, '-', color='black', linewidth=1) 
# -------------------Theoretical spectra
ax1.loglog(n[1:],Ck*((epsilon_f/epsilon_dim)**(2./3.))*(n[1:]**(-5./3.)), '--', color='black', linewidth=2.5)
#ax1.loglog(n[1:],Cz*((beta_sat/beta_dim)**2.)*(n[1:]**(-5.)), '--', color='darkred', linewidth=2.5) 
#----------------------------------------
plt.title('Energy mean z = 8,9,10 each t')
ax1.grid()

ax3 = plt.subplot(323)
for it in range(0,itmax):
	#ax3.loglog(n,np.mean(Enmo[it,0:3,:],1)/E_dim, '-', color='darkred', linewidth=1) 
	ax3.loglog(n,np.mean(Ene[it,6:9,:],0)/E_dim, '-', color='black', linewidth=1) 
# -------------------Theoretical spectra
ax3.loglog(n[1:],Ck*((epsilon_f/epsilon_dim)**(2./3.))*(n[1:]**(-5./3.)), '--', color='black', linewidth=2.5)
#ax3.loglog(n[1:],Cz*((beta_sat/beta_dim)**2.)*(n[1:]**(-5.)), '--', color='darkred', linewidth=2.5) 
#----------------------------------------
plt.title('Energy mean z = 18,19,20 each t')
ax3.grid()

ax4 = plt.subplot(324)
for it in range(0,itmax):
	#ax4.loglog(n,np.mean(Enmo[it,0:3,:],1)/E_dim, '-', color='darkred', linewidth=1) 
	ax4.loglog(n,np.mean(Ene[it,9:12,:],0)/E_dim, '-', color='black', linewidth=1) 
# -------------------Theoretical spectra
ax4.loglog(n[1:],Ck*((epsilon_f/epsilon_dim)**(2./3.))*(n[1:]**(-5./3.)), '--', color='black', linewidth=2.5)
#ax4.loglog(n[1:],Cz*((beta_sat/beta_dim)**2.)*(n[1:]**(-5.)), '--', color='darkred', linewidth=2.5) 
#----------------------------------------
plt.title('Energy mean z = 29,30,31 each t')
ax4.grid()

ax5 = plt.subplot(325)
ax5.loglog(n[0:-en],np.mean(np.mean(Enmo[:,0:-3,0:-en],0),0)/E_dim, '-', color='darkred', linewidth=1.5) 
ax5.loglog(n[0:-en],np.mean(np.mean(Ene[:,0:-3,0:-en],0),0)/E_dim, '-', color='black', linewidth=1.5) 
# -------------------Theoretical spectra
ax5.loglog(n[1:],Ck*((epsilon_f/epsilon_dim)**(2./3.))*(n[1:]**(-5./3.)), '--', color='black', linewidth=2.5)
ax5.loglog(n[1:],Cz*((beta_sat/beta_dim)**2.)*(n[1:]**(-5.)), '--', color='darkred', linewidth=2.5) 
#----------------------------------------
plt.title('Energy mean t & z')
ax5.grid()


# -------------------Typpical scales
trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
ax.axvline(n_beta, color='grey', linestyle='dashed', linewidth=2.5,  alpha=0.7)
ax.axvline(n_R, color='grey', linestyle='dashed', linewidth=2.5,  alpha=0.7)
ax.axvline(n_f, color='red', linestyle='solid', linewidth=1.5)
ax.fill_between(np.arange(n_D,n_D+100,1), 0, 1, facecolor='grey', alpha=0.5, transform=trans)
trans = mtransforms.blended_transform_factory(ax1.transData, ax1.transAxes)
ax1.axvline(n_beta, color='grey', linestyle='dashed', linewidth=2.5,  alpha=0.7)
ax1.axvline(n_R, color='grey', linestyle='dashed', linewidth=2.5,  alpha=0.7)
ax1.axvline(n_f, color='red', linestyle='solid', linewidth=1.5)
ax1.fill_between(np.arange(n_D,n_D+100,1), 0, 1, facecolor='grey', alpha=0.5, transform=trans)
trans = mtransforms.blended_transform_factory(ax3.transData, ax3.transAxes)
ax3.axvline(n_beta, color='grey', linestyle='dashed', linewidth=2.5,  alpha=0.7)
ax3.axvline(n_R, color='grey', linestyle='dashed', linewidth=2.5,  alpha=0.7)
ax3.axvline(n_f, color='red', linestyle='solid', linewidth=1.5)
ax3.fill_between(np.arange(n_D,n_D+100,1), 0, 1, facecolor='grey', alpha=0.5, transform=trans)
trans = mtransforms.blended_transform_factory(ax4.transData, ax4.transAxes)
ax4.axvline(n_beta, color='grey', linestyle='dashed', linewidth=2.5,  alpha=0.7)
ax4.axvline(n_R, color='grey', linestyle='dashed', linewidth=2.5,  alpha=0.7)
ax4.axvline(n_f, color='red', linestyle='solid', linewidth=1.5)
ax4.fill_between(np.arange(n_D,n_D+100,1), 0, 1, facecolor='grey', alpha=0.5, transform=trans)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#								      	  TOTAL ENERGY:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0==0):
	figC = plt.figure()
	plt.plot( dsteps,EZ, '.', color='darkred', linewidth=2.5)# step_days[:-1] 
	plt.plot( dsteps,ER, '.', color='black', linewidth=2.5)
	plt.plot( dsteps,ET, '.', color='navy', linewidth=2.5)
	plt.yscale('log')
	plt.ylabel('Ek')
	plt.xlabel('Time in days')

	_epsilonf_ = np.mean(epsilonf[10:])
	print '< epsilonf > = ',_epsilonf_
	figC = plt.figure()
	plt.plot( dsteps,epsilonf, '.', color='darkred', linewidth=2.5)# step_days[:-1] 
	plt.axhline(_epsilonf_, color='red', linestyle='solid', linewidth=1.5)
	plt.yscale('log')
	plt.ylabel('epsilonf')
	plt.xlabel('Time in days')
	plt.title("< epsilonf > ")
plt.show()
