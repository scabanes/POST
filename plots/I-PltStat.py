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
hn_zm = nc.variables['hn_zm'][:] # [t,Prof,theta,phi] Zonal velocity
fn_zm = nc.variables['fn_zm'][:] # [t,Prof,theta,phi] Zonal velocity
Enmo = nc.variables['Enmo'][:] # [t,Prof,theta,phi] Zonal velocity
Ene = nc.variables['Ene'][:] # [t,Prof,theta,phi] Zonal velocity
hn = nc.variables['hn'][:] # [t,Prof,theta,phi] Zonal velocity
fn = nc.variables['fn'][:] # [t,Prof,theta,phi] Zonal velocity
n = (np.arange(0.,360,1.))
#
itmax  = fn_zm.shape[0]
izmax  = fn.shape[1]
#............LOAD infos
(n,epsilon_dim,epsilon_f,beta_sat,beta_dim,Ck,Cz,E_dim,n_beta,n_R,n_f,n_D) = LoadInfos.LI (  )
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#									 ENERGY FLUXES:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (0==0):
	figA = plt.subplots()
	ax = plt.subplot(321)
	for iz in range(0, 3):
		plt.semilogx(np.mean(fn[:,iz,:],0), '-', color='black', linewidth=1) 
	plt.title('Energy mean t & z = 0,1,2')
	plt.grid()

	ax1 = plt.subplot(322)
	for iz in range(3, 6):
		plt.semilogx(np.mean(fn[:,iz,:],0), '-', color='black', linewidth=1) 
	plt.title('Energy mean t & z = 8,9,10')
	plt.grid()

	ax2 = plt.subplot(323)
	for iz in range(6, 9):
		plt.semilogx(np.mean(fn[:,iz,:],0), '-', color='black', linewidth=1) 
	plt.title('Energy mean t & z = 18,19,20')
	plt.grid()

	ax3 = plt.subplot(324)
	for iz in range(9, 12):
		plt.semilogx(np.mean(fn[:,iz,:],0), '-', color='black', linewidth=1) 
	plt.title('Energy mean t & z = 29,30,31')
	plt.grid()

	ax4 = plt.subplot(325)
	plt.semilogx(np.mean(np.mean(fn[:,:,:],0),0), '-', color='black', linewidth=1.5) 
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
#								      ENSTROPHY FLUXES:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (0==0):
	figA = plt.subplots()
	ax = plt.subplot(321)
	for iz in range(0, 3):
		plt.semilogx(np.mean(hn[:,iz,:],0), '-', color='darkblue', linewidth=1) 
	plt.title('Energy mean t & z = 0,1,2')
	plt.grid()

	ax1 = plt.subplot(322)
	for iz in range(3, 6):
		plt.semilogx(np.mean(hn[:,iz,:],0), '-', color='darkblue', linewidth=1) 
	plt.title('Energy mean t & z = 8,9,10')
	plt.grid()

	ax2 = plt.subplot(323)
	for iz in range(6, 9):
		plt.semilogx(np.mean(hn[:,iz,:],0), '-', color='darkblue', linewidth=1) 
	plt.title('Energy mean t & z = 18,19,20')
	plt.grid()

	ax3 = plt.subplot(324)
	for iz in range(9, 12):
		plt.semilogx(np.mean(hn[:,iz,:],0), '-', color='darkblue', linewidth=1) 
	plt.title('Energy mean t & z = 29,30,31')
	plt.grid()

	ax4 = plt.subplot(325)
	plt.semilogx(np.mean(np.mean(hn[:,:,:],0),0), '-', color='darkblue', linewidth=1.5) 
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
for iz in range(0, 3):
	ax.loglog(n[0:-en],np.mean(Enmo[:,iz,0:-en],0)/E_dim, '-', color='darkred', linewidth=1) 
	ax.loglog(n[0:-en],np.mean(Ene[:,iz,0:-en],0)/E_dim, '-', color='black', linewidth=1)
# -------------------Theoretical spectra
ax.loglog(n[1:],Ck*((epsilon_f/epsilon_dim)**(2./3.))*(n[1:]**(-5./3.)), '--', color='black', linewidth=2.5)
ax.loglog(n[1:],Cz*((beta_sat/beta_dim)**2.)*(n[1:]**(-5.)), '--', color='darkred', linewidth=2.5) 
#----------------------------------------
plt.title('Energy mean t & z = 0,1,2')
ax.grid()

ax1 = plt.subplot(322)
for iz in range(3, 6):
	ax1.loglog(n[0:-en],np.mean(Enmo[:,iz,0:-en],0)/E_dim, '-', color='darkred', linewidth=1) 
	ax1.loglog(n[0:-en],np.mean(Ene[:,iz,0:-en],0)/E_dim, '-', color='black', linewidth=1) 
# -------------------Theoretical spectra
ax1.loglog(n[1:],Ck*((epsilon_f/epsilon_dim)**(2./3.))*(n[1:]**(-5./3.)), '--', color='black', linewidth=2.5)
ax1.loglog(n[1:],Cz*((beta_sat/beta_dim)**2.)*(n[1:]**(-5.)), '--', color='darkred', linewidth=2.5) 
#----------------------------------------
plt.title('Energy mean t & z = 8,9,10')
ax1.grid()

ax3 = plt.subplot(323)
for iz in range(6, 9):
	ax3.loglog(n,np.mean(Enmo[:,iz,:],0)/E_dim, '-', color='darkred', linewidth=1) 
	ax3.loglog(n,np.mean(Ene[:,iz,:],0)/E_dim, '-', color='black', linewidth=1) 
# -------------------Theoretical spectra
ax3.loglog(n[1:],Ck*((epsilon_f/epsilon_dim)**(2./3.))*(n[1:]**(-5./3.)), '--', color='black', linewidth=2.5)
ax3.loglog(n[1:],Cz*((beta_sat/beta_dim)**2.)*(n[1:]**(-5.)), '--', color='darkred', linewidth=2.5) 
#----------------------------------------
plt.title('Energy mean t & z = 18,19,20')
ax3.grid()

ax4 = plt.subplot(324)
for iz in range(9, 12):
	ax4.loglog(n,np.mean(Enmo[:,iz,:],0)/E_dim, '-', color='darkred', linewidth=1) 
	ax4.loglog(n,np.mean(Ene[:,iz,:],0)/E_dim, '-', color='black', linewidth=1) 
# -------------------Theoretical spectra
ax4.loglog(n[1:],Ck*((epsilon_f/epsilon_dim)**(2./3.))*(n[1:]**(-5./3.)), '--', color='black', linewidth=2.5)
ax4.loglog(n[1:],Cz*((beta_sat/beta_dim)**2.)*(n[1:]**(-5.)), '--', color='darkred', linewidth=2.5) 
#----------------------------------------
plt.title('Energy mean t & z = 29,30,31')
ax4.grid()

ax5 = plt.subplot(325)
ax5.loglog(n[0:-en],np.mean(np.mean(Enmo[:,0:-2,0:-en],0),0)/E_dim, '-', color='darkred', linewidth=1.5) 
ax5.loglog(n[0:-en],np.mean(np.mean(Ene[:,0:-2,0:-en],0),0)/E_dim, '-', color='black', linewidth=1.5) 
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
plt.show()
exit()




plt.show()


