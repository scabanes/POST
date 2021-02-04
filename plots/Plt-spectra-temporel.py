#! /usr/bin/env python
from netCDF4 import Dataset as NetCDFFile 
import numpy as np
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt
import matplotlib
import time
import sys
import LoadInfos
############################################################################################
# 		    PTS: "Parameters to set"
############################################################################################
#nbn=7
im = 1 # Il faut indiquer le mode concerne a titre indicatif
ini = 15 #im
stpn = 2 # step on n
#-------axis
xmin=2e-7
xmax=3e-6
# ..........File to load
#zefile = 'TempModalSpectra-im-3-upto-n-30-stratosphere.nc'
zefile = 'TempModalSpectra-im-1-upto-n-35.nc'
#zefile = 'TempModalSpectra-im-3-upto-n-13-troposphere.nc'
nc = NetCDFFile(zefile) # load file
# ..........LOAD data from StatisticalData.nc
EFT = nc.variables['EFT'][:,:,:] # [n,frq,altitude]
omega_RHW = nc.variables['omega_RHW'][:] # [n]
w = nc.variables['frequences'][:] # [frq]
print EFT.shape
print omega_RHW.shape
COLORS = ['magenta','black','darkblue','darkred','darkgreen','red','chartreuse','darkmagenta','navy']
#............LOAD infos
path=0
E=np.array([0,0])
(n,L_dim,T_dim,SatDay_s,omega_sat,epsilon_dim,epsilon_f,epsilon_i,EnstF_dim,beta_sat,beta_dim,Ck,Cz,E_dim,n_beta,n_R,n_f,n_D) = LoadInfos.LI ( path,E  )

ToOMEGA = omega_sat/40.
print omega_sat
#plt.loglog(w,EFT[ni,:,:], '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
ic = -1
#for ni in range(4,nbn+1):
#	ic = ic + 1
ni=ini#im#im
figA, axs = plt.subplots()
plt.title('n = ')
#.................................................Subplots
ax = plt.subplot(331)
Cst = float(ni)*(float(ni)+1.)/4.
plt.loglog(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
plt.axvline(omega_RHW[ni], color=COLORS[3], linestyle='dashed', linewidth=1.)
#plt.axvline(omega_sat*np.pi*float(im)/(float(ni)**2.), color=COLORS[3], linestyle='solid', linewidth=1.)
#plt.axvline(ToOMEGA, color=COLORS[2], linestyle='solid', linewidth=1.,label='$2 \Omega$')
plt.title('n = '+str(ni))
plt.ylabel(r'$<|\Psi^n_m(\omega)|^2 >$', fontsize=13)
#.................................................Subplots
ax1 = plt.subplot(332); ni = ni+stpn
Cst = float(ni)*(float(ni)+1.)/4.
plt.loglog(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
plt.axvline(omega_RHW[ni], color=COLORS[3], linestyle='dashed', linewidth=1.)
#plt.axvline(omega_sat*np.pi*float(im)/(float(ni)**2.), color=COLORS[3], linestyle='solid', linewidth=1.)
plt.title('n = '+str(ni))
#.................................................Subplots
ax1 = plt.subplot(333); ni = ni+stpn
Cst = float(ni)*(float(ni)+1.)/4.
plt.loglog(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
plt.axvline(omega_RHW[ni], color=COLORS[3], linestyle='dashed', linewidth=1.)
#plt.axvline(omega_sat*np.pi*float(im)/(float(ni)**2.), color=COLORS[3], linestyle='solid', linewidth=1.)
plt.title('n = '+str(ni))
#.................................................Subplots
ax1 = plt.subplot(334); ni = ni+stpn
Cst = float(ni)*(float(ni)+1.)/4.
plt.loglog(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
plt.axvline(omega_RHW[ni], color=COLORS[3], linestyle='dashed', linewidth=1.)
#plt.axvline(omega_sat*np.pi*float(im)/(float(ni)**2.), color=COLORS[3], linestyle='solid', linewidth=1.)
plt.title('n = '+str(ni))
plt.ylabel(r'$<|\Psi^n_m(\omega)|^2 >$', fontsize=13)
#.................................................Subplots
ax1 = plt.subplot(335); ni = ni+stpn
Cst = float(ni)*(float(ni)+1.)/4.
plt.loglog(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
plt.axvline(omega_RHW[ni], color=COLORS[3], linestyle='dashed', linewidth=1.)
#plt.axvline(omega_sat*np.pi*float(im)/(float(ni)**2.), color=COLORS[3], linestyle='solid', linewidth=1.)
plt.title('n = '+str(ni))
#.................................................Subplots
ax1 = plt.subplot(336); ni = ni+stpn
Cst = float(ni)*(float(ni)+1.)/4.
plt.loglog(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
plt.axvline(omega_RHW[ni], color=COLORS[3], linestyle='dashed', linewidth=1.)
#plt.axvline(omega_sat*np.pi*float(im)/(float(ni)**2.), color=COLORS[3], linestyle='solid', linewidth=1.)
plt.title('n = '+str(ni))
#.................................................Subplots
ax1 = plt.subplot(337); ni = ni+stpn
Cst = float(ni)*(float(ni)+1.)/4.
plt.loglog(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
plt.axvline(omega_RHW[ni], color=COLORS[3], linestyle='dashed', linewidth=1.)
#plt.axvline(omega_sat*np.pi*float(im)/(float(ni)**2.), color=COLORS[3], linestyle='solid', linewidth=1.)
plt.title('n = '+str(ni))
plt.ylabel(r'$<|\Psi^n_m(\omega)|^2 >$', fontsize=13)
plt.xlabel(r'$\omega$', fontsize=13)
#.................................................Subplots
ax1 = plt.subplot(338); ni = ni+stpn
Cst = float(ni)*(float(ni)+1.)/4.
plt.loglog(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
plt.axvline(omega_RHW[ni], color=COLORS[3], linestyle='dashed', linewidth=1.)
#plt.axvline(omega_sat*np.pi*float(im)/(float(ni)**2.), color=COLORS[3], linestyle='solid', linewidth=1.)
plt.title('n = '+str(ni))
plt.xlabel(r'$\omega$', fontsize=13)
#.................................................Subplots
ax1 = plt.subplot(339); ni = ni+stpn
Cst = float(ni)*(float(ni)+1.)/4.
plt.loglog(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
plt.axvline(omega_RHW[ni], color=COLORS[3], linestyle='dashed', linewidth=1.)
#plt.axvline(omega_sat*np.pi*float(im)/(float(ni)**2.), color=COLORS[3], linestyle='solid', linewidth=1.)
plt.title('n = '+str(ni))
plt.xlabel(r'$\omega$', fontsize=13)

if(0==0):

	ni=ini#im#im
	figA = plt.subplots()
	#.................................................Subplots
	ax = plt.subplot(331)	
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.semilogx(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[3], linestyle='dashed', linewidth=1.)
	plt.axvline(ToOMEGA, color=COLORS[2], linestyle='dashed', linewidth=1.)
	plt.ylabel(r'$<|\Psi^n_m(\omega)|^2 >$', fontsize=13)
	plt.title('n = '+str(ni))
	axes = plt.gca()
	axes.set_xlim([xmin,xmax])
	#.................................................Subplots
	ax1 = plt.subplot(332); ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.semilogx(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni],color=COLORS[3], linestyle='dashed', linewidth=1.)
	plt.title('n = '+str(ni))
	axes = plt.gca()
	axes.set_xlim([xmin,xmax])
	#.................................................Subplots
	ax1 = plt.subplot(333); ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.semilogx(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[3], linestyle='dashed', linewidth=1.)
	plt.title('n = '+str(ni))
	axes = plt.gca()
	axes.set_xlim([xmin,xmax])
	#.................................................Subplots
	ax1 = plt.subplot(334); ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.semilogx(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[3], linestyle='dashed', linewidth=1.)
	plt.ylabel(r'$<|\Psi^n_m(\omega)|^2 >$', fontsize=13)
	plt.title('n = '+str(ni))
	axes = plt.gca()
	axes.set_xlim([xmin,xmax])
	#.................................................Subplots
	ax1 = plt.subplot(335); ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.semilogx(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[3], linestyle='dashed', linewidth=1.)
	plt.title('n = '+str(ni))
	axes = plt.gca()
	axes.set_xlim([xmin,xmax])
	#.................................................Subplots
	ax1 = plt.subplot(336); ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.semilogx(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[3], linestyle='dashed', linewidth=1.)
	plt.title('n = '+str(ni))
	axes = plt.gca()
	axes.set_xlim([xmin,xmax])
	#.................................................Subplots
	ax1 = plt.subplot(337); ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.semilogx(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[3], linestyle='dashed', linewidth=1.)
	plt.ylabel(r'$<|\Psi^n_m(\omega)|^2 >$', fontsize=13)
	plt.title('n = '+str(ni))
	plt.xlabel(r'$\omega$', fontsize=13)
	axes = plt.gca()
	axes.set_xlim([xmin,xmax])
	#.................................................Subplots
	ax1 = plt.subplot(338); ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.semilogx(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[3], linestyle='dashed', linewidth=1.)
	plt.title('n = '+str(ni))
	plt.xlabel(r'$\omega$', fontsize=13)
	axes = plt.gca()
	axes.set_xlim([xmin,xmax])
	#.................................................Subplots
	ax1 = plt.subplot(339); ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.semilogx(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[3], linestyle='dashed', linewidth=1.)
	plt.title('n = '+str(ni))
	plt.xlabel(r'$\omega$', fontsize=13)
	axes = plt.gca()
	axes.set_xlim([xmin,xmax])

#
if(0==0):

	ni=ini#im#im
	figA = plt.subplots()
	#.................................................Subplots
	ax = plt.subplot(111)
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.semilogx(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[0], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[0], linestyle='dashed', linewidth=1.,label='n = '+str(ni))
	ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.semilogx(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[1], linestyle='dashed', linewidth=1.,label='n = '+str(ni))
	ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.semilogx(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[2], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[2], linestyle='dashed', linewidth=1.,label='n = '+str(ni))
	ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.semilogx(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[3], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[3], linestyle='dashed', linewidth=1.,label='n = '+str(ni))
	ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.semilogx(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[4], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[4], linestyle='dashed', linewidth=1.,label='n = '+str(ni))
	ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.semilogx(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[5], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[5], linestyle='dashed', linewidth=1.,label='n = '+str(ni))
	ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.semilogx(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[6], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[6], linestyle='dashed', linewidth=1.,label='n = '+str(ni))
	ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.semilogx(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[7], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[7], linestyle='dashed', linewidth=1.,label='n = '+str(ni))
	ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.semilogx(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[8], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[8], linestyle='dashed', linewidth=1.,label='n = '+str(ni))
#	plt.axvline(omega_RHW[ni], color=COLORS[1], linestyle='dashed', linewidth=1.)
	#plt.axvline(ToOMEGA, color=COLORS[2], linestyle='solid', linewidth=1.,label='$2 \Omega$')
	plt.title('m = '+str(im))
	ax.legend()
	plt.ylabel(r'$<|\Psi^n_m(\omega)|^2 >$', fontsize=13)
	plt.xlabel(r'$\omega$', fontsize=13)
	#axes = plt.gca()
	#axes.set_xlim([xmin,xmax])

if(0==0):

	ni=ini#im#im
	figA = plt.subplots()
	#.................................................Subplots
	ax = plt.subplot(111)
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.loglog(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[0], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[0], linestyle='dashed', linewidth=1.,label='n = '+str(ni))
	ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.loglog(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[1], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[1], linestyle='dashed', linewidth=1.,label='n = '+str(ni))
	ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.loglog(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[2], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[2], linestyle='dashed', linewidth=1.,label='n = '+str(ni))
	ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.loglog(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[3], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[3], linestyle='dashed', linewidth=1.,label='n = '+str(ni))
	ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.loglog(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[4], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[4], linestyle='dashed', linewidth=1.,label='n = '+str(ni))
	ni = ni+stpn
	Cst = float(ni)*(float(ni)+1.)/4.
	plt.loglog(w,Cst*np.mean(EFT[ni,:,:],1), '-', color=COLORS[5], linewidth=1.)#, alpha=ALPHAS	[il])
	plt.axvline(omega_RHW[ni], color=COLORS[5], linestyle='dashed', linewidth=1.,label='n = '+str(ni))
#	plt.axvline(omega_RHW[ni], color=COLORS[1], linestyle='dashed', linewidth=1.)
	#plt.axvline(ToOMEGA, color=COLORS[2], linestyle='solid', linewidth=1.,label='$2 \Omega$')
	plt.title('m = '+str(im))
	ax.legend()
	plt.ylabel(r'$<|\Psi^n_m(\omega)|^2 >$', fontsize=13)
	plt.xlabel(r'$\omega$', fontsize=13)
	#axes = plt.gca()
	#axes.set_xlim([xmin,xmax])
plt.show(figA)	
