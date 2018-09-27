#! /usr/bin/env python
from netCDF4 import Dataset as NetCDFFile 
import numpy as np
import math
import scipy.special as sp
from pylab import *
###########################################################################################################################
#													FUNCTION 1
###########################################################################################################################
def LI ( path ):
####################################################################
# 		    LOAD filePTS ....
####################################################################
	#################################
	if path ==0:
		f = open('filePTS.zono.temp', "r")
	else:
		f = open(path+'filePTS.zono.temp', "r")
	print f.readline()
	print f.readline()
	print f.readline()
	f.readline()
	# ----------------------
	getField= (f.readline().split()[1:])
	# ----------------------
	# Arrays to step in time for averaging. Each array has a given AST:
	# "Average Step Time" which give the sleeping for averaging. One has
	# to concatenate all the arrays in stept_file.
	iAT11,iAT12 = f.readline().split()[1:3]
	iAT21,iAT22 = f.readline().split()[1:3]
	iAT31,iAT32 = f.readline().split()[1:3]
	iAT41,iAT42 = f.readline().split()[1:3]
	#
	arrayT1 = np.arange(int(iAT11),int(iAT12),1)
	arrayT2 = np.arange(int(iAT21),int(iAT22),1)
	arrayT3 = np.arange(int(iAT31),int(iAT32),1)
	arrayT4 = np.arange(int(iAT41),int(iAT42),1)
	#
	AST1 = int(f.readline().split()[1])
	AST2 = int(f.readline().split()[1])
	AST3 = int(f.readline().split()[1])
	AST4 = int(f.readline().split()[1])
	#
	dim = int(f.readline().split()[1])
	#
	print f.readline()
	print f.readline()
	print f.readline()
	omega_sat = float(f.readline().split()[1])
	R_sat = float(f.readline().split()[1])
	SatDay_s = float(f.readline().split()[1])
	Ck = float(f.readline().split()[1])
	epsilon = float(f.readline().split()[1])
	Cz = float(f.readline().split()[1])
	# Brunt Vaisala frequency [NBV] and atmosphere deepness [H]
	H = float(f.readline().split()[1])
	NBV = float(f.readline().split()[1])
	#
	epsilon_f = float(f.readline().split()[1])
	n_f = float(f.readline().split()[1])
	tau_f = float(f.readline().split()[1])
	#
	if(AST1!=0) and (AST2==0):
        	stept_file = arrayT1 - int(iAT11)
        	AST4=AST1
	elif(AST1!=0) and (AST2!=0) and (AST3==0):
        	stept_file = np.concatenate((arrayT1 - iAT11, arrayT2 - iAT21),axis=0)
        	AST4=AST2
	elif(AST1!=0) and (AST2!=0) and (AST3!=0) and (AST4==0):
	        stept_file = np.concatenate((arrayT1 - iAT11, arrayT2 - iAT21, arrayT3 - iAT31),axis=0)
	        AST4=AST3
	elif(AST1!=0) and (AST2!=0) and (AST3!=0) and (AST4!=0):
	        stept_file = np.concatenate((arrayT1 - iAT11, arrayT2 - iAT21, arrayT3 - iAT31, arrayT4 - iAT41),axis=0)
	else:
	        print 'erreur dans filePTS.zono.temp'
	print 'Loaded.....'
	print stept_file
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#===================================================================
	# 			Dimensions
	#===================================================================
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	angle = np.pi/4.
	beta_sat = (2.*omega_sat*np.cos(angle))/R_sat
	ff = 2.*omega_sat*np.sin(angle)
	if dim == 1:
	# Dimensionalisation 1..........
		#k_f = n_f/R_sat# kf = pi/Lf  et Lf = pi*R/nf en m-1
		L_dim = R_sat #rayon en m
		T_dim = (epsilon_f/(R_sat**2.))**(-1./3.) # temp en s
		#
		E_dim = (epsilon_f*R_sat)**(2./3.)
		epsilon_dim = epsilon_f # m^2.s-3
		beta_dim = 1./(L_dim*T_dim) # m-1.s-1
		#
		EnstF_dim = epsilon_f/(R_sat**2.) #s^-3
		#
		#E_dim = (epsilon_f/k_f)**(2./3.)
		#beta_dim = (epsilon_f/E_dim)*k_f # m-1.s-1
		#t_dim = (k_f**2.*epsilon_f)**(-1./3.)
	elif dim == 2:
	# Dimensionalisation 2..........
	# To make quantities non-dimensional we use [omega] = 1/s as unity  
	# for time and [R_sat] = m as length unity.
		L_dim = R_sat #rayon en m
		T_dim = 1./omega_sat # temp en s
		#
		E_dim = ((omega_sat*R_sat)**2.)
		epsilon_dim = (R_sat**2.)*(omega_sat**3.)
		beta_dim = beta_sat
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#===================================================================
	# 			Scaling
	#===================================================================
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	# ---------------------------------------------------------------------
	#	       	Numerical scaling results:
	# Indexes n_R is the Rhines scale; n_beta is the transition scale in
	# the zonostrophic spectra; n_z is the scale at which the flow beco-
	# mes isotropic; n_D is the first Rossby deformation radius. R_beta
	# is the zonostrophic index giving the extension on the zonostrophic
	# inertial range.
	# ---------------------------------------------------------------------	
	NbN = 361
	n = np.arange(0.,NbN,1.)
	URMS = 29.
	print '================================================================'
	print '-------------------- SCALES ------------------------------------'
	n_R = R_sat * (beta_sat/(2*URMS))**0.5 # rad/m
	n_beta = ((Cz/Ck)**(3./10.))*((((beta_sat/beta_dim)**3.)/(epsilon_f/epsilon_dim))**(1./5.))
	n_z = n_beta*(2.*n_beta)**(3./7.)
	L_D = (NBV*H)/(np.pi*ff)
	n_D = (np.pi*R_sat)/L_D
	n_Ro = ((omega_sat*T_dim)/((epsilon/epsilon_dim)**(1./3.)))**(3./2.)
	#n_Ro = ((omega_sat/omega_sat)/((epsilon_f/(epsilon_dim))**(1./3.)))**(3./2.)
	n_S = 2.*np.pi*R_sat/H
	R_beta = n_beta/n_R
	print 'URMS =',URMS
	print 'n_beta =', n_beta 
	print 'n_R = ', n_R
	print 'n_z =', n_z
	print 'n_D =', n_D
	print 'R_beta = ', R_beta
	return (n,epsilon_dim,epsilon_f,beta_sat,beta_dim,Ck,Cz,E_dim,n_beta,n_R,n_f,n_D)		
