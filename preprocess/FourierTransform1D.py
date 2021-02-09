#! /usr/bin/env python
#
def fft1D ( Y, N ):

	import numpy as np
	import math
	import scipy.special as sp
	#import SH
	#import Hanning

#FourierTransform un signal Y de taille N pour lequel on fait une
#transformee de fourier discrete.
# FFT classic
# fft(1)=Reel est la cste +1 du signal y=Amp*sin(x) +1; 
# fft(2:N/2)=imag de la transformee de fourier.
# fft(N/2+1)=Reel est la frq de coupure de lunquist pour un echantillonage N.
# INPUT: le signal Y et sa taille N
# OUTPUT: le signal FFT est de taille N/2+1
# si not un integer alors les indices pour stocker la matrice sont non entier plus bas
	if float(N/2.).is_integer() == False:
		N=N-1
#--------------------------------------------------------------------------
#                             NOTES SUR PSD & FFT
#--------------------------------------------------------------------------
# Note: from ==> http://www.gaussianwaves.com/2013/12/computation-of-power-of-a-signal-in-matlab-simulation-and-verification/
# Note se rapportant a l image EnergySpectrum.jpg au lien C:\Users\simon\Desktop\Num_Jupiter\Nek5000\Routines
# Une fft doit etre normalisee par N dans matlab, ceci n est pas fait
# automatiquement avec la routine fft de matlab. En revanche la ifft
# normalise par N. Conclusion abs(FFT_Y/N) donne une amplitude des modes
# Amp/2 (voir image EnergySpectrum.jpg). Donc je peux multiplier par 2 voir**.
# De meme pour la power spectral density PSD = (1/(N*N)) * abs(xdft).^2 =
# Amp^2/4. Donc je peux multiplier par 4 pour connaitre l amplitude au
# carre voir***: LA BELLE AFFAIRE!!!
#--------------------------------------------------------------------------
#                             NOTES SUR FFT
#--------------------------------------------------------------------------
# Dans FT_Y et FT_Y2 = abs(FT_Y/N) on a: 1/ FT_Y(1) est reel et est la constante du
# signal (non nulle si le signal est de moyenne non nulle) soit le mode
# m=0. 2/ FT_Y(2:N/2) contient tous les modes aux nombre de N/2-1 (m=0 etant a part). N/2 etant la frequence
# max de lunquist contraint par la taille initiale de l echantillonnage
# dans le signal. 3/ FT_Y(N/2+1) est un reel egalement corrspondant a la
# frequence de Lunquist. 4/ FT_Y(N/2+2:end) on retrouve les N/2 modes qui
# sont les conjugues de FT_Y(2:N/2). On a un signal FT_Y de taille L = 2 *
# (N/2-1) + 2.
# La FFT est definit telle que F(k) = Integral f(x) exp[-2ipi x k] dx
# La IFFT est definit comme    f(x) = Integral F(k) exp[ 2ipi x k] dk
# -----Tranformer de fourier par FFT
	FFT_Y = np.fft.fft(Y) # D
	FFT_Y2 = abs(FFT_Y/N)
	FFT_Y1 = FFT_Y2[0:N/2+1] # 
	FFT_Y1[1:-2] = 2*FFT_Y1[1:-2] # asterix **

# -----energy spectra
	xdft = np.fft.fft(Y)
	xdft = xdft[0:N/2+1]
	psdx = (1./(N*N)) * abs(xdft)**2
	psdx[1:-2] = 4.*psdx[1:-2]  # Asterix *** , autre exemple: ==>  http://fr.mathworks.com/help/signal/ug/psd-estimate-using-fft.html
	E_Y = psdx
	m = np.arange(0.,N/2+1) # Attention!!! la gamme des modes demarre a 0 et pas 1, lunquist inclu
	return( FFT_Y1,E_Y,m )
