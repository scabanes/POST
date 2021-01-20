#! /usr/bin/env python
#
def WindowingHanning ( Y,N,L,iP ):
# INPUT: N la taille du signal, Y le signal de N echantillon en temps ou espace.
#	 L longueur du signal dans  l unite souhaitee.
# OUTPUT: Y_Han de taille N 

# % % Y*H(n) % on va perdre ici en amplitude du signal, si AMP=1 dans H
# ==> l'apmplitude est comprise entre 0 (au bord) et 1.
# avec Hanning on va diminuer l amplitude du signal par 2, integral de la
# courbe une fois passee dans la fenetre de hanning en revanche le pic en
# frq est plus marque, le signal etant periodisise.
# ==> En revanche les basses frquencs sont plus energetiques (hanning ajoute
# artificiellement un signal basse frequence) mais les hautes frqce sont
# mieux resolues.
# ==> si le signal est deja periodique on gagne rien avec hanning on perd
# au contraire les resolution basse frq en revanche si non periodique on
# perd bcp dans une fft classique
	import numpy as np
	import math
	import matplotlib.pyplot as plt
	import matplotlib.patches as mpatches

# ==> hanning function
	H=0.5*(1.-np.cos(2.*np.pi*np.arange(0.,N,1.)/(N-1.)))
	r=np.arange(0.,N,1.)/(N-1.)
# ==> windowed signal
	Y_Han = Y*np.transpose(H);
	if(iP==1):
		fig = plt.figure()	
		plt.plot(r*L,Y_Han, color='black')
		plt.plot(r*L,Y, color='darkred')
		plt.xlabel('space/time variable')
		plt.ylabel('Y_windowed')
		plt.grid()
		plt.show(fig)
	return (Y_Han)


def WindowingTukey ( Y,N,L,iP ):
# INPUT: N la taille du signal, Y le signal de N echantillon en temps ou espace.
#	 L longueur du signal dans  l unite souhaitee.
# OUTPUT: Y_Han de taille N 

# % % Y*H(n) % on va perdre ici en amplitude du signal, si AMP=1 dans H
# ==> l'apmplitude est comprise entre 0 (au bord) et 1.
# avec Hanning on va diminuer l amplitude du signal par 2, integral de la
# courbe une fois passee dans la fenetre de hanning en revanche le pic en
# frq est plus marque, le signal etant periodisise.
# ==> En revanche les basses frquencs sont plus energetiques (hanning ajoute
# artificiellement un signal basse frequence) mais les hautes frqce sont
# mieux resolues.
# ==> si le signal est deja periodique on gagne rien avec hanning on perd
# au contraire les resolution basse frq en revanche si non periodique on
# perd bcp dans une fft classique
	import numpy as np
	import math
	import matplotlib.pyplot as plt
	import matplotlib.patches as mpatches

# ==> Tukey function
	T = np.zeros([N])
	alpha = 0.5
	n = np.arange(0.,N,1.)
	n1 = alpha*(N-1)/2.
	n2 = (N-1)*(1-alpha/2.)
	for iN in range(0, len(n)):
		if(n[iN] < n1):
			T[iN] = 0.5*(1.+np.cos(np.pi*(2*n[iN]/(alpha*(N-1.)) - 1.)))
		elif(n[iN] <= n2):
			T[iN] = 1.
		else:
			T[iN] = 0.5*(1.+np.cos(np.pi*(2*n[iN]/(alpha*(N-1.)) - 2./alpha + 1.)))
	r=np.arange(0.,N,1.)/(N-1.)
# ==> windowed signal
	Y_Tukey = Y*np.transpose(T)
	if(iP==1):
		fig = plt.figure()	
		plt.plot(r*L,Y_Tukey, color='black')
		plt.plot(r*L,Y, color='darkred')
		plt.xlabel('space/time variable')
		plt.ylabel('Y_windowed')
		plt.grid()
		plt.show(fig)
	return (Y_Tukey)
