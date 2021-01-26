# POST
Statistical tools to study geophysical flows

```
cd spherepack3.2
make
cd ..
make statistical_analysis_xyzt
```

Steps
1. Preprocess Xhistins (python in preprocess)
2. statistical_analysis_xyzt uvData-istep-0-nstep-40-niz-12.nc
3. e.g. Tevo-PltStat.py StatisticalData.nc

On a:

## A)------------------------------------FOR GENERAL STATISTICAL ANALYSIS

ZUVT_pretraitement_HV.py (->uvData-istep-x-nstep-xxx-niz-xx.nc)
statistical_analysis_xyzt.f90 (->StatisticalData.nc) --> routine de travail de base.
On peut exploiter ces donnees avec les
routines:
+ Tevo-PltStat.py         --> time evolution of the data/spectra/fluxes
+ I-PltStat.py            --> Mean time and each levels
+ comparaison-v2.py	  --> comparaison de toutes les simus
+ comparaison-rotrate.py  --> comparaison des simus ayant differents rotation rates
+ comparaison-kf.py	  --> comparaison des simus ayant different kf de forcage
background routines
+ LoadInfos.py		  --> Load les data from filePTS.zono.temp
+ filePTS.zono.temp	  --> Enter data of the specific simulation

## B)------------------------------------FOR DISSIPATION SPECTRAL STUDIES
ZUVT_pretraitement_dissip.py 
statistical_analysis_dissip.f90 (->StatisticalData.nc) --> il y a egalement une sortie avec des spectres des donnees de dissipation

## C)------------------------------------TEMPORAL SERIES
1. /preprocess/CropData-FT (->Xhistins_X-uvtizX.nc) A bash code to extract selected variable to 
a Xhistins_X.nc file to create a reduced file Xhistins_X-uvtizX.nc with selected 
variables u[:,iz,::],v[:,iz,::],time_counter,lat,lon.
2. /preprocess/UVT_pretraitement_FT.py (->uvData-FullTime-istep-XX-nstep-XXX-iz-X.nc) concatenate
all reduced Xhistins_X-uvtizX.nc into a single file.
3. statistical_analysis_FullTime (->StatisticalData-FullTime.nc) extract u and v to 
compute the streamfunction sf and the associated spherical harmonics coefficients
sf_r and sf_i in long time series
4. /preprocess/spectra-temporel-FT.py  ...

or

ZUVT_pretraitement_FT.py (->uvData-FullTime-istep-219000-nstep-2500-niz-11.nc)
statistical_analysis_FullTime.f90 (>StatisticalData-FullTime.nc) --> Version allegee pour traiter de longue serie temporelle. Enregistre moins de donnees
On peut l exploiter ces donnees avec les
routines:
+ Modal-PltStat.py	  				      --> Plots les spectres modaux pôur des donnees de simus a l equilibre
+ spectra-temporel-FT.py  (->TempModalSpectra-im-Y-in-X.nc)   --> Fait des spectres en frequence pour chaque coeff d HS.
background routines
+ FourierTransform1D.py   --> Fait une transformer de fourrier 1D (FFT)
+ windowing.py		  --> Permet des windowing des signaus selon les finctions d Hanning et Tukey 
