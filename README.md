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
1. <strong>/preprocess/CropData-FT</strong> (needs: ncrcat)   <strong>--></strong> Xhistins_x-uvtizxx.nc <br/>
<em> A bash code that extracts in a Xhistins_x.nc file and creates a reduced file Xhistins_x-uvtizxx.nc with selected variables u[:,iz,::],v[:,iz,::],time_counter,lat,lon at a given level iz.</em>
2. <strong>/preprocess/UVT_pretraitement_FT.py</strong> (needs python)     <strong>--></strong> uvData-FullTime-istep-xx-nstep-xxx-iz-x.nc)  <br/>
<em>Concatenate all reduced Xhistins_x-uvtizx.nc into a single file.</em>
3. <strong>statistical_analysis_FullTime</strong> (needs: Fortran 90 & spherepack)    <strong>--></strong> StatisticalData-FullTime.nc <br/>
<em>Extract u and v to compute the streamfunction sf and the associated spherical harmonics coefficients sf_r and sf_i in long time series.</em>
4. <strong>/preprocess/spectra-temporel-FT.py </strong> (needs: python)  <strong>--></strong> TempModalSpectra-im-0-upto-n-35.nc <br/>
<em>Compute the frequency spectra by applying a fourier transform in time to the harmonics coefficients $sf_r(\omega)$ at given modes $m$ and $n$. A rolling averaged in time is apply to the frequency spectra. </em>
5. <strong>/plots/Plt-spectra-temporel.py</strong> (needs: python, LoadInfos.py & filePTS.zono.temp)   <em>Subplots of the frequency spectra for severl $n$ and a given $m$.</em>

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
