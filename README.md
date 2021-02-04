# POST
Statistical tools to study geophysical flows

```
cd spherepack3.2
make
cd ..
make statistical_analysis_xyzt
```

## A) FOR GENERAL STATISTICAL ANALYSIS

<ins>Steps,</ins>
1. Preprocess Xhistins (python in preprocess)
2. statistical_analysis_xyzt uvData-istep-0-nstep-40-niz-12.nc
3. e.g. Tevo-PltStat.py StatisticalData.nc
ZUVT_pretraitement_HV.py (->uvData-istep-x-nstep-xxx-niz-xx.nc)
statistical_analysis_xyzt.f90 (->StatisticalData.nc) --> routine de travail de base.

<ins>Background functions,</ins>

+ Tevo-PltStat.py         --> time evolution of the data/spectra/fluxes
+ I-PltStat.py            --> Mean time and each levels
+ comparaison-v2.py	  --> comparaison de toutes les simus
+ comparaison-rotrate.py  --> comparaison des simus ayant differents rotation rates
+ comparaison-kf.py	  --> comparaison des simus ayant different kf de forcage
background routines
+ LoadInfos.py		  --> Load les data from filePTS.zono.temp
+ filePTS.zono.temp	  --> Enter data of the specific simulation

## B) TEMPORAL SERIES & SPECTRA

<ins>Steps,</ins>
1. <strong>/preprocess/CropData-FT</strong> (needs: ncrcat)   <strong>--></strong> Xhistins_x-uvtizxx.nc <br/>
<em> A bash code that extracts in a Xhistins_x.nc file and creates a reduced file Xhistins_x-uvtizxx.nc with selected variables u[:,iz,::],v[:,iz,::],time_counter,lat,lon at a given level iz.</em>
2. <strong>/preprocess/UVT_pretraitement_FT.py</strong> (needs python)     <strong>--></strong> uvData-FullTime-istep-xx-nstep-xxx-iz-x.nc)*  <br/>
<em>Concatenate all reduced Xhistins_x-uvtizx.nc into a single file.</em>
3. <strong>statistical_analysis_FullTime</strong> (needs: Fortran 90 & spherepack)    <strong>--></strong> StatisticalData-FullTime.nc* <br/>
<em>Extract u and v to compute the streamfunction sf and the associated spherical harmonics coefficients sf_r and sf_i in long time series.</em>
4. <strong>/preprocess/spectra-temporel-FT.py </strong> (needs: python, FourierTransform1D.py)  <strong>--></strong> TempModalSpectra-im-0-upto-n-35.nc <br/>
<em>Compute the frequency spectra by applying a fourier transform in time to the harmonics coefficients sf_r at given modes m and n. A rolling averaged in time is apply to the frequency spectra. </em>
5. <strong>/plots/Plt-spectra-temporel.py</strong> (needs: python, LoadInfos.py & filePTS.zono.temp)   <em>Subplots of the frequency spectra for severl n and a given m.</em>

<ins>Background functions,</ins>

+ FourierTransform1D.py    --> For 1D Fourier Transform (FFT)
+ windowing.py		   --> For Hanning and Tukey windowing of the signals. 
+ ZUVT_pretraitement_FT.py --> Can extract suited netcdf file to be used with <em>statistical_analysis_FullTime</em>

*Usually uvData-FullTime-istep-xx-nstep-xxx-iz-x.nc & StatisticalData-FullTime.nc cannot exceeds more than ~1500 time steps, then I advise to split in different files and then concatenate StatisticalData-FullTime.nc files in a single one using ncrcat.

## C) FOR DISSIPATION SPECTRAL STUDIES
ZUVT_pretraitement_dissip.py 
statistical_analysis_dissip.f90 (->StatisticalData.nc) --> il y a egalement une sortie avec des spectres des donnees de dissipation
