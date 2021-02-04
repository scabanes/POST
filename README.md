# POST
<strong>Package Of Statistical Tools</strong> to study geophysical flows in spherical geometry

```
cd spherepack3.2
make
cd ..
make statistical_analysis_xyzt
```

## A) FOR GENERAL STATISTICAL ANALYSIS

<ins>Steps</ins>
1. <strong>preprocess/ZUVT_pretraitement_HV.py</strong> (needs: python)    		<strong>--></strong> uvData-istep-xx-nstep-xx-niz-xx.nc <br/>
<em>Select variables u,v,u_SMerid, v_SMerid, and crop dimension [time,altitude,longitude] in files Xhistins_x.nc and concatenate into a single file.</em>
2. <strong>statistical_analysis_xyzt</strong> (needs: Fortran 90 & spherepack)  	<strong>--></strong> StatisticalData.nc <br/>
Arguments:
+ uvData-istep-xx-nstep-xx-niz-xx.nc is the file to read
+ -h    : brief help
+ -mt x : x is the number of step
+ -istp : first temporal step
+ -tstp : final temporal step
3. <strong>plots/Plt-EquiStat.py</strong> (needs: python, LoadInfos.py, filePTS.zono.temp & filePTS.plots.infos)<br/>
<em>Plot average spectral quantity: energy spectra, energy and enstrophy fluxes.</em>

<ins>Background functions</ins>
```
+ LoadInfos.py		  --> Load les data from filePTS.zono.temp
+ filePTS.zono.temp	  --> Input infos of the dataset (from Dev/ > nest in the working directory)
+ filePTS.plots.infos	  --> Input infos to plot (from Dev/ > nest in the working directory)
To be updated:
+ Tevo-PltStat.py         --> time evolution of the data/spectra/fluxes
+ I-PltStat.py            --> Mean time and each levels
+ comparaison-v2.py	  --> comparaison de toutes les simus
+ comparaison-rotrate.py  --> comparaison des simus ayant differents rotation rates
+ comparaison-kf.py	  --> comparaison des simus ayant different kf de forcage
```

## B) TEMPORAL SERIES & SPECTRA

<ins>Steps</ins>
1. <strong>preprocess/CropData-FT</strong> (needs: ncrcat)   <strong>--></strong> Xhistins_x-uvtizxx.nc <br/>
<em> A bash code that extracts in a Xhistins_x.nc file and creates a reduced file Xhistins_x-uvtizxx.nc with selected variables u[:,iz,::],v[:,iz,::],time_counter,lat,lon at a given level iz.</em>
2. <strong>/preprocess/UVT_pretraitement_FT.py</strong> (needs python)     <strong>--></strong> uvData-FullTime-istep-xx-nstep-xxx-iz-x.nc)*  <br/>
<em>Concatenate all reduced Xhistins_x-uvtizx.nc into a single file.</em>
3. <strong>statistical_analysis_FullTime</strong> (needs: Fortran 90 & spherepack)    <strong>--></strong> StatisticalData-FullTime.nc* <br/>
<em>Extract u and v to compute the streamfunction sf and the associated spherical harmonics coefficients sf_r and sf_i in long time series.</em>
4. <strong>preprocess/spectra-temporel-FT.py </strong> (needs: python, FourierTransform1D.py)  <strong>--></strong> TempModalSpectra-im-0-upto-n-35.nc <br/>
<em>Compute the frequency spectra by applying a fourier transform in time to the harmonics coefficients sf_r at given modes m and n. A rolling averaged in time is apply to the frequency spectra. </em>
5. <strong>plots/Plt-spectra-temporel.py</strong> (needs: python, LoadInfos.py & filePTS.zono.temp) <br/>
<em>Subplots of the frequency spectra for several n and a given m.</em>

<ins>Background functions</ins>
```
+ FourierTransform1D.py    --> For 1D Fourier Transform (FFT)
+ windowing.py		   --> For Hanning and Tukey windowing of the signals.
+ LoadInfos.py		   --> Load les data from filePTS.zono.temp
+ filePTS.zono.temp	   --> Input infos of the dataset (> nest in the working directory)
+ ZUVT_pretraitement_FT.py --> Can extract suited netcdf file to be used 
			       with statistical_analysis_FullTime
```

*<em>Usually uvData-FullTime-istep-xx-nstep-xxx-iz-x.nc & StatisticalData-FullTime.nc cannot exceeds more than ~1500 time steps, then I advise to split in different files and then concatenate StatisticalData-FullTime.nc files in a single one using ncrcat.</em>

----- 
## C) FOR DISSIPATION SPECTRAL STUDIES
ZUVT_pretraitement_dissip.py 
statistical_analysis_dissip.f90 (->StatisticalData.nc) --> il y a egalement une sortie avec des spectres des donnees de dissipation
