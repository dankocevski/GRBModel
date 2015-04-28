# GRBModel

GRBModel is an IDL routine that generates synthetic fast rise exponential decay (FRED) GRB light curve pulses.  The code uses the empirical hardness-intensity (Golenetskii 1993) and hardness-fluence (Liang & Kargatis 1996) correlations to model the evolution of a user-supplied GRB spectrum, which is folded through a detector repsonse to produce a count light curve.  The detector responses from the Fermi Gamma-ray Burst Monitor (GBM) and Burst and Transient Telescope (BATSE) instruments are currently supported.  

GRBModel.pro has served as the basis for the population modeling code used in the following papers:

<a href="http://adsabs.harvard.edu/abs/2013ApJ...765..116K">Kocevski & Petrosian 2013</a>  - On the Lack of Time Dilation Signatures in Gamma-Ray Burst Light Curves
<a href="http://adsabs.harvard.edu/abs/2012ApJ...747..146K">Kocevski 2012</a> - On the Origin of High-energy Correlations in Gamma-Ray Bursts

Please refer to [Kocevski (2012)]: http://adsabs.harvard.edu/abs/2012ApJ...747..146K for a more comprehensive description of the GRBModel code.

**Required Arguments**
Epk0_Source								- The initial source frame Epeak (Band et al. 1993)<br>
Luminosity0								- The burst luminosity in photons/cm2/s<br>
alpha									- The low energy power-law slope of the Band function
beta									- The high energy power-law slope of the Band function
Tmax_Source								- The time of the pulse peak

**Optional Arguments**
Time 									- An array containing the time axis data
PhotonFlux_Detector						- An array containing the photon flux light curve in the detector's band bass
EnergyFluence_Estimated 				- The observer frame energy fluence, estimated using the observed duration
EnergyFluence_kCorrected_Estimated		- The observer frame energy fluence, estimated using the observed duration, but k-corrected to a standard band bass
Epk_Observer_Int						- The initial observer frame Epeak (Band et al. 1993) 	
Eiso_Bolometric_True					- The true bolometric isotropic equivelent energy
Eiso_kCorrected_Estimated				- The isotropic equivelent energy estimated using the observed duration * (1+z) and a k-correction to a standard band bass
Epk 									- An array containing the time evolution of Epk in the observer frame
Fpk 									- An array containing the time evolution of the peak flux
t90_count								- The T90 in count space
t90_photon								- The T90 in photon space
t100									- The T100 duration as determined from Bayesian blocks
HR31_Photon								- The hardness ratio between channel 3 and channel 1 in photon space
HR31_Count								- The hardness ratio between channel 3 and channel 1 in count space
HR41_Count								- The hardness ratio between channel 4 and channel 1 in count space
lag31									- The lag between channel 3 and channel 1
flux 									- An array containing the time evolution of flux in the observer frame
fwhm									- The full-width half-max of the pulse
SNRatio 								- The signal to noise of the pulse
SNRatio_Trigger							- The signal to noise of the pulse in the band pass that can trigger the instrument

**Required Keywords**
redshift=value 							- The burst redshift

**Optional Keywords**
rindex=value 							- The r index from the KRL function (kocevski et al. 2005)
dindex=value 							- The d index from the KRL function (kocevski et al. 2005)
POISSON_Median=value					- The background count value
CountsMatrix=CountsMatrix, 
CountsSpectrumCube_Observer=...			- The time-resolved counts spectrum in the observer frame
PhotonSpectrumCube_Observer=...			- The time-resolved photon spectrum in the observer frame
CountRate_Channel1=CountRate_Channel1	- The count light curve for channel 1 
CountRate_Channel2=CountRate_Channel2	- The count light curve for channel 2
CountRate_Channel3=CountRate_Channel3	- The count light curve for channel 3
CountRate_Channel4=CountRate_Channel4	- The count light curve for channel 4
JobNumber=JobNumber						- A user specified job number
Results=Results 						- A structure to house a range of light curve properties
timerange=timerange						- The time range to simulate
timeres=timeres 						- The resolution of the simulation
CountSpectrumNorm=CountSpectrumNorm 	- Renormalizing the count spectrum
yrange_nufnu=yrange_nufnu 				- The y range of the nuFnu plot
/plot 									- Plot the results
/rebin									- Rebin the data to 10 second bins
/ps 									- Save the plot as a postscript
/showbblocks 							- Show the Bayesian block reconstruction
/swift 									- Use the Swift energy range [10,125]
/batse 									- Use the BATSE energy range [20, 1800]
/GBM 									- Use the GBM NaI energy range [8,1000]
/save3d									- Save a plot of the 3D model showing the temporal and spectral evolution 
/MakeFakeBFITS=MakeFakeBFITS
/MakeFakePHA=MakeFakePHA
/lag 									- Calculate the pulse lag
/ShowLagPlot 							- Show the pulse lag calculation
/XSPEC=XSPEC
/CleanUp=CleanUp

**Required Libraryies:
Astrolib - http://idlastro.gsfc.nasa.gov

**Usage Examples: 

Compile the code in IDL
```IDL
    .compile GRBModel
```

Model the triggering pulse from GRB 130427A (from Preece et al. 2013)
```IDL
    z = 1
	GRBModel, 2500*(1+z), 3.20d57, -0.9, -2.66, 0.2, time, photon_flux, /plot, redshift=z, /showbblocks, POISSON_Median=1000, /GBM, xrange=[-1,10], timerange=[-1,10], timeres=0.064, countspectrumnorm = 18, dindex=3
```

**Example Results:
```IDL
Redshift =									1
Tmax_Observer =								0.400000
Epk0_Source =								5000
Epk0_Observer =								2500
Epk_Observer =								NaN
Alpha =										-0.900000
Beta =										-2.66000
Initial Photon Flux =						0.60670793
Peak Photon Flux =							2.3075163
Peak Photon Flux (50-300 keV) =				0.88648505
Peak Energy Flux =							7.5397983e-07
Peak Energy Flux (50-300 keV) =				2.0016942e-07
Peak Photon Luminosity (Input) =			3.2000000e+57
Peak Photon Luminosity =					1.2170687e+58
Peak Energy Luminosity =					3.9767660e+51
Energy Fluence True (Bolometric) =			7.1733277e-06
Energy Fluence True (k-Corrected) =			3.0758651e-06
Energy Fluence Estimated (k-Corrected) =	0.0000000
Energy Fluence Estimated =					0.0000000
Eiso True (Bolometric) 						1.8917380e+52
Eiso Estimated (k-Corrected) =				0.0000000
HR31 Photons =								-NaN
HR31 Counts =								2.3197714
HR41 Counts =								2.8740169
Lag CCF31 =									NaN
Noise STDEV =								32.387506
Peak Counts =								1218.4019
Peak Counts (50-300 keV) =					567.89948
S/N ratio =									7.0470725
S/N ratio (50-300 keV) =					6.7622964
T90 Photon =								2.04800
T90 Count =									-9.53600
T90 Count BBlocks =							9.79200
T100 =										1.02400
```