# GRBModel

**Introduction**

GRBModel is an IDL routine that generates synthetic fast rise exponential decay (FRED) gamma-ray burst light curve pulses.  The code uses the empirical hardness-intensity (<a href="http://adsabs.harvard.edu/abs/1983Natur.306..451G">Golenetskii 1983</a>) and hardness-fluence (<a href="http://adsabs.harvard.edu/abs/1996Natur.381...49L">Liang & Kargatis 1996</a>) correlations to model the evolution of a user-supplied GRB spectrum, which is folded through a detector repsonse to produce a count light curve.  The detector responses from the Fermi Gamma-ray Burst Monitor (GBM) and Burst and Transient Telescope (BATSE) are currently supported.  

GRBModel.pro has served as the basis for the population synthesis code used in the following papers:

<a href="http://adsabs.harvard.edu/abs/2013ApJ...765..116K">Kocevski & Petrosian 2013</a>  - On the Lack of Time Dilation Signatures in Gamma-Ray Burst Light Curves<br>
<a href="http://adsabs.harvard.edu/abs/2012ApJ...747..146K">Kocevski 2012</a> - On the Origin of High-energy Correlations in Gamma-Ray Bursts<br>

Please refer to <a href="http://adsabs.harvard.edu/abs/2012ApJ...747..146K">Kocevski (2012)</a> for a more comprehensive description of the GRBModel code.<br>

**Usage Examples**

Compile the code in IDL
```IDL
    .compile GRBModel
```

Model the triggering pulse of GRB 130427A described in <a href="http://adsabs.harvard.edu/abs/2014Sci...343...51P">Preece et al. (2013)</a>
```IDL
	GRBModel, 2500*(1+0.3399), 3.20d57, -0.9, -2.66, 0.2, time, photon_flux, /plot, redshift=0.3399, /showbblocks, POISSON_Median=1000, /GBM, xrange=[-1,10], timerange=[-1,10], timeres=0.064, countspectrumnorm = 18, dindex=3
```

**Example Results**
```IDL
Redshift =									0.3399
Tmax_Observer =								0.267980
Epk0_Source =								3349.75
Epk0_Observer =								2500.00
Epk_Observer =								NaN
Alpha =										-0.900000
Beta =										-2.66000
Initial Photon Flux =						8.4057639
Peak Photon Flux =							46.170764
Peak Photon Flux (50-300 keV) =				17.611770
Peak Energy Flux =							1.4992391e-05
Peak Energy Flux (50-300 keV) =				3.9771904e-06
Peak Photon Luminosity (Input) =			3.2000000e+57
Peak Photon Luminosity =					1.7576802e+58
Peak Energy Luminosity =					5.7074707e+51
Energy Fluence True (Bolometric) =			4.5280511e-05
Energy Fluence True (k-Corrected) =			3.0258687e-05
Energy Fluence Estimated (k-Corrected) =	0.0000000
Energy Fluence Estimated =					0.0000000
Eiso True (Bolometric) 						1.2865057e+52
Eiso Estimated (k-Corrected) =				0.0000000
HR31 Photons =								-NaN
HR31 Counts =								2.1540972
HR41 Counts =								2.2167317
Lag CCF31 =									NaN
Noise STDEV =								30.771758
Peak Counts =								5114.3130
Peak Counts (50-300 keV) =					2825.0149
S/N ratio =									133.06587
S/N ratio (50-300 keV) =					121.04953
T90 Photon =								1.47200
T90 Count =									1.53600
T90 Count BBlocks =							1.47200
T100 =										1.60000
```
**Example Plots**<br>
Count Light Curve Plot: [LightCurve.png]<br>
Model Diagnostics Plot: [ModelDiagnostics.png]<br>

[LightCurve.png]: https://github.com/dankocevski/GRBModel/blob/master/LightCurve.png "Count Light Curve Plot"
[ModelDiagnostics.png]: https://github.com/dankocevski/GRBModel/blob/master/ModelDiagnostics.png "Model Diagnostics Plot"


**Required Arguments**<br>
* Epk0_Source								- The initial source frame Epeak (Band et al. 1993)<br>
* Luminosity0								- The burst luminosity in photons/cm2/s<br>
* alpha										- The low energy power-law slope of the Band function<br>
* beta										- The high energy power-law slope of the Band function<br>
* Tmax_Source								- The time of the pulse peak<br>

**Optional Arguments**<br>
* Time 										- An array containing the time axis data<br>
* PhotonFlux_Detector						- An array containing the photon flux light curve in the detector's band bass<br>
* EnergyFluence_Estimated 					- The observer frame energy fluence, estimated using the observed duration<br>
* EnergyFluence_kCorrected_Estimated		- The observer frame energy fluence, estimated using the observed duration, but k-corrected to a standard band bass<br>
* Epk_Observer_Int							- The initial observer frame Epeak (Band et al. 1993)<br>
* Eiso_Bolometric_True						- The true bolometric isotropic equivelent energy<br>
* Eiso_kCorrected_Estimated					- The isotropic equivelent energy estimated using the observed duration * (1+z) and a k-correction to a standard band bass<br>
* Epk 										- An array containing the time evolution of Epk in the observer frame<br>
* Fpk 										- An array containing the time evolution of the peak flux<br>
* t90_count									- The T90 in count space<br>
* t90_photon								- The T90 in photon space<br>
* t100										- The T100 duration as determined from Bayesian blocks<br>
* HR31_Photon								- The hardness ratio between channel 3 and channel 1 in photon space<br>
* HR31_Count								- The hardness ratio between channel 3 and channel 1 in count space<br>
* HR41_Count								- The hardness ratio between channel 4 and channel 1 in count space<br>
* lag31										- The lag between channel 3 and channel 1<br>
* flux 										- An array containing the time evolution of flux in the observer frame<br>
* fwhm										- The full-width half-max of the pulse<br>
* SNRatio 									- The signal to noise of the pulse<br>
* SNRatio_Trigger							- The signal to noise of the pulse in the band pass that can trigger the instrument<br>

**Required Keywords**<br>
* redshift=value 								- The burst redshift<br>

**Optional Keywords**
* rindex=value 								- The r index from the KRL function (kocevski et al. 2005)<br>
* dindex=value 								- The d index from the KRL function (kocevski et al. 2005)<br>
* POISSON_Median=value						- The background count value<br>
* CountsMatrix=CountsMatrix 				- The time evolution of the observed counts
* CountsSpectrumCube_Observer=...			- The time-resolved counts spectrum in the observer frame<br>
* PhotonSpectrumCube_Observer=...			- The time-resolved photon spectrum in the observer frame<br>
* CountRate_Channel1=CountRate_Channel1		- The count light curve for channel 1<br>
* CountRate_Channel2=CountRate_Channel2		- The count light curve for channel 2<br>
* CountRate_Channel3=CountRate_Channel3		- The count light curve for channel 3<br>
* CountRate_Channel4=CountRate_Channel4		- The count light curve for channel 4<br>
* JobNumber=value							- A user specified job number<br>
* Results=Results 							- A structure to house a range of light curve properties<br>
* timerange=timerange						- The time range to simulate<br>
* timeres=timeres 							- The resolution of the simulation<br>
* CountSpectrumNorm=CountSpectrumNorm 		- Renormalizing the count spectrum<br>
* yrange_nufnu=yrange_nufnu 				- The y range of the nuFnu plot<br>
* /plot 									- Plot the results<br>
* /rebin									- Rebin the data to 10 second bins<br>
* /ps 										- Save the plot as a postscript<br>
* /showbblocks 								- Show the Bayesian block reconstruction<br>
* /swift 									- Use the Swift energy range [10,125]<br>
* /batse 									- Use the BATSE energy range [20, 1800]<br>
* /GBM 										- Use the GBM NaI energy range [8,1000]<br>
* /save3d									- Save a plot of the 3D model showing the temporal and spectral evolution<br>
* /MakeFakeBFITS
* /MakeFakePHA
* /lag 										- Calculate the pulse lag<br>
* /ShowLagPlot 								- Show the pulse lag calculation<br>
* /XSPEC
* /CleanUp

**Required Libraryies**
Astrolib - http://idlastro.gsfc.nasa.gov