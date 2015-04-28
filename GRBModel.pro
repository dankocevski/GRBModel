;; ==========================================================================================
;; File:  GRBmodel.pro
;; Author:  Daniel Kocevski <kocevski@berkeley.edu>
;; Date Created:  May 21, 2008
;; Last Changed:  -
;; 
;; Purpose:  This program produces a 3-dimensional model of a GRB
;; Dependancies:  None
;; ==========================================================================================


;pro GRBmodel, Epk0, Lum0, alpha, beta, tmax, time, PhotonFlux, redshift=redshift, rindex = r, dindex = d, plot=plot, noise=noise, rebin=rebin, ps=ps, showbblocks=showbblocks
pro GRBmodel, Epk0_Source, Luminosity0, alpha, beta, Tmax_Source, Time, PhotonFlux_Detector, EnergyFluence_Estimated, EnergyFluence_kCorrected_Estimated, Epk_Observer_Int, Eiso_Bolometric_True, Eiso_kCorrected_Estimated, Epk, Fpk, t90_count, t90_photon, t100, HR31_Photon, HR31_Count, HR41_Count, lag31, flux, fwhm, SNRatio, SNRatio_Trigger, redshift=redshift, rindex = r, dindex = d, plot=plot, noise=noise, rebin=rebin, ps=ps, showbblocks=showbblocks, swift=swift, batse=batse, POISSON_Median=POISSON_Median, save3d=save3d, CountsMatrix=CountsMatrix, CountsSpectrumCube_Observer=CountsSpectrumCube_Observer, PhotonSpectrumCube_Observer=PhotonSpectrumCube_Observer, MakeFakeBFITS=MakeFakeBFITS, MakeFakePHA=MakeFakePHA,lag=lag, ShowLagPlot=ShowLagPlot, CountRate_Channel1=CountRate_Channel1, CountRate_Channel2=CountRate_Channel2, CountRate_Channel3=CountRate_Channel3, CountRate_Channel4=CountRate_Channel4, JobNumber=JobNumber, XSPEC=XSPEC, CleanUp=CleanUp, Results=Results, GBM=GBM, _extra=_extra, timerange=timerange, timeres=timeres, CountSpectrumNorM=CountSpectrumNorm, yrange_nufnu=yrange_nufnu

IF (n_params() EQ 0) THEN BEGIN 
	print,'-Usage: GRBmodel, Epk0_Source, Flux0, alpha, beta, tmax, time, PhotonFlux, EnergyFluence, epk_Observer, Eiso, Epk, Fpk, t90, t100, HR31_Photon, HR31_Count, HR41_Count, lag31, flux, fwhm, redshift=redshift, rindex=r, dindex=d, plot=plot, noise=noise, rebin=rebin, swift=swift, batse=batse, lag=lag, ShowLagPlot=ShowLagPlot, JobNumber=JobNumber, XSPEC=XSPEC, CleanUp=CleanUp'
    print, ''
    print, 'Common Usage: grbmodel, 100, 25, -1.0, -2.5, 10, /plot'
    RETURN
ENDIF

;; Setting up a color table
IF KEYWORD_SET(plot) OR KEYWORD_SET(ShowLagPlot) THEN BEGIN
	loadct, 39
	device, decompose=0
ENDIF

IF (KEYWORD_SET(timeres) eq 0) THEN timeres = 1

;; Assign a random job number if one is not specified.  This allows multiple instances of this code to be running at the same time.
IF (n_elements(JobNumber) EQ 0) THEN JobNumber = fix(abs( randomn(seed,1) * 1000 ))
JobNumber = string(JobNumber)

;; Setting the speed of light in cm/s
c = 29979245800D

;; Check that the luminosity is a double.  If not, make it a double
IF (Luminosity0 LT 100) THEN Luminosity0 = 10d^Luminosity0

;; Defining detector bandpasses
SwiftBandPass = [10,125]			;; Swift
SwiftTriggerBandPass = [10,125]		;; Swift

BatseBandPass = [20, 1800]			;; BATSE
BatseTriggerBandPass = [50,300] 	;; BATSE
BatseChannel1BandPass = [25,50]		;; BATSE
BatseChannel2BandPass = [50,100]	;; BATSE
BatseChannel3BandPass = [100,300]	;; BATSE
BatseChannel4BandPass = [300,1800]	;; BATSE

GBMBandPass = [8,1000]				;; GBM NaI
GBMTriggerBandPass = [50,300] 		;; GBM NaI
GBMChannel1BandPass = [25,50]		;; GBM NaI
GBMChannel2BandPass = [50,100]		;; GBM NaI
GBMChannel3BandPass = [100,300]		;; GBM NaI
GBMChannel4BandPass = [300,1800]	;; GBM NaI


;; Use the Swift bandpass by default
IF KEYWORD_SET(batse) THEN BEGIN
	BandPass = BatseBandPass 
	TriggerBandPass = BatseTriggerBandPass
	IF N_ELEMENTS(POISSON_Median) EQ 0 THEN POISSON_Median = 2000.0;
;;	BackgroundCountSpectrumNormalization = 2000.0
	BackgroundCountSpectrumNormalization = POISSON_Median

ENDIF
	
IF KEYWORD_SET(GBM) THEN BEGIN	
	BandPass = GBMBandPass 
	TriggerBandPass = GBMTriggerBandPass
	IF N_ELEMENTS(POISSON_Median) EQ 0 THEN POISSON_Median = 2000.0;
;;	BackgroundCountSpectrumNormalization = 2000.0	
	BackgroundCountSpectrumNormalization = POISSON_Median
ENDIF

IF KEYWORD_SET(Swift) THEN BEGIN	
	BandPass = SwiftBandPass 
	TriggerBandPass = SwiftTriggerBandPass
	IF N_ELEMENTS(POISSON_Median) EQ 0 THEN POISSON_Median = 2000.0;
;;	BackgroundCountSpectrumNormalization = 2000.0	
	BackgroundCountSpectrumNormalization = POISSON_Median
ENDIF

BandPass = BandPass - 1.0

;; Defining a Standard "Bolometric" Bandpass betweeen 10 and 10^4 keV
StandardBandPass = [10,10000]		

;; Definding the time array
IF (N_ELEMENTS(timerange) EQ 0) THEN BEGIN
	time = findgen(300)/1.0
	time_pretrigger = -1 * findgen(50)/1.0
	time_pretrigger = reverse(time_pretrigger(1:*))
ENDIF ELSE BEGIN
	time = findgen(timerange[1]/timeres)*timeres
	time_pretrigger = -1 * findgen(abs(timerange[0])/timeres)*timeres
	time_pretrigger = reverse(time_pretrigger(1:*))

ENDELSE



;; Defining the observer and source frame energy (keV)
;energy = dindgen(1790) + 1.0
Energy_Source = dindgen(20000) + 1.0
;Energy_Source = dindgen(100000) + 1.0
Energy_Observer = Energy_Source / (1 + Redshift)

;; Finding the index of the bandpasses
whereclose, Energy_Observer, bandpass(0), LowerIndex, /silent
whereclose, Energy_Observer, bandpass(1), UpperIndex, /silent
whereclose, Energy_Observer, TriggerBandPass(0), LowerIndex_Trigger, /silent
whereclose, Energy_Observer, TriggerBandPass(1), UpperIndex_Trigger, /silent

;; Find the index of a standardized source frame bandpass for k-correction
whereclose, Energy_Source, StandardBandPass(0), kLowerIndex, /silent
whereclose, Energy_Source, StandardBandPass(1), kUpperIndex, /silent

;; Make sure these indices only have one element
LowerIndex = LowerIndex(0)
UpperIndex = UpperIndex(0)
LowerIndex_Trigger = LowerIndex_Trigger(0)
UpperIndex_Trigger = UpperIndex_Trigger(0)
kLowerIndex = kLowerIndex(0)
kUpperIndex = kUpperIndex(0)

;; Applying a bandpass filter
;energy = energy(bandpass(0):bandpass(1))

;; Making arrays to save the time resolved spectra and bolometric light curves and other things
PhotonSpectrumCube_Source = dblarr(n_elements(time),n_elements(Energy_Source))
Photon_Fnu_SpectrumCube_Source = dblarr(n_elements(time),n_elements(Energy_Source))
Energy_Fnu_SpectrumCube_Source = dblarr(n_elements(time),n_elements(Energy_Source))
Photon_nuFnu_SpectrumCube_Source = dblarr(n_elements(time),n_elements(Energy_Source))

PhotonSpectrumCube_Observer = dblarr(n_elements(time),n_elements(Energy_Observer))
Photon_Fnu_SpectrumCube_Observer = dblarr(n_elements(time),n_elements(Energy_Observer))
Energy_Fnu_SpectrumCube_Observer = dblarr(n_elements(time),n_elements(Energy_Observer))
Photon_nuFnu_SpectrumCube_Observer = dblarr(n_elements(time),n_elements(Energy_Observer))
Photon_Fnu_TimeIntegratedSpectrum_Observer = dblarr(n_elements(Energy_Observer))
Energy_Fnu_TimeIntegratedSpectrum_Observer = dblarr(n_elements(Energy_Observer))
Photon_nuFnu_TimeIntegratedSpectrum_Observer = dblarr(n_elements(Energy_Observer))

PhotonFlux_Detector = dblarr(n_elements(time))
PhotonFlux_DetectorTrigger = dblarr(n_elements(time))
PhotonFlux_kCorrected = dblarr(n_elements(time))
EnergyFlux_Detector = dblarr(n_elements(time))
EnergyFlux_DetectorTrigger = dblarr(n_elements(time))
EnergyFlux_kCorrected = dblarr(n_elements(time))
PhotonLuminosity_Bolometric = dblarr(n_elements(time))
EnergyLuminosity_Bolometric = dblarr(n_elements(time))

channel1 = dblarr(n_elements(time))
channel2 = dblarr(n_elements(time))
channel3 = dblarr(n_elements(time))
channel4 = dblarr(n_elements(time))
FEpk = dblarr(n_elements(time))
CountRate = dblarr(n_elements(time))
CountRate_Channel1 = dblarr(n_elements(time))
CountRate_Channel2 = dblarr(n_elements(time))
CountRate_Channel3 = dblarr(n_elements(time))
CountRate_Channel4 = dblarr(n_elements(time))
CountRate_kCorrected = dblarr(n_elements(time))
CountRate_PreTrigger = dblarr(n_elements(time_pretrigger))
CountRate_PreTrigger_Trigger = dblarr(n_elements(time_pretrigger))
CountRate_Channel1_PreTrigger = dblarr(n_elements(time_pretrigger))
CountRate_Channel2_PreTrigger = dblarr(n_elements(time_pretrigger))
CountRate_Channel3_PreTrigger = dblarr(n_elements(time_pretrigger)) 
CountRate_Channel4_PreTrigger = dblarr(n_elements(time_pretrigger)) 
CountRate_Trigger = dblarr(n_elements(time))


;; Defining the angular timescale.  Typical values: Gamma = 100, Radius = 1E13, RoG = 1E12
;tau = Radius / (2.0 * Gamma * c)

;; Time dilating and redshifting if necessary
IF (N_ELEMENTS(redshift) NE 0) THEN BEGIN
	dlum, redshift, dl, /silent, /cm
	Epk0_Observer = Epk0_Source / (1 + redshift)
	Tmax_Observer = Tmax_Source * (1 + redshift)
	Flux0 = Luminosity0 / (4*!PI*dl^2.0)
ENDIF

;; Defining the Epk evolution given some initial Epk0
;Epk = Epk0 / (1 + (time/tau))

;; Defining the Flux evolution given some initial F0
;Flux = ( Flux0 / ((1 + (time/tau))^2.0) ) * (time) / (1 + redshift) 
;Flux = ( Flux0 / ((1 + (time/tau))^2.0) )

;; Using The KRL Function (This seems to be the only way to consistantly include the rise phase)
IF (N_ELEMENTS(r) EQ 0) THEN r = 1.0
IF (N_ELEMENTS(d) EQ 0) THEN d = 2.2

;; Caculating Source Frame Luminosity(t) and Epk_Source(t)
Luminosity = Luminosity0 * ((time/Tmax_Source)^r) * ( (d/(d+r)) + ((r/(d+r))*((time/Tmax_Source)^(r+1))) )^(-1*(r+d)/(r+1))		;; photons/s (depends on the units of Luminosity0)
Epk_Source = Epk0_Source * ( (d/(d+r)) + (r/(d+r))*((time/Tmax_Source)^(r+1)) )^(-1*(d-1)/(r+1))								;; keV
	
;; Caculating Observer Frame Flux(t) and Epk_Observer(t)
Flux = Flux0 * ((time/Tmax_Observer)^r) * ( (d/(d+r)) + ((r/(d+r))*((time/Tmax_Observer)^(r+1))) )^(-1*(r+d)/(r+1))		;; photons/s (depends on the units of Luminosity0)
Epk_Observer = Epk0_Observer * ( (d/(d+r)) + (r/(d+r))*((time/Tmax_Observer)^(r+1)) )^(-1*(d-1)/(r+1))								;; keV



;; Generating a Band model spectrum using the Epk and Flux values calculated above
print, "Generating event..."
FOR i=0D, N_ELEMENTS(time)-1 DO BEGIN
	
	;; The flux term includes a correction in order to keep the amplitude constant for different Epk values.  This ensures that flux decrease is soley do to doppler effects.
;	band, Energy_Source, [Flux(i) * (Epk0/Epk(i)), Epk(i),alpha,beta], Photons
;	band, Energy_Source, [Flux(i), Epk(i),alpha,beta], Photons


	; ==== Source Frame ===== ;

	;; Calculating the source frame photon spectra (Band function)
	band, Energy_Source, [1, Epk_Source(i),alpha,beta], PhotonSpectrum_Source		;; This returns a photon spectrum in photons cm^-2 s^-1 keV^-1
		
	;; Renormalizing the photon spectrum to 1 and then setting the Luminosity @ Epk_Source to the model Luminosity(t)
	IF (beta LE -2) THEN BEGIN
		PhotonSpectrumCube_Source(i,*) = ( PhotonSpectrum_Source / max(PhotonSpectrum_Source) ) * Luminosity(i)
	ENDIF ELSE BEGIN
		PhotonSpectrumCube_Source(i,*) = ( PhotonSpectrum_Source / max(PhotonSpectrum_Source) ) * Luminosity(i)
	ENDELSE
	
	;; Defining the energy and nuFnu spectra in the source frame
	Photon_Fnu_SpectrumCube_Source(i,*) = PhotonSpectrumCube_Source(i,*) * (Energy_Source)									;; photons
	Energy_Fnu_SpectrumCube_Source(i,*) = PhotonSpectrumCube_Source(i,*) * (Energy_Source) * 1000 * 1.60217646D-12			;; ergs	
	Photon_nuFnu_SpectrumCube_Source(i,*) = PhotonSpectrumCube_Source(i,*) * (Energy_Source^2.0)							;; photons keV

	;; Finding the source frame bolometric energy and photon flux light curves
	PhotonLuminosity_Bolometric(i) = int_tabulated(Energy_Source,PhotonSpectrumCube_Source(i,*))
	EnergyLuminosity_Bolometric(i) = int_tabulated(Energy_Source,Energy_Fnu_SpectrumCube_Source(i,*))


	; ==== Observer Frame ===== ;
	
	;; Calculating the observer frame photon spectra (Band function)
	band, Energy_Observer, [1, Epk_Observer(i),alpha,beta], PhotonSpectrum_Observer		;; This returns a photon spectrum in photons cm^-2 s^-1 keV^-1
		
	;; Renormalizing the photon spectrum to 1 and then setting the Flux @ Epk to the model Flux(t)
	IF (beta LE -2) THEN BEGIN
		PhotonSpectrumCube_Observer(i,*) = ( PhotonSpectrum_Observer / max(PhotonSpectrum_Observer) ) * Flux(i)
	ENDIF ELSE BEGIN
		PhotonSpectrumCube_Observer(i,*) = ( PhotonSpectrum_Observer / max(PhotonSpectrum_Observer) ) * Flux(i)
	ENDELSE	
	
	;; Defining the energy and nuFnu spectra in the source frame
	Photon_Fnu_SpectrumCube_Observer(i,*) = PhotonSpectrumCube_Observer(i,*) * (Energy_Observer)									;; photons cm^-2 s^-1
	Energy_Fnu_SpectrumCube_Observer(i,*) = PhotonSpectrumCube_Observer(i,*) * (Energy_Observer) * 1000 * 1.60217646D-12			;; ergs cm^-2 s^-1	
	Photon_nuFnu_SpectrumCube_Observer(i,*) = PhotonSpectrumCube_Observer(i,*) * (Energy_Observer^2.0)								;; photons keV cm^-2 s^-1
	
	;; Finding the observer frame energy and photon flux light curve in the detector band pass
	PhotonFlux_Detector(i) = int_tabulated(Energy_Observer(LowerIndex:UpperIndex),PhotonSpectrumCube_Observer(i,LowerIndex:UpperIndex))
	EnergyFlux_Detector(i) = int_tabulated(Energy_Observer(LowerIndex:UpperIndex),Energy_Fnu_SpectrumCube_Observer(i,LowerIndex:UpperIndex))
	
	;; Finding the observer frame energy and photon flux light curve in band pass over which the instrument triggers
	PhotonFlux_DetectorTrigger(i) = int_tabulated(Energy_Observer(LowerIndex_Trigger:UpperIndex_Trigger),PhotonSpectrumCube_Observer(i,LowerIndex_Trigger:UpperIndex_Trigger))
	EnergyFlux_DetectorTrigger(i) = int_tabulated(Energy_Observer(LowerIndex_Trigger:UpperIndex_Trigger),Energy_Fnu_SpectrumCube_Observer(i,LowerIndex_Trigger:UpperIndex_Trigger))
	
	;; Finding the k-corrected energy and photon flux light curve in a standard band pass (1 to 10000 keV)
	PhotonFlux_kCorrected(i) = int_tabulated(Energy_Source(kLowerIndex:kUpperIndex),PhotonSpectrumCube_Observer(i,kLowerIndex:kUpperIndex))
	EnergyFlux_kCorrected(i) = int_tabulated(Energy_Source(kLowerIndex:kUpperIndex),Energy_Fnu_SpectrumCube_Observer(i,kLowerIndex:kUpperIndex))

ENDFOR

;; Finding the total source frame bolometric photon and energy fluence (Eiso) over thr true duration
Siso_Bolometric_True = int_tabulated(time, PhotonLuminosity_Bolometric) 
Eiso_Bolometric_True = int_tabulated(time, EnergyLuminosity_Bolometric) 

;; Finding the tru bolometric energy fluence from the source frame quantities 
EnergyFluence_Bolometric_True  = Eiso_Bolometric_True / ((4*!PI*(dl^2.0))/(1+redshift))

;; Finding the total observer frame photon and energy fluence over the true duration
PhotonFluence_Detector_True = int_tabulated(time, PhotonFlux_Detector)
EnergyFluence_Detector_True = int_tabulated(time, EnergyFlux_Detector)

;; Finding the total observer frame k-corrected photon and energy fluence over the true duration
PhotonFluence_kCorrected_True = int_tabulated(time, PhotonFlux_kCorrected) 
EnergyFluence_kCorrected_True = int_tabulated(time, EnergyFlux_kCorrected) 



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Create a count rate light curve using a detector data response matrix
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

IF KEYWORD_SET(BATSE) THEN BEGIN

	;; Load the data response matrix
	DRMFilename = './Resources/her_drm_1_907.fits'
	DRM_PHA1Filename = './Resources/her_drm_1_907.rmf'		
	BFITSFilename = './Resources/her_bfits_1_907.fits'	
	BFITSFilenameNew = './Resources/her_bfits_1_Sim' + JobNumber + '.fits'
	PHA1Filename = './Resources/her_bfits_1_Sim_FromRMFit.pha'
	PHA1FilenameNew = './Resources/her_bfits_1_Sim' + JobNumber + '.pha'
	BackgroundPHA1FilenameNew = './Resources/her_bfits_1_Sim' + JobNumber + '.bak'	
	
	;; Make sure the new names don't have any spaces in them
	BFITSFilenameNew = strcompress(BFITSFilenameNew, /remove_all)
	PHA1FilenameNew = strcompress(PHA1FilenameNew, /remove_all)
	BackgroundPHA1FilenameNew = strcompress(BackgroundPHA1FilenameNew, /remove_all)
	
	;; Read in the data response matrix 
	print, 'Using DRM File: ', drmfilename
	makebatsedrm, drmfilename, resp_matrix, e_edges, c_edges, n_zeros

ENDIF

IF KEYWORD_SET(GBM) THEN BEGIN

	;; Load the data response matrix  (140 input energies, 128 output energies)
	DRMFilename = './Resources/glg_cspec_n9_bn130427324_v00.rsp2'
	BFITSFilename = './Resources/glg_cspec_n9_bn130504978_v01.pha'
	
	data1 = mrdfits(DRMFilename,1,hdr1)
	data2 = mrdfits(DRMFilename,2,hdr2)
	resp_matrix = data2.matrix

	;; Define the input energy boundries (128)
	ENERG_LO = data2.ENERG_LO
	ENERG_HI = data2.ENERG_HI
	E_EDGES	= [ENERG_LO,ENERG_HI(n_elements(ENERG_HI)-1)]
		
	;; Define the output energy boundries (140)
	E_MIN = data1.E_MIN
	E_MAX = data1.E_MAX
	C_EDGES = [E_MIN,E_MAX(n_elements(E_MAX)-1)]


	
ENDIF 

	
	;; Create an array containing the centers and widths of the channel energy bins.  Note that the definition of e_edges and c_edges are reversed in bfit files
	ChannelEnergy = ( C_EDGES(0:n_elements(C_EDGES)-2) + C_EDGES(1:*)) / 2.0
	ChannelEnergyBinWidth = C_EDGES(1:*)  - C_EDGES(0:n_elements(C_EDGES)-2)
;	PhotonEnergyBinWidth = E_EDGES(1:*) - E_EDGES(0:n_elements(E_EDGES)-2)
	PhotonEnergyBinWidth = C_EDGES(1:*) - C_EDGES(0:n_elements(C_EDGES)-2)
			
	;; Find the channel number that corresponds to the lower and upper detector bandpass
	whereclose, ChannelEnergy, BandPass(0), Clower, /silent
	whereclose, ChannelEnergy, BandPass(1), Cupper, /silent
	
	;; Find the channel number that corresponds to the lower and upper detector bandpass over which the detector triggers!
	whereclose, ChannelEnergy, TriggerBandPass(0), Clower_Trigger, /silent
	whereclose, ChannelEnergy, TriggerBandPass(1), Cupper_Trigger, /silent
	
	;; Find the channel number that corresponds to the lower and upper boundries of BATSE broad energy "channels" 1-4
	whereclose, ChannelEnergy, BatseChannel1BandPass(0), C1lower, /silent
	whereclose, ChannelEnergy, BatseChannel1BandPass(1), C1upper, /silent
	whereclose, ChannelEnergy, BatseChannel2BandPass(0), C2lower, /silent
	whereclose, ChannelEnergy, BatseChannel2BandPass(1), C2upper, /silent	
	whereclose, ChannelEnergy, BatseChannel3BandPass(0), C3lower, /silent
	whereclose, ChannelEnergy, BatseChannel3BandPass(1), C3upper, /silent	
	whereclose, ChannelEnergy, BatseChannel4BandPass(0), C4lower, /silent
	whereclose, ChannelEnergy, BatseChannel4BandPass(1), C4upper, /silent
	
	;; Create an empty CountSpectrumCube
	CountsSpectrumCube = dblarr(n_elements(time),n_elements(ChannelEnergy))
	CountsSpectrumCube_UnNormalized = dblarr(n_elements(time),n_elements(ChannelEnergy))
	
	;; Create an empty time integrated count spectrum
	CountsSpectrum_TimeIntegrated = dblarr(n_elements(ChannelEnergy))
	CountsSpectrum_TimeIntegratedError = dblarr(n_elements(ChannelEnergy))
	
	;; Create an empty time integrated background count spectrum
	BackgroundCountSpectrum_PreTrigger_TimeIntegrated = dblarr(n_elements(ChannelEnergy))
	BackgroundCountSpectrum_PreTrigger_TimeIntegratedError = dblarr(n_elements(ChannelEnergy))

	;; Extracting a real background count spectrum
	BFITSData = mrdfits(BFITSFilename, 2, hdr2)

	IF KEYWORD_SET(BATSE) THEN BEGIN
		Rates = BFITSData.Rates
	ENDIF 
	IF KEYWORD_SET(GBM) THEN BEGIN
		Counts = BFITSData.Counts		
		Exposure = BFITSData.Exposure
		Rates = Counts
		for i=0, n_elements(Rates(0,*))-1 do Rates(*,i) = Rates(*,i)/Exposure(i)		
	ENDIF
	
	BackgroundCountSpectrum_Template = fltarr(n_elements(Rates(*,0)))
	BackgroundCountSpectrumCube = fltarr(n_elements(time), n_elements(Rates(*,0))) 
	BackgroundCountSpectrumCube_PreTrigger = fltarr(n_elements(time_pretrigger), n_elements(Rates(*,0))) 
	
	;; Average over ten time bins to get an time averaged background count spectrum
	BackgroundStart = 0
	BackgroundStop = 20
	FOR i=0D, n_elements(BackgroundCountSpectrum_Template)-1 DO BackgroundCountSpectrum_Template(i) = median(Rates(i,BackgroundStart:BackgroundStop))
	
	;; Normalize the background by the total background counts 
	BackgroundCountSpectrum_TemplateNormalized = BackgroundCountSpectrum_Template/total(BackgroundCountSpectrum_Template)
	
	;; Set a new integrated background level 
	BackgroundCountSpectrum_TemplateAmplified = BackgroundCountSpectrum_TemplateNormalized * BackgroundCountSpectrumNormalization
	
	;; Create a random Poisson distributed background count value for each channel, centered on the counts level in the BackgroundCountSpectrum, for each time bin
	FOR i=1, n_elements(BackgroundCountSpectrum_Template)-1 DO BEGIN
		IF BackgroundCountSpectrum_TemplateAmplified(i) EQ 0 THEN BEGIN
			BackgroundCountSpectrumCube(*,i) = 0.0
		ENDIF ELSE BEGIN
			BackgroundCountSpectrumCube(*,i) = randomn(seed,n_elements(time), POISSON = BackgroundCountSpectrum_TemplateAmplified(i))
			BackgroundCountSpectrumCube_PreTrigger(*,i) = randomn(seed,n_elements(time_pretrigger), POISSON = BackgroundCountSpectrum_TemplateAmplified(i))
		ENDELSE
	ENDFOR

	print, 'Simulating Detector Response...'
	FOR i=0D, N_ELEMENTS(time)-1 DO BEGIN
	
		;; Rebin the simulated photon spectrum g to a 128 energy bin photon spectrum of the same units
;		rebinspectrum, Energy, PhotonSpectrumCube(i,*), Energy_Bin, PhotonSpectrum_bin, E_EDGES(0:n_elements(E_EDGES)-2), E_EDGES(1:*)
;		rebinspectrum, Energy, PhotonSpectrumCube(i,*), Energy_Bin, PhotonSpectrum_bin, C_EDGES(0:n_elements(C_EDGES)-2), C_EDGES(1:*)
		
		;; Rebining the simulated photon spectrum * energy (cm^-2 s^-1) and then converting it back into a photon spectrum (cm^-2 s^-1 keV^-1) after the rebin
		rebinspectrum, Energy_Observer, PhotonSpectrumCube_Observer(i,*) * Energy_Observer, Energy_Observer_Bin, PhotonSpectrum_Observer_Bin, E_EDGES(0:n_elements(E_EDGES)-2), E_EDGES(1:*)
		PhotonSpectrum_Observer_Bin = PhotonSpectrum_Observer_Bin / Energy_Observer_Bin
	
		;; Convolve the data response matrix with the binned photon spectrum to create a counts matrix and a time integrated count rate spectrum (counts s-1 kev-1)
;		MakeCountsDataCube, resp_matrix, PhotonSpectrum_bin, PhotonEnergyBinWidth, ChannelEnergyBinWidth, CountsMatrix, CountSpectrum, /batse

		IF KEYWORD_SET(batse) THEN BEGIN
			MakeCountsDataCube, resp_matrix, PhotonSpectrum_Observer_Bin, Energy_Observer_Bin, ChannelEnergy, CountsMatrix, CountSpectrum, /batse
		ENDIF
		IF KEYWORD_SET(GBM) THEN BEGIN
			MakeCountsDataCube, resp_matrix, PhotonSpectrum_Observer_Bin, Energy_Observer_Bin, ChannelEnergy, CountsMatrix, CountSpectrum
		ENDIF	
		
		;; Renormalizing the count spectrum to make the results reasonable.  **This is a temporary fix!**
		CountSpectrum_UnNormalized = CountSpectrum
		IF (N_ELEMENTS(CountSpectrumNorm) EQ 0) THEN BEGIN
			CountSpectrumNorm = 12
		ENDIF 
		CountSpectrum = CountSpectrum / CountSpectrumNorm

		;; Add a real background count spectrum to the CountSpectrum.
		CountSpectrum = CountSpectrum + BackgroundCountSpectrumCube(i,*)
		
		;; Fill the CountSpectrumCube
		CountsSpectrumCube_UnNormalized(i,*) = CountSpectrum_UnNormalized	
		CountsSpectrumCube(i,*) = CountSpectrum
				
		;; Create a count rate light curve, both over the entire detector bandpass and over just the triggerable bandpass
		CountRate(i) = total(CountSpectrum(Clower:Cupper))
		CountRate_Trigger(i) = total(CountSpectrum(Clower_Trigger:Cupper_Trigger))
		
		;; Creating a count rate light curve for each individual channels
		CountRate_Channel1(i) = total(CountSpectrum(C1lower:C1upper))
		CountRate_Channel2(i) = total(CountSpectrum(C2lower:C2upper))
		CountRate_Channel3(i) = total(CountSpectrum(C3lower:C3upper))
		CountRate_Channel4(i) = total(CountSpectrum(C4lower:C4upper))
		
		;; Saving a bogus k-corrected coutrate light curve
		;CountRate_kCorrected(i) = total(CountSpectrum(Clower:Cupper))
		
	ENDFOR
		
	;; Create a 50s pretrigger background
	FOR i=0D, n_elements(CountRate_PreTrigger)-1 DO BEGIN
		CountRate_PreTrigger(i) = total(BackgroundCountSpectrumCube_PreTrigger(i,Clower:Cupper))
		CountRate_PreTrigger_Trigger(i) = total(BackgroundCountSpectrumCube_PreTrigger(i,Clower_Trigger:Cupper_Trigger))
		CountRate_Channel1_PreTrigger(i) = total(BackgroundCountSpectrumCube_PreTrigger(i,C1lower:C1upper))
		CountRate_Channel2_PreTrigger(i) = total(BackgroundCountSpectrumCube_PreTrigger(i,C2lower:C2upper))
		CountRate_Channel3_PreTrigger(i) = total(BackgroundCountSpectrumCube_PreTrigger(i,C3lower:C3upper))
		CountRate_Channel4_PreTrigger(i) = total(BackgroundCountSpectrumCube_PreTrigger(i,C4lower:C4upper))	
	ENDFOR
	
	

;ENDIF ELSE BEGIN

	;; Create a count rate light curve
;	photon_background = POISSON_Median / count_rsp	;photons
;	CountRate = (PhotonFlux_Detector * count_rsp)
;	CountRate_kCorrected = (PhotonFlux_kCorrected * count_rsp)
	
;ENDELSE


;; Calculating the standard deviation of the background noise 
STDEV_Noise = STDEV(CountRate_PreTrigger)


;; Adding background signal prior to the trigger
;time_pretrig = reverse(time(1:50) * (-1.0))
time = [time_pretrigger,time]
CountRate = [CountRate_PreTrigger,CountRate]
CountRate_Trigger = [CountRate_PreTrigger_Trigger,CountRate_Trigger]
CountRate_Channel1 = [CountRate_Channel1_PreTrigger,CountRate_Channel1]
CountRate_Channel2 = [CountRate_Channel2_PreTrigger,CountRate_Channel2]
CountRate_Channel3 = [CountRate_Channel3_PreTrigger,CountRate_Channel3]
CountRate_Channel4 = [CountRate_Channel4_PreTrigger,CountRate_Channel4]

;; Adding empty arrays to the photon and energy flux light curves (i.e. no "real" signal is coming from the source pre-trigger)
PhotonFlux_Detector = [fltarr(n_elements(TIME_PRETRIGGER)), PhotonFlux_Detector]
EnergyFlux_Detector = [fltarr(n_elements(TIME_PRETRIGGER)), EnergyFlux_Detector]

PhotonFlux_DetectorTrigger = [fltarr(n_elements(TIME_PRETRIGGER)), PhotonFlux_DetectorTrigger]
EnergyFlux_DetectorTrigger = [fltarr(n_elements(TIME_PRETRIGGER)), EnergyFlux_DetectorTrigger]

PhotonFlux_kCorrected = [fltarr(n_elements(TIME_PRETRIGGER)), PhotonFlux_kCorrected]
EnergyFlux_kCorrected = [fltarr(n_elements(TIME_PRETRIGGER)), EnergyFlux_kCorrected]	

PhotonLuminosity_Bolometric = [fltarr(n_elements(TIME_PRETRIGGER)), PhotonLuminosity_Bolometric]
EnergyLuminosity_Bolometric = [fltarr(n_elements(TIME_PRETRIGGER)), EnergyLuminosity_Bolometric]	


;; rebinning if needed
IF KEYWORD_SET(rebin) THEN BEGIN
	rebin, time, timebin, 10
	rebin, PhotonFlux_Detector, PhotonFlux_Detector_bin, 10
	rebin, PhotonFlux_DetectorTrigger, PhotonFlux_DetectorTrigger_bin, 10
	rebin, PhotonFlux_kCorrected, PhotonFlux_kCorrected_bin, 10
	rebin, EnergyFlux_Detector, EnergyFlux_Detector_bin, 10
	rebin, EnergyFlux_DetectorTrigger, EnergyFlux_DetectorTrigger_bin, 10	
	rebin, EnergyFlux_kCorrected, EnergyFlux_kCorrected_bin, 10
	rebin, Epk, Epkbin, 10
	rebin, Flux, Fluxbin, 10
	rebin, Fepk, Fepkbin, 10
	rebin, noise_pretrig, noise_pretrig_bin, 10
	rebin, noise, noise_bin, 10
	rebin, CountRate, CountRate_bin, 10
	rebin, CountRate_Trigger, CountRate_Trigger_bin, 10
	
	time_full = time		
	time = timebin
	PhotonFlux_Detector = PhotonFlux_Detector_bin * 2.0
	PhotonFlux_trigger = PhotonFlux_trigger * 2.0
	PhotonFlux_kcorrected = PhotonFlux_kcorrected_bin * 2.0	
	EnergyFlux = EnergyFlux_bin
	EnergyFlux_kCorrected = EnergyFlux_kCorrected_bin	
	Epk = Epkbin
	flux = fluxbin
	Fepk = Fepkbin
	noise_pretrig = noise_pretrig_bin * 2.0
	CountRate = CountRate_bin * 2.0
	CountRate_Trigger = CountRate_Trigger_bin * 2.0

ENDIF


;; Estimating the background subtracted signal to noise ratio
Signal = max(CountRate)
Background = median(CountRate_PreTrigger)
Signal_Trigger = max(CountRate_Trigger)
Background_Trigger = median(CountRate_PreTrigger_Trigger)

;; Original Definition of signal to noise used in Kocevski 2011
;SNRatio = Signal/STDEV_Noise
;SNRatio_Trigger = Signal_Trigger/STDEV_Noise

;; Gruber et al 20011 definition of signal to noise (Changed 07/02/2013)
SNRatio = (Signal-Background)/SQRT(Background)
SNRatio_Trigger = (Signal_Trigger-Background_Trigger)/SQRT(Background_Trigger)

;; Calculating the burst duration using the running, background substracted, count fluence
t90, time, CountRate-median(CountRate_PreTrigger), t90_count, t90min_count, t90max_count, CountFluence_Running, MedianMaxCount

;; Calculating the burst duration using the running photon fluence
t90, time, PhotonFlux_Detector, t90_photon, t90min_photon, t90max_photon, PhotonFluence_Detector_Running, MedianMaxPhotonFluence
t90, time, PhotonFlux_kCorrected, t90_photon_kCorrected, t90min_photon_kCorrected, t90max_photon_kCorrected, PhotonFluence_kCorrected_Running

;; Calculating the burst duration using the running photon fluence
t90, time, EnergyFlux_Detector, t90_Energy, t90min_Energy, t90max_Energy, EnergyFluence_Detector_Running, MedianMaxEnergyFluence
t90, time, EnergyFlux_kCorrected, t90_Energy_kCorrected, t90min_Energy_kCorrected, t90max_Energy_kCorrected, EnergyFluence_kCorrected_Running


;; Measure the fwhm

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Calculating the t90 using a Baysian Block analysis
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

IF (!version.os EQ 'darwin') THEN BEGIN
	CountRateFilename = strcompress('CountRate_' + JobNumber + '.txt', /remove_all)
	BayesianBlockFilename = strcompress('BayesianBlock_' + JobNumber + '.txt', /remove_all)
	CountRateFilename = CountRateFilename[0]
	BayesianBlockFilename = BayesianBlockFilename[0]
ENDIF ELSE BEGIN
	CountRateFilename = strcompress('PhotonFlux_' + JobNumber + '.txt', /remove_all)
	BayesianBlockFilename = strcompress('BayesianBlock_' + JobNumber + '.txt', /remove_all)
	CountRateFilename = CountRateFilename[0]
	BayesianBlockFilename = BayesianBlockFilename[0]
	
ENDELSE


openw, 1, CountRateFilename
openw, 2, BayesianBlockFilename

FOR i=0D, n_elements(time)-1 DO BEGIN
	;; Converting from a photon light curve count light curve
;	printf, 1, ((PhotonFlux(i) + photon_background) * count_rsp)
	printf, 1, CountRate(i)
	printf, 2, CountRate(i)
	
ENDFOR

close, 1
close, 2

print, "Calling BayesianBlocks.py..."
;spawn, 'python ~/bin/custom/BayesianBlocks_python_idl.py'
;readcol, '/Users/Kocevski/Temp/idl_tmp_bb.txt', time_bb, counts_bb, /silent

;; Calling the BayesianBlocks_python.py.  The script reads in the count rate light curve file and then overwrites this file with the bb light curve
pythoncmd = string('python BayesianBlocks.py' + ' ' + BayesianBlockFilename)

print, pythoncmd
spawn, pythoncmd
	
;; Reading in the resulting Bayesian block reconstruction
readcol, BayesianBlockFilename, time_bb, counts_bb, /silent

;; Accounting for the time resolution and the duration of the pretrigger background
time_bb = (time_bb * timeres) - abs(min(time))

;; Delete the CountRateFilename
rmcmd = string('rm ' + BayesianBlockFilename)
;;spawn, rmcmd

;; Extract the signal that is both at least 5% of the background subtracted max signal and above the median poisson noise plus sqrt(N) 
NoiseCutoff = counts_bb(0) + ( SQRT(counts_bb(0)) * 1.0 )
SignalCutoff = max(counts_bb-counts_bb(0)) * 0.05
above_noise = where(counts_bb-counts_bb(0) GT SignalCutoff AND counts_bb GT NoiseCutoff)
below_noise = where(counts_bb-counts_bb(0) LE SignalCutoff AND counts_bb GT NoiseCutoff)			
			
;; Define T100
IF (N_ELEMENTS(above_noise) GT 1) THEN BEGIN

	;; Define the source interval
	FirstRealChangePoint = above_noise(0)
	LastRealChangePoint = above_noise(n_elements(above_noise)-1)
	T100_min = time_bb(FirstRealChangePoint)
	T100_max = time_bb(LastRealChangePoint)

	IF (T100_max EQ 350.00) THEN BEGIN
		T100 = 0
		T100_min = -1
		T100_max = -1
	ENDIF ELSE BEGIN
		T100 = time_bb(above_noise(N_ELEMENTS(above_noise)-1)) - time_bb(above_noise(0))
		;	print, 'Noise STDEV BB =', STDEV_Noise_pretrig_bb
	ENDELSE
	
ENDIF ELSE BEGIN

	T100 = 0
	T100_min = -1
	T100_max = -1
	
ENDELSE


;; Calculating the burst duration using the running photon fluence
t90, time_bb, counts_bb, t90_bblocks, t90min_bblocks, t90max_bblocks, /bblocks

;; Finding the index of the start and stop values associated with the T100 duration
whereclose, time, T100_min - 100.0 , BurstStartIndex, /silent
whereclose, time, T100_max - 101.0 , BurstStopIndex, /silent
IF (BurstStopIndex EQ 200) THEN BurstStopIndex = 199

;; Calculating the time integrated spectra over the T100 duration
FOR i=0D, n_elements(Energy_Observer)-1 DO BEGIN

	Photon_Fnu_TimeIntegratedSpectrum_Observer(i) = mean(Photon_nuFnu_SpectrumCube_Observer(BurstStartIndex:BurstStopIndex,i))
	Photon_nuFnu_TimeIntegratedSpectrum_Observer(i) = mean(Photon_nuFnu_SpectrumCube_Observer(BurstStartIndex:BurstStopIndex,i))
	Energy_Fnu_TimeIntegratedSpectrum_Observer(i) = mean(Energy_Fnu_SpectrumCube_Observer(BurstStartIndex:BurstStopIndex,i))
	
ENDFOR

;; Calculating the observer frame photon and energy flux from the time integrated band spectrum, but using the source frame energy to do the k-correction
PhotonFlux_TimeIntegrated_Observer = int_tabulated(Energy_Observer(LowerIndex:UpperIndex),Photon_Fnu_TimeIntegratedSpectrum_Observer(LowerIndex:UpperIndex)) 								;; photons cm^-2 s^-1
EnergyFlux_TimeIntegrated_Observer = int_tabulated(Energy_Observer(LowerIndex:UpperIndex),Energy_Fnu_TimeIntegratedSpectrum_Observer(LowerIndex:UpperIndex))								;; ergs cm^-2 s^-1	
PhotonFlux_TimeIntegrated_kCorrected_Observer = int_tabulated(Energy_Source(kLowerIndex:kUpperIndex),Photon_Fnu_TimeIntegratedSpectrum_Observer(kLowerIndex:kUpperIndex)) 		;; photons cm^-2 s^-1
EnergyFlux_TimeIntegrated_kCorrected_Observer = int_tabulated(Energy_Source(kLowerIndex:kUpperIndex),Energy_Fnu_TimeIntegratedSpectrum_Observer(kLowerIndex:kUpperIndex))		;; ergs cm^-2 s^-1	

;; Finding the estimated (i.e. using T100) k-corrected, photon and energy fluence
PhotonFluence_Estimated = PhotonFlux_TimeIntegrated_Observer * T100
EnergyFluence_Estimated = EnergyFlux_TimeIntegrated_Observer * T100
PhotonFluence_kCorrected_Estimated = PhotonFlux_TimeIntegrated_kCorrected_Observer * T100
EnergyFluence_kCorrected_Estimated = EnergyFlux_TimeIntegrated_kCorrected_Observer * T100

;; Finding the estimated Eiso
Eiso_kCorrected_Estimated = EnergyFluence_kCorrected_Estimated * ((4*!PI*(dl^2.0))/(1+redshift))

;; Finding the indices representing the Channel 1 and 3 energy boundries
Channel1LowerEnergy = 20
Channel1UpperEnergy = 50
Channel2LowerEnergy = 50
Channel2UpperEnergy = 100
Channel3LowerEnergy = 100
Channel3UpperEnergy = 300
Channel4LowerEnergy = 300
Channel4UpperEnergy = BandPass[n_elements(BandPass)-1]

whereclose, Energy_Observer, Channel1LowerEnergy, Channel1LowerIndex, /silent
whereclose, Energy_Observer, Channel1UpperEnergy, Channel1UpperIndex, /silent
whereclose, Energy_Observer, Channel3LowerEnergy, Channel3LowerIndex, /silent
whereclose, Energy_Observer, Channel3UpperEnergy, Channel3UpperIndex, /silent	
whereclose, Energy_Observer, Channel4LowerEnergy, Channel4LowerIndex, /silent
whereclose, Energy_Observer, Channel4UpperEnergy, Channel4UpperIndex, /silent

;; Calculating the photon flux in Channels 1 and 3.  This should really be done using the total photons (or counts) in each channel and not the mean photon spectra
PhotonFlux_Channel1TimeIntegrated_Observer = int_tabulated(Energy_Observer(Channel1LowerIndex:Channel1UpperIndex),Photon_Fnu_TimeIntegratedSpectrum_Observer(Channel1LowerIndex:Channel1UpperIndex)) 								;; photons cm^-2 s^-1
PhotonFlux_Channel3TimeIntegrated_Observer = int_tabulated(Energy_Observer(Channel3LowerIndex:Channel3UpperIndex),Photon_Fnu_TimeIntegratedSpectrum_Observer(Channel3LowerIndex:Channel3UpperIndex)) 								;; photons cm^-2 s^-1

;; Calculating the hardness ratio
HR31_Photon = PhotonFlux_Channel3TimeIntegrated_Observer/PhotonFlux_Channel1TimeIntegrated_Observer
HR31_Count = total(CountRate_Channel3) / total(CountRate_Channel1)
HR41_Count = total(CountRate_Channel4) / total(CountRate_Channel1)



;; Calculating the spectroscopic lag
IF KEYWORD_SET(lag) THEN BEGIN
		
	;; Create a lag array
	lag = findgen(10000)/100.
	lag = lag - 50.0
	
	;; Find the cross correlation Pxy(L) of two light curves as a function of the lag
	CCF = c_correlate(CountRate_Channel3, CountRate_Channel1, lag)

	;; Find the peak of the cross correlation function
	PeakCCF = where(ccf eq max(ccf))
	PeakCCF = PeakCCF(0)
	
	BracketIndex = 600
	
	; Performing polynomial fit to determine the exact maximum of the cross correlation function
	xlag = lag(PeakCCF-BracketIndex:PeakCCF+BracketIndex)
	yccf = ccf(PeakCCF-BracketIndex:PeakCCF+BracketIndex)
	fit = POLY_FIT(xlag, yccf, 3)
	fitcurve = fit(0) + fit(1)*xlag + fit(2)*xlag^2.0 + fit(3)*xlag^3.0
	maxCCFIndex = where(fitcurve eq max(fitcurve))
	maxCCFIndex = maxCCFIndex + (PeakCCF-BracketIndex)
	Lag31 = lag(maxCCFIndex)

	IF KEYWORD_SET(ShowLagPlot) THEN BEGIN
		window, 5
		plot, xlag, yccf, xtitle='Lag (seconds)', ytitle = 'CCF'
		oplot, xlag, fitcurve, color=150, linestyle=1
	ENDIF
	
ENDIF ELSE BEGIN

	Lag31 = float('NaN')

ENDELSE


;; Recording the time integrated Epk only if beta is less than -2.0
IF (beta LT -2.0) THEN BEGIN
	
	;; Find the peak of the time integrated nuFnu spectrum in the observer frame
	Epk_Observer_Int = Energy_Observer(wheremax(Photon_nuFnu_TimeIntegratedSpectrum_Observer))
	
	IF (N_ELEMENTS(Epk_Observer_Int) GT 1) THEN BEGIN
		Epk_Observer_Int = float('NaN')
	ENDIF	
ENDIF ELSE BEGIN
	Epk_Observer_Int = float('NaN')
ENDELSE

;; Write out a fake time resolved BFITS file for use in rmfit	
IF KEYWORD_SET(MakeFakeBFITS) THEN BEGIN

	print, 'Template BFITS Filen: ', BFITSFilename
	print, 'New BFITS File: ', BFITSFilenameNew	
	
	;; Read in the data
	Data0 = mrdfits(BFITSFilename,0,hdr0)
	Data1 = mrdfits(BFITSFilename,1,hdr1)
	Data2 = mrdfits(BFITSFilename,2,hdr2)
	
	;; Create a BATSE time array to go with the CountsSpectrumCube
	times = fltarr(2,n_elements(CountsSpectrumCube(*,0)))
	times(0,*) = findgen(n_elements(CountsSpectrumCube(*,0)))
	times(1,*) = times(0,*) + 1.0
	
	;; Creating and filling a structure to hold all the parameter distributions
	Data2New = {times:fltarr(2,n_elements(CountsSpectrumCube(*,0))), rates:fltarr(n_elements(CountsSpectrumCube(0,*)),n_elements(CountsSpectrumCube(*,0))), errors:fltarr(n_elements(CountsSpectrumCube(0,*)),n_elements(CountsSpectrumCube(*,0))) }
	
	Data2New.times = times
	Data2New.rates = transpose(CountsSpectrumCube)
;		Data2New.errors = (SQRT(transpose(CountsSpectrumCube_UnNormalized))/transpose(CountsSpectrumCube_UnNormalized)) * transpose(CountsSpectrumCube)
	Data2New.errors = (transpose(SQRT(CountsSpectrumCube)))		
;		Data2New.errors(*) = 1.0
	
	mwrfits, Data0, BFITSFilenameNew, hdr0, /Create
	mwrfits, Data1, BFITSFilenameNew, hdr1
	mwrfits, Data2New, BFITSFilenameNew, hdr2

		
ENDIF
	
	
;; Write out a fake PHA1 file for use in xspec	
IF KEYWORD_SET(MakeFakePHA) THEN BEGIN
	
	print, 'Template PHA File: ', PHA1Filename
	print, 'New PHA File: ', PHA1FilenameNew	
	
	;; Read in the data
	Data0 = mrdfits(PHA1Filename,0,hdr0)
	Data1 = mrdfits(PHA1Filename,1,hdr1)
	Data2 = mrdfits(PHA1Filename,2,hdr2)
	Data3 = mrdfits(PHA1Filename,3,hdr3)

	;; Make a time integrated count spectrum for the entire burst
;		FOR i=0, n_elements(CountsSpectrumCube(*,0))-1 DO BEGIN
;			CountsSpectrum_TimeIntegrated = CountsSpectrum_TimeIntegrated + CountsSpectrumCube(i,*)
;		ENDFOR

	;; Using the Bayesian Block determined duration as the time over which to integrate the source spectra
	IntegrationTime = BurstStopIndex		;; Sec
	
	;; Producing a time integrated count spectrum
	FOR i=0D, n_elements(CountsSpectrumCube(0,*))-1 DO BEGIN
	
		;; Note that taking the mean is equivilent to finding the rate because we are using 1 second resolution spectra.  Otherwise we would have to multiply by the exposure of each bin, then add, then devide by the total exposure
		CountsSpectrum_TimeIntegrated(i) = mean(CountsSpectrumCube(0:IntegrationTime,i))

		;; Finding the error on the integrated spectrum.  Doing this by taking the SQRT of the total and then finding the fractional error with the eventual rate (i.e. the mean as defined above)
		CountsSpectrum_TimeIntegratedError(i) = ( SQRT(total(CountsSpectrumCube(0:IntegrationTime,i))) / total(CountsSpectrumCube(0:IntegrationTime,i)) ) * CountsSpectrum_TimeIntegrated(i)

;			CountsSpectrum_TimeIntegrated_SumError = sqrt( total( (sqrt( CountsSpectrumCube(*,i) )/ CountsSpectrumCube(*,i) )^2.0)) * total(CountsSpectrumCube(*,i))			
;			CountsSpectrum_TimeIntegratedError(i) = ( CountsSpectrum_TimeIntegrated_SumError /  total(CountsSpectrumCube(*,i)) ) * CountsSpectrum_TimeIntegrated(i)
;			CountsSpectrum_TimeIntegratedError(i) = CountsSpectrum_TimeIntegratedError(i) / 2.0
	ENDFOR
	
	
;		CountsSpectrum_TimeIntegrated(i) = mean(Photon_nuFnu_SpectrumCube(BurstStartIndex:BurstStopIndex,i))

	;; Creating a new pha data structure for extension 1
	Data1New = {CHANNEL:intarr(n_elements(Data1.Channel)), E_MIN:fltarr(n_elements(Data1.E_MIN)), E_MAX:fltarr(n_elements(Data1.E_MAX)) }
	Data1New.CHANNEL = Data1.Channel
	Data1New.E_MIN = C_EDGES(0:127)
	Data1New.E_MAX = C_EDGES(1:128)	
	
	;; Creating a new pha data structure for extension 2
	Data2New = {CHANNEL:intarr(n_elements(Data2.Channel)), RATE:fltarr(n_elements(CountsSpectrum_TimeIntegrated)), STAT_ERR:fltarr(n_elements(CountsSpectrum_TimeIntegrated)) }
	Data2New.CHANNEL = Data2.Channel
	Data2New.RATE = CountsSpectrum_TimeIntegrated
	Data2New.STAT_ERR = CountsSpectrum_TimeIntegratedError	

	sxaddpar, hdr0, 'TRIGTIME', 0.0
	sxaddpar, hdr0, 'TSTART', 0.0
	sxaddpar, hdr0, 'TSTOP', n_elements(CountsSpectrumCube(0:IntegrationTime,0))

	sxaddpar, hdr1, 'TRIGTIME', 0.0
	sxaddpar, hdr1, 'TSTART', 0.0
	sxaddpar, hdr1, 'TSTOP', n_elements(CountsSpectrumCube(0:IntegrationTime,0))
	
	sxaddpar, hdr2, 'TRIGTIME', 0.0
	sxaddpar, hdr2, 'TSTART', 0.0
	sxaddpar, hdr2, 'TSTOP', n_elements(CountsSpectrumCube(0:IntegrationTime,0))
	sxaddpar, hdr2, 'TTYPE3', 'POISSON'
	
	sxaddpar, hdr3, 'TRIGTIME', 0.0
	sxaddpar, hdr3, 'TSTART', 0.0
	sxaddpar, hdr3, 'TSTOP', n_elements(CountsSpectrumCube(0:IntegrationTime,0))

	sxaddpar, hdr3, 'EXPOSURE', n_elements(CountsSpectrumCube(0:IntegrationTime,0))
	sxaddpar, hdr3, 'TZERO1', 0.0
	sxaddpar, hdr3, 'TZERO2', n_elements(CountsSpectrumCube(0:IntegrationTime,0))
	
	mwrfits, Data0, PHA1FilenameNew, hdr0, /Create
	mwrfits, Data1New, PHA1FilenameNew, hdr1
	mwrfits, Data2New, PHA1FilenameNew, hdr2
	mwrfits, Data3, PHA1FilenameNew, hdr3
	
	
	;; Make a time integrated count spectrum just for the background
	FOR i=0D, n_elements(BackgroundCountSpectrumCube_PreTrigger(0,*))-1 DO BEGIN
		BackgroundCountSpectrum_PreTrigger_TimeIntegrated(i) = mean(BackgroundCountSpectrumCube_PreTrigger(*,i))
		BackgroundCountSpectrum_PreTrigger_TimeIntegratedError(i) = ( SQRT(total(BackgroundCountSpectrumCube_PreTrigger(*,i))) / total(BackgroundCountSpectrumCube_PreTrigger(*,i)) ) * BackgroundCountSpectrum_PreTrigger_TimeIntegrated(i)

		;; XSPEC doesn't like -NaN values, so setting those values equal to an artifically small number.  It is recommended that these channels be ignored
		if (BackgroundCountSpectrum_PreTrigger_TimeIntegrated(i) EQ 0) then BackgroundCountSpectrum_PreTrigger_TimeIntegrated(i) = 1e-10
		if (finite(BackgroundCountSpectrum_PreTrigger_TimeIntegratedError(i)) EQ 0) then BackgroundCountSpectrum_PreTrigger_TimeIntegratedError(i) = 1e-10

	ENDFOR

	;; Creating a new pha data structure for extension 1
	BackgroundData1New = {CHANNEL:intarr(n_elements(Data1.Channel)), E_MIN:fltarr(n_elements(Data1.E_MIN)), E_MAX:fltarr(n_elements(Data1.E_MAX)) }
	BackgroundData1New.CHANNEL = Data1.Channel
	BackgroundData1New.E_MIN = C_EDGES(0:127)
	BackgroundData1New.E_MAX = C_EDGES(1:128)	
	
	;; Creating a new pha data structure for extension 2
	BackgroundData2New = {CHANNEL:intarr(n_elements(Data2.Channel)), RATE:fltarr(n_elements(BackgroundCountSpectrum_PreTrigger_TimeIntegrated)), STAT_ERR:fltarr(n_elements(BackgroundCountSpectrum_PreTrigger_TimeIntegrated)) }
	BackgroundData2New.CHANNEL = Data2.Channel
	BackgroundData2New.RATE = BackgroundCountSpectrum_PreTrigger_TimeIntegrated
	BackgroundData2New.STAT_ERR = BackgroundCountSpectrum_PreTrigger_TimeIntegratedError	

	sxaddpar, hdr0, 'TRIGTIME', 0.0
	sxaddpar, hdr0, 'TSTART', 0.0
;;		sxaddpar, hdr0, 'TSTOP', n_elements(BackgroundCountSpectrumCube_PreTrigger(*,0))
	sxaddpar, hdr0, 'TSTOP', n_elements(BackgroundCountSpectrumCube_PreTrigger(*,0))

	sxaddpar, hdr1, 'TRIGTIME', 0.0
	sxaddpar, hdr1, 'TSTART', 0.0
;;		sxaddpar, hdr1, 'TSTOP',  n_elements(BackgroundCountSpectrumCube_PreTrigger(*,0))
	sxaddpar, hdr1, 'TSTOP',  n_elements(BackgroundCountSpectrumCube_PreTrigger(*,0))
	
	sxaddpar, hdr2, 'TRIGTIME', 0.0
	sxaddpar, hdr2, 'TSTART', 0.0
;;		sxaddpar, hdr2, 'TSTOP',  n_elements(BackgroundCountSpectrumCube_PreTrigger(*,0))
	sxaddpar, hdr2, 'TSTOP',  n_elements(BackgroundCountSpectrumCube_PreTrigger(*,0))
	
	sxaddpar, hdr3, 'TRIGTIME', 0.0
	sxaddpar, hdr3, 'TSTART', 0.0
;;		sxaddpar, hdr3, 'TSTOP',  n_elements(BackgroundCountSpectrumCube_PreTrigger(*,0))
	sxaddpar, hdr3, 'TSTOP',  n_elements(BackgroundCountSpectrumCube_PreTrigger(*,0))

;;		sxaddpar, hdr3, 'EXPOSURE',  n_elements(BackgroundCountSpectrumCube_PreTrigger(*,0))
	sxaddpar, hdr3, 'EXPOSURE',  n_elements(BackgroundCountSpectrumCube_PreTrigger(*,0))
	sxaddpar, hdr3, 'TZERO1', 0.0
;;		sxaddpar, hdr3, 'TZERO2', n_elements(BackgroundCountSpectrumCube_PreTrigger(*,0))
	sxaddpar, hdr3, 'TZERO2', n_elements(BackgroundCountSpectrumCube_PreTrigger(*,0))
	

	mwrfits, Data0, BackgroundPHA1FilenameNew, hdr0, /Create
	mwrfits, BackgroundData1New, BackgroundPHA1FilenameNew, hdr1
	mwrfits, BackgroundData2New, BackgroundPHA1FilenameNew, hdr2
	mwrfits, Data3, BackgroundPHA1FilenameNew, hdr3
		
ENDIF
	


;; Printing diagnostics
print, ''
print, "Redshift = ", redshift
print, "Tmax_Observer = ", Tmax_Observer
print, "Epk0_Source =", Epk0_Source
print, "Epk0_Observer =", Epk0_Source / (1+redshift)
print, "Epk_Observer =", Epk_Observer_Int
print, "Alpha = ", alpha
print, "Beta = ", beta
print, "Initial Photon Flux =", Flux0
print, "Peak Photon Flux =", max(PhotonFlux_Detector)
print, "Peak Photon Flux (50-300 keV) =", max(PhotonFlux_DetectorTrigger)
print, "Peak Energy Flux =", max(EnergyFlux_Detector)
print, "Peak Energy Flux (50-300 keV) =", max(EnergyFlux_DetectorTrigger)
print, "Peak Photon Luminosity (Input) =", Luminosity0
print, "Peak Photon Luminosity =", max(PhotonFlux_Detector) * (4*!PI*dl^2.0)
print, "Peak Energy Luminosity =", max(EnergyFlux_Detector) * (4*!PI*dl^2.0)
print, "Energy Fluence True (Bolometric) =", EnergyFluence_Bolometric_True
print, "Energy Fluence True (k-Corrected) =", EnergyFluence_kCorrected_True
print, "Energy Fluence Estimated (k-Corrected) =", EnergyFluence_kCorrected_Estimated
print, "Energy Fluence Estimated =", EnergyFluence_Estimated
print, "Eiso True (Bolometric) =", Eiso_Bolometric_True
;print, "Eiso True (k-Corrected) =", Eiso_kCorrected_True
print, "Eiso Estimated (k-Corrected) =", Eiso_kCorrected_Estimated
print, "HR31 Photons = ", HR31_Photon
print, "HR31 Counts = ", HR31_Count
print, "HR41 Counts = ", HR41_Count
print, 'Lag CCF31 = ', Lag31
print, "Noise STDEV = ", STDEV_Noise
print, 'Peak Counts = ', max(CountRate)
print, 'Peak Counts (50-300 keV) = ', max(CountRate_Trigger)
print, "S/N ratio =", SNRatio
print, "S/N ratio (50-300 keV) =", SNRatio_Trigger
print, "T90 Photon =", t90_photon
print, "T90 Count =", t90_count
print, "T90 Count BBlocks =", t90_bblocks
print, "T100 =", T100
print, ''





;; Plotting the results
IF (KEYWORD_SET(plot) or KEYWORD_SET(ps)) THEN BEGIN

	;; Plotting the Bayesian Blocks Reconstruction
	IF KEYWORD_SET(ps) THEN BEGIN
		psopen, 'CountLightCurve'
	ENDIF ELSE BEGIN
		window, 2
		loadct, 39, /silent
	ENDELSE 
	IF N_ELEMENTS(_extra) EQ 0 THEN BEGIN
		plot,  time, CountRate , psym=10, yrange=[min(CountRate),max(CountRate)], xtitle='Time (sec)', ytitle='Counts', xrange=[-50,300], xstyle=1
	ENDIF ELSE BEGIN
		plot,  time, CountRate , psym=10, xtitle='Time (sec)', ytitle='Counts', xstyle=1, _extra=_extra
	ENDELSE

	oplot, time_bb, counts_bb, color=150
	oplot, [T100_min,T100_min], [-1e10,1e10], color=150, linestyle=1
	oplot, [T100_max,T100_max], [-1e10,1e10], color=150, linestyle=1
	
	IF KEYWORD_SET(ps) THEN ps_close
	
	;; Plot the photon lightcurve
	IF KEYWORD_SET(ps) THEN BEGIN
		psopen, 'PhotonFluxVsTime'
	ENDIF ELSE BEGIN
		window, 0
		!p.multi=[0,2,2]
		device, decompose = 0
		loadct, 39, /silent
	ENDELSE
	
		IF N_ELEMENTS(_extra) EQ 0 THEN BEGIN
			plot, time, PhotonFlux_Detector, xtitle='Time (sec)', ytitle=textoidl('Photon Flux (photons cm^{-2} s^{-1})'), psym=10, yrange=[-0.1,max(PhotonFlux_Detector)], xrange=[-50,300], xstyle=1, ystyle=1
		ENDIF ELSE BEGIN
			plot, time, PhotonFlux_Detector, xtitle='Time (sec)', ytitle=textoidl('Photon Flux (photons cm^{-2} s^{-1})'), psym=10, yrange=[-0.1,max(PhotonFlux_Detector)], xstyle=1, ystyle=1
		ENDELSE		
	;	oplot, time, (PhotonFlux_kcorrected/max(PhotonFlux_kcorrected))*max(PhotonFlux), psym=10, color=75, linestyle=1

		; Plotting the count fluence derived T90
		oplot, [t90min_count,t90min_count],[-1e10,1e10], linestyle=1, color=75
		oplot, [t90max_count,t90max_count],[-1e10,1e10], linestyle=1, color=75
		
		; Plotting the photon fluence derived T90
		oplot, [t90min_photon,t90min_photon],[-1e10,1e10], linestyle=1, color=190
		oplot, [t90max_photon,t90max_photon],[-1e10,1e10], linestyle=1, color=190
		
		oplot, [T100_min,T100_min]-50, [-1e10,1e10], color=150, linestyle=1
		oplot, [T100_max,T100_max]-50, [-1e10,1e10], color=150, linestyle=1
	
	IF KEYWORD_SET(ps) THEN ps_close
	
	
	;; Plot the Epk vs Time	
	IF KEYWORD_SET(ps) THEN BEGIN
		psopen, 'EpkVsTime'
	ENDIF	
	
		whereclose, time, 0, i, /silent
		plot, time(i:*), Epk_Observer, xtitle='Time (sec)', ytitle='Epk (keV)', /ylog, /xlog, xrange=[0.1,timerange[1]]
;		oplot, [-1E15,1E15], [BandPass(1)+1,BandPass(1)+1], linestyle=1, color=150
;		oplot, [-1E15,1E15], [BandPass(0)+1,BandPass(0)+1], linestyle=1, color=150
	
	IF KEYWORD_SET(ps) THEN ps_close

	
	IF KEYWORD_SET(ps) THEN BEGIN
		psopen, 'nuFnuVsTime'
	ENDIF
	
		if (n_elements(yrange_nufnu) eq 0) then begin
			yrange_nufnu = [min(Photon_nuFnu_SpectrumCube_Observer(fix(Tmax_Observer*1),*)),max(Photon_nuFnu_SpectrumCube_Observer)]
		endif
	
		;; Plot the evolution of the spectral model
		plot, Energy_Observer, Photon_nuFnu_SpectrumCube_Observer(fix(Tmax_Observer*1),*), /ylog, /xlog, xtitle='Energy (keV)', ytitle=textoidl('\nu F \nu (photons keV cm^{-2} s^{-1})'), yrange=yrange_nufnu
			oplot, Energy_Observer, Photon_nuFnu_SpectrumCube_Observer(20,*)
			oplot, Energy_Observer, Photon_nuFnu_SpectrumCube_Observer(30,*)
			oplot, Energy_Observer, Photon_nuFnu_SpectrumCube_Observer(40,*)	
			oplot, Energy_Observer, Photon_nuFnu_SpectrumCube_Observer(50,*)	
			oplot, Energy_Observer, Photon_nuFnu_SpectrumCube_Observer(60,*)	
			oplot, Energy_Observer, Photon_nuFnu_SpectrumCube_Observer(70,*)		
			oplot, Energy_Observer, Photon_nuFnu_SpectrumCube_Observer(80,*)
			oplot, Energy_Observer, Photon_nuFnu_SpectrumCube_Observer(90,*)
			oplot, Energy_Observer, Photon_nuFnu_SpectrumCube_Observer(99,*)
			oplot, Energy_Observer, (Photon_nuFnu_TimeIntegratedSpectrum_Observer/max(Photon_nuFnu_TimeIntegratedSpectrum_Observer)) * max(Photon_nuFnu_TimeIntegratedSpectrum_Observer) * 2, color=254, linestyle=1
			oplot, [BandPass(0)+1,BandPass(0)+1], [1E-5,1E10], linestyle=1, color=150
			oplot, [BandPass(1)+1,BandPass(1)+1], [1E-5,1E10], linestyle=1, color=150


	IF KEYWORD_SET(ps) THEN ps_close
	
	IF KEYWORD_SET(ps) THEN BEGIN
		psopen, 'EpkVsFlux'
	ENDIF
	
		;; Plotting Epk vs Energy Flux @ Epk
		plot, Epk_Observer, Flux, /xlog, /ylog, xrange = [1,10000], yrange=[min(Flux(1:*)),max(Flux)], ystyle=1, ytitle='Energy Flux @ Epk', xtitle='Epk (keV)'
		oplot, [BandPass(1)+1,BandPass(1)+1], [0.001,1E15], linestyle=1, color=150
		oplot, [BandPass(0)+1,BandPass(0)+1], [0.001,1E15], linestyle=1, color=150


	IF KEYWORD_SET(ps) THEN ps_close

	IF KEYWORD_SET(ps) THEN BEGIN
		psopen, 'PhotonFluenceVsTime'
	ENDIF ELSE BEGIN
		window, 1
	ENDELSE
	
		;; finding the index of t = 0
		whereclose, time, 0, zerotime, /silent
	
		;; Plotting time vs energy flux @ Epk
	;	plot, time, EnergyFluence_Detector_Running, ytitle=textoidl('Energy Fluence (ergs cm^{-2}'), xtitle='Time (sec)'
		plot, time, PhotonFluence_Detector_Running, ytitle=textoidl('Photon Fluence (photons cm^{-2}'), xtitle='Time (sec)', yrange=[-1,max(PhotonFluence_Detector_Running)]
		
		; Plotting the count fluence derived T90
		oplot, [T90Min_count,T90Min_count],[-1e10,1e10], linestyle=1, color=75
		oplot, [T90Max_count,T90Max_count],[-1e10,1e10], linestyle=1, color=75
		
		; Plotting the photon fluence derived T90
		oplot, [T90Min_Photon,T90Min_Photon], [-1D10,1D10], color=190, linestyle=1
		oplot, [T90Max_Photon,T90Max_Photon], [-1D10,1D10], color=190, linestyle=1
		
		; Plotting the Bayesian Block derived T100
		oplot, [T100_Min - 50, T100_Min - 50], [-1D10,1D10], color=150, linestyle=1
		oplot, [T100_Max - 50, T100_Max - 50], [-1D10,1D10], color=150, linestyle=1 
		oplot, [-1000,1000],[MedianMaxPhotonFluence,MedianMaxPhotonFluence], linestyle=1, color=254

	IF KEYWORD_SET(ps) THEN ps_close


	IF KEYWORD_SET(ps) THEN BEGIN
		psopen, 'CountFluenceVsTime'
	ENDIF ELSE BEGIN
	ENDELSE
	
		;; finding the index of t = 0
		whereclose, time, 0, zerotime, /silent
	
		;; Plotting time vs energy flux @ Epk
	;	plot, time, energy_fluence, ytitle=textoidl('Energy Fluence (ergs cm^{-2}'), xtitle='Time (sec)'
	;	plot, time, PhotonFluence_Running, ytitle=textoidl('Photon Fluence (photons cm^{-2}'), xtitle='Time (sec)', yrange=[-1,max(PhotonFluence_Running)]
		plot, time, CountFluence_Running, ytitle=textoidl('Count Fluence (photons cm^{-2}'), xtitle='Time (sec)', yrange=[-1,max(CountFluence_Running)]
		
		; Plotting the count fluence derived T90
		oplot, [T90Min_Count,T90Min_Count],[-1e10,1e10], linestyle=1, color=75
		oplot, [T90Max_Count,T90Max_Count],[-1e10,1e10], linestyle=1, color=75
		
		; Plotting the photon fluence derived T90
		oplot, [T90Min_Photon,T90Min_Photon], [-1D10,1D10], color=190, linestyle=1
		oplot, [T90Max_Photon,t90Max_Photon], [-1D10,1D10], color=190, linestyle=1
		
		; Plotting the Bayesian Block derived T100
		oplot, [T100_Min - 50, T100_Min - 50], [-1D10,1D10], color=150, linestyle=1
		oplot, [T100_Max - 50, T100_Max - 50], [-1D10,1D10], color=150, linestyle=1 
		oplot, [-1000,1000],[MedianMaxCount,MedianMaxCount], linestyle=1, color=254

	IF KEYWORD_SET(ps) THEN ps_close
	
	
;	IF KEYWORD_SET(ps) THEN BEGIN
;		psopen, 'EnergyFluxVsTime'
;	ENDIF 
;	
;		;; Plotting time vs energy flux @ Epk
;		plot, time(zerotime:*), flux, ytitle=textoidl('Photon Spectrum Amplitude'), xtitle='Time (sec)', /ylog, /xlog, xrange=[min(time),max(time)], yrange=[min(flux), max(flux)]
;	
;	IF KEYWORD_SET(ps) THEN ps_close

	IF KEYWORD_SET(ps) THEN BEGIN
		psopen, 'PhotonFluenceVsEpk'
	ENDIF 
	
		;; Plotting Energy Fluence vs. Epk
		plot, PhotonFluence_Detector_Running(zerotime:*), Epk_Observer, xrange=[10,max(PhotonFluence_Detector_Running)], xtitle='Photon Fluence', ytitle='Epk (keV)', /ylog
		oplot, PhotonFluence_kCorrected_Running(zerotime:*), Epk_Observer, color=75, linestyle=1
	;	oplot, [-1E15,1E15], [BandPass(1)+1,BandPass(1)+1], linestyle=1, color=150
		oplot, [-1E15,1E15], [100,100], linestyle=1, color=150
		oplot, [-1E15,1E15], [BandPass(0)+1,BandPass(0)+1], linestyle=1, color=150

	IF KEYWORD_SET(ps) THEN ps_close

	IF KEYWORD_SET(ps) THEN BEGIN
		psopen, '3DModel'
	ENDIF 
	
		;; Plot a 3D representation of the GRB model
		;shade_surf, model, az=-75, /ylog, yrange=[1,10000],/zlog, zrange=[min(Fpk(1:*)),max(Fpk)]
		;show3d, model, E_SURFACE={az:-75}, E_CONTOUR={nlevels:8} 	
		shade_surf, Photon_nuFnu_SpectrumCube_Observer(0:long(t90_photon),0:1000), az=-65
		
	
	IF KEYWORD_SET(ps) THEN ps_close

ENDIF

IF KEYWORD_SET(save3d) THEN BEGIN

		!p.multi=[0,0,1]
		loadct, 39, /silent
	;	psopen, 'GRBmodel3D'
		ps_open, 'GRBmodel3D', thick=3, /color, /ps_font

		;; Plot a 3D shaded surface
		shade_surf, Photon_nuFnu_SpectrumCube_Observer(0:50,0:10000), az=-65

		ps_close

ENDIF

;; Analyze the results with xspec
IF KEYWORD_SET(XSPEC) THEN BEGIN
	
	;; Make the xspec script file
	IF (!version.os EQ 'darwin') THEN BEGIN
		XSPECScript = strcompress('GRBModel' + JobNumber + '.tcl', /remove_all) 
	ENDIF ELSE BEGIN
		XSPECScript = strcompress('/nfs/slac/g/ki/ki08/kocevski/Temp/GRBModel' + JobNumber + '.tcl', /remove_all) 
	ENDELSE 

	Openw, 1, XSPECScript
	
	printf, 1, 'data ' + PHA1FilenameNew + '{1}'
	printf, 1, 'backgrnd ' + BackgroundPHA1FilenameNew + '{1}'
	printf, 1, 'response ' + DRM_PHA1Filename
	printf, 1, 'model grbm'
	printf, 1, '-1'
	printf, 1, '-2'
	printf, 1, '150'
	printf, 1, '1'
	printf, 1, 'ignore 1:1-6 1:121-127'
	printf, 1, 'fit'
	printf, 1, 'exit'
	close, 1

	IF (!version.os EQ 'darwin') THEN BEGIN
		SpawnCommand = 'xspec - ' + XSPECScript
	ENDIF ELSE BEGIN
		SpawnCommand = 'xspec.csh ' + XSPECScript
	ENDELSE
	
	print, ''
	print, 'Calling XSPEC:'
	print, SpawnCommand
	spawn, SpawnCommand

ENDIF

;; Cleanup generated data products
IF KEYWORD_SET(CleanUp) THEN BEGIN

	cmd1 = 'rm ' + PHA1FilenameNew
	cmd2 = 'rm ' + BackgroundPHA1FilenameNew
	cmd3 = 'rm ' + XSPECScript
	
	print, ''
	print, 'Cleaning up...'
	spawn, cmd1
	spawn, cmd2
	spawn, cmd3
	print, 'Done.'
	
ENDIF

IF KEYWORD_SET(Results) THEN BEGIN

	;; Define a new structure to house the data
	Results = {PhotonLuminosity0:dblarr(n_elements(PhotonLuminosity0)), PhotonLuminosity:dblarr(n_elements(PhotonLuminosity)), EnergyLuminosity:dblarr(n_elements(EnergyLuminosity)), Epk0_Source:fltarr(n_elements(Epk0_Source)), Epk0_Observer:fltarr(n_elements(Epk0_Observer)), Epk_Observer:fltarr(n_elements(Epk_Observer)), Eiso_Bolometric_True:dblarr(n_elements(Eiso_Bolometric_True)), Eiso_kCorrected_Estimated:dblarr(n_elements(Eiso_kCorrected_Estimated)), T90Photon:fltarr(n_elements(T90Photon)), T90Count:fltarr(n_elements(T90Count)), T100:fltarr(n_elements(T100)), PeakPhotonFlux:fltarr(n_elements(PeakPhotonFlux)), PeakPhotonFluxTrigger:fltarr(n_elements(PeakPhotonFluxTrigger)), PeakEnergyFlux:fltarr(n_elements(PeakEnergyFlux)), PeakEnergyFluxTrigger:fltarr(n_elements(PeakEnergyFluxTrigger)), HR31_photon:fltarr(n_elements(HR31_photon)), HR31_count:fltarr(n_elements(HR31_count)), HR41_count:fltarr(n_elements(HR41_count)), Lag31:fltarr(n_elements(Lag31)), EnergyFluence_Bolometric_True:fltarr(n_elements(EnergyFluence_Bolometric_True)), EnergyFluence_kCorrected_True:fltarr(n_elements(EnergyFluence_kCorrected_True)), EnergyFluence_kCorrected_Estimated:fltarr(n_elements(EnergyFluence_kCorrected_Estimated)), SNRatio:fltarr(n_elements(SNRatio)), SNRatioTrigger:fltarr(n_elements(SNRatioTrigger))   }

	;; Filling the structure
	Results.PhotonLuminosity0 = PhotonLuminosity0
	Results.PhotonLuminosity = PhotonLuminosity
	Results.EnergyLuminosity = EnergyLuminosity
	Results.Epk0_Source = Epk0_Source
	Results.Epk0_Observer = Epk0_Observer 
	Results.Epk_Observer = Epk_Observer
	Results.Eiso_Bolometric_True = Eiso_Bolometric_True
	Results.Eiso_kCorrected_Estimated = Eiso_kCorrected_Estimated
	Results.T90Photon = T90Photon  
	Results.T90Count = T90Count  
	Results.T100 = T100    
	Results.PeakPhotonFlux = PeakPhotonFlux
	Results.PeakPhotonFluxTrigger = PeakPhotonFluxTrigger
	Results.PeakEnergyFlux = PeakEnergyFlux
	Results.PeakEnergyFluxTrigger = PeakEnergyFluxTrigger
	Results.HR31_photon = HR31_photon
	Results.HR31_count = HR31_count
	Results.HR41_count = HR41_count	
	Results.Lag31 = Lag31 	 
	Results.EnergyFluence_Bolometric_True = EnergyFluence_Bolometric_True
	Results.EnergyFluence_kCorrected_True = EnergyFluence_kCorrected_True
	Results.EnergyFluence_kCorrected_Estimated = EnergyFluence_kCorrected_Estimated	
	Results.SNRatio = SNRatio
	Results.SNRatioTrigger = SNRatioTrigger
	Results.Eiso_Bolometric_True = Eiso_Bolometric_True
	Results.Eiso_kCorrected_Estimated = Eiso_kCorrected_Estimated
	
	
ENDIF


!p.multi=[0,0,1]
	
	
RETURN
END