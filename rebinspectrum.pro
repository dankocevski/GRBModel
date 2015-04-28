;; ==========================================================================================
;; File:  rebinspectrum.pro
;; Author:  Daniel Kocevski <kocevski@berkeley.edu>
;; Date Created:  circa 2002
;; Last Changed:  -
;; 
;; Purpose:  This program rebins a spectrum based on user defined energy boundries
;; Dependancies:  None
;; ==========================================================================================


pro rebinspectrum, energy, spectrum, energy_bin, newspectrum, minboundry, maxboundry

IF (n_params() EQ 0) THEN BEGIN 
	print,'-Usage: rebinspectrum, energy, spectrum, newspectrum, minboundry, maxboundry
    RETURN
ENDIF 

;; Using the center of the energy bin to create a binned energy array
energy_bin = ( minboundry + maxboundry ) / 2.0
;energy_bin = SQRT (( minboundry * maxboundry ))

newspectrum = fltarr(n_elements(minboundry))

FOR i=0, N_ELEMENTS(newspectrum)-1 DO BEGIN

	whereclose, energy, minboundry(i), j, /silent
	whereclose, energy, maxboundry(i), k, /silent

;	print, energy(j),energy(k), mean(spectrum(j:k))
	if (j EQ k) then begin
		newspectrum(i) = spectrum(j)
	endif else begin
		newspectrum(i) = mean(spectrum(j:k))
	endelse 
	
ENDFOR


RETURN
END
