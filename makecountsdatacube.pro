;; ==========================================================================================
;; File:  MakeCountsDataCube
;; Author:  Daniel Kocevski <kocevski@berkeley.edu>
;; Date Created:  circa 2002
;; Last Changed:  -
;; 
;; Purpose:  This program takes a drm matrix and a photon spectrum and creates a count rate data cube
;; Dependancies:  None
;; ==========================================================================================


pro MakeCountsDataCube, matrix, photonspectrum, photonenergywidth, channelenergywidth, countsmatrix, countspectrum, batse=batse

IF (n_params() EQ 0) THEN BEGIN 
	print,'-Usage: MakeCountsDataCube, maxtrix, photonspectrum, photonenergywidth, channelenergywidth, countsmatrix, countspectrum, batse=batse
    RETURN
ENDIF 

IF KEYWORD_SET(batse) THEN BEGIN
	countsmatrix = transpose(matrix)
ENDIF ELSE BEGIN
	countsmatrix = matrix
ENDELSE

countspectrum = countsmatrix(*,0)
countspectrum(*) = 0

;; Convolve the data response matrix with the photon spectrum
FOR i=0, N_ELEMENTS(countsmatrix(*,0))-1 DO BEGIN

	countsmatrix(i,*) = countsmatrix(i,*) * photonspectrum 
	countsmatrix(i,*) = countsmatrix(i,*) * photonenergywidth

ENDFOR

;; Divide by the channel energy widths to get the proper units (counts s-1 kev-1)
FOR j=0, N_ELEMENTS(countsmatrix(0,*))-1 DO BEGIN

	;; This has the units of counts s-1 kev-1
;;	countsmatrix(*,j) = countsmatrix(*,j)/channelenergywidth

ENDFOR

;; Collapse the counts matrix to form a 128 channel count rate spectrum
FOR j=0, N_ELEMENTS(countsmatrix(*,0))-1 DO BEGIN

	;; This has the units of counts s-1 kev-1
	countspectrum(j) = total( countsmatrix(j,*) ) 

ENDFOR



RETURN
END
