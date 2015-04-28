pro t90, Time_User, Flux_User, t90, tmin, tmax, fluence, MedianMaxFluence, plot=plot, bblocks=bblocks

IF (n_params() EQ 0) THEN BEGIN 
	print,'-Usage: t90, time, flux, t90, tmin, tmax, fluence, plot=plot
    RETURN
ENDIF 

;; Making copies of the input arrays so that they're returned inspoiled
time = Time_User
flux = Flux_User
	
IF KEYWORD_SET(bblocks) THEN BEGIN

	time = time(uniq(time))
	flux = flux(uniq(time))
	time = time(1:*)
	flux = flux(1:*)
END


fluence = dblarr(n_elements(time))

IF (N_ELEMENTS(time) EQ 0) THEN time = findgen(n_elements(flux))

FOR i=1, n_elements(time)-1 DO BEGIN

	fluence(i) = int_tabulated(time(0:i),flux(0:i))

ENDFOR

;; Finding the median of the fluence after the burst has returned 
TailEndIndex= long(n_elements(fluence)*0.65)
MedianMaxFluence = median( fluence( TailEndIndex:* ) )

;tmin = time(max(where(fluence LE max(fluence)*0.05)))
;tmax = time(min(where(fluence GE max(fluence)*0.95)))
tmin = time(max(where(fluence LE MedianMaxFluence*0.05)))
tmax = time(min(where(fluence GE MedianMaxFluence*0.95)))


t90 = tmax - tmin

IF KEYWORD_SET(plot) THEN BEGIN

	plot, time, flux
	oplot, [tmin,tmin],[1D-5,1D55], linestyle=1
	oplot, [tmax,tmax],[1D-5,1D55], linestyle=1

ENDIF


return
end



