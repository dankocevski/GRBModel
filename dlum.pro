;This program calculates the luminosity distance of a user defined redshift for a set cosmology.
;
;Dan Kocevski  4-9-02
;
;with code borrowed from Jay Norris

pro dlum, z, d_lum, d_ang, d_motion, co_volume, silent=silent, cm=cm, Ho=Ho, Mpc=Mpc

    IF (n_params() EQ 0) THEN BEGIN 
      print,'-Usage: dlum, z, d_lum, d_ang, d_motion, co_volume, silent=silent, cm=cm, Ho=Ho, Mpc=Mpc'
      return
  ENDIF 
    
    
    ; Setting the speed of light
       c=3.0d5

       ; Setting the default Hubble Constant
       IF (n_elements(Ho) eq 0) THEN Ho=71.d0
       rh = c/Ho
       c20 = 2.d0
       c10 = 1.d0

;
; ENTER COSMOLOGY
;

       lambda=0.73d0
       omega = 0.27d0

       IF keyword_set(silent) THEN BEGIN

       ENDIF ELSE BEGIN
       		   print, 'Ho     = ', Ho
               print, 'Lambda = ',lambda
               print, 'Omegao = ',omega
       ENDELSE


; Lumosity Distance in Gpc
       d_lum_sub=dang_std(rh,lambda,omega,0.d0,z)
       d_lum=((1.d0+z)^2)*d_lum_sub*0.001d0


;  Angular dist
       d_ang =  dang_std(rh,lambda,omega,0.d0,z) * 0.001d0
;      print,'Ang=',d_ang


;  Proper distance
       dproper = D_proper(rh,lambda,omega,0.d0,z)*0.001d0
;      print*,'Proper=',dproper


; Proper motion distance
        d_motion = (1.d0+z) * dang_std(rh,lambda,omega,0.d0,z) * 0.001d0
;       print*,'Motion=',d_motion


; Andular diamater distance segun Kochanek, ApJ, 384, 1, 1992.
       x2 = 1.+z
       x1 = 1.
       g1 = 1.
       g2 = sqrt(1.+ omega*z)
       d_kochanek = 2.d0*rh*(c10-omega-g1*g2)*(g1-g2)/(omega*omega*x1*x2*x2) * 0.001d0


; Proper distance differential
       dDolp = d_olp(rh,lambda,omega,z,0.1d0)
;      print,'dDolp = ',dDolp


;Comoving Volume (Hogg astro-ph9905116 v4)
co_volume = ((4.0/3.0)*!PI)*(d_motion^3.0)

IF KEYWORD_SET(cm) THEN BEGIN

    d_motion = d_motion*3.08568D27
    d_ang = d_ang*3.08568D27
    d_lum = d_lum*3.08568D27
    co_volume = ((4.0/3.0)*!PI)*(d_motion^3.0)

ENDIF

IF KEYWORD_SET(Mpc) THEN BEGIN

    d_motion = d_motion*1E3
    d_ang = d_ang*1E3
    d_lum = d_lum*1E3
    co_volume = ((4.0/3.0)*!PI)*(d_motion^3.0)

ENDIF

IF keyword_set(silent) THEN BEGIN

ENDIF ELSE BEGIN

    print, 'Redshift of GRB =', z
    print, 'Angular Distance = ', d_ang
    print, 'Comoving Distance =', d_motion
    print, 'Luminosity Distance =', d_lum
    print, 'Comoving Volume =', co_volume

ENDELSE


RETURN
END


;-------------------------------------------------------
       function dang_std, rh,lambda,omega,z1,z2
;--------------------------------------------------------
; Filled beam or smooth Universe, alpha = 1
;--------------------------------------------------------
;       implicit real*8(a-h,o-z)
;       parameter(epsilon = 1.d-5)
       epsilon = 1.d-5


;       real*8 int1,Integral,lambda
       c10 = 1.d0
       c20 = 2.d0
       c30 = 3.d0

       x1 = c10 + z1
       x2 = c10 + z2

       aux = lambda + omega
;      if(abs(lambda) gt epsilon) then
          int1=Integral(omega,lambda,z1,z2)
          aux2 = sqrt(abs(omega + lambda - c10))
          factor = aux2*int1
          if (aux gt c10) then begin        ;open universe,  k=-1
            Dang = rh/(x2*aux2)*sin(factor)
          endif
          if (aux lt c10) then begin    ;closed universe,k=1
            Dang = rh/(x2*aux2)*sinh(factor)
          endif else begin
            Dang = rh/x2*int1         ;flat universe,k=0
          endelse
;      else
;         Dang = c20*rh*(c10 - omega - sqrt(x1*x2))*(sqrt(x1)-sqrt(x2))/
;    &                 (omega*omega*x1*x2*x2)
;        Dang=c20*rh*((c20-omega + omega*z2)*sqrt(c10+omega*z1) -
;    &    (c20 - omega + omega*z1)*sqrt(c10 + omega*z2))/(omega*omega*x1*x2*x2)
;      endif
       dang_std = dang
       return, dang_std
       end

;--------------------------------------------------------
        function D_proper, rh,lambda,omega,z1,z2
;--------------------------------------------------------
; Filled beam or smooth Universe, alpha = 1
;--------------------------------------------------------
;       implicit real*8(a-h,o-z)
       epsilon = 1.d-5

;       real*8 int1,rIntegral_proper,lambda

       c10 = 1.d0
       c20 = 2.d0
       c30 = 3.d0

       x1 = c10 + z1
       x2 = c10 + z2

       aux = lambda + omega
;      if(abs(lambda).gt.epsilon) then
          int1=rIntegral_proper(rh,lambda,omega,z1,z2)
          aux2 = sqrt(abs(omega + lambda - c10))
          factor = aux2*int1
          if (aux gt c10) then begin         ;open universe,  k=-1
            Dang = rh/(x2*aux2)*sin(factor)
          endif
          if(aux lt c10) then begin    ;closed universe,k=1
            Dang = rh/(x2*aux2)*sinh(factor)
          endif else begin
            Dang = rh/x2*int1         ;flat universe,k=0
          endelse
;      else
;         Dang = c20*rh*(c10 - omega - sqrt(x1*x2))*(sqrt(x1)-sqrt(x2))/
;    &                 (omega*omega*x1*x2*x2)
;        Dang=c20*rh*((c20-omega + omega*z2)*sqrt(c10+omega*z1) -
;    &    (c20 - omega + omega*z1)*sqrt(c10 + omega*z2))/(omega*omega*x1*x2*x2)
;      endif
       D_proper = Dang
       return, D_proper
       end




;--------------------------------------------------------
       function Integral, omega,lambda,x1,x2
;--------------------------------------------------------
;       implicit real*8(a-h,o-z)
       nsteps= 1000
       epsilon = 1.d-5

;       real*8 lambda

       c00 = 0.d0
       c10 = 1.d0
       c20 = 2.d0
       c30 = 3.d0
       c1000 = 1000.d0
       c05 = 0.5d0


;l WE HAVE TO SOLVE THIS: INTEGRAL TRAPEZOIDAL METHOD
       du = (x2-x1)/double(nsteps)
       u = x1
       f = c10/sqrt(omega*(c10+u)*(c10+u)*(c10+u) + (c10 - omega - lambda)*(c10+u)*(c10+u) + lambda)*du*c05
       for i=1,nsteps-1 do begin
         u = u + du
         f = f + c10/sqrt(omega*(c10+u)*(c10+u)*(c10+u) + (c10 - omega - lambda)*(c10+u)*(c10+u) + lambda)*du
       endfor
       u = x2
       f = f + c10/sqrt(omega*(c10+u)*(c10+u)*(c10+u) + (c10 - omega - lambda)*(c10+u)*(c10+u) + lambda)*du*c05
       Integral= f
       return, integral
       end

;-----------------------------------------------------
       function rIntegral_proper, rh,lambda,omega,z1,z2
;-----------------------------------------------------
;      implicit real*8(a-h,o-z)
       epsilon = 1.d-5
       nsteps = 1000

;      real*8 lambda

       c00 = 0.d0
       c10 = 1.d0
       c20 = 2.d0
       c05 = 0.5d0

       dz = (z2-z1)/double(nsteps)
       z = z1
       f= c10/((c10+z)*sqrt((c10 + omega*z)*(c10+z)*(c10+z) - z*(c20+z)*lambda))*c05
       for i=1,nsteps-1 do begin
         z = z + dz
         f = f + c10/((c10+z)*sqrt((c10 + omega*z)*(c10+z)*(c10+z) - z*(c20+z)*lambda))
       endfor
       z = z2
       f = f + c10/((c10+z)*sqrt((c10 + omega*z)*(c10+z)*(c10+z) - z*(c20+z)*lambda))*c05
       rIntegral_proper=f*dz
       return, rIntegral_proper
       end
;-----------------------------------------------------------
       function d_olp, rh,lambda,omega,zgrb,dzgrb
;-----------------------------------------------------------
; differential element of proper distance
;-----------------------------------------------------------
;       implicit real*8(a-h,o-z)

;      real*8 lambda
       c00 = 0.d0
       c10 = 1.d0
       c20 = 2.d0
       c30 = 3.d0
       O_k = c10 - (omega + lambda)
       c10z = c10 + zgrb
       d_olp =rh*dzgrb/c10z/sqrt(omega*c10z^3 + O_k*c10z^2+ lambda)




       return,  d_olp
       end






