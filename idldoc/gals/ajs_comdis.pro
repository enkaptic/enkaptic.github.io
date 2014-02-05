; docformat = 'rst'
;+
; Return comoving distance to z, in Mpc
;-


;+
; Integrand for comoving distance integral
;-
FUNCTION ajs_comdis_er, z
  compile_opt idl2, hidden

  COMMON cosmology, omega_m, omega_l
  return, 1.0 / sqrt(omega_m * (1 + z) ^ 3 $
                     + (1 - omega_m - omega_l) * (1 + z) ^ 2 $
                     + omega_l)
END


;+
; Integrand for lookback time integral
;-
FUNCTION ajs_lookback_er, z
  compile_opt idl2, hidden

  return, ajs_comdis_er(z) / (1 + z)
END


;+
; Return comoving distance to z, in Mpc
; :Params:
;    z, in, required
;       Redshift (or array of values)
; :Keywords:
;    omega_m : in, optional
;       Matter density (default 0.3)
;    omega_l : in, optional
;       Vacuum energy density (default 0.7)
;    h0 : in, optional
;       Hubble constant (default 100)
;    c : in, optional
;       Speed of light (default 299792.458 km/s)
;    slow : in, optional
;       Set /slow to force calculation for each value (defaults to
;       interpolation if input array has more then 10,000 elements)
;    lookback_time : in, optional
;       Set /lookback_time to return lookback time in seconds.  See
;       ajs_lookback_time and ajs_zage.
; :History:
;    13 Mar 2008 Written, Anthony Smith
;
;    10 Oct 2008 Added lookback_time, AJS
;
;    16 Jan 2009 Corrected bug, which was affecting large input arrays
;-
FUNCTION ajs_comdis, z, omega_m=omega_m, omega_l=omega_l, h0=h0, c=c, $
                     slow=slow, lookback_time=lookback_time
  compile_opt idl2
  COMMON cosmology, o_m, o_l
  IF n_elements(omega_m) EQ 0 THEN $
     o_m = 0.3 $
  ELSE $
     o_m = omega_m
  IF n_elements(omega_l) EQ 0 THEN $
     o_l = 0.7 $
  ELSE $
     o_l = omega_l
  IF n_elements(h0) EQ 0 THEN $
     h0 = 100
  IF n_elements(c) EQ 0 THEN $ 
     c = 299792.458
  
  z = double(z)
  IF n_elements(z) GT 1e4 AND NOT keyword_set(slow) THEN BEGIN
      z_arr = ajs_linspace(min(z), max(z), 5001)
      d_arr = ajs_comdis(z_arr, omega_m=o_m, omega_l=o_l, h0=h0, c=c, $
                         lookback_time=lookback_time)
      d = interpol(d_arr, z_arr, z, /quadratic)
  ENDIF ELSE BEGIN
      IF keyword_set(lookback_time) THEN $
         d = qromb('ajs_lookback_er', 0, z, /double) / h0 * 3.085678e19 $
      ELSE $
         d = qromb('ajs_comdis_er', 0, z, /double) * c / h0
  ENDELSE
  return, d
END
