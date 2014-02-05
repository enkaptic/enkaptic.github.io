; docformat = 'rst'
;+
; Convert redshift to age of Universe, in years
;
; See ajs_comdis for keyword parameters
; :Params:
;    z, in, required
;       Redshift (or array of values)
; :Keywords:
;    h0, in, optional
;       Hubble constant (default 70)
; :History:
;    10 Oct 2008 Written, Anthony Smith
;-
FUNCTION ajs_zage, z, omega_m=omega_m, omega_l=omega_l, h0=h0, _REF_EXTRA=e
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
     h0 = 70

  ;; Analytic form from Liddle (2003) p. 59
  IF o_m + o_l EQ 1 THEN $
     age =  2. / 3 * 1 / sqrt(1 - o_m) $
            * alog((1 + sqrt(1 - o_m)) / sqrt(o_m)) $
            * 3.085678e19 / h0 / 60. / 60. / 24. / 365.25 $
  ELSE $
     age = ajs_lookback_time(1e4, omega_m=o_m, omega_l=o_l, h0=h0, $
                             _STRICT_EXTRA=e)

  return, age - ajs_lookback_time(z, omega_m=o_m, omega_l=o_l, h0=h0, $
                                  _STRICT_EXTRA=e)
END
