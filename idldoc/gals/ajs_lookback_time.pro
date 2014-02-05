; docformat = 'rst'
;+
; Return lookback time to z, in years
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
FUNCTION ajs_lookback_time, z, h0=h0, _REF_EXTRA=e
  compile_opt idl2
  IF n_elements(h0) EQ 0 THEN $
     h0 = 70

  return, ajs_comdis(z, /lookback_time, h0=h0, _strict_extra=e) $
          / 60. / 60. / 24. / 365.25
END
