; docformat = 'rst'
;+
; Return distance modulus to z
;
; See ajs_comdis for keyword parameters
; :Params:
;    z, in, required
;       Redshift (or array of values)
;-
FUNCTION ajs_distmod, z, _REF_EXTRA=e
  compile_opt idl2

  return, 5 * alog10(ajs_comdis(z, _strict_extra=e) * (1 + z) * 1e6 / 10)
END
