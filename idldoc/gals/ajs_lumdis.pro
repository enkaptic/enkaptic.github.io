; docformat = 'rst'
;+
; Return Luminosity distance to z
;
; See ajs_comdis for keyword parameters
; :Params:
;    z, in, required
;       Redshift (or array of values)
;-
FUNCTION ajs_lumdis, z, _REF_EXTRA=e
  compile_opt idl2

  return, ajs_comdis(z, _strict_extra=e) * (1 + z)
END
