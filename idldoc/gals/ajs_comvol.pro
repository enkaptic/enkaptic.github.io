; docformat = 'rst'
;+
; Return comoving volume out to redshift z (flat universe only)
;
; See ajs_comdis for other keyword parameters
; :Params:
;    z : in, required
;       Redshift (or array of values)
; :Keywords:
;    area : in, optional
;       Area in square degrees. Defaults to whole sky.
;-
FUNCTION ajs_comvol, z, area=area, _REF_EXTRA=e
  compile_opt idl2

  IF n_elements(area) GT 0 THEN $
     frac_sky = area / (360. ^ 2 / !PI) $
  ELSE $
     frac_sky = 1

  return, frac_sky * 4 * !PI / 3 * ajs_comdis(z, _strict_extra=e) ^ 3
END
