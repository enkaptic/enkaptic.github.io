; docformat = 'rst'
;+
; This procedure reads LF (luminosity function) data from a binary
; file creaed using ajs_lf_write.
;
; :History: 11 Jan 2008 Created, Anthony Smith
;-
PRO ajs_lf_read, datfile, bincentres, philf, philferr
  compile_opt idl2

  ;; Read LF dat file

  openr, unit, datfile, /get_lun

  ;; Find dimensions of inputs and read
  bincentres_type = 0L
  bincentres_n_elements = 0L
  readu, unit, bincentres_type, bincentres_n_elements
  CASE bincentres_type OF
      4: bincentres = fltarr(bincentres_n_elements)
      5: bincentres = dblarr(bincentres_n_elements)
  ENDCASE
  readu, unit, bincentres

  philf_type = 0L
  philf_n_elements = 0L
  readu, unit, philf_type, philf_n_elements
  CASE philf_type OF
      4: philf = fltarr(philf_n_elements)
      5: philf = dblarr(philf_n_elements)
  ENDCASE
  readu, unit, philf

  IF eof(unit) EQ 0 THEN BEGIN
      philferr_type = 0L
      philferr_n_elements = 0L
      readu, unit, philferr_type, philferr_n_elements
      IF philferr_type GT 0 THEN BEGIN
          CASE philferr_type OF
              4: philferr = fltarr(philferr_n_elements)
              5: philferr = dblarr(philferr_n_elements)
          ENDCASE
          readu, unit, philferr
      ENDIF
  ENDIF

  free_lun, unit
END 
