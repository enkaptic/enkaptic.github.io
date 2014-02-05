; docformat = 'rst'
;+
; This procedure takes arrays describing the LF (luminosity function)
; and writes the data either to a binary file or to a text file.
;
; :History: 11 Jan 2008 Created, Anthony Smith
;-
PRO ajs_lf_write, bincentres, philf, philferr, $
                  datfile=datfile, txtfile=txtfile
  compile_opt idl2

  ;; Write LF to file

  IF keyword_set(datfile) THEN BEGIN
      ;; Binary file
      openw, unit, datfile, /get_lun
      IF n_params() EQ 3 THEN BEGIN
          writeu, unit, size(bincentres, /type), $
                  size(bincentres, /n_elements), $
                  bincentres, $
                  size(philf, /type), size(philf, /n_elements), philf, $
                  size(philferr, /type), size(philferr, /n_elements), philferr
      ENDIF ELSE IF n_params() EQ 2 THEN BEGIN
          writeu, unit, size(bincentres, /type), $
                  size(bincentres, /n_elements), $
                  bincentres, $
                  size(philf, /type), size(philf, /n_elements), philf
      ENDIF
      free_lun, unit     
  ENDIF
  IF keyword_set(txtfile) THEN BEGIN
      ;; Text file
      openw, unit, txtfile, /get_lun
      FOR i = 0, n_elements(bincentres) - 1 DO BEGIN
          IF n_params() EQ 3 THEN $
             printf, unit, bincentres[i], philf[i], philferr[i] $
          ELSE IF n_params() EQ 2 THEN $
             printf, unit, bincentres[i], philf[i]
      ENDFOR 
      free_lun, unit
  ENDIF
END  
