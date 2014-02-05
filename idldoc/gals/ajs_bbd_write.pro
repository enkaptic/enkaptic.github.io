; docformat = 'rst'
;+
; This procedure takes arrays describing the BBD (bivariate brightness
; distribution) and writes the data either to a binary file or to a
; text file.
;
; :Params:
;    bincentres1 : in, required, type="fltarr or dblarr(nbins1)"
;    bincentres2 : in, required, type="fltarr or dblarr(nbins2)"
;    phibbd : in, required, type="fltarr or dblarr(nbins1, nbins2)"
; :Keywords:
;    datfile : in, optional, type=string
;       Name of binary file for output (read using ajs_bbd_read)
;    txtfile : in, optional, type=string
;       Name of text file for output
; :History:
;    11 Jan 2008 Created, Anthony Smith
;-
PRO ajs_bbd_write, bincentres1, bincentres2, phibbd, phibbderr, $
                   datfile=datfile, txtfile=txtfile
  compile_opt idl2

  ;; Write BBD to file

  IF n_elements(datfile) GT 0 THEN BEGIN
      ;; Binary file
      openw, unit, datfile, /get_lun
      writeu, unit, size(bincentres1, /type), size(bincentres1, /n_elements), $
              bincentres1, $
              size(bincentres2, /type), size(bincentres2, /n_elements), $ 
              bincentres2, $
              size(phibbd, /type), size(phibbd, /dimensions), phibbd, $
              size(phibbderr, /type), size(phibbderr, /dimensions), phibbderr
      free_lun, unit     
  ENDIF
  IF n_elements(txtfile) GT 0 THEN BEGIN
      ;; Text file
      openw, unit, txtfile, /get_lun
      FOR i = 0, n_elements(bincentres1) - 1 DO BEGIN
          FOR j = 0, n_elements(bincentres2) - 1  DO BEGIN
              printf, unit, bincentres1[i], bincentres2[j], phibbd[i,j], $
                      phibbderr[i,j]
          ENDFOR
      ENDFOR 
      free_lun, unit
  ENDIF
END  
