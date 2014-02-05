; docformat = 'idl'

;+
; NAME:
;	ajs_bbd_read
;
;
; PURPOSE:
;	This procedure reads BBD (bivariate brightness distribution)
;	data from a binary file created using ajs_bbd_write.
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;
;
;
; INPUTS:
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;	Written by: Anthony Smith
;	July 1994 Describe modifications
;	11 Jan 2008 Created, Anthony Smith
;-
PRO ajs_bbd_read, datfile, bincentres1, bincentres2, phibbd, phibbderr
  compile_opt idl2

  ;; Read BBD dat file

  openr, unit, datfile, /get_lun

  ;; Find dimensions of inputs and read
  bincentres1_type = 0L
  bincentres1_n_elements = 0L
  readu, unit, bincentres1_type, bincentres1_n_elements
  CASE bincentres1_type OF
      4: bincentres1 = fltarr(bincentres1_n_elements)
      5: bincentres1 = dblarr(bincentres1_n_elements)
  ENDCASE
  readu, unit, bincentres1

  bincentres2_type = 0L
  bincentres2_n_elements = 0L
  readu, unit, bincentres2_type, bincentres2_n_elements
  CASE bincentres2_type OF
      4: bincentres2 = fltarr(bincentres2_n_elements)
      5: bincentres2 = dblarr(bincentres2_n_elements)
  ENDCASE
  readu, unit, bincentres2

  phibbd_type = 0L
  phibbd_dimensions = lonarr(2)
  readu, unit, phibbd_type, phibbd_dimensions
  CASE phibbd_type OF
      4: phibbd = fltarr(phibbd_dimensions)
      5: phibbd = dblarr(phibbd_dimensions)
  ENDCASE
  readu, unit, phibbd

  phibbderr_type = 0L
  phibbderr_dimensions = lonarr(2)
  readu, unit, phibbderr_type, phibbderr_dimensions
  CASE phibbderr_type OF
      4: phibbderr = fltarr(phibbderr_dimensions)
      5: phibbderr = dblarr(phibbderr_dimensions)
  ENDCASE
  readu, unit, phibbderr

  free_lun, unit
END 
