; docformat = 'idl'
;+
; NAME:
;	ajs_linspace
;
;
; PURPOSE:
;	This function returns an array of points equally (linearly)
;	spaced between two extremes (as Python's numpy.linspace or Matlab).
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;	Result = ajs_linspace(start, stop)
;
;
; INPUTS:
;	start	First point in output array
;	stop	Last point in output array
;
; OPTIONAL INPUTS:
;	num	Number of points in output array (default 51, or 50 bincentres)
;
;
; KEYWORD PARAMETERS:
;	bincentres	Set /bincentres to return centres of bins
;			bounded by start and stop
;
;
; OUTPUTS:
;	Returns array of floating point numbers, or double-precision if
;	either start or stop is double precision
;
;
; OPTIONAL OUTPUTS:
;	step	Step size used
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
;	Result = ajs_linspace(10, 20)
;
;
; MODIFICATION HISTORY:
;	17 Jan 2008 Created, Anthony Smith
;-
FUNCTION ajs_linspace, start, stop, num, step=step, bincentres=bincentres
  compile_opt idl2

  IF n_elements(num) EQ 0 THEN $
     IF keyword_set(bincentres) THEN $
        num = 50 $
     ELSE $
        num = 51

  IF (size(start, /type) EQ 5) OR (size(stop, /type) EQ 5) THEN BEGIN
      integers = dindgen(num)   ; Double precision
      num = double(num)
  ENDIF ELSE BEGIN
      integers = findgen(num)   ; Floating point 
      num = float(num)
  ENDELSE

  IF keyword_set(bincentres) THEN BEGIN
      step = (stop - start) / num
      start = start + step / 2
      stop = stop - step / 2
  ENDIF ELSE $
     step = (stop - start) / (num - 1.)

  linspace_array = start + integers * step

  return, linspace_array

END
