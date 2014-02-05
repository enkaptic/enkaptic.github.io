; docformat = 'idl'
;+
; NAME:
;	ajs_logspace
;
;
; PURPOSE:
;	This function returns an array of points equally (logarithmically)
;	spaced between two extremes (as Python's numpy.logspace or Matlab).
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;	Result = ajs_logspace(start, stop)
;
;
; INPUTS:
;	start	Logarithm of first point in output array
;	stop	Logarithm of last point in output array
;
; OPTIONAL INPUTS:
;	num	Number of points in output array (default 51)
;
; KEYWORD PARAMETERS:
;	base	Base for logarithm (default 10)
;
; OUTPUTS:
;	Returns array of floating point numbers, or double-precision if
;	either start or stop is double precision
;
;
; OPTIONAL OUTPUTS:
;	step	Step size (logarithmic) used
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
;	Result = ajs_logspace(-10, 15)
;
;
; MODIFICATION HISTORY:
;	17 Jan 2008 Created, Anthony Smith
;-
FUNCTION ajs_logspace, start, stop, num, base=base, step=step
  compile_opt idl2

  IF n_params() EQ 2 THEN $
     num = 51
  IF NOT keyword_set(base) THEN $
     base = 10

  IF (size(start, /type) EQ 5) OR (size(stop, /type) EQ 5) THEN $
     integers = dindgen(num) $  ; Double precision
  ELSE $
     integers = findgen(num)    ; Floating point 

  step = (stop - start) / (num - 1.)
  linspace_array = start + integers * step
  logspace_array = 10. ^ (linspace_array * alog10(base))

  return, logspace_array

END
