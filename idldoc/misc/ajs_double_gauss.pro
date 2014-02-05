; docformat = 'rst'

;+
; This function returns the sum of two Gaussians, using gauss1
; (mpfit). Suitable for curve fitting.
;
; :Params:
;    x : in, required
;       Array of X values
;    p : in, required
;       [mean1, sigma1, area1, mean2, sigma2, area2] parameters for
;       the two Gaussian components
; :Keywords:
;    _ref_extra
;       Extra parameters to be passed to gauss1
; :History: 11 Jan 2008	Created, Anthony Smith
;-
FUNCTION ajs_double_gauss, x, p, _REF_EXTRA=e
  compile_opt idl2

  return, gauss1(x, p[0:2], _strict_extra=e) $
          + gauss1(x, p[3:5], _strict_extra=e)
END

