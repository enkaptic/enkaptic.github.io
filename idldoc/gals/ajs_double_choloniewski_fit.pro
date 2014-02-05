; docformat = 'idl'
;+
; PURPOSE:
;       (NB: unlikely to be useful!) This function finds the
;       best-fitting double Choloniewski function for input arrays of
;       absolute magnitudes, surface brightness, phi and phi_err.
;-

;+
; NAME:
;	ajs_double_choloniewski_fit
;
;
; PURPOSE:
;	This function finds the best-fitting double Choloniewski function for
;	input arrays of absolute magnitudes, surface brightness, phi
;	and phi_err.
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
;	absmag
;	sbr
;	phi
;	phierr
;
; OPTIONAL INPUTS:
;	params_init	Initial guess for 
;			[mstar,alpha,phistar,sbstar,sigmasb,beta]
;			twice (or default)
;
; KEYWORD PARAMETERS:
;	txtfile		Output text file
;
;
; OUTPUTS:
;	Result	[mstar,alpha,phistar,sbstar,sigmasb,beta]
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
;	Uses mpfit2d (idlutils)
;	ajs_choloniewski
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
;	17 Jan 2008 Created, Anthony Smith
;-
FUNCTION ajs_double_choloniewski_fit, absmag, sb, phi, phierr, params_init, $
                                      txtfile=txtfile
  compile_opt idl2

  IF n_elements(params_init) EQ 0 THEN $
     params_init=[median(absmag) - 1, -0.5, 0.01, $
                  median(sb) - 1, 0.5, 0.2, $
                  median(absmag) + 1, -1.5, 0.01, $
                  median(sb) + 1, 1.5, 0.6]
  params_init = double(params_init)

  absmag_arr = absmag # (sb * 0 + 1)
  sb_arr = (absmag * 0 + 1) # sb

  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
                       limits:[0.D,0]}, 12)
  parinfo[1].limited = [1, 1]
  parinfo[1].limits  = double([1e-3, 1e-1])
  parinfo[7].limited = [1, 1]
  parinfo[7].limits  = [1e-3, 1e-1]

  params=mpfit2dfun('ajs_double_choloniewski', $
                    absmag_arr, $
                    sb_arr, $
                    phi, $
                    phierr, $
                    params_init, /quiet, covar=covar, perror=perror, $
                    parinfo=parinfo)
;;   print, 'PError:', float(perror)
;;   print, 'Covariance:'
;;   print, float(covar)
;;   PCOR = covar * 0
;;   n = 6
;;   FOR i = 0, n-1 DO FOR j = 0, n-1 DO $
;;      PCOR(i,j) = covar(i,j)/sqrt(covar(i,i)*covar(j,j))
;;   print, 'Correlation:'
;;   print, float(pcor)

  IF keyword_set(txtfile) THEN BEGIN
      ;; Text file
      openw, unit, txtfile, /get_lun
      printf, unit, params
      free_lun, unit
  ENDIF

  RETURN,params

END

PRO ajs_choloniewski_fit_test

  ;; This works:
;;   params_true = [-24, 0, 1e-2, 17, 0.5, 0, $
;;                  -21, -1, 1e-2, 18, 0.9, 0.5]                 
  params_true = [-24, 0, 1e-2, 17, 0.5, 0, $
                 -22, -1, 1e-2, 18, 0.9, 0.5]                 
  m = ajs_linspace(-16, -26, 21)
  sb = ajs_linspace(22, 14, 21)
  print,'Starting with'
  print, params_true[0:5]
  print, params_true[6:11]
  chol = ajs_double_choloniewski(m, sb, params_true)
  cholerr = chol / 10
  ajs_bbd_plot, m, sb, chol, /show_plot, plotfile='/tmp/bbd.eps'

  ;; Fit
;  params_init = params_true * 1.1
  params = ajs_double_choloniewski_fit(m, sb, chol, cholerr);, params_init)
  print,'Fitted parameters:'
  print, params[0:5]
  print, params[6:11]
  cholfit = ajs_double_choloniewski(m, sb, params)
  ajs_bbd_plot, m, sb, cholfit, /show_plot, plotfile='/tmp/bbd2.eps'


END
