; docformat = 'idl'
;+
; PURPOSE:
;	This function finds the best-fitting Schechter function for a
;	given array of absolute magnitudes, phi and phi_err
;-

;+
; NAME:
;	ajs_schechter_fit
;
;
; PURPOSE:
;	This function finds the best-fitting Schechter function for a
;	given array of absolute magnitudes, phi and phi_err
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
;	phi
;	phierr
;
; OPTIONAL INPUTS:
;	params_init	Initial guess for [mstar,alpha,phistar] (or default)
;	txtfile		Output text file
;
; KEYWORD PARAMETERS:
;	loglum	Set to input log10 luminosity instead of absolute magnitudes
;
; OUTPUTS:
;	Result	[mstar,alpha,phistar]
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
;	Uses mpfit (idlutils)
;	ajs_schechter
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
;	5 Sep 2007 Created, Anthony Smith
;	17 Jan 2008 Added loglum keyword and text file output, AJS
;-
FUNCTION ajs_schechter_fit,absmag, phi, phierr, params_init, txtfile=txtfile, $
                           loglum=loglum

  IF n_elements(params_init) EQ 0 THEN $
     params_init=[median(absmag), -0.5, 0.005]

  IF NOT keyword_set(loglum) THEN $
     loglum = 0

  IF (where(phierr EQ 0))[0] NE -1 THEN $
     print, 'ajs_schechter_fit warning: ignoring bins with zero error.'
  params=mpfitfun('ajs_schechter', $
                  absmag[where(phi GT 0)], $
                  phi[where(phi GT 0)], $
                  phierr[where(phi GT 0)], $
                  params_init,/quiet, covar=covar, perror=perror, $
                  functargs=create_struct('loglum', loglum))
;;   print, 'PError:', float(perror)
;;   print, 'Covariance:'
;;   print, float(covar)
;;   PCOR = covar * 0
;;   n = 3
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

;+
; PURPOSE:
;	Test ajs_schechter_fit
;-
PRO ajs_schechter_fit_test

  params_true = [-25, -1.3, 0.02123]
  m = ajs_linspace(-20, -26, 51)
  print,'Starting with', params_true
  phi = ajs_schechter(m, params_true)
  phierr = phi / 10
  plot, m, phi, /ylog

  ;; Fit
  params = ajs_schechter_fit(m, phi, phierr)
  print,'Fitted parameters:', params
  phifit = ajs_schechter(m, params)
  oplot, m, phifit, color=1, linestyle=5


END
