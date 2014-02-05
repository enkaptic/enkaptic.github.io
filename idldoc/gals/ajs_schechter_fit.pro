; docformat = 'rst'
;+
; This function finds the best-fitting Schechter function for a
; given array of absolute magnitudes, phi and phi_err
;-

;+
; This function finds the best-fitting Schechter function for a
; given array of absolute magnitudes, phi and phi_err
; :Returns:
;    [mstar, alpha, phistar]
; :Params:
;    absmag : in, required
;    phi : in, required
;    phierr : in, required
;    params_init : in, optional
;       Initial guess for [mstar,alpha,phistar] (or default)
; :Keywords:
;    txtfile : in, optional
;       Output text file for Schechter parameters
;    loglum : in, optional
;       Set to input log10 luminosity instead of absolute magnitudes
;    cov_mat : out, optional
;       Covariance matrix
;    range : in, optional, type=fltarr(2)
;       [min, max] absmag for fitting
;    lum_dens : out, optional
;       Luminosity density.  Return
;       value is j = ... x 10^(0.4 M_sun) h L_sun Mpc^-3
;    err_lum_dens : out, optional
;       Error in luminosity density
;    _REF_EXTRA : in, optional
;       Extra keywords for mpfitfun
; :Uses:
;    mpfit (idlutils)
;
;    ajs_schechter
; :History:
;    5 Sep 2007 Created, Anthony Smith
;
;    17 Jan 2008 Added loglum keyword and text file output, AJS
;
;    4 Apr 2008 Added range keyword
;
;    15 Apr 2008 Added lum_dens
;-
FUNCTION ajs_schechter_fit,absmag, phi, phierr, params_init, txtfile=txtfile, $
                           loglum=loglum, cov_mat=cov_mat, range=range, $
                           lum_dens=lum_dens, err_lum_dens=err_lum_dens, $
                           _REF_EXTRA=e
  compile_opt idl2
  debug = ajs_debug()
  IF debug GE 1 THEN message, 'Finding Schechter fit', /inf

  IF n_elements(params_init) EQ 0 THEN $
     params_init=[median(absmag), -0.5, 0.005]

  IF NOT keyword_set(loglum) THEN $
     loglum = 0

  IF n_elements(range) GT 1 THEN BEGIN
      IF range[0] GT range[1] THEN $
         range = reverse(range) ; [min, max] rather than [max, min]
      fit = where(phierr GT 0 AND absmag GE range[0] AND absmag LE range[1], $
                  ncomplement=n_no_error_or_outside_range)
      res = where(absmag GE range[0] AND absmag LE range[1], $
                  ncomplement=n_outside_range)
      n_no_error = n_no_error_or_outside_range - n_outside_range
  ENDIF ELSE $
     fit = where(phierr GT 0, ncomplement=n_no_error)

  ;; Undefined errors?
  IF n_no_error GT 0 THEN $
     message, 'Ignoring bins with zero/undefined error.', /inf

  ;; Perform fit
  params=mpfitfun('ajs_schechter', $
                  absmag[fit], $
                  phi[fit], $
                  phierr[fit], $
                  params_init, /quiet, covar=cov_mat, perror=perror, $
                  functargs=create_struct('loglum', loglum), _STRICT_EXTRA=e)

  ;; Text file
  IF n_elements(txtfile) GT 0 THEN BEGIN
      openw, unit, txtfile, /get_lun
      printf, unit, params
      free_lun, unit
  ENDIF

  ;; Luminosity density
  IF keyword_set(loglum) THEN $
     lum_dens = params[2] * 10. ^ (params[0]) * gamma(params[1] + 2) $
  ELSE $
     lum_dens = params[2] * 10. ^ ((- params[0]) / 2.5) * gamma(params[1] + 2)
  ;; Variance: sigma_g^2 = sum_i,j ( dg/dx_i dg/dx_j cov(x_i, x_j) )
  ;; (partial derivatives)
  IF arg_present(err_lum_dens) THEN BEGIN
      IF keyword_set(loglum) THEN BEGIN 
          partial_derivs = [alog(10) * lum_dens, $
                            ajs_digamma(params[1] + 2) * lum_dens, $
                            lum_dens / params[2]]
          var_lum_dens = partial_derivs ## cov_mat ## (1 # partial_derivs)
          err_lum_dens = sqrt(var_lum_dens)
      ENDIF ELSE BEGIN
          partial_derivs = [- alog(10) / 2.5 * lum_dens, $
                            ajs_digamma(params[1] + 2) * lum_dens, $
                            lum_dens / params[2]]
          var_lum_dens = partial_derivs ## cov_mat ## (1 # partial_derivs)
          err_lum_dens = sqrt(var_lum_dens)
      ENDELSE
  ENDIF

  RETURN, params
END


;+
; Test ajs_schechter_fit
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
