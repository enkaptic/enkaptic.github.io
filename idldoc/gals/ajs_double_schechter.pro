; docformat = 'rst'
;+
; This function generates a double Schechter function (single value of
; M-star, but two values for alpha and phi-star) along an input array of
; absolute magnitudes.
;
; Compatible with mpfit
;
; :Params:
;    absmag : in, required
;       Array of absolute magnitudes
;    params : in, optional
;       [mstar, alpha_1, phistar_1, alpha_2, phistar_2] OR input separately:
; :Keywords:
;    mstar : in, optional
;    alpha_1 : in, optional
;    phistar_1 : in, optional
;    alpha_2 : in, optional
;    phistar_2 : in, optional
;    loglum : in, optional
;       Set to /loglum to input log10 luminosity instead of absolute magnitudes
; :Returns:
;    Array of phi for corresponding absmag
; :History:
;    21 April 2008 Written, Anthony Smith
;-
FUNCTION ajs_double_schechter, $
   absmag, params, mstar=mstar, alpha_1=alpha_1, $
   phistar_1=phistar_1, alpha_2=alpha_2, $
   phistar_2=phistar_2, loglum=loglum
  compile_opt idl2

  IF n_params() EQ 2 THEN BEGIN
      mstar = params[0]
      alpha_1 = params[1]
      phistar_1 = params[2]
      alpha_2 = params[3]
      phistar_2 = params[4]
  ENDIF 

  absmag = double(absmag)
  IF keyword_set(loglum) THEN BEGIN
      phi = alog(10) $
            * (phistar_1 * 10 ^ (- (alpha_1 + 1) * (mstar - absmag)) $
               + phistar_2 * 10 ^ (- (alpha_2 + 1) * (mstar - absmag))) $
            * exp(-10 ^ (- (mstar - absmag)))
  ENDIF ELSE BEGIN
      phi = 0.4 * alog(10) $
            * (phistar_1 * 10 ^ (0.4 * (alpha_1 + 1) * (mstar - absmag)) $
               + phistar_2 * 10 ^ (0.4 * (alpha_2 + 1) * (mstar - absmag))) $
            * exp(-10 ^ (0.4 * (mstar - absmag)))
  ENDELSE

  RETURN, phi
END
