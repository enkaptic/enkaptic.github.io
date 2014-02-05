; docformat = 'rst'
;+
; This function generates a Schechter function along an input array of
; absolute magnitudes.
;
; Compatible with mpfit
;
; :Params:
;    absmag : in, required
;       Array of absolute magnitudes
;    params : in, optional
;       [mstar,alpha,phistar] OR input separately:
; :Keywords:
;    mstar : in, optional
;    alpha : in, optional
;    phistar : in, optional
;    loglum : in, optional
;       Set to /loglum to input log10 luminosity instead of absolute magnitudes
; :Returns:
;    Array of phi for corresponding absmag
; :History:
;    5 Sep 2007 Created, Anthony Smith
;
;    17 Jan 2008 Added loglum keyword, AJS
;-
FUNCTION ajs_schechter,absmag,params,mstar=mstar,alpha=alpha,phistar=phistar, $
                       loglum=loglum
  compile_opt idl2

  IF n_params() EQ 2 THEN BEGIN
      mstar=params[0]
      alpha=params[1]
      phistar=params[2]
  ENDIF 

  absmag=double(absmag)
  IF keyword_set(loglum) THEN BEGIN
      phi = alog(10) * phistar $
            * 10 ^ (- (alpha + 1) * (mstar - absmag)) $
            * exp(-10 ^ (- (mstar - absmag)))
  ENDIF ELSE BEGIN
      phi = 0.4 * alog(10) * phistar $
            * 10 ^ (0.4 * (alpha + 1) * (mstar - absmag)) $
            * exp(-10 ^ (0.4 * (mstar - absmag)))
  ENDELSE

  RETURN,phi
END
