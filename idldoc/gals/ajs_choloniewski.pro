; docformat = 'idl'

;+
; NAME:
;	ajs_choloniewski
;
;
; PURPOSE:
;	This function generates a Choloniewski function along an input
;	array of absolute magnitudes and surface brightnesses.
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;	result = ajs_schechter(absmag,sb,params)
;	result = ajs_schechter(absmag,sb, $
;                              sbstar=sbstar, sigmasb=sigmasb, beta=beta, $
;			       mstar=mstar,alpha=alpha,phistar=phistar)
;
; INPUTS:
;	absmag
;	sb
;
;
; OPTIONAL INPUTS:
;	params:	[mstar,alpha,phistar,sbstar,sigmasb,beta] OR input separately:
;
; KEYWORD PARAMETERS:
;	mstar
;	alpha
;	phistar
;	sbstar
;	sigmasb
;	beta
;
; OUTPUTS:
;	This function returns an rray of phi for corresponding absmag, sbr
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
;	Compatible with mpfit2d.
;
;	Magnitudes only
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
FUNCTION ajs_choloniewski, absmag, sb, params, $
                           mstar=mstar, alpha=alpha, phistar=phistar, $
                           sbstar=sbstar, sigmasb=sigmasb, beta=beta
  compile_opt idl2

  IF n_params() EQ 3 THEN BEGIN
      mstar = params[0]
      alpha = params[1]
      phistar = params[2]
      sbstar = params[3]
      sigmasb = params[4]
      beta = params[5]
  ENDIF

  ;; Convert 1-d inputs to 2-d if necessary
  IF size(absmag, /n_dimensions) EQ 1 THEN $
     absmag_arr = double(absmag # (sb * 0 + 1)) $
  ELSE $
     absmag_arr = double(absmag)
  IF size(sb, /n_dimensions) EQ 1 THEN $
     sb_arr = double((absmag * 0 + 1) # sb) $
  ELSE $
     sb_arr = double(sb)
  
  phi = 0.4 * alog(10) / sqrt(2 * !PI) / sigmasb * phistar $
        * 10 ^ (0.4 * (alpha + 1) * (mstar - absmag_arr)) $
        * exp(-10 ^ (0.4 * (mstar - absmag_arr))) $
        * exp(-0.5 * ((sb_arr - sbstar - beta * (absmag_arr - mstar)) $
                      / sigmasb) ^ 2)

  RETURN, phi
END
