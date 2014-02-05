; docformat = 'idl'

;+
; NAME:
;	ajs_double_choloniewski
;
;
; PURPOSE:
;	This function generates a double Choloniewski function along an input
;	array of absolute magnitudes and surface brightnesses.
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;	result = ajs_schechter(absmag,sb,params)
;
; INPUTS:
;	absmag
;	sb
;
;
; OPTIONAL INPUTS:
;	params:	[mstar,alpha,phistar,sbstar,sigmasb,beta] twice
;
; KEYWORD PARAMETERS:
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
;	15 Feb 2008 Created, Anthony Smith
;-
FUNCTION ajs_double_choloniewski, absmag, sb, params
  compile_opt idl2

  params1 = params[0:5]
  params2 = params[6:11]

  phi = ajs_choloniewski(absmag, sb, params1) $
        + ajs_choloniewski(absmag, sb, params2)

  RETURN, phi
END
