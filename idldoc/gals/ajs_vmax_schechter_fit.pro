; docformat = 'idl'
;+
; NAME:
;	ajs_vmax_schechter_fit
;
; PURPOSE:
;	This function calculates 1/Vmax LF and finds best-fit
;	Schechter function for input arrays of absolute magnitude and Vmax.
;
; CATEGORY:
;
; CALLING SEQUENCE:
;	Result = ajs_vmax_schechter_fit(absmag, vmax, bincentres,
;		meanbin, phi, phierr, phierr2, ninbin, absmagfaint, absmagbright
;
; INPUTS:
;	absmag -> ninbin: see ajs_vmax_lf.pro
;	absmagfaint	Faint absolute magnitude limit for Schechter fit
;	absmagbright	Bright " "
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	Result: [mstar,alpha,phistar]
;	See ajs_vmax_lf for more
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;	Calls ajs_vmax_lf
;	Calls ajs_schechter_fit
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;	5 Sep 2007: Created, Anthony Smith
;-
FUNCTION ajs_vmax_schechter_fit,absmag,vmax,bincentres,meanbin,phi,phierr,phierr2,ninbin,absmagfaint,absmagbright
  compile_opt idl2

  ajs_vmax_lf,absmag,vmax,bincentres,meanbin,phi,phierr,phierr2,ninbin

  ;; Best fit Schechter function: observed (within absmag range)
  good=where(meanbin GE absmagbright $
             AND meanbin LE absmagfaint $
             AND finite(phi))

  params=ajs_schechter_fit(meanbin[good], phi[good], phierr2[good])

  RETURN,params
END
