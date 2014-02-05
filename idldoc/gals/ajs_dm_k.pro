; docformat = 'rst'
;+
; This function returns some combination of the distance modulus,
; K-correction and evolution correction for a galaxy at a redshift z.
;-

;+
; Includes K-correction if band keyword is set, and evolution correction if q0
; and/or q1 is set. K(z) and/or E(z) without the distance modulus if
; /nodm is set.
;
; Uses Michael Blanton's K-correct code, with coefficients of a
; typical galaxy for the K-corrections.
;
; :Returns: DM(z) + K(z) - E(z), where E(z) = q0 * (1 + q1 * z) * z
;
; :Params:
;    z : in, required, type="fltarr or dblarr"
; :Keywords:
;    band : in, optional, type=string
;       Band for K-corrections: choose from ugrizYJHK
;    nodm : in, optional, type=bool
;       Set /nodm to return K(z) - E(z) only (no distance modulus)
;    q0 : in, optional, type=float
;    q1 : in, optional, type=float
;       q0 and q1 define the evolution correction, E(z) = q0 * (1 + q1
;       * z) * z
;    coeffs : in, optional
;       kcorrect coefficients. Default to typical coefficients for a
;       galaxy. Array of 5 elements to give same coefficients to all
;       galaxies, or array of [5, n_elements(z)] to specify
;       coefficients for each galaxy
;    kc : out, optional
;       k (+ e) corrections, without distance modulus
;    vmatrix : in, out, optional
;       supply a variable name here to speed up multiple executions
;    lambda : in, out, optional
;       supply a variable name here to speed up multiple executions
;    rmatrix : in, out, optional
;       for repeated execution with same z, same band, but with
;       different coeffs, supply a variable name here
;    zvals : in, out, optional
;       for repeated execution with same z, same band, but with
;       different coeffs, supply a variable name here
; :Author: Anthony Smith
; :History: 25 Jan 2008: Added evolution corrections
;
;    5 Mar 2008: Added coeffs, kc, vmatrix, lambda, rmatrix, zvals keywords
;-
FUNCTION ajs_dm_k, z, band=band, nodm=nodm, q0=q0, q1=q1, coeffs=coeffs, $
                   kc=kc, vmatrix=vmatrix, lambda=lambda, rmatrix=rmatrix, $
                   zvals=zvals
  compile_opt idl2
  debug = ajs_debug()
  IF debug GE 3 THEN BEGIN
      message, 'Finding DM(z)/K-corrections/E-corrections', /inf
      message, ajs_kw_string(q0=q0, q1=q1, nodm=nodm, $
                             band=band), /inf
  ENDIF

  ;; Number of K-correct templates
  n_templates = 5

  ;; Check dimensions of coeffs
  IF n_elements(coeffs) GT 0 THEN BEGIN
      IF n_elements(coeffs) NE n_templates THEN BEGIN          
          IF product(size(coeffs, /dimensions) $
                     EQ [n_templates, n_elements(z)]) EQ 0 THEN $
                      message, 'Incorrect dimensions for coeffs' $
          ELSE BEGIN
              IF debug GE 3 THEN message, 'Coeffs for each galaxy', /inf
              many_coeffs = 1
          ENDELSE
      ENDIF ELSE BEGIN
          IF debug GE 3 THEN message, 'Coeffs for 1 galaxy: duplicating', /inf
          many_coeffs = 0
      ENDELSE
  ENDIF ELSE BEGIN
      IF debug GE 3 THEN message, 'Coeffs not specified: default', /inf
      ;; Typical coefficients: median(table.coeffs, dimension=2)
      many_coeffs = 0
      coeffs = [1.89821e-05, 2.37252e-09, 1.97939e-06, 7.54209e-08, $
                1.93300e-08]
;;       coeffs = [6.585688e-06, 8.208907e-10, 1.208443e-06, 3.311702e-05, $
;;                 2.723333e-07]
  ENDELSE

  ;; Choose filter
  IF n_elements(band) EQ 0 THEN BEGIN
      kc = 0
  ENDIF ELSE BEGIN
      IF band EQ 'u' THEN filterlist = ['sdss_u0.par']
      IF band EQ 'g' THEN filterlist = ['sdss_g0.par']
      IF band EQ 'r' THEN filterlist = ['sdss_r0.par']
      IF band EQ 'i' THEN filterlist = ['sdss_i0.par']
      IF band EQ 'z' THEN filterlist = ['sdss_z0.par']
      IF band EQ 'Y' THEN filterlist = ['wfcam_Y.par']
      IF band EQ 'J' THEN filterlist = ['wfcam_J.par']
      IF band EQ 'H' THEN filterlist = ['wfcam_H.par']
      IF band EQ 'K' THEN filterlist = ['wfcam_K.par']

      ;; Prepare for applying K-corrections
      IF n_elements(rmatrix) EQ 0 OR n_elements(zvals) EQ 0 THEN BEGIN  
          IF n_elements(vmatrix) EQ 0 OR n_elements(lambda) EQ 0 THEN $ 
             k_load_vmatrix, vmatrix, lambda
          ;zmin = min(z)
          ;zmax = max(z)
          ;nz = n_elements(z) + 10000
          k_projection_table, rmatrix, vmatrix, lambda, zvals, filterlist, $
                              zmin=zmin, zmax=zmax, nz=nz, /silent
      ENDIF

      ;; Apply K-corrections
      IF debug GE 3 THEN message, 'Starting reconstruct_maggies', /inf
      IF many_coeffs THEN BEGIN
          k_reconstruct_maggies, reform(coeffs, n_elements(coeffs), 1), $
                                 replicate(0, n_elements(z)), $
                                 flux0, $
                                 rmatrix=rmatrix, zvals=zvals
          k_reconstruct_maggies, reform(coeffs, n_elements(coeffs), 1), $
                                 z, $
                                 fluxz, $
                                 rmatrix=rmatrix, zvals=zvals
      ENDIF ELSE BEGIN
          k_reconstruct_maggies, [[coeffs], $
                                  [coeffs # (dblarr(n_elements(z)) + 1)]], $
                                 [0, z], $
                                 flux, $
                                 rmatrix=rmatrix, zvals=zvals
          flux0 = flux[0]
          fluxz = flux[1:*]
      ENDELSE

      ;; K-corrections
      kc = 2.5 * alog10(flux0 / fluxz)
      kc = reform(kc)      
  ENDELSE

  ;; Evolution corrections
  IF n_elements(q0) GT 0 THEN BEGIN
      IF n_elements(q1) GT 0 THEN BEGIN
          kc = kc - q0 * (1 + q1 * z) * z
      ENDIF ELSE BEGIN
          kc = kc - q0 * z
      ENDELSE
  ENDIF

  IF keyword_set(nodm) THEN $
     dm_k = kc $
  ELSE $                        
     dm_k = ajs_distmod(z) + kc      ; Distance modulus + K(z) = m - M

  RETURN, dm_k

END 

;+
; Test ajs_dm_k
;-
PRO ajs_dm_k_test

  z = ajs_linspace(0, 0.5, 201)
  dm_k = ajs_dm_k(z, band='K', /nodm)
  plot, z, dm_k, yrange=[-2, 0]
  dm_k = ajs_dm_k(z, band='K', /nodm, q0=1)
  oplot, z, dm_k, color=1
  dm_k = ajs_dm_k(z, band='K', /nodm, q0=1, q1=1)
  oplot, z, dm_k, color=2

  z = [0.05, 0.1, 0.15]
  coeffs = [1.585688e-06, 1.208907e-10, 0, 0, 0]
  print, 'should be    -0.0816033    -0.147924    -0.199019:'
  print, ajs_dm_k(z, band='K', /nodm, coeffs=coeffs)
  z = [0.1, 0.1, 0.1]
  coeffs = [[1.585688e-06, 1.208907e-10, 0, 0, 0], $
            [1.585688e-06, 1.208907e-10, 0, 0, 0], $
            [1.585688e-06, 1.208907e-10, 0, 0, 0]]
  print, 'should be  -0.147924    -0.147924    -0.147924'
  print, ajs_dm_k(z, band='K', /nodm, coeffs=coeffs)
  coeffs = [[6.585688e-06, 8.208907e-10, 1.208443e-06, 3.311702e-05, $
             2.723333e-07], $
            [5.585688e-06, 7.208907e-10, 1.308443e-06, 3.311702e-05, $
             2.723333e-07], $
            [1.585688e-06, 1.208907e-10, 0, 0, 0]]
  print, 'should be -0.198890    -0.202550    -0.147924'
  print, ajs_dm_k(z, band='K', /nodm, coeffs=coeffs)

END
