; docformat = 'rst'
;+
; This procedure calculates the 1/Vmax bivariate brightness
; distribution from input arrays of two quantities (e.g., absolute
; magnitude and surface brightness) and weights (e.g., 1/Vmax)
;-


;+
; This procedure calculates the 1/Vmax bivariate brightness
; distribution from input arrays of two quantities (e.g., absolute
; magnitude and surface brightness) and weights (e.g., 1/Vmax)
;
; Bins in each dimension must have equal width
;
; :Returns: fltarr/dblarr
;    Returns phi, value of space density
; :Params:
;    m1 : in, required
;    m2 : in, required
;    weights : in, optional
;       E.g., 1./Vmax. Set to single value to give all galaxies the
;       same weight, or leave empty to give all galaxies a weight of 1.
; :Keywords:
;    bincentres1 : in, out, optional
;       Centre of each bin in m1
;    bincentres2 : in, out, optional
;       Centre of each bin in m2
;    err_phi : out, optional
;       Poisson errors = phi / sqrt(n)
;    ninbin : out, optional
;       Number of galaxies in each bin
;    nbins1 : in, optional
;       Number of bins in m1 (ignored if bincentres set)
;    nbins2 : in, optional
;       Number of bins in m2 (ignored if bincentres set)
;    xrange : in, optional
;       Minimum and maximum m1 (ignored if bincentres1 set)
;    yrange : in, optional
;       Minimum and maximum m2 (ignored if bincentres2 set)
;    jackknife : in, optional
;       Integer for each galaxy giving the jackknife region in which
;       the galaxy lies
;    choloniewski : in, out, optional
;       Choloniewski fit
;    chol_xrange : in, optional, type=fltarr(2)
;       xrange for Choloniewski fit (ajs_choloniewski_fit)
;    chol_yrange : in, optional, type=fltarr(2)
;       yrange for Choloniewski fit (ajs_choloniewski_fit)
;    cov_mat : out, optional
;       Covariance matrix of Choloniewski fit
;    vars : in, optional
;       See ajs_choloniewski_fit
;    zmin : in, optional
;       See ajs_choloniewski_fit
;    zmax : in, optional
;       See ajs_choloniewski_fit
;    area : in, optional
;       See ajs_choloniewski_fit
;    q0 : in, optional
;       See ajs_choloniewski_fit
;    q1 : in, optional
;       See ajs_choloniewski_fit
; :History:
;    6 Mar 2008 Created (Anthony Smith)
;
;    18 Mar 2008 Added jackknife
;
;    7 Apr 2008 Added jackknife Choloniewski fit
;-
FUNCTION ajs_vmax_bbd, m1, m2, weights, bincentres1=bincentres1, $
                       bincentres2=bincentres2, $
                       err_phi=err_phi, ninbin=ninbin, nbins1=nbins1, $
                       nbins2=nbins2, xrange=xrange, yrange=yrange, $
                       jackknife=jackknife, $
                       choloniewski=choloniewski, chol_xrange=chol_xrange, $
                       chol_yrange=chol_yrange, cov_mat=cov_mat, $
                       vars=vars, zmin=zmin, zmax=zmax, area=area, q0=q0, q1=q1

  compile_opt idl2
  debug = ajs_debug()
  IF debug GE 1 THEN BEGIN
      message, 'Estimating 1/Vmax BBD', /inf
      message, ajs_kw_string(nbins1=nbins1, nbins2=nbins2, xrange=xrange, $
                             yrange=yrange, choloniewski=choloniewski, $
                             chol_xrange=chol_xrange, $
                             chol_yrange=chol_yrange, $
                             zmin=zmin, zmax=zmax, area=area, q0=q0, q1=q1), $
               /inf
  ENDIF

  ;; Preliminaries
  IF n_elements(weights) EQ 1 THEN $
     weights = replicate(weights, n_elements(m1)) $ 
  ELSE IF n_elements(weights) EQ 0 THEN $
     weights = replicate(1, n_elements(m1))
  IF n_elements(bincentres1) EQ 0 THEN BEGIN
      IF n_elements(nbins1) EQ 0 THEN $
         nbins1 = 24
      IF n_elements(xrange) GT 0 THEN $
;;          bincentres1 = ajs_linspace(max([min(xrange), min(m1)]), $
;;                                     min([max(xrange), max(m1) + 1d-5]), $
         bincentres1 = ajs_linspace(min(xrange), max(xrange), $
                                    nbins1, $
                                    /bincentres) $
      ELSE $
         bincentres1 = ajs_linspace(min(m1), max(m1) + 1e-6, nbins1, $
                                    /bincentres)
  ENDIF
  IF n_elements(bincentres2) EQ 0 THEN BEGIN
      IF n_elements(nbins2) EQ 0 THEN $
         nbins2 = 24
      IF n_elements(yrange) GT 0 THEN $
;;          bincentres2 = ajs_linspace(max([min(yrange), min(m2)]), $
;;                                     min([max(yrange), max(m2) + 1d-5]), $
         bincentres2 = ajs_linspace(min(yrange), max(yrange), $
                                    nbins2, $
                                    /bincentres) $
      ELSE $
         bincentres2 = ajs_linspace(min(m2), max(m2) + 1d-5, nbins2, $
                                    /bincentres)
  ENDIF
  binsize1 = (max(bincentres1) - min(bincentres1)) $
             / (n_elements(bincentres1) - 1)
  m1_min = min(bincentres1) - binsize1 / 2.
  m1_max = max(bincentres1) + binsize1 / 2. - 1d-5
  binsize2 = (max(bincentres2) - min(bincentres2)) $
             / (n_elements(bincentres2) - 1)
  m2_min = min(bincentres2) - binsize2 / 2.
  m2_max = max(bincentres2) + binsize2 / 2. - 1d-5

  ;; Do the whole thing whether or not jackknife is going to be done
  incl = where(m1 GE m1_min AND m1 LE m1_max $
               AND m2 GT m2_min AND m2 LE m2_max)
  phi = hist2d(m1[incl], m2[incl], weights[incl], $
               min1=m1_min, max1=m1_max, $
               binsize1=binsize1, obin1=obin1, omin1=omin1, omax1=omax1, $
               min2=m2_min, max2=m2_max, $
               binsize2=binsize2, obin2=obin2, omin2=omin2, omax2=omax2, $
               density=ninbin) / binsize1 / binsize2
  err_phi = phi / sqrt(ninbin)
  IF arg_present(choloniewski) THEN BEGIN
      choloniewski = ajs_choloniewski_fit(bincentres1, bincentres2, $
                                          phi, err_phi, $
                                          xrange=chol_xrange, $
                                          yrange=chol_yrange, vars=vars, $
                                          zmin=zmin, zmax=zmax, area=area, $
                                          q0=q0, q1=q1, cov_mat=cov_mat)
  ENDIF

  ;; Check okay
  IF n_elements(obin1) NE n_elements(bincentres1) $
     OR n_elements(obin2) NE n_elements(bincentres2) THEN $
     message, 'Conflicting bins'
  IF total(ninbin) NE n_elements(incl) THEN $
     message, 'Different total number after making histogram'

  IF n_elements(jackknife) GT 0 THEN BEGIN
      IF debug GE 1 THEN message, 'Using jackknife resampling', /inf

      ;; Original (full) sample
      phi_orig = phi
      err_phi_orig = err_phi
      ninbin_orig = ninbin
      IF arg_present(choloniewski) THEN $
         choloniewski_orig = choloniewski
      
      ;; Prepare for jackknife resampling
      min_jack = min(jackknife)
      n_jack = max(jackknife) - min_jack + 1
      phi_arr = dblarr(n_elements(bincentres1), n_elements(bincentres2), $
                       n_jack)
      ninbin_arr = dblarr(n_elements(bincentres1), n_elements(bincentres2), $
                          n_jack)
      IF arg_present(choloniewski) THEN $
         choloniewski_arr = dblarr(6, n_jack)

      ;; Results excluding each jackknife sample in turn
      FOR i = 0, n_jack - 1 DO BEGIN
          IF debug GE 2 THEN message, 'Jackknife ' + string(i), /inf
          incl = where(m1 GE m1_min AND m1 LE m1_max $
                       AND m2 GT m2_min AND m2 LE m2_max $
                       AND jackknife NE i + min_jack, n_current)
          phi_arr[*, *, i] = $
             hist2d(m1[incl], m2[incl], weights[incl], $
                    min1=m1_min, max1=m1_max, $
                    binsize1=binsize1, obin1=obin1, omin1=omin1, omax1=omax1, $
                    min2=m2_min, max2=m2_max, $
                    binsize2=binsize2, obin2=obin2, omin2=omin2, omax2=omax2, $
                    density=ninbin) / binsize1 / binsize2 $
             ;; NB weights adjusted according to no. of galaxies
             ;; If equal area: * n_jack / (n_jack - 1)
             * n_elements(m1) / n_current
          ;; NB original err_phi may have been altered by chol. fit to
          ;; supply values for empty bins: don't want to do
          ;; that every time!
          err_phi_tmp = phi_arr[*, *, i] / sqrt(ninbin)
          no_err = where(err_phi_tmp LE 0 OR finite(err_phi_tmp, /nan), $
                         n_no_err)
          IF n_no_err GT 0 THEN $
             err_phi_tmp[no_err] = err_phi_orig[no_err] $
                                   * sqrt(total(ninbin_orig) / total(ninbin))
          ninbin_arr[*, *, i] = ninbin
          IF arg_present(choloniewski) THEN BEGIN
              choloniewski_arr[*, i] = $
                 ajs_choloniewski_fit(bincentres1, bincentres2, $
                                      phi_arr[*, *, i], err_phi_tmp, $
                                      xrange=chol_xrange, $
                                      yrange=chol_yrange, vars=vars, $
                                      zmin=zmin, zmax=zmax, area=area, $
                                      q0=q0, q1=q1)
          ENDIF
      ENDFOR

      ;; Jackknife quantities
      phi = dblarr(n_elements(bincentres1), n_elements(bincentres2))
      err_phi = dblarr(n_elements(bincentres1), n_elements(bincentres2))
      FOR i = 0, n_elements(bincentres1) - 1 DO BEGIN
          FOR j = 0, n_elements(bincentres2) - 1 DO BEGIN
              err_phi[i, j] = ajs_jackknife(phi_arr[i, j, *], $
                                            bias_correction=bc_tmp, $
                                            original_estimate=phi_orig[i, j])
              phi[i, j] = phi_orig[i, j] + bc_tmp
;;               err_phi[i, j] = ajs_jackknife(phi_arr[i, j, *], $
;;                                             jackknife_mean=phi_tmp)
;;               phi[i, j] = phi_tmp
          ENDFOR
          ;; No zero errors
          zeros = where(err_phi EQ 0, n_zeros)
          IF n_zeros GT 0 THEN $
             err_phi[zeros] = err_phi_orig[zeros]
      ENDFOR

      IF arg_present(choloniewski) THEN BEGIN
          choloniewski = dblarr(6)
          cov_mat = ajs_jackknife(choloniewski_arr, bias_correction=bc_tmp, $
                                  original_estimate=choloniewski_orig)
          choloniewski = choloniewski_orig + bc_tmp
      ENDIF
      
      ;; Number in each bin (no jackknife)
      ninbin = ninbin_orig
;;       incl = where(m1 GE m1_min AND m1 LE m1_max $
;;                    AND m2 GT m2_min AND m2 LE m2_max)
;;       ninbin = $
;;          hist2d(m1[incl], m2[incl], $
;;                 min1=m1_min, max1=m1_max, $
;;                 binsize1=binsize1, obin1=obin1, omin1=omin1, omax1=omax1, $
;;                 min2=m2_min, max2=m2_max, $
;;                 binsize2=binsize2, obin2=obin2, omin2=omin2, omax2=omax2)
  ENDIF ;ELSE BEGIN
;;       incl = where(m1 GE m1_min AND m1 LE m1_max $
;;                    AND m2 GT m2_min AND m2 LE m2_max)
;;       phi = hist2d(m1[incl], m2[incl], weights[incl], $
;;                    min1=m1_min, max1=m1_max, $
;;                    binsize1=binsize1, obin1=obin1, omin1=omin1, omax1=omax1, $
;;                    min2=m2_min, max2=m2_max, $
;;                    binsize2=binsize2, obin2=obin2, omin2=omin2, omax2=omax2, $
;;                    density=ninbin) / binsize1 / binsize2
;;       err_phi = phi / sqrt(ninbin)
;;   ENDELSE

  IF debug GE 1 THEN $
     message, 'Total number of galaxies: ' + strtrim(total(ninbin), 2), /inf

;;   empty_bins = where(ninbin EQ 0, nempty)
;;   IF nempty GT 0 THEN BEGIN
;;       phi[empty_bins] = 0
;;       err_phi[empty_bins] = !values.f_nan
;;   ENDIF

  return, phi
END


;+
; Test ajs_vmax_lf
;-
PRO ajs_vmax_bbd_test
  compile_opt idl2

  m1 = ajs_linspace(-26, -16, 50)
  m2 = ajs_linspace(14, 24, 50)
  weights = intarr(n_elements(m1)) + 1e-6 
  phi = ajs_vmax_bbd(m1, m2, weights, nbins1=30, nbins2=30, ninbin=ninbin)
  print, phi
  print, total(ninbin)

END
