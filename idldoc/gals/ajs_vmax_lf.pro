; docformat = 'rst'
;+
; This procedure calculates the 1/Vmax luminosity function from input
; arrays of absolute magnitude and weights (e.g., 1/Vmax)
;-


;+
; This procedure calculates the 1/Vmax luminosity function from input
; arrays of absolute magnitude and weights (e.g., 1/Vmax)
;
; Bins must have equal width
;
; :Returns: fltarr/dblarr
;    Returns phi, value of luminosity function
; :Params:
;    absmag : in, required
;    weights : in, optional
;       E.g., 1./Vmax. Set to a single value to give all galaxies that
;       weight, or leave blank to give all galaxies a weight of 1.
; :Keywords:
;    bincentres : in, out, optional
;       Centre of each bin in absolute magnitude
;    err_phi : out, optional
;       Poisson errors = phi / sqrt(n)
;    ninbin : out, optional
;       Number of galaxies in each bin
;    nbins : in, optional
;       Number of bins (ignored if bincentres set)
;    xrange : in, optional
;       Minimum and maximum absmag (ignored if bincentres set)
;    jackknife : in, optional
;       Integer for each galaxy giving the jackknife region in which
;       the galaxy lies
;    lum_dens : out, optional
;       Luminosity density.  Return value is sum(lum * weight),
;       j = ... x 10^(0.4 M_sun) h L_sun Mpc^-3
;    err_lum_dens : out, optional
;       Error on luminosity density, calculated from errors on phi, or
;       from jackknife resampling
;    loglum : in, optional
;       For Schechter fits: /loglum indicates log(luminosity) inputs
;    schechter : out, optional
;       Returns best-fitting Schechter function (using jackknife if present)
;    cov_mat : out, optional
;       Returns covariance matrix of Schechter function fit (using
;       jackknife if present)
;    sch_range : in, optional
;       Range of absolute magnitudes for (double) Schechter function fit
;    sch_lum_dens : out, optional
;       As lum_dens, but for the Schechter function fit
;    sch_err_lum_dens : out, optional
;       Error on Schechter luminosity density
;    double_schechter : out, optional
;       Returns best-fitting double Schechter function (using
;       jackknife if present)
;    double_cov_mat : out, optional
;       Returns covariance matrix of double Schechter function fit (using
;       jackknife if present)
;    double_sch_lum_dens : out, optional
;       As lum_dens, but for the double Schechter function fit
;    double_sch_err_lum_dens : out, optional
;       Error on double Schechter luminosity density
; :History:
;    5 Sep 2007 Created (Anthony Smith)
;
;    6 Mar 2008 Re-written
;
;    18 Mar 2008 Added jackknife
;
;    7 Apr 2008 Jackknife estimation of Schechter function fit
;
;    15 Apr 2008 Added luminosity density
;
;    21 Apr 2008 Added double Schechter function
;-
FUNCTION ajs_vmax_lf, $
   absmag, weights, bincentres=bincentres, err_phi=err_phi, ninbin=ninbin, $
   nbins=nbins, xrange=xrange, jackknife=jackknife, lum_dens=lum_dens, $
   err_lum_dens=err_lum_dens, loglum=loglum, $
   schechter=schechter, cov_mat=cov_mat, $
   sch_range=sch_range, sch_lum_dens=sch_lum_dens, $
   sch_err_lum_dens=sch_err_lum_dens, $
   double_schechter=double_schechter, double_cov_mat=double_cov_mat, $
   double_sch_lum_dens=double_sch_lum_dens, $
   double_sch_err_lum_dens=double_sch_err_lum_dens
  compile_opt idl2
  debug = ajs_debug()
  IF debug GE 1 THEN message, 'Estimating 1/Vmax LF', /inf

  IF n_elements(weights) EQ 1 THEN $
     weights = replicate(weights, n_elements(absmag))
  IF n_elements(bincentres) EQ 0 THEN BEGIN
      IF n_elements(nbins) EQ 0 THEN $
         nbins = 24
      IF n_elements(xrange) GT 0 THEN $
         bincentres = ajs_linspace(min(xrange), max(xrange), nbins, $
                                   /bincentres) $
      ELSE $
         bincentres = ajs_linspace(min(absmag), max(absmag) + 1e-6, nbins, $
                                   /bincentres)
  ENDIF
  binsize = (max(bincentres) - min(bincentres)) / (n_elements(bincentres) - 1)
  m_min = min(bincentres) - binsize / 2.

  ;; Do the whole thing whether or not jackknife is going to be done
  phi = hist1d(absmag, weights, min=m_min, nbins=n_elements(bincentres), $
               binsize=binsize, obin=obin, omin=omin, omax=omax, $
               density=ninbin) / binsize
  err_phi = phi / sqrt(ninbin)
  IF arg_present(schechter) THEN BEGIN
      schechter = ajs_schechter_fit(bincentres, phi, err_phi, $
                                    cov_mat=cov_mat, range=sch_range, $
                                    lum_dens=sch_lum_dens, $
                                    err_lum_dens=sch_err_lum_dens, $
                                    loglum=loglum)
  ENDIF
  IF arg_present(double_schechter) THEN BEGIN
      double_schechter = ajs_double_schechter_fit( $
                         bincentres, phi, err_phi, $
                         cov_mat=double_cov_mat, range=sch_range, $
                         lum_dens=double_sch_lum_dens, $
                         err_lum_dens=double_sch_err_lum_dens, $
                         loglum=loglum)
  ENDIF
  IF keyword_set(loglum) THEN BEGIN 
;;       lum_dens = total(phi * 10. ^ bincentres) * binsize
      lum_dens = total(10. ^ absmag * weights)
      err_lum_dens = sqrt(total(((err_phi * 10. ^ bincentres) $
                                 * binsize) ^ 2))
  ENDIF ELSE BEGIN
;;       lum_dens = total(phi * 10. ^ ((-bincentres) / 2.5)) * binsize
      lum_dens = total(10. ^ ((-absmag) / 2.5) * weights)
      err_lum_dens = sqrt(total(((err_phi * 10. ^ ((-bincentres) / 2.5)) $
                                 * binsize) ^ 2))
  ENDELSE

  ;; Jackknife
  IF n_elements(jackknife) GT 0 THEN BEGIN
      IF debug GE 1 THEN message, 'Using jackknife resampling', /inf

      ;; Original (full) sample
      phi_orig = phi
      ninbin_orig = ninbin
      lum_dens_orig = lum_dens
      IF arg_present(schechter) THEN BEGIN
          schechter_orig = schechter
          sch_lum_dens_orig = sch_lum_dens
      ENDIF
      IF arg_present(double_schechter) THEN BEGIN
          double_schechter_orig = double_schechter
          double_sch_lum_dens_orig = double_sch_lum_dens
      ENDIF

      ;; Prepare for jackknife resampling
      min_jack = min(jackknife, max=max_jack)
      jacks = min_jack + $
              (indgen(max_jack - min_jack + 1))[where(histogram(jackknife) $
                                                      GT 0, n_jack)]
      phi_arr = dblarr(n_elements(bincentres), n_jack)
      ninbin_arr = lonarr(n_elements(bincentres), n_jack)
      lum_dens_arr = dblarr(n_jack) 
      IF arg_present(schechter) THEN BEGIN
          schechter_arr = dblarr(3, n_jack)
          sch_lum_dens_arr = dblarr(n_jack)
      ENDIF
      IF arg_present(double_schechter) THEN BEGIN
          double_schechter_arr = dblarr(5, n_jack)
          double_sch_lum_dens_arr = dblarr(n_jack)
      ENDIF

      ;; Results excluding each jackknife sample
      FOR i = 0, n_jack - 1 DO BEGIN
          current_sample = where(jackknife NE jacks[i], n_current)
          phi_arr[*, i] = hist1d(absmag[current_sample], $
                                 weights[current_sample], $
                                 min=m_min, nbins=n_elements(bincentres), $
                                 binsize=binsize, obin=obin, omin=omin, $
                                 omax=omax, density=ninbin) / binsize $
                          ;; NB weights adjusted according to no. of galaxies
                          ;; If equal area: * n_jack / (n_jack - 1)
                          * n_elements(absmag) / n_current
          err_phi_tmp = phi_arr[*, i] / sqrt(ninbin)
          ninbin_arr[*, i] = ninbin
          IF keyword_set(loglum) THEN $
             lum_dens_arr[i] = total(10. ^ absmag[current_sample] $
                                     * weights[current_sample]) $
                               * n_elements(absmag) / n_current $
;;              lum_dens_arr[i] = total(phi_arr[*, i] * 10. ^ bincentres) $
;;                                * binsize $
          ELSE $
             lum_dens_arr[i] = total(10. ^ ((-absmag[current_sample]) / 2.5) $
                                     * weights[current_sample]) $
                               * n_elements(absmag) / n_current
;;              lum_dens_arr[i] = total(phi_arr[*, i] $
;;                                      * 10. ^ ((-bincentres) / 2.5)) * binsize
          IF arg_present(schechter) THEN BEGIN
              schechter_arr[*, i] = $
                 ajs_schechter_fit(bincentres, phi_arr[*, i], err_phi_tmp, $
                                   range=sch_range, lum_dens=sch_lum_dens, $
                                   loglum=loglum)
              sch_lum_dens_arr[i] = sch_lum_dens
          ENDIF
          IF arg_present(double_schechter) THEN BEGIN
              double_schechter_arr[*, i] = $
                 ajs_double_schechter_fit( $
                 bincentres, phi_arr[*, i], err_phi_tmp, $
                 range=sch_range, lum_dens=double_sch_lum_dens, $
                 loglum=loglum)
              double_sch_lum_dens_arr[i] = double_sch_lum_dens
          ENDIF
      ENDFOR

      ;; Jackknife quantities
      phi = dblarr(n_elements(bincentres))
      err_phi = dblarr(n_elements(bincentres))
      FOR i = 0, n_elements(bincentres) - 1 DO BEGIN
          err_phi[i] = ajs_jackknife(phi_arr[i, *], bias_correction=bc_tmp, $
                                     original_estimate=phi_orig[i])
          phi[i] = phi_orig[i] + bc_tmp
      ENDFOR
      err_lum_dens = ajs_jackknife(lum_dens_arr, bias_correction=bc_tmp, $
                                   original_estimate=lum_dens_orig)
      lum_dens = lum_dens_orig + bc_tmp
      IF arg_present(schechter) THEN BEGIN
          schechter = dblarr(3)
          cov_mat = ajs_jackknife(schechter_arr, bias_correction=bc_tmp, $
                                  original_estimate=schechter_orig)
          schechter = schechter_orig + bc_tmp
          sch_err_lum_dens = ajs_jackknife(sch_lum_dens_arr, $
                                           bias_correction=bc_tmp, $
                                           original_estimate=sch_lum_dens_orig)
          sch_lum_dens = sch_lum_dens_orig + bc_tmp
      ENDIF
      IF arg_present(double_schechter) THEN BEGIN
          double_schechter = dblarr(5)
          double_cov_mat = ajs_jackknife(double_schechter_arr, $
                                  bias_correction=bc_tmp, $
                                  original_estimate=double_schechter_orig)
          double_schechter = double_schechter_orig + bc_tmp
          double_sch_err_lum_dens = ajs_jackknife( $
                                    double_sch_lum_dens_arr, $
                                    bias_correction=bc_tmp, $
                                    original_estimate=double_sch_lum_dens_orig)
          double_sch_lum_dens = double_sch_lum_dens_orig + bc_tmp
      ENDIF

      ;; Number in each bin (no jackknife)
      ninbin = ninbin_orig
  ENDIF

  ;; Empty bins: NaN
  empty_bins = where(ninbin EQ 0, nempty)
  IF nempty GT 0 THEN BEGIN
      phi[empty_bins] = !values.f_nan
      err_phi[empty_bins] = !values.f_nan
  ENDIF

  return, phi
END


;+
; Test ajs_vmax_lf
;-
PRO ajs_vmax_lf_test
  compile_opt idl2

  absmag = ajs_linspace(-26, -16, 50)
  weights = intarr(n_elements(absmag)) + 1e-6 
  phi = ajs_vmax_lf(absmag, weights, nbins=30, ninbin=ninbin)
  print, phi
  print, total(ninbin)

END
