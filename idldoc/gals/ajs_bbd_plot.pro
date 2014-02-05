; docformat = 'rst'

;+
; This procedure produces a surface/contour plot of the
; bivariate brightness distribution (number of galaxies per unit
; volume as a function of absolute magnitude and surface brightness)
;-

;+
; Main procedure.
;
; :Params:
;    bincentres1 : in, out, optional
;       Centres of absolute magnitude bins (x-axis). Optional
;       if calculating BBD from m1, m2 and weights
;    bincentres2 : in, out, optional
;       Centres of surface brightness bins (y-axis). Optional
;       if calculating BBD from m1, m2 and weights
;    phibbd : in, optional
;       Value of phi at bincentres1[i], bincentres2[j]. Not required
;       if calculating BBD from m1, m2 and weights.
;
; :Keywords:
;    phibbderr : in, optional
;       1-sigma errors in phibbd: converted to 1-sigma errors in
;       corresponding chi^2, then plotted as dotted/dashed contours. Set
;       /phibbderr, or supply variable, with m1, m2 to calculate then
;       plot errors
;    error_contours : in, optional
;       Set /error_contours to force plotting of errors (switches off
;       if Choloniewski function is plot)
;    ct : in, optional
;       Color table to use (default 1)
;    plotfile : in, optional
;       Destination for output EPS plot
;    show_plot : in, optional
;       Open EPS plot for display
;    band : in, optional, type=string
;       Name of waveband (for axis titles)
;    title : in, optional
;    xtitle : in, optional
;    ytitle : in, optional
;    colorbartitle : in, optional
;    min_log_phi : in, optional
;    max_log_phi : in, optional
;    choloniewski : in, out, optional
;       If a Choloniewski fit is performed, error contours will not be
;       displayed, unless /error_contours is set
;    chol_xrange : in, optional, type=fltarr(2)
;       xrange for Choloniewski fit (ajs_choloniewski_fit)
;    chol_yrange : in, optional, type=fltarr(2)
;       yrange for Choloniewski fit (ajs_choloniewski_fit)
;    cov_mat : out, optional
;       Covariance matrix of Choloniewski fit
;    m1 : in, optional
;       Values of first magnitude for each galaxy (if phibbd empty)
;    m2 : in, optional
;       Values of second magnitude for each galaxy (if phibbd empty)
;    weights : in, optional
;       Weight (e.g., 1/Vmax) for each galaxy (if phibbd empty). See
;       ajs_vmax_bbd for more details.
;    jackknife : in, optional
;       Integer for each galaxy giving the jackknife region in which
;       the galaxy lies
;    xrange : in, optional
;       Range of plot
;    yrange : in, optional
;       Range of plot
;    nbins1 : in, optional
;       Number of bins for BBD calculation from m2 and weights
;    nbins2 : in, optional
;       Number of bins for BBD calculation from m2 and weights
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
; :Examples:
;    ajs_bbd_plot, bincentres1, bincentres2, phibbd
;
;    ajs_bbd_plot, m1=abs_mag, m2=surface_brightness, weights=weights
; :History:
;    11 Jan 2008 Created, Anthony Smith
;
;    18 Mar 2008 Added jackknife
;
;    25 Mar 2008 Added keywords for Choloniewski fit
;
;    7 Apr 2008 Choloniewski included in ajs_vmax_bbd call
;-
PRO ajs_bbd_plot, bincentres1, bincentres2, phibbd, phibbderr=phibbderr, $
                  error_contours=error_contours, $
                  plotfile=plotfile, ct=ct, show_plot=show_plot, $
                  band=band, title=title, xtitle=xtitle, ytitle=ytitle, $
                  colorbartitle=colorbartitle, min_log_phi=min_log_phi, $
                  max_log_phi=max_log_phi, $
                  choloniewski=choloniewski, chol_xrange=chol_xrange, $
                  chol_yrange=chol_yrange, cov_mat=cov_mat, $
                  m1=m1, m2=m2, weights=weights, jackknife=jackknife, $
                  xrange=xrange, yrange=yrange, $
                  nbins1=nbins1, nbins2=nbins2, $
                  vars=vars, zmin=zmin, zmax=zmax, area=area, q0=q0, q1=q1
  compile_opt idl2

  debug = ajs_debug()
  IF debug GE 1 THEN BEGIN
      message, 'BBD plot', /inf
      message, ajs_kw_string(zmin=zmin, zmax=zmax, area=area, q0=q0, q1=q1), $
               /inf
  ENDIF

  ;; Deal with parameters & keywords
  IF n_elements(error_contours) EQ 0 THEN $ 
     ;; Show error contours if phibbderr set but not if Choloniewski fit
     error_contours = (n_elements(phibbderr) GT 0 OR arg_present(phibbderr)) $
                      AND ~ ((n_elements(choloniewski) LE 1 $
                              AND (arg_present(choloniewski) $
                                   OR keyword_set(choloniewski))))
  IF n_elements(ct) EQ 0 THEN $
     ct = 1                     ; Colour table
  IF n_elements(band) EQ 0 THEN $
     band = ''
  IF n_elements(xtitle) EQ 0 THEN $ 
     xtitle = 'M!D' + band + '!N-5log h'
  IF n_elements(ytitle) EQ 0 THEN $ 
     ytitle = '!4l!3!D' + band + '!N'
  IF n_elements(colorbartitle) EQ 0 THEN $  
     colorbartitle = 'log !4u!3'

  ;; Calculate BBD from raw data?
  IF n_elements(phibbd) EQ 0 THEN BEGIN
      IF n_elements(m1) EQ n_elements(m2) THEN BEGIN
          ;; Include Choloniewski fit?
          IF (n_elements(choloniewski) LE 1 $
              AND (arg_present(choloniewski) $
                   OR keyword_set(choloniewski))) THEN BEGIN
              IF debug GE 2 THEN message, '1/Vmax with Choloniewski', /inf
              phibbd = $
                 ajs_vmax_bbd(m1, m2, weights, bincentres1=bincentres1, $
                              bincentres2=bincentres2, err_phi=phibbderr, $
                              xrange=xrange, yrange=yrange, $
                              nbins1=nbins1, nbins2=nbins2, $
                              jackknife=jackknife, choloniewski=choloniewski, $
                              chol_xrange=chol_xrange, $
                              chol_yrange=chol_yrange, vars=vars, $
                              zmin=zmin, zmax=zmax, area=area, $
                              q0=q0, q1=q1, cov_mat=cov_mat)
          ENDIF ELSE BEGIN
              IF debug GE 2 THEN message, '1/Vmax (no Choloniewski)', /inf
              phibbd = $
                 ajs_vmax_bbd(m1, m2, weights, bincentres1=bincentres1, $
                              bincentres2=bincentres2, err_phi=phibbderr, $
                              xrange=xrange, yrange=yrange, $
                              nbins1=nbins1, nbins2=nbins2, $
                              jackknife=jackknife)
          ENDELSE
      ENDIF
  ENDIF

  IF n_elements(min_log_phi) EQ 0 THEN BEGIN
      min_log_phi = floor(alog10(min(phibbd[where(phibbd GT 0)])))
      min_log_phi = max([min_log_phi, -10])
  ENDIF
  IF n_elements(max_log_phi) EQ 0 THEN $
     max_log_phi = ceil(alog10(max(phibbd[where(finite(phibbd))])))
  levels_phi_range_log = max_log_phi - min_log_phi
  
  ;; Contour levels (NB if the first is '' it substitutes numbers for all)
  lev_all = [1e-10, 10.^(-9.5), 1e-9, 10.^(-8.5), 1e-8, 10.^(-7.5), 1e-7, $
             10.^(-6.5), $
             1e-6, 10.^(-5.5), 1e-5, 10.^(-4.5), 1e-4, 10.^(-3.5), 1e-3, $
             10.^(-2.5), 1e-2, $
             10.^(-1.5), 1e-1, 10.^(-0.5), 1, 10.^(0.5), 1e1, 10.^(1.5), $
             1e2, 10.^(2.5), $
             1e3, 10.^(3.5), 1e4, 10.^(4.5), 1e5, 10.^(5.5), 1e6, 10.^(6.5), $
             1e7, 10.^(7.5), $
             1e8, 10.^(8.5), 1e9]
  c_annot_all = (['1e-10', '', '1e-9', '', '1e-8', '', '1e-7', '', '1e-6', $
                   '', '1e-5', '', '1e-4', '', '1e-3', '', '1e-2', '', $
                   '1e-1', '', '1e0', '', '1e1', '', '1e2', '', '1e3', '', $
                   '1e4', '', '1e5', '', '1e6', '', '1e7', '', '1e8', '', $
                   '1e9'])
  c_annotation = c_annot_all[ $
                 where(lev_all GE 10. ^ min_log_phi $
                       AND lev_all LE 10. ^ max_log_phi)]
  levels = lev_all[where(lev_all GE  10. ^ min_log_phi $
                         AND lev_all LE 10. ^ max_log_phi)]

  ;; Continuous contour levels
  levels_phi = 10 ^ (findgen(256) / 255 * levels_phi_range_log $
                     + min_log_phi)
  c_colors_phi = indgen(256)

  ;; Ranges, backwards (faint -> bright)
  binsize1 = (max(bincentres1) - min(bincentres1)) / n_elements(bincentres1) 
  binsize2 = (max(bincentres2) - min(bincentres2)) / n_elements(bincentres2) 
  IF n_elements(xrange) EQ 0 THEN $ 
     xrange = [max(bincentres1) + binsize1 / 2., $
               min(bincentres1) - binsize1 / 2.]
  IF n_elements(yrange) EQ 0 THEN $
     yrange = [max(bincentres2) + binsize2 / 2., $
               min(bincentres2) - binsize2 / 2.]

  ;; Plot
  eps = keyword_set(open) OR keyword_set(show_plot)
  ajs_plot_start, plotfile=plotfile, eps=eps, $
                  close=close, bits=8 ; 8 bits per colour = 24 bits

  ;; Axes
  set_basic_colours
  plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
        xtitle=xtitle, ytitle=ytitle, /xstyle, /ystyle, $
        position=[0.1275, 0.11, 0.862, 0.945], title=title

  ;; Shaded (filled) contours
  loadct, ct, /silent
  ;; Smooth contours
;;   xout = findgen(201) / 200 * (max(bincentres1) - min(bincentres1)) $
;;          + min(bincentres1)
;;   yout = findgen(201) / 200 * (max(bincentres2) - min(bincentres2)) $
;;          + min(bincentres2)
;;   contour, $
;;      min_curve_surf(phibbd, xvalues=bincentres1, yvalues=bincentres2, $
;;                     xout=xout, yout=yout), $
;;      xout, yout, /overplot, c_charsize=!P.charsize, color=1, $
;;      levels=levels_phi, $
;;      c_color=c_colors_phi, /fill, position=[0.1275, 0.11, 0.862, 0.945]
  contour, $
     phibbd, $
     bincentres1, bincentres2, /overplot, $;, color=1, $
     levels=levels_phi, $
     c_color=c_colors_phi, /fill, position=[0.1275, 0.11, 0.862, 0.945]

  ;; Colour bar
  colorbar, ncolors=n_elements(c_colors_phi), $
            minrange=alog10(min(levels_phi)), $
            maxrange=alog10(max(levels_phi)), $
            /vertical, position=[0.912, 0.11, 0.962, 0.945], $
            title=colorbartitle, AnnotateColor='black', Color=255, $
            divisions=levels_phi_range_log

  ;; Labelled contours
  set_basic_colours
  ;; Smooth contours
;;   xout = findgen(201) / 200 * (max(bincentres1) - min(bincentres1)) $
;;          + min(bincentres1)
;;   yout = findgen(201) / 200 * (max(bincentres2) - min(bincentres2)) $
;;          + min(bincentres2)
;;   contour, $
;;      min_curve_surf(phibbd, xvalues=bincentres1, yvalues=bincentres2, $
;;                     xout=xout, yout=yout), $
;;      xout, yout, /overplot, $
;;      levels=levels, $
;;      c_annotation=c_annotation, c_charsize=!P.charsize, $
;;      position=[0.1275, 0.11, 0.862, 0.945], c_charthick=2, $
;;      thick=3
  contour, $
     phibbd, $
     bincentres1, bincentres2, /overplot, $
     levels=levels, $
     c_annotation=c_annotation[*], $
     position=[0.1275, 0.11, 0.862, 0.945]
     thick=1.5 * (1 > !P.thick)

  ;; Choloniewski fit (only if parameters not supplied, and if not done above)
  IF (n_elements(choloniewski) LE 1 $
      AND (arg_present(choloniewski) $
           OR keyword_set(choloniewski))) THEN BEGIN
      ;; Calculate best-fitting Choloniewski function
      ;; NB phibbderr modified by this
      choloniewski = ajs_choloniewski_fit(bincentres1, bincentres2, $
                                          phibbd, phibbderr, $
                                          xrange=chol_xrange, $
                                          yrange=chol_yrange, vars=vars, $
                                          zmin=zmin, zmax=zmax, area=area, $
                                          q0=q0, q1=q1, cov_mat=cov_mat)
  ENDIF

  IF error_contours THEN BEGIN
      IF debug GE 1 THEN message, 'Doing error contours', /inf
      ;; Convert input errors to 1-sigma errors in chi^2, by scaling
      ;; so that chi^2 in each bin is equal, and taking into account
      ;; the number of degrees of freedom
      IF min(phibbderr) LE 0 THEN $
         message, 'Some bins have zero error: bit unrealistic?', /inf
      k = n_elements(phibbderr)
      chi2_bin = sqrt(2. / k)
      phibbderr_chi2 = sqrt(chi2_bin) * phibbderr
      contour, $ 
         phibbd + phibbderr_chi2, $
         bincentres1, bincentres2, /overplot, $
         levels=levels[1:*:2], c_linestyle=1, $
         position=[0.1275, 0.11, 0.862, 0.945], thick=1.5 * (1 > !P.thick)
      contour, $
         phibbd - phibbderr_chi2, $
         bincentres1, bincentres2, /overplot, $
         levels=levels[1:*:2], c_linestyle=2, $
         position=[0.1275, 0.11, 0.862, 0.945], thick=1.5 * (1 > !P.thick)

  ENDIF

  IF n_elements(choloniewski) GT 1 THEN BEGIN
      ;; Plot Choloniewski fit contours
      chol_m1 = ajs_linspace(min(xrange), max(xrange), 101)
      chol_m2 = ajs_linspace(min(yrange), max(yrange), 101)
      chol_bbd = ajs_choloniewski(chol_m1, chol_m2, $
                                  choloniewski)
      contour, $
         chol_bbd, $
         chol_m1, chol_m2, /overplot, $
         levels=lev_all, $
         c_annotation=c_annot_all, $
         position=[0.1275, 0.11, 0.862, 0.945], $
         thick=1.5 * (1 > !P.thick), c_linestyle=1, c_colors=1
      ;; Display Choloniewski parameters on the plot
;;       xyouts, xrange[0] + (xrange[1] - xrange[0]) / 25., $
;;               yrange[1] - (yrange[1] - yrange[0]) / 15., $
;;               string(choloniewski[0], format='(F0.2)') + ', ' $
;;               + string(choloniewski[1], format='(F0.2)') + ', ' $
;;               + string(choloniewski[2], format='(e0.1)') + ', ' $
;;               + string(choloniewski[3], format='(F0.2)') + ', ' $
;;               + string(choloniewski[4], format='(F0.2)') + ', ' $
;;               + string(choloniewski[5], format='(F0.2)')
;;       xyouts, xrange[0] + (xrange[1] - xrange[0]) / 25., $
;;               yrange[1] - 1.5 * (yrange[1] - yrange[0]) / 15., $
;;               string(choloniewski[3], format='(F0.2)') + ', ' $
;;               + string(choloniewski[4], format='(F0.2)') + ', ' $
;;               + string(choloniewski[5], format='(F0.2)')
  ENDIF  

  ajs_plot_stop, plotfile=plotfile, show_plot=show_plot, open=open, close=close

END

;+
; Test ajs_bbd_plot
;-
PRO ajs_bbd_plot_test, show_plot, _REF_EXTRA=e
  compile_opt idl2

  IF n_elements(show_plot) EQ 0 THEN $
     show_plot = 0

;;   default_dir = '/Users/anthonys/Documents/UKIDSS/DR3plus/'
;;   ajs_bbd_read, default_dir + '/multiLum/default/multiLum/bbd_vmax.dat', $
;;                 bincentres1, bincentres2, phibbd, phibbderr

  bincentres1 = ajs_linspace(-16, -26, 24)
  bincentres2 = ajs_linspace(21, 14, 24)
  phibbd = ajs_choloniewski(bincentres1, bincentres2, $
                            [-23, -1, 2e-2, 17, 0.5, 0.3])
  phibbderr = sqrt(phibbd)

  ;; Choloniewski fit
  ;; Cannot have zero errors: where err is zero, substitue 1e-4 (fudge)
  chol_phierr = phibbderr
  IF min(phibbderr) LE 0 THEN $
     chol_phierr[where(phibbderr EQ 0)] = 1e-4
  
  ;; Fit Choloniewski function
;;   chol_params = ajs_choloniewski_fit(bincentres1, bincentres2, $
;;                                      phibbd, chol_phierr, $
;;                                      txtfile=txtfile)

  ajs_bbd_plot, bincentres1, bincentres2, phibbd, phibbderr=chol_phierr, $
                show_plot=show_plot, $
                _STRICT_EXTRA=e, choloniewski=choloniewski

  print, choloniewski

  IF show_plot THEN $
     wait, 2 $
  ELSE $
     window, !d.window + 1
  
  ajs_bbd_plot, bincentres1, bincentres2, phibbd, $
                show_plot=show_plot, $
                _STRICT_EXTRA=e, choloniewski=choloniewski


END
