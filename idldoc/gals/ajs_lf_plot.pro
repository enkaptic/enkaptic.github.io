; docformat = 'rst'
;+
; This procedure produces a plot of the luminosity function.
;-

;+
; Choose and set xrange
;-
FUNCTION ajs_lf_plot_xrange, bincentres=bincentres, loglum=loglum, $
                             schechter=schechter, $
                             double_schechter=double_schechter
  compile_opt idl2, hidden

  ;; M-star?
  IF n_elements(schechter) GT 1 THEN $
     mstar = schechter[0] $
  ELSE IF n_elements(double_schechter) GT 1 THEN $
     mstar = double_schechter[0]

  ;; Set xrange limits
  IF n_elements(bincentres) GT 0 THEN BEGIN
      ;; ... from data points
      binwidth = (max(bincentres) - min(bincentres)) $
                 / (n_elements(bincentres) - 1)
      IF keyword_set(loglum) THEN $
         xrange = [min(bincentres) - binwidth / 2., $
                   max(bincentres) + binwidth / 2.] $
      ELSE $
         xrange = [max(bincentres) + binwidth / 2., $
                   min(bincentres) - binwidth / 2.]
  ENDIF ELSE IF n_elements(mstar) GT 0 THEN BEGIN
      ;; ... or from Schechter function
      IF keyword_set(loglum) THEN $
         xrange = [mstar - 6, $
                   mstar + 4] $
      ELSE $
         xrange = [mstar + 6, $
                   mstar - 4]
  ENDIF

  return, xrange
END


;+
; Choose and set yrange
;-
FUNCTION ajs_lf_plot_yrange, phi=phi, input_logphi=input_logphi, $
                             plot_logphi=plot_logphi, schechter=schechter, $
                             double_schechter=double_schechter
  compile_opt idl2, hidden

  ;; Phi-star?
  IF n_elements(schechter) GT 1 THEN $
     phistar = schechter[2] $
  ELSE IF n_elements(double_schechter) GT 1 THEN $
     phistar = double_schechter[2] > double_schechter[4]
  
  ;; Set yrange limits
  IF n_elements(phi) GT 0 THEN BEGIN
      ;; ... from data points
      IF keyword_set(input_logphi) THEN $
         yrange = 10. ^ [min(phi[where(finite(phi))]) - 10, $
                         max(phi[where(finite(phi))]) + 2] $
      ELSE $
         yrange = [min(phi[where(finite(phi))]) / 10, $
                   max(phi[where(finite(phi))]) * 2]
  ENDIF ELSE BEGIN
      ;; ... or from Schechter function
      IF keyword_set(input_logphi) THEN $ 
         yrange = 10. ^ [phistar - 4, phistar + 2] $
      ELSE $
         yrange = [phistar / 1e4, phistar * 1e2]
  ENDELSE

  ;; Y-axis log(phi) rather than phi?
  IF keyword_set(plot_logphi) THEN $
     yrange = alog10(yrange)

  return, yrange
END


;+
; Add legend text to the plot
;-
PRO ajs_lf_plot_legend, legend_text=legend_text, legend_pos=legend_pos, $
                        arg_legend_pos=arg_legend_pos, $
;;                         xrange=xrange, yrange=yrange, loglum=loglum, $
;;                         plot_logphi=plot_logphi, $
                        schechter=schechter, $
                        double_schechter=double_schechter, $
                        color=color, psym=psym, _REF_EXTRA=e
  compile_opt idl2, hidden
  debug = ajs_debug()

  ;; Create legend text?
  IF n_elements(legend_text) GT 0 THEN BEGIN 
      IF n_elements(legend_pos) LE 1 THEN BEGIN
          ;; Position for legend item
;;           IF n_elements(xrange) LE 1 OR n_elements(yrange) LE 1 THEN BEGIN
;;               message, 'Cannot position legend text: supply xrange and yrange'
;;           ENDIF ELSE BEGIN
;;               xpos = (xrange[1] - xrange[0]) / 10. + xrange[0]
;;               IF keyword_set(plot_logphi) THEN $
;;                  ypos = (yrange[1] - yrange[0]) / 10. + yrange[0] $
;;               ELSE $
;;                  ypos = (yrange[1] / yrange[0]) ^ 0.1 * yrange[0]
;;           ENDELSE
          xpos = 0.1
          ypos = 0.1
      ENDIF ELSE BEGIN
          xpos = legend_pos[0]
          ypos = legend_pos[1]
      ENDELSE 
          
      ;; Symbol/line (to left on plot = + if x-axis reversed)
;;       xpos_sym = xpos + 0.12 * (!x.crange[1] GT !x.crange[0] ? -1 : 1)
;;       xpos_line = [xpos + 0.21 * (!x.crange[1] GT !x.crange[0] ? -1 : 1), $
;;                    xpos + 0.03 * (!x.crange[1] GT !x.crange[0] ? -1 : 1)]
;;       IF keyword_set(plot_logphi) THEN $
;;          ypos_sym = ypos + 0.05 $
;;       ELSE $
;;          ypos_sym = ypos * 1.15
;;       IF n_elements(psym) GT 0 THEN $
;;          oplot, [xpos_sym], [ypos_sym], psym=psym, color=color
;;       IF n_elements(schechter) GT 1 OR n_elements(double_schechter) GT 1 $
;;          OR psym LE 0 THEN $
;;          oplot, xpos_line, [ypos_sym, ypos_sym], $
;;                 color=color, _STRICT_EXTRA=e

;;       xpos_sym = ajs_plot_pos(xpos - 0.05)
;;       xpos_line = ajs_plot_pos([xpos - 0.08, xpos - 0.02])
;;       ypos_sym = ajs_plot_pos(y_in=ypos + 0.01)
      IF n_elements(psym) GT 0 THEN $
         oplot, [ajs_plot_pos(xpos - 0.04)], $
                [ajs_plot_pos(y_in=ypos + 0.01)], psym=psym, color=color
      IF n_elements(schechter) GT 1 OR n_elements(double_schechter) GT 1 $
         OR psym LE 0 THEN $
            oplot, ajs_plot_pos([xpos - 0.06, xpos - 0.02]), $
                   ajs_plot_pos(y_in=[ypos + 0.01, ypos + 0.01]), $
                   color=color, _STRICT_EXTRA=e
      
      ;; Text
      IF debug GE 2 THEN $
         message, 'xpos, ypos = ' + string(xpos) + string(ypos), /inf
;;       xyouts, xpos, ypos, legend_text, color=color
      xyouts, ajs_plot_pos(xpos), ajs_plot_pos(y_in=ypos), legend_text, $
              color=color
      
      ;; Position for next legend item
      IF arg_legend_pos THEN BEGIN
;;          IF n_elements(yrange) GT 1 THEN BEGIN
;;               IF keyword_set(plot_logphi) THEN BEGIN 
;;                   ydelta = (yrange[1] - yrange[0]) / 20.
;;                   ypos = ypos + ydelta              
;;               ENDIF ELSE BEGIN
;;                   ydelta = (yrange[1] / yrange[0]) ^ (1. / 20)
;;                   ypos = ypos * ydelta
;;               ENDELSE
;;               legend_pos = [xpos, ypos]
          legend_pos = [xpos, ypos + 0.05]
;;           ENDIF ELSE BEGIN
;;               message, 'Supply yrange for new legend_pos values', /inf
;;               legend_pos = 0
;;           ENDELSE
      ENDIF ELSE BEGIN
          message, 'Supply variable name for legend_pos for new values', /inf
      ENDELSE          
  ENDIF
END


;+
; This procedure produces a plot of the luminosity function.
;
; Data may be entered in four forms:
;
; (1) Values at particular absolute magnitudes (bincentres, phi[, err_phi])
;
; (2) (double) Schechter function parameters (always plotted if supplied)
;
; (3) Array of galaxy absolute magnitudes and weights (e.g., 1/Vmax) -
; ignored if (1) set
;
; (4) From a file, to be read using ajs_lf_read - ignored if (1) or (3) set
;
; In addition, a (double) Schechter function fit can be performed.
;
; :Keywords:
;    bincentres : in, optional
;    phi : in, optional
;    err_phi : in, optional
;    neg_err_phi : in, optional
;       Negative error, for asymmetric
;    pos_err_phi : in, optional
;       Positive error, for asymmetric
;    lower_limit_where : in, optional
;       Indices of values of phi where the upper error bar is a lower
;       limit (i.e., plot arrowhead and empty circle, if psym=16 set)
;    schechter : in, optional
;       Plot Schechter function. Set to undefined variable to
;       calculate best-fitting Schechter function and return the
;       parameters. Or set /schechter to calculate and print the
;       parameters
;    sch_range : in, optional
;       Range of absolute magnitudes for (double) Schechter fit
;       (ajs_schechter_fit)
;    cov_mat : out, optional
;       Covariance matrix of Schechter function fit
;    double_schechter: in, out, optional
;       Plot double Schechter function. Set to undefined variable to
;       calculate best-fitting double Schechter function and return the
;       parameters. Or set /double_schechter to calculate and print the
;       parameters
;    double_cov_mat : out, optional
;       Covariance matrix of double Schechter function fit
;    abs_mag : in, optional
;       Array of absolute magnitudes (use with weights)
;    weights : in, optional
;       Array of weights (e.g., 1/Vmax) for estimating luminosity
;       function
;    jackknife : in, optional
;       Used with absmag, this is an integer for each galaxy
;       giving the jackknife region in which the galaxy lies
;    nbins : in, optional
;       Number of bins for abs_mag & weights (use with xrange or bincentres)
;    datfile : in, optional
;       Input data: binary format file created by ajs_lf_write
;    overplot : in, optional
;       Set /overplot to overplot
;    plotfile : in, optional
;       EPS file name
;    show_plot : in, optional
;       Set /show_plot to open the EPS file for viewing
;    open : in, optional
;       Set /open to open the EPS file and leave it for overplotting
;    close : in, optional
;       Set /close to close the EPS plotfile
;    xrange : in, out, optional
;    yrange : in, out, optional
;    input_logphi : in, optional
;       Set /input_logphi to indicate data values are log(phi) rather
;       than phi
;    plot_logphi : in, optional
;       Set /plot_logphi to plot log(phi) rather than phi (redundant
;       if overplotting)
;    loglum : in, optional
;       Set /loglum for log luminosities rather than magnitudes
;    ylog : in, optional, default=1-plot_logphi
;    psym : in, optional
;       Uses D. Fanning's symcat (e.g., 16 = filled circle) for
;       non-standard psym values
;    symcat_thick : in, optional
;       Thick keyword for symcat
;    band_name : in, optional
;       For the x-axis: M_<band_name> - 5 log h
;    color : in, optional
;    xstyle : in, optional
;       Default /xstyle
;    ystyle : in, optional
;       Default /ystyle
;    title : in, optional
;       Title for plot
;    xtitle : in, optional
;       Title for x-axis
;    ytitle : in, optional
;       Title for y-axis
;    legend_pos : in, out, optional, type=fltarr(2)
;       [xpos, ypos] for legend text. Output is position for next
;       legend label, for overplotting. varies from [0, 0] = [left,
;       bottom] to [1, 1] = [right, top].
;    legend_text : in, optional
;       Legend text
;    _REF_EXTRA : in, optional
;       Extra keywords for plots
;
; :Examples:
;    bincentres = [-20, -22, -24]
;
;    phi = [1e-2, 1e-2, 1e-3]
;
;    err_phi = [1e-2, 1e-3, 1e-4]
;
;    ajs_lf_plot, bincentres=bincentres, phi=phi
;
;    ajs_lf_plot, bincentres=bincentres, phi=phi, err_phi=err_phi
;
;    ajs_lf_plot, bincentres=bincentres, phi=phi, err_phi=err_phi, psym=16
;
;    ajs_lf_plot, schechter=[-24,-1,1e-2]
;
;    ajs_lf_plot, schechter=[-24,-1,1e-2], color=2, bincentres=bincentres, phi=phi, err_phi=err_phi, xrange=[-16,-26], psym=16, plotfile='~/lf.eps', /show_plot
;
; :History:
;    13 Feb 2008 Created, Anthony Smith
;
;    5 Mar 2008 err_phi no longer required
;
;    6 Mar 2008 Added abs_mag, weights and nbins
;
;    18 Mar 2008 Added jackknife
;
;    25 Mar 2008 Added legend_text and legend_pos
;
;    4 Apr 2008 Added sch_range (ajs_schechter_fit range)
;
;    7 Apr 2008 Schechter function info included in call to ajs_vmax_lf
;
;    21 Apr 2008 double Schechter function input
;
;    23 Apr 2008 psym now uses D. Fanning's symcat routine;
;    added lower_limit_where keyword
;
;    1 May 2008 Uses ajs_plot_pos for legend position
;-
PRO ajs_lf_plot, bincentres=bincentres, $
                 phi=phi, $
                 err_phi=err_phi, $
                 neg_err_phi=neg_err_phi, $
                 pos_err_phi=pos_err_phi, $
                 lower_limit_where=lower_limit_where, $
                 schechter=schechter, $
                 sch_range=sch_range, $
                 double_schechter=double_schechter, $
                 double_cov_mat=double_cov_mat, $
                 abs_mag=abs_mag, $
                 weights=weights, $
                 jackknife=jackknife, $
                 nbins=nbins, $
                 datfile=datfile, $
                 overplot=overplot, $
                 plotfile=plotfile, $
                 show_plot=show_plot, $
                 open=open, $
                 close=close, $
                 xrange=xrange, $
                 yrange=yrange, $
                 input_logphi=input_logphi, $
                 plot_logphi=plot_logphi, $
                 loglum=loglum, $
                 ylog=ylog, $
                 psym=psym, $
                 symcat_thick=symcat_thick, $
                 band_name=band_name, $
                 color=color, $
                 xstyle=xstyle, $
                 ystyle=ystyle, $
                 title=title, $
                 xtitle=xtitle, $
                 ytitle=ytitle, $
                 legend_pos=legend_pos, $
                 legend_text=legend_text, $
                 cov_mat=cov_mat, $
                 _REF_EXTRA=e

  compile_opt idl2
  debug = ajs_debug()
  IF debug GE 1 THEN BEGIN 
      message, 'Plotting LF', /inf
      message, ajs_kw_string(overplot=overplot, schechter=schechter, $
                             plotfile=plotfile, show_plot=show_plot, $
                             open=open, close=close, $
                             double_schechter=double_schechter), /inf
  ENDIF 

  ;; Keywords
  ;; Won't go wrong if symcat not available
  IF n_elements(psym) GT 0 THEN BEGIN
      IF abs(psym) EQ 9 OR abs(psym) GT 10 THEN $
         psym_use = symcat(abs(psym), thick=symcat_thick) $
      ELSE $
         psym_use = psym
      psym_use = abs(psym_use) * (psym GT 0 ? 1 : -1)
  ENDIF ELSE $
     psym_use = 0
     
  ;; Plot something?
  IF NOT keyword_set(close) THEN BEGIN
      ;; Any data to plot?
      IF n_elements(datfile) EQ 0 AND n_elements(phi) EQ 0 $
         AND NOT keyword_set(open) EQ 0 $
         AND n_elements(schechter) LE 1 $
         AND n_elements(double_schechter) LE 1 $ 
         AND n_elements(absmag) EQ 0 THEN $
            message, 'No data!'

      ;; Setup plot
      eps = keyword_set(open) OR keyword_set(show_plot)
      ajs_plot_start, plotfile=plotfile, eps=eps

      ;; Choose source of phi data: (1) phi, (2) absmag or (3) datfile
      IF n_elements(phi) EQ 0 THEN BEGIN
          IF n_elements(abs_mag) GT 0 THEN BEGIN
              ;; Data from absolute magnitudes & weights
              ;; Calculate luminosity function
              IF debug GE 1 THEN message, 'Calculating from absmags', /inf
              ;; Include Schechter fit?
              IF (n_elements(schechter) LE 1 $
                  AND (arg_present(schechter) $
                       OR keyword_set(schechter))) THEN BEGIN
                  IF (n_elements(double_schechter) LE 1 $
                      AND (arg_present(double_schechter) $
                           OR keyword_set(double_schechter))) THEN BEGIN
                      IF debug GE 2 THEN message, $
                         '1/Vmax with Schechter and double Schechter', /inf
                      phi = ajs_vmax_lf(abs_mag, weights, $
                                        bincentres=bincentres, $
                                        err_phi=err_phi, nbins=nbins, $
                                        xrange=xrange, $
                                        jackknife=jackknife, $
                                        schechter=schechter, $
                                        cov_mat=cov_mat, sch_range=sch_range, $
                                        double_schechter=double_schechter, $
                                        double_cov_mat=double_cov_mat, $
                                        loglum=loglum)
                  ENDIF ELSE BEGIN
                      IF debug GE 2 THEN message, $
                         '1/Vmax with Schechter', /inf
                      phi = ajs_vmax_lf(abs_mag, weights, $
                                        bincentres=bincentres, $
                                        err_phi=err_phi, nbins=nbins, $
                                        xrange=xrange, $
                                        jackknife=jackknife, $
                                        schechter=schechter, $
                                        cov_mat=cov_mat, sch_range=sch_range, $
                                        loglum=loglum)
                  ENDELSE
              ENDIF ELSE BEGIN
                  IF (n_elements(double_schechter) LE 1 $
                      AND (arg_present(double_schechter) $
                           OR keyword_set(double_schechter))) THEN BEGIN
                      IF debug GE 2 THEN message, $
                         '1/Vmax with double Schechter', /inf
                      phi = ajs_vmax_lf(abs_mag, weights, $
                                        bincentres=bincentres, $
                                        err_phi=err_phi, nbins=nbins, $
                                        xrange=xrange, $
                                        jackknife=jackknife, $
                                        sch_range=sch_range, $
                                        double_schechter=double_schechter, $
                                        double_cov_mat=double_cov_mat, $
                                        loglum=loglum)
                  ENDIF ELSE BEGIN
                      IF debug GE 2 THEN message, '1/Vmax (no Schechter)', /inf
                      phi = ajs_vmax_lf(abs_mag, weights, $
                                        bincentres=bincentres, $
                                        err_phi=err_phi, nbins=nbins, $
                                        xrange=xrange, $
                                        jackknife=jackknife, $
                                        loglum=loglum)
                  ENDELSE
              ENDELSE 
          ENDIF ELSE IF n_elements(datfile) GT 0 THEN BEGIN
              ;; Data from file
              IF debug GE 1 THEN message, 'Reading data from file', /inf
              ajs_lf_read, datfile, bincentres, phi, err_phi
          ENDIF
      ENDIF ELSE $
         IF debug GE 1 THEN message, 'Using supplied values of phi', /inf

      ;; Create axes (if not overplot)
      IF NOT keyword_set(overplot) THEN BEGIN
          ;; Titles not set
          IF n_elements(ytitle) EQ 0 THEN BEGIN
              IF keyword_set(loglum) THEN BEGIN
                  IF keyword_set(plot_logphi) THEN $
                     ytitle = 'log (!4U!3 / (h!U3!N Mpc!U-3!N dex!U-1!N) )' $
                  ELSE $
                     ytitle = '!4U!3 / (h!U3!N Mpc!U-3!N dex!U-1!N)'
              ENDIF ELSE BEGIN
                  IF keyword_set(plot_logphi) THEN $
                     ytitle = 'log (!4U!3 / (h!U3!N Mpc!U-3!N mag!U-1!N) )' $
                  ELSE $
                     ytitle = '!4U!3 / (h!U3!N Mpc!U-3!N mag!U-1!N)'
              ENDELSE
          ENDIF
          IF n_elements(xtitle) EQ 0 THEN BEGIN
              IF keyword_set(loglum) THEN BEGIN
                  IF n_elements(band_name) GT 0 THEN $
                     xtitle = 'log L!D' + band_name + '!N' $
                  ELSE $
                     xtitle = 'log L'
              ENDIF ELSE BEGIN
                  IF n_elements(band_name) GT 0 THEN $
                     xtitle = 'M!D' + band_name + '!N-5log h' $
                  ELSE $
                     xtitle = 'M-5log h'
              ENDELSE
          ENDIF

          ;; No xrange/yrange set
          IF n_elements(xrange) LE 1 THEN $
             xrange = ajs_lf_plot_xrange(bincentres=bincentres, $
                                         loglum=loglum, schechter=schechter, $
                                         double_schechter=double_schechter)
          IF n_elements(yrange) LE 1 THEN $
             yrange = ajs_lf_plot_yrange(phi=phi, input_logphi=input_logphi, $
                                         plot_logphi=plot_logphi, $
                                         schechter=schechter, $
                                         double_schechter=double_schechter)

          ;; ylog if not plot_logphi
          IF n_elements(ylog) EQ 0 AND NOT keyword_set(plot_logphi) THEN $
             ylog = 1

          ;; Set /xstyle & /ystyle by default
          IF n_elements(xstyle) EQ 0 THEN $
             xstyle=1
          IF n_elements(ystyle) EQ 0 THEN $
             ystyle=1

          ;; Plot empty axes
          plot,[0],[1], xrange=xrange, yrange=yrange, /nodata, $
               title=title, ytitle=ytitle, xtitle=xtitle, ylog=ylog, $
               xstyle=xstyle, ystyle=ystyle, _STRICT_EXTRA=e
      ENDIF
      plot_logphi = 1 - !y.type

      ;; Define filled circle (IDL Help): psym=8
      ;; Obselete: now uses D. Fanning's symcat routine
;;       A = FINDGEN(17) * (!PI * 2 / 16.)  
;;       USERSYM, COS(A), SIN(A), /FILL  

      ;; Plot data points?
      IF n_elements(phi) GT 0 THEN BEGIN 
          ;; Error bars
          IF n_elements(err_phi) GT 0 THEN BEGIN
              neg_err_phi = err_phi
              pos_err_phi = err_phi
          ENDIF
          IF n_elements(neg_err_phi) GT 0 THEN BEGIN 
              phi_err_max = phi + pos_err_phi
              phi_err_min = phi - neg_err_phi
              IF NOT keyword_set(input_logphi) THEN BEGIN 
                  ;; Infinite plot range?
                  too_small = where(phi_err_min LE 0)
                  IF too_small[0] NE -1 THEN $
                     phi_err_min[too_small] = 1e-30 ; log(phi) > 0
              ENDIF
              IF keyword_set(input_logphi) $
                 AND NOT keyword_set(plot_logphi) THEN BEGIN
                  phi_err_max = 10. ^ phi_err_max
                  phi_err_min = 10. ^ phi_err_min                 
              ENDIF ELSE IF NOT keyword_set(input_logphi) $
                 AND keyword_set(plot_logphi) THEN BEGIN
                  phi_err_max = alog10(phi_err_max)
                  phi_err_min = alog10(phi_err_min)
              ENDIF
          ENDIF

          ;; log(phi) ? (NB don't change phi itself!)
          IF keyword_set(input_logphi) $
             AND NOT keyword_set(plot_logphi) THEN $
                phi_plot = 10. ^ phi $
          ELSE IF NOT keyword_set(input_logphi) $
             AND keyword_set(plot_logphi) THEN $
                phi_plot = alog10(phi) $
          ELSE $
             phi_plot = phi

          ;; Plot the data
          IF n_elements(neg_err_phi) GT 0 THEN BEGIN
              ;; With error bars
              IF n_elements(lower_limit_where) GT 0 THEN BEGIN
                  ;; ... and mark out bins which are lower limits
                  all_indices = indgen(n_elements(phi_plot))
                  not_lower_limit_bool = intarr(n_elements(phi_plot)) + 1
                  not_lower_limit_bool[lower_limit_where] = 0
                  not_lower_limit_where = $
                     all_indices[where(not_lower_limit_bool)]
                  oploterror, bincentres[not_lower_limit_where], $
                              phi_plot[not_lower_limit_where], $
                              phi_plot[not_lower_limit_where] $
                              - phi_err_min[not_lower_limit_where], $
                              /lobar, psym=psym_use, $
                              color=color, errcolor=color, _STRICT_EXTRA=e
                  oploterror, bincentres[not_lower_limit_where], $
                              phi_plot[not_lower_limit_where], $
                              phi_err_max[not_lower_limit_where] $
                              - phi_plot[not_lower_limit_where], $
                              /hibar, psym=psym_use, $
                              color=color, errcolor=color, _STRICT_EXTRA=e
                  ;; Open circles?
                  IF n_elements(psym) GT 0 THEN $
                     IF abs(psym) EQ 16 THEN $
                        psym_use = symcat(9 * abs(psym) / 16, $
                                          thick=symcat_thick)
                  oploterror, bincentres[lower_limit_where], $
                              phi_plot[lower_limit_where], $
                              phi_plot[lower_limit_where] $
                              - phi_err_min[lower_limit_where], $
                              /lobar, psym=psym_use, $
                              color=color, errcolor=color, _STRICT_EXTRA=e
                  oploterror, bincentres[lower_limit_where], $
                              phi_plot[lower_limit_where], $
                              phi_err_max[lower_limit_where] $
                              - phi_plot[lower_limit_where], $
                              /hibar, psym=psym_use, /nohat, $
                              color=color, errcolor=color, _STRICT_EXTRA=e
                  ;; Define arrowhead
                  usersym, [-1, 0, 1], [-1, 0, -1]
                  ;; Plot arrows for lower limits
                  oplot, bincentres[lower_limit_where], $
                         phi_err_max[lower_limit_where], psym=8, color=color
                  ;; Restore psym
                  IF psym_use EQ 8 THEN $
                     IF psym EQ 9 OR psym GT 10 THEN $
                        psym_use = symcat(psym, thick=symcat_thick)
              ENDIF ELSE BEGIN
                  ;; Nothing special for lower-limit bins
                  oploterror, bincentres, phi_plot, $
                              phi_plot - phi_err_min, $
                              /lobar, psym=psym_use, $
                              color=color, errcolor=color, _STRICT_EXTRA=e
                  oploterror, bincentres, phi_plot, $
                              phi_err_max - phi_plot, $
                              /hibar, psym=psym_use, $
                              color=color, errcolor=color, _STRICT_EXTRA=e
              ENDELSE
          ENDIF ELSE BEGIN
              ;; Without error bars
              oplot, bincentres, phi_plot, $
                     psym=psym_use, $
                     color=color, _STRICT_EXTRA=e
          ENDELSE
      ENDIF

      ;; Calcluate Schechter function (if not done with LF calculation above)?
      IF (n_elements(schechter) LE 1 $
          AND (arg_present(schechter) $
               OR keyword_set(schechter))) THEN BEGIN
          ;; Calculate best-fitting Schechter function
          schechter = ajs_schechter_fit(bincentres, phi, err_phi, $
                                        cov_mat=cov_mat, range=sch_range)
      ENDIF

      ;; Calcluate double Schechter function (if not done above)?
      IF (n_elements(double_schechter) LE 1 $
          AND (arg_present(double_schechter) $
               OR keyword_set(double_schechter))) THEN BEGIN
          ;; Calculate best-fitting double Schechter function
          double_schechter = $
             ajs_double_schechter_fit(bincentres, phi, err_phi, $
                                      cov_mat=double_cov_mat, range=sch_range)
      ENDIF

      ;; Plot Schechter function?
      IF n_elements(schechter) GT 0 THEN BEGIN
          IF debug GE 1 THEN message, 'Plotting Schechter function', /inf
          mstar = schechter[0]
          alpha = schechter[1]
          IF keyword_set(input_logphi) THEN $ 
             phistar = 10. ^ schechter[2] $
          ELSE $
             phistar = schechter[2]

          IF n_elements(xrange) GT 0 THEN $ 
             m = ajs_linspace(min(xrange), max(xrange), 101) $
          ELSE $
             m = ajs_linspace(mstar + 20, mstar - 10, 501)
          phi_sch = ajs_schechter(m, [mstar, $
                                      alpha, $
                                      phistar], loglum=loglum)

          IF keyword_set(plot_logphi) THEN $
             phi_sch = alog10(phi_sch)
          
          oplot, m, phi_sch, color=color, _STRICT_EXTRA=e
      ENDIF

      ;; Plot double Schechter function?
      IF n_elements(double_schechter) GT 0 THEN BEGIN
          IF debug GE 1 THEN message, 'Plotting double Schechter fn', /inf
          mstar = double_schechter[0]
          alpha_1 = double_schechter[1]
          alpha_2 = double_schechter[3]
          IF keyword_set(input_logphi) THEN BEGIN
              phistar_1 = 10. ^ double_schechter[2]
              phistar_2 = 10. ^ double_schechter[4]
          ENDIF ELSE BEGIN
              phistar_1 = double_schechter[2]
              phistar_2 = double_schechter[4]
          ENDELSE
              
          IF n_elements(xrange) GT 0 THEN $ 
             m = ajs_linspace(min(xrange), max(xrange), 101) $
          ELSE $
             m = ajs_linspace(mstar + 20, mstar - 10, 501)
          phi_double_sch = ajs_double_schechter( $
                           m, [mstar, $
                               alpha_1, phistar_1, alpha_2, $
                               phistar_2], loglum=loglum)
          
          IF keyword_set(plot_logphi) THEN $
             phi_double_sch = alog10(phi_double_sch)
          
          oplot, m, phi_double_sch, color=color, _STRICT_EXTRA=e
      ENDIF

      ;; Add legend?
      ajs_lf_plot_legend, legend_text=legend_text, legend_pos=legend_pos, $
                          arg_legend_pos=arg_present(legend_pos), $
;;                           xrange=xrange, yrange=yrange, loglum=loglum, $
;;                           plot_logphi=plot_logphi, $
                          schechter=schechter, $
                          double_schechter=double_schechter, $
                          color=color, psym=psym_use, $
                          _EXTRA=e
  ENDIF

  ajs_plot_stop, plotfile=plotfile, show_plot=show_plot, open=open, close=close
      
  RETURN
END
