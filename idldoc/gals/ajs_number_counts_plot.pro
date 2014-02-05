; docformat = 'rst'
;+
; This procedure produces a plot of the number counts.
;-


;+
; Choose and set xrange
;-
FUNCTION ajs_number_counts_plot_xrange, $
   bincentres=bincentres, logflux=logflux
  compile_opt idl2, hidden

  ;; Set xrange limits
  IF n_elements(bincentres) GT 0 THEN BEGIN
      ;; ... from data points
      binwidth = (max(bincentres) - min(bincentres)) $
                 / (n_elements(bincentres) - 1)
      xrange = [min(bincentres) - binwidth / 2., $
                max(bincentres) + binwidth / 2.]
  ENDIF ELSE $
     message, 'Cannot determine xrange'

  return, xrange
END


;+
; Choose and set yrange
;-
FUNCTION ajs_number_counts_plot_yrange, $
   ngals=ngals, input_logngals=input_logngals, $
   plot_logngals=plot_logngals, err_ngals=err_ngals, $
   pos_err_ngals=pos_err_ngals, neg_err_ngals=neg_err_ngals
  compile_opt idl2, hidden

  ;; Set yrange limits
  ;; Todo: include error bars
  IF n_elements(ngals) GT 0 THEN BEGIN
      ;; ... from data points
      IF keyword_set(input_logngals) THEN $
         yrange = 10. ^ [min(ngals[where(finite(ngals))]) - 10, $
                         max(ngals[where(finite(ngals))]) + 2] $
      ELSE $
         yrange = [min(ngals[where(finite(ngals))]) / 10, $
                   max(ngals[where(finite(ngals))]) * 2]
  ENDIF ELSE $
     message, 'Cannot determine yrange'

  ;; Y-axis log(ngals) rather than ngals?
  IF keyword_set(plot_logngals) THEN $
     yrange = alog10(yrange)

  return, yrange
END


;+
; Add legend text to the plot
;-
PRO ajs_number_counts_plot_legend, $
   legend_text=legend_text, legend_pos=legend_pos, $
   arg_legend_pos=arg_legend_pos, $
   color=color, psym=psym, _REF_EXTRA=e
  compile_opt idl2, hidden
  debug = ajs_debug()

  ;; Create legend text?
  IF n_elements(legend_text) GT 0 THEN BEGIN 
      IF n_elements(legend_pos) LE 1 THEN BEGIN
          xpos = 0.1
          ypos = 0.9
      ENDIF ELSE BEGIN
          xpos = legend_pos[0]
          ypos = legend_pos[1]
      ENDELSE 

      ;; Symbol/line (to left on plot = + if x-axis reversed)
      IF n_elements(psym) GT 0 THEN $
         oplot, [ajs_plot_pos(xpos - 0.04)], $
                [ajs_plot_pos(y_in=ypos + 0.01)], psym=psym, color=color
      IF psym LE 0 THEN $
            oplot, ajs_plot_pos([xpos - 0.06, xpos - 0.02]), $
                   ajs_plot_pos(y_in=[ypos + 0.01, ypos + 0.01]), $
                   color=color, _STRICT_EXTRA=e
      
      ;; Text
      IF debug GE 2 THEN $
         message, 'xpos, ypos = ' + string(xpos) + string(ypos), /inf
      xyouts, ajs_plot_pos(xpos), ajs_plot_pos(y_in=ypos), legend_text, $
              color=color
      
      ;; Position for next legend item
      IF arg_legend_pos THEN BEGIN
          legend_pos = [xpos, ypos - 0.05]
      ENDIF ELSE BEGIN
          message, 'Supply variable name for legend_pos for new values', /inf
      ENDELSE          
  ENDIF
END


;+
; This procedure produces a plot of the number counts.
;
; :Keywords:
;    bincentres : in, optional
;       x-position (e.g., magnitude, flux, log-flux)
;    ngals : in, optional
;       Number density, per square degree, per unit (x_unit)
;    err_ngals : in, optional
;       Uncertainty in ngals
;    neg_err_ngals : in, optional
;       Negative error in ngals, for asymmetric
;    pos_err_ngals : in, optional
;       Positive error in ngals, for asymmetric
;    mag : in, optional
;       Array of magnitudes (or log(flux) with /logflux)
;    area : in, optional
;       Area in square degrees (for use with mag)
;    nbins : in, optional
;       Number of bins for mag (use with xrange or bincentres)
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
;    input_logngals : in, optional
;       Set /input_logngals to indicate data values are log(n) rather
;       than n
;    plot_logngals : in, optional
;       Set /plot_logngals to plot log(n) rather than n (redundant
;       if overplotting)
;    logflux : in, optional
;       Set /logflux for log flux rather than magnitudes
;    ylog : in, optional, default=1-plot_logn
;    psym : in, optional
;       Uses D. Fanning's symcat (e.g., 16 = filled circle) for
;       non-standard psym values
;    symcat_thick : in, optional
;       Thick keyword for symcat
;    color : in, optional
;    xstyle : in, optional
;       Default xstyle=1
;    ystyle : in, optional
;       Default ystyle=1
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
;    euclid : in, optional
;       Set /euclid to subtract Euclidean slope, or euclid=x to then
;       divide by x (specify axis ranges and titles explicitly)
;    _REF_EXTRA : in, optional
;       Extra keywords for plots
;
; :Examples:
;    bincentres = [10, 11, 12]
;
;    ngals = [10, 100, 1000]
;
;    err_ngals = [1, 10, 100]
;
;    ajs_number_counts_plot, bincentres=bincentres, ngals=ngals,
;    err_ngals=err_ngals
;
; :History:
;    29 Oct 2008 Created, Anthony Smith
;
;    16 Jan 2009 Added euclid keyword, AJS
;-
PRO ajs_number_counts_plot, $
   bincentres=bincentres, $
   ngals=ngals, $
   err_ngals=err_ngals, $
   neg_err_ngals=neg_err_ngals, $
   pos_err_ngals=pos_err_ngals, $
   mag=mag, $
   area=area, $
   nbins=nbins, $
   overplot=overplot, $
   plotfile=plotfile, $
   show_plot=show_plot, $
   open=open, $
   close=close, $
   xrange=xrange, $
   yrange=yrange, $
   input_logngals=input_logngals, $
   plot_logngals=plot_logngals, $
   logflux=logflux, $
   ylog=ylog, $
   psym=psym, $
   symcat_thick=symcat_thick, $
   color=color, $
   xstyle=xstyle, $
   ystyle=ystyle, $
   title=title, $
   xtitle=xtitle, $
   ytitle=ytitle, $
   legend_pos=legend_pos, $
   legend_text=legend_text, $
   euclid=euclid, $
   _REF_EXTRA=e

  compile_opt idl2
  debug = ajs_debug()
  IF debug GE 1 THEN BEGIN 
      message, 'Plotting number counts', /inf
      message, ajs_kw_string(overplot=overplot, $
                             plotfile=plotfile, show_plot=show_plot, $
                             open=open, close=close), /inf
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
      IF n_elements(ngals) EQ 0 AND n_elements(mag) EQ 0 $ 
         AND NOT keyword_set(open) EQ 0 THEN $
            message, 'No data!'

      ;; Setup plot
      eps = keyword_set(open) OR keyword_set(show_plot)
      ajs_plot_start, plotfile=plotfile, eps=eps


      ;; Choose source of phi data: (1) ngals or (2) mag
      IF n_elements(ngals) EQ 0 THEN BEGIN
          IF n_elements(mag) GT 0 THEN BEGIN
              ;; Data from magnitudes (or log fluxes)
              ;; Calculate number counts
              IF debug GE 1 THEN message, 'Calculating from mags', /inf
              ngals = ajs_number_counts(mag, area=area, $
                                        bincentres=bincentres, $
                                        err_ngals=err_ngals, nbins=nbins, $
                                        xrange=xrange)
              input_logngals = 0
          ENDIF
      ENDIF ELSE $
         IF debug GE 1 THEN message, 'Using supplied values of phi', /inf

      ;; Create axes (if not overplot)
      IF NOT keyword_set(overplot) THEN BEGIN
          ;; Titles not set
          IF n_elements(ytitle) EQ 0 THEN BEGIN
              IF keyword_set(logflux) THEN BEGIN
                  IF keyword_set(plot_logngals) THEN $
                     ytitle = 'log (N / (dex!U-1!N deg!U-2!N) )' $
                  ELSE $
                     ytitle = 'N / (dex!U-1!N deg!U-2!N)'
              ENDIF ELSE BEGIN
                  IF keyword_set(plot_logngals) THEN $
                     ytitle = 'log (N / (mag!U-1!N deg!U-2!N) )' $
                  ELSE $
                     ytitle = 'N / (mag!U-1!N deg!U-2!N)'
              ENDELSE
          ENDIF
          IF n_elements(xtitle) EQ 0 THEN BEGIN
              IF keyword_set(logflux) THEN BEGIN
                  xtitle = 'log flux'
              ENDIF ELSE BEGIN
                  xtitle = 'magnitude'
              ENDELSE
          ENDIF

          ;; No xrange/yrange set
          IF n_elements(xrange) LE 1 THEN $
             xrange = ajs_number_counts_plot_xrange(bincentres=bincentres, $
                                                    logflux=logflux)
          IF n_elements(yrange) LE 1 THEN $
             yrange = ajs_number_counts_plot_yrange( $
                      ngals=ngals, input_logngals=input_logngals, $
                      plot_logngals=plot_logngals, err_ngals=err_ngals, $
                      pos_err_ngals=pos_err_ngals, neg_err_ngals=neg_err_ngals)

          ;; ylog if not plot_logn
          IF n_elements(ylog) EQ 0 AND NOT keyword_set(plot_logngals) THEN $
             ylog = 1

          ;; Set xstyle & ystyle = 1 by default
          IF n_elements(xstyle) EQ 0 THEN $
             xstyle=1
          IF n_elements(ystyle) EQ 0 THEN $
             ystyle=1

          ;; Plot empty axes
          plot,[0],[1], xrange=xrange, yrange=yrange, /nodata, $
               title=title, ytitle=ytitle, xtitle=xtitle, ylog=ylog, $
               xstyle=xstyle, ystyle=ystyle, _STRICT_EXTRA=e
      ENDIF
      plot_logngals = 1 - !y.type

      ;; Define filled circle (IDL Help): psym=8
      ;; Obselete: now uses D. Fanning's symcat routine
;;       A = FINDGEN(17) * (!PI * 2 / 16.)  
;;       USERSYM, COS(A), SIN(A), /FILL  

      ;; Plot data points?
      IF n_elements(ngals) GT 0 THEN BEGIN 
          ;; Error bars
          IF n_elements(err_ngals) GT 0 THEN BEGIN
              neg_err_ngals = err_ngals
              pos_err_ngals = err_ngals
          ENDIF
          IF n_elements(neg_err_ngals) GT 0 THEN BEGIN 
              ngals_err_max = ngals + pos_err_ngals
              ngals_err_min = ngals - neg_err_ngals
              IF NOT keyword_set(input_logngals) THEN BEGIN 
                  ;; Infinite plot range?
                  too_small = where(ngals_err_min LE 0)
                  IF too_small[0] NE -1 THEN $
                     ngals_err_min[too_small] = 1e-30 ; log(ngals) > 0
              ENDIF
              IF keyword_set(input_logngals) $
                 AND NOT keyword_set(plot_logngals) THEN BEGIN
                  ngals_err_max = 10. ^ ngals_err_max
                  ngals_err_min = 10. ^ ngals_err_min                 
              ENDIF ELSE IF NOT keyword_set(input_logngals) $
                 AND keyword_set(plot_logngals) THEN BEGIN
                  ngals_err_max = alog10(ngals_err_max)
                  ngals_err_min = alog10(ngals_err_min)
              ENDIF
          ENDIF

          ;; log(ngals) ? (NB don't change ngals itself!)
          IF keyword_set(input_logngals) $
             AND NOT keyword_set(plot_logngals) THEN $
                ngals_plot = 10. ^ ngals $
          ELSE IF NOT keyword_set(input_logngals) $
             AND keyword_set(plot_logngals) THEN $
                ngals_plot = alog10(ngals) $
          ELSE $
             ngals_plot = ngals

          ;; Plot the data
          IF n_elements(neg_err_ngals) GT 0 THEN BEGIN
              IF n_elements(euclid) GT 0 THEN BEGIN
                  ;; With error bars
                  oploterror, bincentres, $
                              ngals_plot / 10 ^ (0.6 * bincentres) / euclid, $
                              (ngals_plot - ngals_err_min) $
                               / 10 ^ (0.6 * bincentres) / euclid, $
                              /lobar, psym=psym_use, $
                              color=color, errcolor=color, _STRICT_EXTRA=e
                  oploterror, bincentres, $
                              ngals_plot / 10 ^ (0.6 * bincentres) / euclid, $
                              (ngals_err_max - ngals_plot) $
                               / 10 ^ (0.6 * bincentres) / euclid, $
                              /hibar, psym=psym_use, $
                              color=color, errcolor=color, _STRICT_EXTRA=e
                  
              ENDIF ELSE BEGIN
                  ;; With error bars
                  oploterror, bincentres, ngals_plot, $
                              ngals_plot - ngals_err_min, $
                              /lobar, psym=psym_use, $
                              color=color, errcolor=color, _STRICT_EXTRA=e
                  oploterror, bincentres, ngals_plot, $
                              ngals_err_max - ngals_plot, $
                              /hibar, psym=psym_use, $
                              color=color, errcolor=color, _STRICT_EXTRA=e
              ENDELSE
          ENDIF ELSE BEGIN
              IF n_elements(euclid) GT 0 THEN BEGIN
                  ;; Without error bars
                  oplot, bincentres, $
                         ngals_plot / 10 ^ (0.6 * bincentres) / euclid, $
                         psym=psym_use, $
                         color=color, _STRICT_EXTRA=e                  
              ENDIF ELSE BEGIN
                  ;; Without error bars
                  oplot, bincentres, ngals_plot, $
                         psym=psym_use, $
                         color=color, _STRICT_EXTRA=e
              ENDELSE
          ENDELSE
      ENDIF

      ;; Add legend?
      ajs_number_counts_plot_legend, $
         legend_text=legend_text, legend_pos=legend_pos, $
         arg_legend_pos=arg_present(legend_pos), $
;;                           xrange=xrange, yrange=yrange, logflux=logflux, $
;;                           plot_logn=plot_logn, $
         color=color, psym=psym_use, $
         _EXTRA=e
  ENDIF

  ajs_plot_stop, plotfile=plotfile, show_plot=show_plot, open=open, close=close
      
  RETURN
END
