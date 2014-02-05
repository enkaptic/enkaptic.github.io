; docformat = 'rst'
;+
; Plot ra and dec on equal-area plot with reversed x-axis
;
; Wrapper on plot/oplot
; :Params:
;    x : in, required
;       ra, degrees
;    y : in, required
;       dec, degrees
; :Keywords:
;    overplot : in, optional
;       Set /overplot to overplot
;    xrange : in, optional
;       As for plot (reversed if not specified)
;    yrange : in, optional
;       As for plot
;    xstyle : in, optional
;       As for plot
;    ystyle : in, optional
;       As for plot
;    _REF_EXTRA : in, optional
;       Extra keywords for [o]plot
; :History:
;    17 Apr 2008 Written, Anthony Smith
;-
PRO ajs_plot_radec, x, y, overplot=overplot, xrange=xrange, yrange=yrange, $
                    xstyle=xstyle, ystyle=ystyle, _REF_EXTRA=e
  compile_opt idl2

  IF n_elements(y) EQ 0 THEN $
     message, 'Input values for x and y'
  IF max(y) GT 90 OR min(y) LT -90 THEN $
     message, 'Values for y must be between -90 and 90'

  IF ~ keyword_set(overplot) THEN BEGIN 
      ;; Axis ranges
      IF n_elements(xrange) GT 0 THEN $
         xrange_plot = xrange $
      ELSE $
         xrange_plot = [max(x), min(x)] ; Reversed by default
      IF n_elements(yrange) GT 0 THEN $
         yrange_plot = yrange $
      ELSE $
         yrange_plot = [min(y), max(y)]

      ;; Exact axis range?
      IF n_elements(xstyle) GT 0 THEN $
         xstyle_plot = xstyle $
      ELSE $
         xstyle_plot = 0
      IF n_elements(ystyle) GT 0 THEN $
         ystyle_plot = ystyle $
      ELSE $
         ystyle_plot = 0

      ;; Use default range if ticks would be in 10s or less
      xspan = abs(xrange_plot[1] - xrange_plot[0])
      yspan = abs(yrange_plot[1] - yrange_plot[0])
      xdefault = (xstyle_plot AND xspan LE 70) $
                 OR (ceil(max(xrange_plot) / 10.) * 10 $
                     - floor(min(xrange_plot) / 10.) * 10 LE 60)
      ydefault = (ystyle_plot AND yspan LE 70) $
                 OR (ceil(max(yrange_plot) / 10.) * 10 $
                     - floor(min(yrange_plot) / 10.) * 10 LE 60)
      IF xdefault OR ydefault THEN BEGIN
          plot, x, y, _STRICT_EXTRA=e, xtick_get=xtickv, ytick_get=ytickv, $
                xrange=xrange_plot, yrange=yrange_plot, $
                xstyle=xstyle_plot, ystyle=ystyle_plot
          IF xdefault THEN $
             xrange_plot = !x.crange
          IF ydefault THEN $
             yrange_plot = !y.crange
      ENDIF

      ;; x : Non-default ranges (large, in 15s or 30s)
      IF ~ xdefault THEN BEGIN
          ;; Try in 15s
          IF NOT xstyle_plot THEN $
             xrange_tmp = [floor(min(xrange_plot) / 15.) * 15, $
                           ceil(max(xrange_plot) / 15.) * 15] $
          ELSE $
             xrange_tmp = [min(xrange_plot), max(xrange_plot)]
          n_xticks = floor(xrange_tmp[1] / 15.) $
                     - ceil(xrange_tmp[0] / 15.) + 1
          xtickv = (indgen(n_xticks) + ceil(xrange_tmp[0] / 15.)) * 15

          ;; Too many ticks?  In 30s
          IF n_elements(xtickv) GT 7 THEN BEGIN
              IF NOT xstyle_plot THEN $
                 xrange_tmp = [floor(min(xrange_plot) / 30.) * 30, $
                               ceil(max(xrange_plot) / 30.) * 30]
              n_xticks = floor(xrange_tmp[1] / 30.) $
                         - ceil(xrange_tmp[0] / 30.) + 1
              xtickv = (indgen(n_xticks) + ceil(xrange_tmp[0] / 30.)) * 30
          ENDIF

          ;; Too many ticks?  In 45s
          IF n_elements(xtickv) GT 7 THEN BEGIN
              IF NOT xstyle_plot THEN $
                 xrange_tmp = [floor(min(xrange_plot) / 45.) * 45, $
                               ceil(max(xrange_plot) / 45.) * 45]
              n_xticks = floor(xrange_tmp[1] / 45.) $
                         - ceil(xrange_tmp[0] / 45.) + 1
              xtickv = (indgen(n_xticks) + ceil(xrange_tmp[0] / 45.)) * 45
          ENDIF

          ;; Too many ticks?  In 60s
          IF n_elements(xtickv) GT 7 THEN BEGIN
              IF NOT xstyle_plot THEN $
                 xrange_tmp = [floor(min(xrange_plot) / 60.) * 60, $
                               ceil(max(xrange_plot) / 60.) * 60]
              n_xticks = floor(xrange_tmp[1] / 60.) $
                         - ceil(xrange_tmp[0] / 60.) + 1
              xtickv = (indgen(n_xticks) + ceil(xrange_tmp[0] / 60.)) * 60
          ENDIF

          IF xrange_plot[1] LT xrange_plot[0] THEN $
             xrange_plot = reverse(xrange_tmp) $
          ELSE $
             xrange_plot = xrange_tmp
      ENDIF

      ;; y : Non-default ranges (large, in 15s or 30s)
      IF ~ ydefault THEN BEGIN
          ;; In 15s
          IF NOT ystyle_plot THEN $
             yrange_tmp = [floor(min(yrange_plot) / 15.) * 15, $
                           ceil(max(yrange_plot) / 15.) * 15] $
          ELSE $
             yrange_tmp = [min(yrange_plot), max(yrange_plot)]

          n_yticks = floor(max(yrange_tmp) / 15.) $
                     - ceil(min(yrange_tmp) / 15.) + 1
          ytickv = (indgen(n_yticks) + ceil(min(yrange_tmp) / 15.)) * 15

          IF yrange_plot[1] LT yrange_plot[0] THEN $
             yrange_plot = reverse(yrange_tmp) $
          ELSE $
             yrange_plot = yrange_tmp          
      ENDIF

      ;; Strict ranges for the plots
      IF NOT xstyle_plot THEN $
         xstyle_plot += 1
      IF NOT ystyle_plot THEN $
         ystyle_plot += 1

      ;; Minor ticks for x-axis (not for y-axis: need to be irregular)
      xminor = n_elements(xtickv) GT 6 ? 5 : 10
      yminor = 1

      ;; Labels for y-axis: don't label the 75s if large range
      IF ytickv[1] - ytickv[0] GE 1 THEN $
         ytickname = string(ytickv, format='(I0)') $
      ELSE BEGIN
          sf = 0
          REPEAT BEGIN
              sf += 1
              format = '(F0.' + string(sf, format='(I0)') + ')'
              ytickname = string(ytickv, format=format)
          ENDREP UNTIL total(ytickname[1:*] $
                             EQ ytickname[0:n_elements(ytickname) - 2]) EQ 0
      ENDELSE
      IF max(yrange_plot) EQ 90 AND min(yrange_plot) LT 0 THEN $
         ytickname[where(ytickv EQ 75)] = ' '
      IF min(yrange_plot) EQ - 90 AND max(yrange_plot) GT 0 THEN $
         ytickname[where(ytickv EQ - 75)] = ' '

      plot, x, sin(y / !RADEG), _STRICT_EXTRA=e, $
            xrange=xrange_plot, yrange=sin(yrange_plot / !RADEG), $
            xticks=n_elements(xtickv) - 1, xtickv=xtickv, xminor=xminor, $
            xstyle=xstyle_plot, $
            yticks=n_elements(ytickv) - 1, ytickv=sin(ytickv / !RADEG), $
            yminor=yminor, ytickname=ytickname, $
            ystyle=ystyle_plot
  ENDIF ELSE $
     oplot, x, sin(y / !RADEG), _STRICT_EXTRA=e

END
