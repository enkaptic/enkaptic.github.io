; docformat = 'rst'
;+
; Return the location on the plot in data co-ordinates given input
; value(s) between 0 and 1.
;
; Like using xyouts or plots with /normal but between the axis limits
; rather than the whole extent of the plot area
; :Returns:
;    Position in data units corresponding to input x, y, [x, y], [x],
;    [y] or [[x], [y]]
; :Params:
;    x : in, optional
;       x-position, between 0 (left-hand side) and 1 (right-hand
;       side). Or [N, 2] values for [x, y]
;    y : in, optional
;       y-position, between 0 (bottom) and 1 (top). Ignored if x is
;       [N, 2]
; :Keywords:
;    x_in : in, optional
;       x-position (as x parameter; ignored if x present)
;    y_in : in, optional
;       y-position (as y parameter; ignored if y present)
;    xrange : in, optional
;       Overrides !x.crange
;    yrange : in, optional
;       Overrides !y.crange
;    xlog : in, optional
;       Use with xrange to specify logarithmic axis. If xrange is not
;       specified, reads from !x.type
;    ylog : in, optional
;       Use with yrange to specify logarithmic axis. If yrange is not
;       specified, reads from !y.type
; :Examples:
;    plot, [0], xrange=[-10, 20], yrange=[10, 15]
;
;    print, ajs_plot_pos(0.5)
;
;           5.0000000
;
;    print, ajs_plot_pos(0.5, 0.8)
;
;           5.0000000       14.000000
;
;    print, ajs_plot_pos(y_in=0.8)
;
;           14.000000
;
;    print, ajs_plot_pos([0.3, 0.5, 0.7], [0.2, 0.4, 0.6])
;
;         -0.99999964       5.0000000       11.000000
;
;           11.000000       12.000000       13.000000
;
;    print, ajs_plot_pos(y_in=0.8, yrange=[0.1, 0.2])
;
;         0.180000
;
;    print, ajs_plot_pos(y_in=0.5, yrange=[1, 100], /ylog)
;
;          10.0000
; :History:
;    30 Apr 2008 Written, Anthony Smith
;-
FUNCTION ajs_plot_pos, x, y, x_in=x_in, y_in=y_in, $
                       xrange=xrange, yrange=yrange, xlog=xlog, ylog=ylog
  compile_opt idl2

  ;; Input values
  IF n_elements(x_in) GT 0 THEN $
     xx = x_in
  IF n_elements(y_in) GT 0 THEN $
     yy = y_in
  IF n_elements(x) GT 0 THEN BEGIN
      IF size(x, /n_dim) EQ 2 THEN BEGIN
          xx = x[*, 0]
          yy = y[*, 1]
      ENDIF ELSE BEGIN
          xx = x
          IF n_elements(y) GT 0 THEN $ ; y present only if x present (params)
             yy = y
      ENDELSE
  ENDIF
  IF n_elements(xx) EQ 0 AND n_elements(yy) EQ 0 THEN $
     message, 'No input values for x or y!'

  ;; Output values
  xtype = keyword_set(xlog)     ; linear or logarithmic
  ytype = keyword_set(ylog) 
  IF n_elements(xx) GT 0 THEN BEGIN
      IF n_elements(xrange) EQ 0 THEN BEGIN
          xr = !x.crange
          xtype = !x.type       ; Overrides input xlog
      ENDIF ELSE IF xtype THEN $
         xr = alog10(xrange)
      IF xtype THEN $
         x_out = 10 ^ (xx * (xr[1] - xr[0]) + xr[0]) $
      ELSE $
         x_out = xx * (xr[1] - xr[0]) + xr[0]
  ENDIF
  IF n_elements(yy) GT 0 THEN BEGIN
      IF n_elements(yrange) EQ 0 THEN BEGIN
          yr = !y.crange
          ytype = !y.type
      ENDIF ELSE IF ytype THEN $
         yr = alog10(yrange)
      IF ytype THEN $
         y_out = 10 ^ (yy * (yr[1] - yr[0]) + yr[0]) $
      ELSE $
         y_out = yy * (yr[1] - yr[0]) + yr[0]
  ENDIF

  ;; Return values
  IF n_elements(x_out) GT 0 AND n_elements(y_out) GT 0 THEN BEGIN
      IF n_elements(x_out) + n_elements(y_out) GT 2 THEN $
         return, [[x_out], [y_out]] $
      ELSE $
         return, [x_out, y_out]
  ENDIF ELSE IF n_elements(x_out) GT 0 THEN $
     return, x_out $
  ELSE $
     return, y_out
END
