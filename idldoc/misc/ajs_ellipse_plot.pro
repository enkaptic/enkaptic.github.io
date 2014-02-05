; docformat = 'rst'
;+
; Plot an ellipse
;
; Wrapper on plot or oplot
; :Params:
;    a : in, required
;       Semimajor/minor x-axis
;    b : in, required
;       Semimajor/minor y-axis
;    xc : in, optional, default=0
;       x-co-ordinate of centre of ellipse
;    yc : in, optional, default=0
;       y-co-ordinate of centre of ellipse
; :Keywords:
;    ang : in, optional
;       Positional angle (degrees anticlockwise)
;    rho : in, optional
;       Correlation between x and y, with a and b now interpreted as
;       standard deviation in x and y, to give <n_sigma>-sigma contours
;
;       Equivalently, plot an ellipse within a rectangle bounded by xc
;       +/- a and yc +/- b, where the ellipse touches the rectangle at
;       +/- (xc + rho * a, yc + b) and +/- (xc + a, yc + rho * b)
;    n_sigma : in, optional, default=1
;       With rho, plot <n_sigma>-sigma contour(s) (equivalent to
;       multiplier on a and b). Can be an array of values
;    overplot : in, optional
;       Set /overplot to overplot (uses oplot rather than plot)
;    plotfile : in, optional
;       See ajs_plot_start
;    show_plot : in, optional
;       See ajs_plot_start
;    open : in, optional
;       See ajs_plot_start
;    close : in, optional
;       See ajs_plot_start
;    _REF_EXTRA
;       Keywords for plot/oplot
; :Examples:
;    ajs_ellipse_plot, 1, 1
;
;    ajs_ellipse_plot, 2, 3, 10, 20
;
;    ajs_ellipse_plot, 4, 1, ang=30
;
;    ajs_ellipse_plot, 1, 1, rho=0.8, /overplot
;
;    ajs_ellipse_plot, 3, 1, rho=0.5, n_sigma=2
;
;    ajs_ellipse_plot, 3, 1, rho=-0.5, n_sigma=[1,3,2]
; :History:
;    28 Mar 2008 Written, Anthony Smith
;-
PRO ajs_ellipse_plot, a, b, xc, yc, ang=ang, rho=rho, overplot=overplot, $
                      n_sigma=n_sigma, plotfile=plotfile, $
                      show_plot=show_plot, open=open, close=close, _REF_EXTRA=e
  compile_opt idl2
  debug = ajs_debug()
  IF debug GE 2 THEN BEGIN 
      message, 'Plotting ellipse', /inf
      message, ajs_kw_string(a=a, b=b, xc=xc, xy=xy, ang=ang, rho=rho, $
                             overplot=overplot, n_sigma=n_sigma, $
                             plotfile=plotfile, show_plot=show_plot, $
                             open=open, close=close), /inf
  ENDIF 

  ;; Setup plot
  eps = keyword_set(open) OR keyword_set(show_plot)
  ajs_plot_start, plotfile=plotfile, eps=eps

  IF NOT keyword_set(close) AND n_elements(a) GT 0 THEN BEGIN 
      ;; Plot <n_sigma>-sigma contours?
      IF n_elements(n_sigma) EQ 0 THEN BEGIN
          ;; Change names so input arrays do not get altered
          c = a
          d = b
      ENDIF ELSE IF n_elements(n_sigma) EQ 1 THEN BEGIN
          c = a * n_sigma
          d = b * n_sigma
      ENDIF ELSE BEGIN
          ;; Loop over different values of n_sigma and then return
          n_sigma_sorted = n_sigma[reverse(sort(n_sigma))]
          
          ;; Optional overplot for first (largest) contour
          IF n_elements(plotfile) GT 0 OR keyword_set(show_plot) THEN $
             open_tmp = 1 $
          ELSE $
             open_tmp = keyword_set(open)
          ajs_ellipse_plot, a, b, xc, yc, ang=ang, rho=rho, $
                            overplot=overplot, $
                            n_sigma=n_sigma_sorted[0], plotfile=plotfile, $
                            open=open_tmp, show_plot=show_plot, _STRICT_EXTRA=e

          ;; Loop over the rest
          FOR i = 1, n_elements(n_sigma) - 1 DO BEGIN
              ajs_ellipse_plot, a, b, xc, yc, ang=ang, rho=rho, /overplot, $
                                n_sigma=n_sigma_sorted[i], open=open_tmp, $
                                _STRICT_EXTRA=e
          ENDFOR
          
          IF n_elements(plotfile) GT 0 AND NOT keyword_set(open) THEN $ 
             ajs_ellipse_plot, plotfile=plotfile, /close, show_plot=show_plot

          ;; Done
          return
      ENDELSE

      ;; Make a circle
      theta = ajs_linspace(0, 360, 181) / !RADEG ; Radians
      x = cos(theta)
      y = sin(theta)
      
      ;; Bounded by a rectangle?
      IF n_elements(rho) GT 0 THEN BEGIN
          ;; Centred on origin, equation for bivariate Gaussian gives
          ;; 1-sigma contour as (ellipse)
          ;; 1 = (1 / (1 - rho^2)) * [ x^2 / sigma_x^2 + y^2 / sigma_y^2 
          ;;                           - 2 rho x y / (sigma_x sigma_y) ]
          ;; Quadratic function => semimajor/minor axes given by eigenvalues
          ;; 1 = e_1 X^2 + e_2 Y^2
          ;; where ellipse is given by   e_1 = 1 / c^2
          ;;                             e_2 = 1 / d^2
          
          ;; E-vals of symm. matrix of above fn, with c = sigma_x, d = sigma_y
          evals = eigenql([[1. / (c ^ 2. * (1. - rho ^ 2.)), $
                            - rho / (c * d * (1. - rho ^ 2.))], $
                           [- rho / (c * d * (1. - rho ^ 2.)), $
                            1. / (d ^ 2. * (1. - rho ^ 2.))]], $
                          eigenvectors=evecs, /double)
          
          ;; Axes of unrotated ellipse (to be rotated later on)
          c = 1. / sqrt(evals[0])
          d = 1. / sqrt(evals[1])
          
          ;; Matrix of eigenvectors will be the rotation matrix => find angle
          ang = atan(evecs[0, 1] / evecs[0, 0]) * !RADEG
          
          ;; Eigenvalues are in descending order, so may need to reverse rot.
          ;; (Seems to work, not sure why!)
          IF rho GT 0 THEN $
             ang = -abs(ang) $
          ELSE $
             ang = abs(ang)
      ENDIF
      
      ;; Change its size
      x = x * c
      y = y * d
      
      ;; Rotate it
      IF n_elements(ang) GT 0 THEN BEGIN 
          x0 = x
          y0 = y
          cosang = cos(ang / !RADEG)
          sinang = sin(ang / !RADEG)
          x = x0 * cosang - y0 * sinang
          y = x0 * sinang + y0 * cosang
      ENDIF 
      
      ;; Shift it
      IF n_elements(xc) GT 0 THEN $
         x = x + xc
      IF n_elements(yc) GT 0 THEN $
         y = y + yc

      ;; Plot it  
      IF keyword_set(overplot) THEN $
         oplot, x, y, _EXTRA=e $
      ELSE $
         plot, x, y, _STRICT_EXTRA=e
  ENDIF

  ajs_plot_stop, plotfile=plotfile, show_plot=show_plot, open=open, close=close
END
