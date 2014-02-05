; docformat = 'rst'
;+
; Plot 1- and 2-sigma contours for covariance matrix
;-

;+
; Plot 1- and 2-sigma contours for covariance matrix
; :Params:
;    cov_mat : in, required
;       [n, n] covariance matrix (symmetric)
; :Keywords:
;    means : in, optional, default=0
;       [n] mean of each parameter
;    labels : in, optional, type=strarr
;       [n] label for plot axes
;    reverse : in, optional
;       [n] Reverse parameters: 1 for yes, 0 for no
;    true_means : in, optional
;       [n] "True" values to show on plots
;    plotfile : in, optional
;       See ajs_plot_start
;    show_plot : in, optional
;       See ajs_plot_start
; :History:
;    31 Mar 2008 Written, Anthony Smith
;
;    14 May 2008 Added true_means
;-
PRO ajs_cov_mat_plot, cov_mat, means=means, labels=labels, reverse=reverse, $
                      true_means=true_means, plotfile=plotfile, $
                      show_plot=show_plot
  compile_opt idl2
  
  ;; Check dimensions of covariance matrix
  dims_mat = size(cov_mat, /dimensions)
  IF n_elements(dims_mat) NE 2 THEN $
     message, 'Wrong number of dimensions for covariance matrix [n, n]' $
  ELSE IF dims_mat[0] NE dims_mat[1] THEN $
     message, 'Inconsistent dimensions for covariance matrix [n, n]'
  n_vars = dims_mat[0]

  ;; Means and labels
  IF n_elements(means) GT 0 AND n_elements(means) NE n_vars THEN $
     message, 'Incorrect number of elements for means' $
  ELSE IF n_elements(means) EQ 0 THEN $
     means = intarr(n_vars)
  IF n_elements(labels) GT 0 AND n_elements(labels) NE n_vars THEN $
     message, 'Incorrect number of elements for labels' $
  ELSE IF n_elements(labels) EQ 0 THEN $
     labels = strarr(n_vars)
  IF n_elements(reverse) GT 0 AND n_elements(reverse) NE n_vars THEN $ 
     message, 'Incorrect number of elements for reverse' $
  ELSE IF n_elements(reverse) EQ 0 THEN $
     reverse = intarr(n_vars)
  IF NOT array_equal(cov_mat, transpose(cov_mat)) THEN $
     message, 'Covariance matrix not symmetrical'

  ;; Clear plot
  IF !d.name EQ 'X' THEN $
     erase
  IF n_vars GT 3 THEN $
     csize = 0.3 $
  ELSE $
     csize = 0.5
  ygap = csize / 30.

  ;; Start plot
  IF n_elements(plotfile) GT 0 OR keyword_set(show_plot) THEN BEGIN
      open = 1
      ajs_ellipse_plot, plotfile=plotfile, open=open
      thick = 2
  ENDIF ELSE BEGIN
     open = 0
     thick = 1
     csize = csize * 1.5
 ENDELSE

  ;; Disable Gaussian plots (all identical - unnecessary!)
  FOR i = 1, n_vars - 1 DO BEGIN ; Row, counting from top
      FOR j = 0, n_vars - 2 DO BEGIN ; Column, counting from left

          IF (i GT 0 OR j LT n_vars - 1) AND (i - j LT 2) THEN BEGIN
              ;; Plot something (else leave blank space)
              
              ;; Axis labels
              IF i - j EQ 1 THEN BEGIN
                  xtitle = labels[j + 1]
                  ytitle = labels[i - 1]
                  doxaxis = 1
                  doyaxis = 1
              ENDIF ELSE BEGIN
                  xtitle = ''
                  ytitle = ''
                  doxaxis = 0
                  doyaxis = 0
              ENDELSE

              ;; Next multiplot subplot (goes left then down)
              IF i EQ 1 AND j EQ 0 THEN $
                 multiplot, [n_vars - 1, n_vars - 1], doxaxis=doxaxis, $
                            doyaxis=doyaxis, ygap=ygap $
              ELSE $
                 multiplot, doxaxis=doxaxis, doyaxis=doyaxis
                 
              
              IF i EQ 0 OR j EQ n_vars - 1 THEN BEGIN
                  ;; Gaussian (disabled: check FOR loops)
;;                   IF i EQ 0 THEN BEGIN
;;                       ;; Upright Gaussian
;;                       mean = means[j + 1]
;;                       sigma = sqrt(cov_mat[j + 1, j + 1])
;;                       x = ajs_linspace(mean - 3 * sigma, mean + 3 * sigma, 101)
;;                       gauss = exp(- (x - mean) ^ 2 / 2 / sigma ^ 2) $
;;                           / sigma / sqrt(2 * !PI)
;;                       plot, x, gauss / max(gauss), /xstyle, /ystyle, $
;;                             xtitle=xtitle, ytitle=ytitle
;;                   ENDIF ELSE BEGIN
;;                       ;; Sideways Gaussian
;;                       mean = means[i - 1]
;;                       sigma = sqrt(cov_mat[i - 1, i - 1])
;;                       x = ajs_linspace(mean - 3 * sigma, mean + 3 * sigma, 101)
;;                       gauss = exp(- (x - mean) ^ 2 / 2 / sigma ^ 2) $
;;                           / sigma / sqrt(2 * !PI)
;;                       plot, gauss / max(gauss), x, /xstyle, /ystyle, $
;;                             xtitle=xtitle, ytitle=ytitle
;;                   ENDELSE
              ENDIF ELSE BEGIN
                  ;; 1- and 2-sigma contour plot
                  ;; Correlation
                  sigma_x = sqrt(cov_mat[j + 1, j + 1])
                  sigma_y = sqrt(cov_mat[i - 1, i - 1])
                  rho = cov_mat[i - 1, j + 1] / sigma_x / sigma_y
                  xrange = [means[j + 1] - 2 * sigma_x, $
                            means[j + 1] + 2 * sigma_x]
                  yrange = [means[i - 1] - 2 * sigma_y, $
                            means[i - 1] + 2 * sigma_y]
                  IF n_elements(true_means) GT 1 THEN BEGIN
                      ;; Plot '+' symbol for "true" values
                      xrange[0] = min([xrange[0], true_means[j + 1] $
                                       - abs(true_means[j + 1]) * 1e-6])
                      xrange[1] = max([xrange[1], true_means[j + 1] $
                                       + abs(true_means[j + 1]) * 1e-6])
                      yrange[0] = min([yrange[0], true_means[i - 1] $
                                       - abs(true_means[i - 1]) * 1e-6]) 
                      yrange[1] = max([yrange[1], true_means[i - 1] $
                                       + abs(true_means[i - 1]) * 1e-6])
                  ENDIF
                  IF reverse[j + 1] THEN $
                     xrange = reverse(xrange)
                  IF reverse[i - 1] THEN $
                     yrange = reverse(yrange)
                  ajs_ellipse_plot, sigma_x, sigma_y, $
                                    means[j + 1], means[i - 1], $
                                    rho=rho, n_sigma=[1, 2], $
                                    xrange=xrange, yrange=yrange, $
                                    ;;/xstyle, /ystyle, $
                                    xtitle=xtitle, ytitle=ytitle, $
                                    charsize=csize, xthick=thick, $
                                    ythick=thick, charthick=thick
                  IF n_elements(true_means) GT 1 THEN $
                     ;; Plot '+' symbol for "true" values
                     oplot, [true_means[j + 1]], [true_means[i - 1]], psym=1
              ENDELSE
          ENDIF ELSE BEGIN 
              ;; Next subplot (leave blank)
              IF i EQ 0 AND j EQ 0 THEN $
                 multiplot, [n_vars, n_vars], gap=0 $
              ELSE $
                 multiplot
          ENDELSE
          
      ENDFOR
  ENDFOR
  
  ;; Close plot (if necessary)
  IF open EQ 1 THEN $
     ajs_ellipse_plot, plotfile=plotfile, /close, show_plot=show_plot

  ;; Reset multiplot options
;;   IF n_vars GT 2 THEN $
;;      !P.charsize = charsizetemp
  multiplot, /default
END


;+
; Test ajs_cov_mat_plot
;-
PRO ajs_cov_mat_plot_test
  means = [-23.042634, -0.38052661, 0.017801707, 17.324397, 0.58152665, 0.24535320]
  reverse = [1, 0, 0, 1, 0, 0]
  labels = ['M*', '!4a!3', '!4u!3*', '!4l!3*', '!4r!Dl!U!3', '!4b!3']
  cov_mat = [[0.00011233944, 0.00012717791, 1.5205410e-06, 2.5459600e-05, 4.4394161e-06, 5.2019502e-06], $
             [0.00012717791, 0.00020450268, 1.5074665e-06, 2.4230836e-05, 3.2560968e-07, -4.1112576e-06], $
             [1.5205410e-06, 1.5074665e-06, 3.4861713e-08, 4.4584887e-07, 6.8802010e-08, 1.0825542e-07], $
             [2.5459600e-05, 2.4230836e-05, 4.4584887e-07, 2.5795079e-05, 2.0363024e-06, 1.1370081e-05], $
             [4.4394161e-06, 3.2560968e-07, 6.8802010e-08, 2.0363024e-06, 4.6717352e-06, 1.4206954e-06], $
             [5.2019502e-06, -4.1112576e-06, 1.0825542e-07, 1.1370081e-05, 1.4206954e-06, 1.5966081e-05]]
  ajs_cov_mat_plot, cov_mat, means=means, reverse=reverse, labels=labels, $
                    /show_plot
END
