; docformat = 'rst'
;+
; This procedure takes a two-dimensional array and creates a plot displaying
; the value of the array at each point of the array.
;
; Would be nice to plot shaded rectangles, but I can't figure
; that out. (See ajs_pixel_plot.)
;
; :Todo:
;    max_n : maximum value for n (if different from n itself)
;
;    label_max : display integer value of n if n < label_max
;
; :Params:
;    x : in, required, type="fltarr or dblarr(nx)"
;       Bin centres, evenly spaced
;    y : in, required, type="fltarr or dblarr(ny)"
;       Bin centres, evenly spaced
;    n : in, required, type="fltarr or dblarr(nx, ny) or (nx[=ny])"
;       Value in each bin
; :Keywords:
;    min_n : in, optional, type=int, default=1
;       Minimum value to display text for n
;    overplot : in, optional, type=bool, default=0
;       Set this keyword to overplot on current axes
;    _REF_EXTRA : in, optional
;       Additional keywords for plot procedure
;
; :Author: Anthony Smith
;-
PRO ajs_density_plot, x, y, n, min_n=min_n, overplot=overplot, _REF_EXTRA=e
  compile_opt idl2
  
  IF NOT keyword_set(min_n) THEN min_n = 1 
  n_x = n_elements(x)
  n_y = n_elements(y) 
  ;; Check dimensions of n: either [n_x,n_y] or [nx] == [ny]
  IF ((size(n))[0] EQ 2 AND (size(n))[1] EQ n_x AND (size(n))[2] EQ n_y) $
  THEN BEGIN
      n_print = strarr(n_x, n_y)
      dims = 2
  ENDIF ELSE IF ((size(n))[0] EQ 1 AND n_elements(n) EQ n_x AND n_x EQ n_y) $
  THEN BEGIN
      n_print = strarr(n_x)
      dims = 1
  ENDIF ELSE $
     message, 'Incorrect dimensions'
  x_binsize = min(x[where(x GT min(x))]) - min(x)
  y_binsize = min(y[where(y GT min(y))]) - min(y)
  IF NOT keyword_set(max_n) THEN max_n = max(n)

  ;; Convert numbers to strings (n>1000: '1e3' etc)
  big_numbers = where(abs(n) GE 1000, n_big_numbers, $
                      complement=small_numbers, ncomplement=n_small_numbers)
  IF n_big_numbers GT 0 THEN BEGIN
      pow = floor(alog10(abs(n[big_numbers])))
      n_print[big_numbers] = $
         string(round(n[big_numbers] / 10.^pow), format='(I0)') $
         + 'e' + string(pow, format='(I0)')
  ENDIF
  IF n_small_numbers GT 0 THEN BEGIN
      n_print[small_numbers] = string(round(n[small_numbers]), format='(I0)')
  ENDIF

  IF NOT keyword_set(overplot) THEN BEGIN
      plot, [x[0] - x_binsize / 2., x[n_x - 1] + x_binsize / 2.], $
            [y[0] - y_binsize / 2., y[n_y - 1] + y_binsize / 2.], $
            /nodata, /xstyle, /ystyle, _strict_extra=e
  ENDIF
  
  IF dims EQ 1 THEN BEGIN
      FOR i = 0, n_x - 1 DO BEGIN
          IF n[i] GE min_n THEN $
             xyouts, x[i], y[i], n_print[i], alignment=0.5
      ENDFOR
  ENDIF ELSE IF dims EQ 2 THEN BEGIN
      FOR i = 0, n_x - 1 DO BEGIN
          FOR j = 0, n_y - 1  DO BEGIN
              IF n[i, j] GE min_n THEN $
                 xyouts, x[i], y[j], n_print[i, j], alignment=0.5
          ENDFOR
      ENDFOR
  ENDIF

END

