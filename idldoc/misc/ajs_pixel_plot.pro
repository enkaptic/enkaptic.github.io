; docformat = 'rst'
;+
; Plot (large number of) points as pixels (rectangles), with colour
; (white->black) representing density 
;
; :Params:
;    x : in, required
;       x-coordinate of points
;    y : in, required
;       y-coordinate of points
; :Keywords:
;    binsize1 : in, optional
;       Size of bin in x-direction
;    binsize2 : in, optional
;       Size of bin in y-direction
;    nbins1 : in, optional
;       Number of bins in x-direction (ignored if xrange and binsize1
;       set; default 50)
;    nbins2 : in, optional
;       Number of bins in y-direction (ignored if yrange and binsize2
;       set; default 50)
;    xrange : in, optional
;       Range for plot (if not overplot) and range for histogram
;    yrange : in, optional
;       Range for plot (if not overplot) and range for histogram
;    overplot : in, optional
;       Set /overplot to overplot
;    ct : in, optional
;       Colour table: -1 = use current, -2 = fixed colour (color keyword),
;       -3 (or not specified) = default (white -> black)
;       If color keyword set, and ct != -2, generate colour table from color
;    color : in, optional
;       Colour to use instead of grayscale/black
;    _REF_EXTRA : in, optional
;       Extra keywords for plot
; :Examples:
;    ajs_plot_radec, ra, dec, /nodata, xtitle='RA / degrees',
;    ytitle='dec / degrees'
;
;    ajs_pixel_plot, ra, sin(dec / !RADEG), nbins1=200, nbins2=200,
;    /overplot, ct=40
;
;    Sample output (click to enlarge):
;
;    <a href="http://astronomy.sussex.ac.uk/~anthonys/images/sdss_radec.png"><img src="http://astronomy.sussex.ac.uk/~anthonys/images/sdss_radec.png" width="200"></a>
; :History:
;    18 Apr 2008 Written, Anthony Smith
;-
PRO ajs_pixel_plot, x, y, binsize1=binsize1, binsize2=binsize2, $
                    nbins1=nbins1, nbins2=nbins2, $
                    xrange=xrange, yrange=yrange, overplot=overplot, $
                    ct=ct, color=color, _REF_EXTRA=e
  compile_opt idl2

  ;; No nbins? (Redundant if binsizes are set)
  IF n_elements(nbins1) EQ 0 THEN $
     nbins1 = 50
  IF n_elements(nbins2) EQ 0 THEN $
     nbins2 = 50

  ;; Overplot: get xrange & yrange from plot?
  IF keyword_set(overplot) THEN BEGIN
      IF n_elements(xrange) EQ 0 THEN $
         xrange = !x.crange
      IF n_elements(yrange) EQ 0 THEN $
         yrange = !y.crange
  ENDIF 

  ;; No bin sizes?
  IF n_elements(binsize1) EQ 0 THEN $
     IF n_elements(xrange) EQ 0 THEN $ 
        binsize1 = (max(x) - min(x)) / float(nbins1) $
     ELSE $
        binsize1 = (max(xrange) - min(xrange)) / float(nbins1)
  IF n_elements(binsize2) EQ 0 THEN $
     IF n_elements(yrange) EQ 0 THEN $ 
        binsize2 = (max(y) - min(y)) / float(nbins2) $
     ELSE $
        binsize2 = (max(yrange) - min(yrange)) / float(nbins2)

  ;; Range for plot (if not set)
  IF n_elements(xrange) EQ 0 THEN $ 
     xrange = [min(x), $
               (fix((max(x) - min(x)) / binsize1) + 1) * binsize1 + min(x)]
  IF n_elements(yrange) EQ 0 THEN $
     yrange = [min(y), $
               (fix((max(y) - min(y)) / binsize2) + 1) * binsize2 + min(y)]

  ;; Make histogram
  min = [min(xrange), min(yrange)]
  max = [max(xrange), max(yrange)]
  hist = hist_nd([1#x, 1#y], [binsize1, binsize2], $
                 min=min, max=max)
  dims = size(hist, /dimensions)

  ;; Edges of the bins (left/bottom)
  obin1 = dindgen(dims[0]) * binsize1 + min[0]
  obin2 = dindgen(dims[1]) * binsize2 + min[1]

  ;; Convert obin1 and obin2 to array with same dimensions as hist
  xbin = rebin(obin1, n_elements(obin1), n_elements(obin2))
  ybin = rebin(transpose(obin2), n_elements(obin1), n_elements(obin2))

  ;; Axes
  IF ~ keyword_set(overplot) THEN $ 
     plot, [0], xrange=xrange, yrange=yrange, _STRICT_EXTRA=e, /nodata

  ;; Plot pixels
  boxx = [0, 1, 1, 0, 0] * binsize1 * 1.05 ; No gaps between pixels
  boxy = [0, 0, 1, 1, 0] * binsize2 * 1.05

  nonzero = where(hist GT 0, n_nonzero)  
  IF n_nonzero GT 0 THEN BEGIN
      ;; Get current colour table (restore later)
      tvlct, r_orig, g_orig, b_orig, /get

      IF n_elements(ct) EQ 0 THEN $
         ct = -3

      IF n_elements(color) GT 0 AND ct NE -2 THEN BEGIN
          ;; Generate colour table 
          ;; From white (0) to the input colour (255)
          color_convert, r_orig, g_orig, b_orig, $
                         h_orig, l_orig, s_orig, /rgb_hls
          hls = double([h_orig[color], l_orig[color], s_orig[color]])
          h_curr = replicate(hls[0], 256)
          l_curr = 1. - findgen(256) * (1. - hls[1]) / 255
          s_curr = replicate(hls[2], 256)
          tvlct, h_curr, l_curr, s_curr, /hls
          colors = (hist[nonzero] * 255. / max(hist))
      ENDIF ELSE IF ct EQ -2 THEN BEGIN
          ;; Fixed colour
          IF n_elements(color) GT 0 THEN $ 
             colors = replicate(color, n_nonzero) $
          ELSE $
             colors = replicate(0, n_nonzero)
      ENDIF ELSE IF ct GE 0 THEN BEGIN 
          ;; Load colour table
          loadct, ct
          colors = (hist[nonzero] * 255. / max(hist))          
      ENDIF ELSE IF ct EQ -3 THEN BEGIN 
          ;; Default behaviour (white -> black)
          loadct, 0
          colors = 255 - (hist[nonzero] * 255. / max(hist))
      ENDIF ELSE IF ct EQ -1 THEN BEGIN
          ;; Use current colour table
          colors = (hist[nonzero] * 255. / max(hist))              
      ENDIF ELSE $
         message, 'Don''t know which colours to use!'

      ;; Plot
      FOR i = 0, n_nonzero - 1 DO BEGIN
          polyfill, xbin[nonzero[i]] + boxx, ybin[nonzero[i]] + boxy, $
                    color=colors[i]
      ENDFOR 

      ;; Restore colour table
      tvlct, r_orig, g_orig, b_orig
  ENDIF

END
