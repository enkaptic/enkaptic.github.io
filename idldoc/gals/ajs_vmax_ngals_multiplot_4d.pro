; docformat = 'idl'
;+
; NAME:
;	ajs_vmax_ngals_multiplot_4d
;
; PURPOSE:
;	This procedure provides a visualisation of a four-dimensional
;	space density (number of galaxies per unit volume per unit A, B, C,
;	D - absolute magnitudes, etc.). For a four-dimensional parameter
;	space (luminosities, SBs,
;	etc), fixing the value of two of these parameters, plot
;	contours of constant Vmax along with the number of objects in
;	each bin, for the other two parameters.  Do this for all
;	combinations of parameters and all bins (six plots in all).
;
; CATEGORY:
;
; CALLING SEQUENCE:
;	ajs_vmax_ngals_multiplot_4d, $
; 		plotfile_base, axis_labels, varmin, varmax, $
;   		nbins, $
;   		bincentres1, bincentres2, bincentres3, bincentres4
;
; INPUTS:
;	plotfile_base : E.g., '/Users/.../vmax_ngalx_multiplot_'
;			Appends '01.eps' etc.
;	axis_labels :	['M_K', 'K SB', 'log10(R_K)', 'M_r']
;	varmin :	[4] Minimum values of the four parameters
;	varmax :	[4] Maximum values of the four parameters
;	nbins :		[4] Number of bins for ngalsall for each parameter
;	bincentres<n> :	[nbins[<n-1>]] Value of parameter at bin centre
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;	phi_arr :	[nbins[0] , ... , nbins[3]] Phi in each bin
;	ngalsall :	[nbins[0] , ... , nbins[3]] number in each bin
;	vmax_arr :	[nbinsvmax[0] , ...] Vmax in each bin
;	nbinsvmax :	[4] Number of bins for vmax_arr for each parameter
;			(Must be an integer multiple of nbins)
;	bincentresvmax<n> : [nbinsvmax[<n-1>]] Value of parameter at bincentre
;	show_plots :	Open plot files for viewing
;
; OUTPUTS:
;	Six plots (EPS files)
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;	Resets multiplot settings (multiplot,/default)
;
; RESTRICTIONS:
;	Parameters 1, 2 & 4 reversed to give faint -> bright
;	Parameter 3 not reversed (radius)
;
; PROCEDURE:
;
; EXAMPLE:
;	Use calling sequence verbatim
;
; MODIFICATION HISTORY:
;	23 Nov 2007 : Created, Anthony Smith
;-
PRO ajs_vmax_ngals_multiplot_4d, $
   plotfile_base, axis_labels, varmin, varmax, $
   nbins, $
   bincentres1, bincentres2, bincentres3, bincentres4, $
   phi_arr=phi_arr, ngalsall=ngalsall, $
   vmax_arr=vmax_arr, nbinsvmax=nbinsvmax, $
   bincentresvmax1=bincentresvmax1, bincentresvmax2=bincentresvmax2, $
   bincentresvmax3=bincentresvmax3, bincentresvmax4=bincentresvmax4, $
   show_plots=show_plots

  compile_opt idl2

  plot_ngals = keyword_set(ngalsall)   ; Plot ngals?
  plot_vmax = keyword_set(vmax_arr)   ; Plot vmax?
  plot_phi = keyword_set(phi_arr)     ; Plot phi?

  IF plot_vmax THEN BEGIN
      fine = nbinsvmax[0] / nbins[0]
      levels = [1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9]
      c_annotation = ['1e0', '1e1', '1e2', '1e3', '1e4', '1e5', $
                      '1e6', '1e7', '1e8', '1e9']    
  ENDIF
  IF plot_phi THEN BEGIN
      levels_phi_min_log = alog10(min(phi_arr[where(phi_arr GT 0)]))
      levels_phi_range_log = alog10(max(phi_arr)) - levels_phi_min_log
      
      levels_phi = 10 ^ (findgen(256) / 255 * levels_phi_range_log $
                           + levels_phi_min_log)
      c_colors_phi = indgen(256)
  ENDIF

  set_plot,'ps'                 ; device for graphics output
  charsizetemp = !P.charsize
  ;; Loop over six massive plots
  FOR px = 0, 2 DO BEGIN
      FOR py = px + 1, 3 DO BEGIN
          ;; Put four parameters in correct order
          p_others = ([0,1,2,3])[where([0,1,2,3] NE px AND [0,1,2,3] NE py)]
          r = execute('bc1 = bincentres' + string(px + 1, format='(I0)'))
          r = execute('bc2 = bincentres' + string(py + 1, format='(I0)'))
          r = execute('bc3 = bincentres' + string(p_others[0] + 1, $
                                                  format='(I0)'))
          r = execute('bc4 = bincentres' + string(p_others[1] + 1, $
                                                  format='(I0)'))
          IF plot_vmax THEN BEGIN
              r = execute('bcvm1 = bincentresvmax' + string(px + 1, $
                                                            format='(I0)'))
              r = execute('bcvm2 = bincentresvmax' + string(py + 1, $
                                                            format='(I0)'))
              r = execute('bcvm3 = bincentresvmax' + string(p_others[0] + 1, $
                                                            format='(I0)'))
              r = execute('bcvm4 = bincentresvmax' + string(p_others[1] + 1, $
                                                            format='(I0)'))
          ENDIF

          ;; Setup postscript file
          multiplotfile = $
             plotfile_base + string(px, format='(I0)') $
             + string(py, format='(I0)') $
             + '.eps'
          device,filename=multiplotfile,/col,/encapsulated

          !P.charsize = 0.5
          
          ;; Start new multiplot, label axes/title
          erase
          multiplot, $
             [nbins[p_others[0]],nbins[p_others[1]]], gap=0, $
             mXtitle=axis_labels[p_others[0]], $
             mYtitle=axis_labels[p_others[1]], $
             mTitle= $
             'Volume probed (h!U-3!N Mpc!U3!N) and number of galaxies, ' $
             + 'varying ' + axis_labels[p_others[0]] + ' AND ' $
             + axis_labels[p_others[1]], $
             mTitSize=0.5, mxTitSize=0.5, myTitSize=0.5, $
             mTitOffset=-0.7, mxTitOffset=-1.7, myTitOffset=-2
          
          !P.charsize=0.2
          
          ;; Loop over each of the subplots
          FOR i = 0, nbins[p_others[0]] * nbins[p_others[1]] - 1 DO BEGIN
              ;; Which plot is this? Convert i to (j, k)
              IF p_others[0] EQ 2 THEN $ 
                 ;; Left to right
                 j = i - i / nbins[p_others[0]] * nbins[p_others[0]] $
              ELSE $
                 ;; Right to left (reverse)
                 j = nbins[p_others[0]] - 1 $
                 - (i - i / nbins[p_others[0]] * nbins[p_others[0]])
              IF p_others[1] EQ 2 THEN $
                 ;; Bottom to top (reverse)
                 k = nbins[p_others[1]] - 1 - (i / nbins[p_others[0]]) $ 
              ELSE $
                 ;; Top to bottom
                 k = (i / nbins[p_others[0]])

              ;; Make arrays of the two parameters to be varied in subplot
              IF plot_ngals THEN BEGIN
                  CASE px*10 + py OF
                      01: ngalsall_temp = ngalsall[*, *, j, k]
                      02: ngalsall_temp = ngalsall[*, j, *, k]
                      03: ngalsall_temp = ngalsall[*, j, k, *]
                      12: ngalsall_temp = ngalsall[j, *, *, k]
                      13: ngalsall_temp = ngalsall[j, *, k, *]
                      23: ngalsall_temp = ngalsall[j, k, *, *]
                  ENDCASE
                  ngalsall_temp = reform(ngalsall_temp, nbins[px], nbins[py])
              ENDIF
              IF plot_phi THEN BEGIN
                  CASE px*10 + py OF
                      01: phi_arr_temp = phi_arr[*, *, j, k]
                      02: phi_arr_temp = phi_arr[*, j, *, k]
                      03: phi_arr_temp = phi_arr[*, j, k, *]
                      12: phi_arr_temp = phi_arr[j, *, *, k]
                      13: phi_arr_temp = phi_arr[j, *, k, *]
                      23: phi_arr_temp = phi_arr[j, k, *, *]
                  ENDCASE
                  phi_arr_temp = reform(phi_arr_temp, nbins[px], nbins[py])
              ENDIF
              IF plot_vmax THEN BEGIN
                  CASE px*10 + py OF
                      01: vmax_arr_temp = vmax_arr[*, *, j*fine + (fine-1)/2., $
                                                   k*fine + (fine-1)/2.]
                      02: vmax_arr_temp = vmax_arr[*, j*fine + (fine-1)/2., *, $
                                                   k*fine + (fine-1)/2.]
                      03: vmax_arr_temp = vmax_arr[*, j*fine + (fine-1)/2., $
                                                   k*fine + (fine-1)/2., *]
                      12: vmax_arr_temp = vmax_arr[j*fine + (fine-1)/2., *, *, $
                                                   k*fine + (fine-1)/2.]
                      13: vmax_arr_temp = vmax_arr[j*fine + (fine-1)/2., *, $
                                                   k*fine + (fine-1)/2., *]
                      23: vmax_arr_temp = vmax_arr[j*fine + (fine-1)/2., $
                                                   k*fine + (fine-1)/2., *, *]
                  ENDCASE
                  vmax_arr_temp = reform(vmax_arr_temp, $
                                         nbinsvmax[px], nbinsvmax[py])
              ENDIF
              
              ;; Range for the subplot axes
              IF px EQ 2 THEN $
                 xrange = [varmin[px],varmax[px]] $
              ELSE $
                 xrange = [varmax[px],varmin[px]]
              IF py EQ 2 THEN $
                 yrange = [varmin[py],varmax[py]] $
              ELSE $
                 yrange = [varmax[py],varmin[py]]
              
              ;; Titles for the subplot axes (only at the edges)
              IF i MOD nbins[p_others[0]] EQ 0 THEN $
                 ytitle = axis_labels[py] $
              ELSE $
                 ytitle = ''
              IF (nbins[p_others[0]] * nbins[p_others[1]] - i) $
                 LE nbins[p_others[0]] THEN $
                    xtitle = axis_labels[px] $
              ELSE $
                 xtitle = ''

              ;; Plot data on subplot
              plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
                    xtitle=xtitle, ytitle=ytitle, /xstyle, /ystyle
;;               IF plot_ngals THEN BEGIN
;;                   ajs_density_plot, bc1, bc2, $
;;                                     ngalsall_temp, $
;;                                     xrange=xrange, yrange=yrange, $
;;                                     xtitle=xtitle, ytitle=ytitle
;;               ENDIF ELSE BEGIN
;;                   plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
;;                         xtitle=xtitle, ytitle=ytitle, /xstyle, /ystyle
;;               ENDELSE 
              IF plot_phi THEN BEGIN
                  loadct, 33, /silent
                  contour, $
                     phi_arr_temp, $
                     bc1, bc2, /overplot, c_charsize=!P.charsize, color=1, $
                     levels=levels_phi, $ ;c_annotation=c_annotation_phi, $
                     c_color=c_colors_phi, /fill
                  set_basic_colours
              ENDIF
              IF plot_ngals THEN BEGIN
                  ajs_density_plot, bc1, bc2, $
                                    ngalsall_temp, /overplot
              ENDIF
              IF plot_vmax THEN BEGIN
                  contour, $
                     vmax_arr_temp, $
                     bcvm1, bcvm2, /overplot, $
                     levels=levels, $
                     c_annotation=c_annotation, c_charsize=!P.charsize, color=2
              ENDIF

              legend, $
                 [axis_labels[p_others[0]] + '=' $
                  + string(bc3[j], $
                           format='(G0.3)'), $
                  axis_labels[p_others[1]] + '=' $
                  + string(bc4[k], $
                           format='(G0.3)')], $
                  ;'Max V=' + string(max(vmax_arr_temp), format='(G0.3)')], $
                 /bottom, /right
              multiplot
          endfor
          device,/close
          
          multiplot,/default
          
          IF keyword_set(show_plots) THEN $
             ajs_open_file, multiplotfile
      ENDFOR 
  ENDFOR 
  
  !P.charsize = charsizetemp
  set_plot,'X'                  ; back to normal
  
END
