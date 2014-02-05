; docformat = 'rst'
;+
; Start plot (device/EPS settings, if required)
;
; Double thickness for all lines
; :Keywords:
;    plotfile : in, optional
;       Name of EPS file. If EPS file required but not specified,
;       set /eps and it defaults to '/tmp/idl_plot.eps'
;    eps : in, optional
;       Set /eps for EPS file
;    _REF_EXTRA : in, optional
;       Keywords for device
; :History:
;    1 Apr 2008 Written, Anthony Smith
;
;    17 Apr 2008 Replaced open and show_plot keywords with eps
;-
PRO ajs_plot_start, plotfile=plotfile, eps=eps, _REF_EXTRA=e
  compile_opt idl2
  debug = ajs_debug()
  IF debug GE 2 THEN BEGIN 
      message, 'Starting plot', /inf
      message, ajs_kw_string(plotfile=plotfile, show_plot=show_plot, $
                             open=open), /inf
  ENDIF 

  ;; Setup plot (if not already open)
  IF !D.unit EQ 0 AND (keyword_set(eps) OR n_elements(plotfile) GT 0) $
  THEN BEGIN
      ;; Default filename for EPS plot (if required)
      IF n_elements(plotfile) EQ 0 THEN $
         plotfile = '/tmp/idl_plot.eps'         

      ;; Open Postscript file
      set_plot, 'ps'            ; device for graphics output
      device, filename=plotfile, /col, /encapsulated, _STRICT_EXTRA=e

      ;; Change thickness
      !P.thick = 2
      !P.charthick = 2
      !X.thick = 2
      !Y.thick = 2
      !Z.thick = 2
  ENDIF ELSE $
     IF debug GE 2 THEN message, 'File already open!', /inf
  
END
