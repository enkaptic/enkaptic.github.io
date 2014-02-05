; docformat = 'rst'
;+
; Stop plot (device/EPS settings, if required)
;
; :Keywords:
;    plotfile : in, optional
;       Name of EPS file. If EPS file required but not specified,
;       defaults to '/tmp/idl_plot.eps'
;    show_plot : in, optional
;       Close EPS file (for writing) and display the file (for viewing)
;    open : in, optional
;       Leave the plot open, unless /close is set
;    close : in, optional
;       Set /close to close an already open EPS file.
;-
PRO ajs_plot_stop, plotfile=plotfile, show_plot=show_plot, open=open, $
                   close=close
  compile_opt idl2
  debug = ajs_debug()
  IF debug GE 2 THEN BEGIN 
      message, 'Stopping plot', /inf
      message, ajs_kw_string(plotfile=plotfile, show_plot=show_plot, $
                             open=open, close=close), /inf
  ENDIF 

  ;; Close plot (if open)
  IF keyword_set(close) $
     OR ((n_elements(plotfile) GT 0 OR keyword_set(show_plot)) $
         AND NOT keyword_set(open)) $
;;           OR (NOT keyword_set(close) AND NOT keyword_set(plotfile) $
;;               AND NOT keyword_set(open) AND NOT keyword_set(show_plot))) $
  THEN BEGIN
      ;; Close Postscript file
      device, /close
      set_plot, 'X'

      ;; Change thickness
      !P.thick = 0
      !P.charthick = 0
      !X.thick = 0
      !Y.thick = 0
      !Z.thick = 0

      ;; Show plot
      IF keyword_set(show_plot) THEN BEGIN
          IF n_elements(plotfile) EQ 0 THEN $
             plotfile = '/tmp/idl_plot.eps'
          ajs_open_file, plotfile
      ENDIF
  ENDIF 

END
