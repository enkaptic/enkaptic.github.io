; docformat = 'rst'
;+
; Return an integer specifying how verbose the debugging messages
; should be.
;
; 0 none, 1 basics, 2 verbose, 3 ridiculous
; :Returns:
;    The returned value is stored in a common block (unless debug is specified)
; :Params:
;    debug : in, optional
;       If debug is specified, this overrides (but does not change)
;       the shared value 
; :Keywords:
;    set_global : in, optional
;       New global debug value
; :History:
;    3 Apr 2008 Written, Anthony Smith
;-
FUNCTION ajs_debug, debug, set_global=set_global
  compile_opt idl2

  COMMON debug_block, debug_shared
  IF n_elements(debug_shared) EQ 0 THEN $
     debug_shared = 0

  ;; Set global value
  IF n_elements(set_global) GT 0 THEN $
     debug_shared = set_global

  ;; Override with input value
  IF n_elements(debug) EQ 0 THEN $
      debug = debug_shared

  return, debug
END
