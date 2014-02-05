; docformat = 'rst'

;+
; Open a file using the default system application.
;
; Works on Macs only!
;
; :Params:
;    file : in, required, type=string
;       Path to the file to be opened.
;-
PRO ajs_open_file, file
  compile_opt idl2
  on_error, 2

  IF !version.os_name EQ 'Mac OS X' THEN BEGIN
      spawn, 'open ' + file
  ENDIF ELSE BEGIN
      message, 'Don''t know how to launch file: ' + file
  ENDELSE

END
