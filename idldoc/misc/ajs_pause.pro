; docformat = 'rst'
;+
; This procedure invites the user to press return before continuing.
;-
PRO ajs_pause
  COMPILE_OPT idl2
  
  ans=''
  read,ans,prompt='Return to continue: '
  
END
