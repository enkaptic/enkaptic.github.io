; docformat = 'idl'

;+
; PURPOSE:
;	This function displays a widget asking for one or more options
;	to be selected.
;-

;+
; PURPOSE:
;	This procedure handles the widget event.
;-
PRO ajs_widget_choice_event, ev
  compile_opt idl2

  COMMON selected_values, button_names
  IF ev.value EQ 'Done' THEN BEGIN
      IF ev.SELECT THEN WIDGET_CONTROL, ev.TOP, /DESTROY
  ENDIF ELSE BEGIN
      IF ev.SELECT THEN WIDGET_CONTROL, ev.TOP
      IF where(button_names EQ ev.value) EQ -1 THEN $
         button_names = [button_names, ev.value] $
      ELSE $
         button_names = button_names[where(button_names NE ev.value)]
  ENDELSE
END

;+
; NAME:
;	ajs_widget_choice
;
; PURPOSE:
;	This function displays a widget asking for one or more options
;	to be selected.
;
; CATEGORY:
;	Widgets
;
; CALLING SEQUENCE:
;	Result = AJS_WIDGET_CHOICE(Values)
;
; INPUTS:
;	Values:	Array of strings
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	This function returns the selected items from Values (array of
;	strings).
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; EXAMPLE:
;	result = AJS_WIDGET_CHOICE(['Hello', 'Boo'])
;
; MODIFICATION HISTORY:
;	22 Nov 2007: Created, Anthony Smith
;-
FUNCTION ajs_widget_choice, values
  compile_opt idl2

  COMMON selected_values, button_names
  button_names = strarr(1)

  base = WIDGET_BASE(/column)
  bgroup1 = CW_BGROUP(base, values, /nonexclusive, /return_name)
  bgroup2 = CW_BGROUP(base, ['Done'], /return_name)
  
  WIDGET_CONTROL, base, /REALIZE 
  XMANAGER, 'ajs_widget_choice', base

  return, button_names
END
