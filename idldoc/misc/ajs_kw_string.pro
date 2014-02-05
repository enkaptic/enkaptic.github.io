; docformat = 'rst'
;+
; This function returns a string specifying values of keywords
;-

;+
; Return the string
;-
FUNCTION ajs_kw_string_values, empty_string, names, _EXTRA=e
  compile_opt idl2

  IF n_elements(e) EQ 0 THEN $
     e = {thisstructureisempty:'empty'}
  str_list = ['']
  FOR i = 0, n_elements(names) - 1 DO BEGIN
      tag_value = where(tag_names(e) EQ names[i])
      IF tag_value NE -1 THEN $
         str_list = [str_list, names[i] + '=' $
                     + strjoin(strtrim(e.(tag_value), 2), ' ')] $
      ELSE $
         str_list = [str_list, names[i] + '=' + empty_string]
  ENDFOR
  return, strjoin(str_list[1:*], ' ')
END

;+
; This function returns a string specifying values of keywords
; :Params:
;    empty_string, in, optional
;       String to return for undefined keywords (default "empty")
; :Examples:
;    print, ajs_kw_string(a=2, b=b)
;
;    A=2 B=empty
; :History:
;    12 Mar 2008 Written, Anthony Smith
;-
FUNCTION ajs_kw_string, empty_string, _REF_EXTRA=re
  compile_opt idl2

  IF n_elements(empty_string) EQ 0 THEN $
     empty_string = 'empty'
  return, ajs_kw_string_values(empty_string, re, _EXTRA=re)
END
