; docformat = 'rst'

;+
; Add new tags to an existing structure
;
; :Params:
;    struct : in, out, required
;       Array of structures with values
;    new_tag_struct : in, required
;       Single structure, giving the names and data types of the new
;       tags
; :Examples:
;    s = [{a:1, b:2}, {a:3, b:4}]
;
;    ajs_add_struct_tags, s, {c:0.0, d:''}
;
;    help, s, /structure
;
;    print, s.c
;-
PRO ajs_add_struct_tags, struct, new_tag_struct
  compile_opt idl2

  struct_tmp = struct
  struct = replicate(create_struct(struct_tmp[0], new_tag_struct), $
                     n_elements(struct_tmp)) 
  struct_assign, struct_tmp, struct

END
