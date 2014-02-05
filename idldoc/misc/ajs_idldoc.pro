; docformat = 'rst'
;+
; This procedure generates HTML documentation for all routines
; in Anthony Smith's IDL library, using IDLdoc
;
; Creates .tar.gz of all .pro files contained within subdirectories only.
;
; ; docformat = 'rst' (for example) should be specified at the
; top of each file, or uses idldoc default.
;
; :Examples: ajs_idldoc
;
; :History: 28 Jan 2008 Created
;
; :Author: Anthony Smith
;-
PRO ajs_idldoc
  compile_opt idl2

  ;; Public routines
  root = '/Users/anthonys/Documents/idl/ajs/'
  output = '/Users/anthonys/Sites/idldoc/'
  title = 'Anthony Smith''s IDL routines'
  subtitle = 'Galaxies etc.'
  overview = '/Users/anthonys/Documents/idl/ajs/overview.txt'
  footer = '/Users/anthonys/Documents/idl/ajs/footer.txt'

  idldoc, root=root, $
          output=output, $
          title=title, $
          subtitle=subtitle, $
          overview=overview, $
          footer=footer

  spawn, 'cd ' + file_dirname(root) + '; tar czf ' $
         + output + file_basename(root) + '.tar.gz ' $
         + file_basename(root) + '/*/*.pro'
  mg_open_url, 'http://localhost/~anthonys/idldoc/'

  ;; Private galaxy routines
  root = '/Users/anthonys/Documents/idl/gals/'
  output = '/Users/anthonys/Sites/idldoc_gals/'
  title = 'Anthony Smith''s IDL routines (galaxies)'
  subtitle = 'Galaxies etc.'
  overview = '/Users/anthonys/Documents/idl/gals/overview.txt'

  idldoc, root=root, $
          output=output, $
          title=title, $
          subtitle=subtitle, $
          overview=overview

  mg_open_url, 'http://localhost/~anthonys/idldoc_gals/'

  ;; Private misc routines
  root = '/Users/anthonys/Documents/idl/misc/'
  output = '/Users/anthonys/Sites/idldoc_misc/'
  title = 'Anthony Smith''s IDL routines (misc)'
  subtitle = 'Galaxies etc.'
  overview = '/Users/anthonys/Documents/idl/misc/overview.txt'

  idldoc, root=root, $
          output=output, $
          title=title, $
          subtitle=subtitle, $
          overview=overview

  mg_open_url, 'http://localhost/~anthonys/idldoc_misc/'


END
