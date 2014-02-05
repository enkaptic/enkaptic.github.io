; docformat = 'rst'
;+
; Return 1 if the array is monotonic, 0 if not
;-

;+
; Return 1 if the array is monotonic, 0 if not
;
; :Params:
;    array, in, required, type=array
; :History:
;    12 Mar 2008 Created, Anthony Smith
; :Examples:
;    print, ajs_is_monotonic([1,2,3])
;
;       1
;
;    print, ajs_is_monotonic([1,1,3])
;
;       0
;
;    print, ajs_is_monotonic([3,2,1])
;
;       1
;-
FUNCTION ajs_is_monotonic, array
  compile_opt idl2

  ;; Change from decreasing to increasing
  increasing = array[n_elements(array) - 1] GT array[0]
  IF NOT increasing THEN $
     array = -array

  ;; Difference between consecutive elements of array
  array_difference = array[1:*] - array[0:n_elements(array) - 1]

  ;; Any neighbouring elements equal or decreasing?
  res = where(array_difference LE 0, n_decreasing)
  is_monotonic = n_decreasing EQ 0

  return, is_monotonic
END

;+
; Test ajs_is_monotonic
;-
PRO ajs_is_monotonic_test
  compile_opt idl2

  IF ajs_is_monotonic([1,2,3]) EQ 1 $
     AND ajs_is_monotonic([1,1,3]) EQ 0 $
     AND ajs_is_monotonic([3,2,1]) EQ 1 $
     AND ajs_is_monotonic([0,0,0]) EQ 0 THEN $
        print, 'Passed tests'
END
