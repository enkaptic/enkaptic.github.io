; docformat = 'rst'
;+
; Return random indices, or sort an array into random order
;
; ajs_random_sort(n) returns a random permutation of the integers 0 to
; (n - 1)
;
; ajs_random_sort(n, min=min, max=max) returns n (different) random
; integers between min and (max - 1)
;
; ajs_random_sort(array) returns the elements of array in a random order
;
; :Params:
;    input : in, required
;       Either (1) A single integer giving the number of array indices
;       to return or (2) an input array to be sorted into random order
; :Keywords:
;    min : in, optional, default=0
;       Used with a single integer input, this is the minimum array
;       index that may be returned.
;    max : in, optional, default=input
;       Used with a single integer input, this number ** - 1 ** is the
;       maximum array index that may be returned.
;    seed : in, optional
;       Seed for randomu
; :Examples:
;    print, ajs_random_sort([1, 2, 3])
;
;           1           2           3
;
;    print, ajs_random_sort(4)
;
;           3           1           0           2
;
;    print, ajs_random_sort(4, min=10, max=15)
;
;          12          10          13          14
; :History:
;    22 Apr 2008 Written, Anthony Smith
;-
FUNCTION ajs_random_sort, input, min=min, max=max, seed=seed
  compile_opt idl2

  ;; How many elements should the return array have?
  n_out = n_elements(input) 
  sort_input = n_out GT 1        ; Sort the original input?
  IF n_out EQ 1 THEN $
     n_out = input

  ;; Max & min
  IF n_elements(max) EQ 0 THEN $
     max = n_out
  IF n_elements(min) EQ 0 THEN $
     min = 0

  ;; Generate indices
  indices_unsorted = lindgen(max - min) + min
  rand = randomu(seed, max - min)
  indices_rand = indices_unsorted[sort(rand)]

  ;; Return array
  IF sort_input THEN $
     return, input[indices_rand] $
  ELSE $
     return, indices_rand[0 : n_out - 1]
END
