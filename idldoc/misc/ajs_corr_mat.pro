; docformat = 'rst'
;+
; Convert covariance matrix to correlation matrix
; :Returns:
;    Correlation matrix
; :Params:
;    cov_mat : in, required
;       Covariance matrix
; :History:
;    7 Apr 2008 Written, Anthony Smith
;-
FUNCTION ajs_corr_mat, cov_mat  
  compile_opt idl2

  corr_mat = cov_mat / float(cov_mat) ; Size, type, diagonal 1's
  FOR i = 0, sqrt(n_elements(cov_mat)) - 1 DO BEGIN
      FOR j = i, sqrt(n_elements(cov_mat)) - 1 DO BEGIN
          corr_mat[i, j] = cov_mat[i, j] / sqrt(cov_mat[i, i] * cov_mat[j, j])
          corr_mat[j, i] = corr_mat[i, j]
      ENDFOR 
  ENDFOR

  return, corr_mat
END
