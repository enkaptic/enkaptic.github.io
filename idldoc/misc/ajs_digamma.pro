; docformat = 'rst'
;+
; Return digamma function (derivative of gamma function divided by
; gamma function)
;-


;+
; Return digamma function (derivative of gamma function divided by
; gamma function)
; :Params:
;    x, in, required
;       Argument for digamma function
;    eps, in, optional, default=1e-12
;       Precision for result
; :History:
;    15 Apr 2008 Written, Anthony Smith
;-
FUNCTION ajs_digamma, x, eps
  compile_opt idl2

  x = double(x)
  IF n_elements(eps) EQ 0 THEN $
     eps = 1e-12
  
  ;; Euler-Mascheroni constant
  gamma = 0.57721566490153286060651209008240243104215933593992d

  psi = - gamma
  n = 1
  REPEAT BEGIN
      delta_psi = (x - 1) / (n * (n + x - 1))
      psi += delta_psi
      n += 1
  ENDREP UNTIL (abs(delta_psi) LT eps)

  return, psi
END


;+
; Test ajs_digamma
;-
PRO ajs_digamma_test
  gamma = 0.57721566490153286060651209008240243104215933593992d
  print, 'Zeros:'
  print, ajs_digamma(1) - (- gamma)
  print, ajs_digamma(0.5) - (- gamma - 2 * alog(2))
  print, ajs_digamma(1./3) - (- !PI / 2 / sqrt(3) - 3. / 2 * alog(3) - gamma)
END
