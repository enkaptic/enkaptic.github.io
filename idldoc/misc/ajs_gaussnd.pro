; docformat = 'rst'
;+
; This function returns the probability density at a given point
; _x_ for a multivariate normal distribution with specified
; mean, covariance matrix and amplitude.
;
; Compatible with mpfit
; :Params:
;    x : in, required
;       [n,N] Array of n-dimensional input values
;    params : in, optional
;       [1 + n + (n*(n+1))/2] = [amp,mean,cov]
;
;	        WHERE	amp	amplitude (scalar)
;
;		mean	[n] Mean of each dimension
;
;		cov	[(n*(n+1))/2] Unique components of symmetric
;			covariance matrix, e.g., 00, 01, 02, 11, 12, 22 for 3-d
; :Keywords:
;    amp : in, optional
;       amplitude (optional; default=1)
;    mean : in, optional
;       [n] Mean of each dimension
;    cov : in, optional
;       [n,n] Covariance matrix
; :Returns:
;    [N] probability density of specified normal distribution
;    at _x_
; :Examples:
;    Identical result:
;
;    y=ajs_gaussnd([1,2,3],mean=[0,0,0],cov=[[1,0,0],[0,1,0],[0,0,1]])
;
;    y=ajs_gaussnd([1,2,3],amp=1,mean=[0,0,0],cov=[[1,0,0],[0,1,0],[0,0,1]])
;
;    y=ajs_gaussnd([1,2,3],[1,0,0,0,1,0,0,1,0,1])
; :History:
;    7 Sep 2007: Created, Anthony Smith
;-
FUNCTION ajs_gaussnd,x,params,amp=amp,mean=mean,cov=cov
  compile_opt idl2

  IF n_params() EQ 2 THEN BEGIN
      n=(-3 + sqrt(9+8*(n_elements(params)-1)))/2
      amp = params[0]
      mean = params[1:n]
      cov=dblarr(n,n)
      k=0
      FOR i=0,n-1 DO BEGIN
          FOR j=i,n-1 DO BEGIN
              cov[i,j]=params[n+1+k]
              cov[j,i]=params[n+1+k]
              k=k+1
          ENDFOR
      ENDFOR
  ENDIF ELSE BEGIN
      n = n_elements(mean)
      IF n_elements(amp) EQ 0 THEN amp=1
  ENDELSE

  ;; For each element in the input array, x, calculate y=f(x)
  y=dblarr(n_elements(x)/n) 
  a = amp * 1./((2*!PI)^(n/2.)*sqrt(abs(determ(cov,/double))))
  invertcov = invert(cov,/double)
  FOR i = 0,n_elements(x)/n-1 DO BEGIN
      y[i] = a * exp( -1./2 * (x[*,i]-mean) ## invertcov ## (1#(x[*,i]-mean)) )
  ENDFOR 

  RETURN,y
END
