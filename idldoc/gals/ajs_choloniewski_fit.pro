; docformat = 'rst'
;+
; This function finds the best-fitting Choloniewski function for
; input arrays of absolute magnitudes, surface brightness, phi
; and phi_err.
;-


;+
; Estimate phierr for bins with zero phierr
;-
PRO ajs_choloniewski_fit_phierr, absmag, sb, phi, phierr, $
                                 vars=vars, zmin=zmin, zmax=zmax, area=area, $
                                 q0=q0, q1=q1
  compile_opt idl2

  debug = ajs_debug()
  IF debug GE 1 THEN BEGIN
      message, 'Estimating phierr for bins with zero phierr', /inf
      message, ajs_kw_string(zmin=zmin, zmax=zmax, area=area, q0=q0, q1=q1), $
               /inf
  ENDIF

  IF n_elements(vars) EQ 0 $
     OR n_elements(zmin) EQ 0 $
     OR n_elements(zmax) EQ 0 $
     OR n_elements(area) EQ 0 THEN $
        message, 'Zero phierr: supply vars, zmin, zmax and area' $
  ELSE BEGIN
      IF debug GE 2 THEN print, alog10(phierr), format='(' $
                                + strtrim(n_elements(absmag), 2) + 'I2)'
      ;; Multiplier: phierr = 1. / vmax * err_multip
      ;; Taking into account the number of degrees of freedom
      ;; so that 1-sigma errors (chi^2) will give phierr (why? fudge again??)
      ;; See ajs_bbd_plot
      k = n_elements(phierr)
      chi2_bin = sqrt(2. / k)
      err_multip = 1. / sqrt(chi2_bin)
      ;; Scale according to bin widths (larger error for smaller bins)
      err_multip = err_multip / abs(absmag[1] - absmag[0]) / abs(sb[1] - sb[0])
      err_multip *= 10 ; fudge
      IF debug GE 2 THEN $
         message, 'Errors multiplied by ' + strtrim(err_multip, 2), /inf

      ;; Which bins need attention?
      zeros = where(phierr LE 0)
      
      ;; How many extra variables are we dealing with?
      n_other_vars = n_elements(vars) - 2

      ;; Default coeffs (see ajs_dm_k.pro)
      coeffs = [1.89821e-05, 2.37252e-09, 1.97939e-06, 7.54209e-08, $
                1.93300e-08]

      ;; Create structure for table
      table = create_struct('coeffs', fltarr(5))
      FOR i = 0, n_elements(vars) - 1 DO BEGIN 
          table = create_struct(vars[i].abs_name, 0.0, table)
      ENDFOR 

      IF n_other_vars EQ -1 THEN BEGIN
          ;; Consider only absmag
          IF debug GE 2 THEN message, 'Using only absmag', /inf

          coeffs = coeffs # (intarr(n_elements(zeros)) + 1) 

          ;; Create table of galaxies at centre of each zero bin
          table = replicate(table, n_elements(zeros))

          ;; Assign data to table
          table.coeffs = coeffs
          table.(where(tag_names(table) EQ strupcase(vars[0].abs_name))) $
             = absmag[reform((array_indices(phi, zeros))[0, *])]

          ;; Find vmax and set error = 1. / vmax
          vmax = ajs_vmax(table, vars, zmin=zmin, zmax=zmax, area=area, $
                          q0=q0, q1=q1)
          phierr[zeros] = 1. / vmax * err_multip
      ENDIF ELSE IF n_other_vars EQ 0 THEN BEGIN
          ;; Consider only absmag and sb
          IF debug GE 2 THEN message, 'Using only absmag & sb', /inf

          coeffs = coeffs # (intarr(n_elements(zeros)) + 1) 

          ;; Create table of galaxies at centre of each zero bin
          table = replicate(table, n_elements(zeros))

          ;; Assign data to table
          table.coeffs = coeffs
          table.(where(tag_names(table) EQ strupcase(vars[0].abs_name))) $
             = absmag[reform((array_indices(phi, zeros))[0, *])]
          table.(where(tag_names(table) EQ strupcase(vars[1].abs_name))) $
             = sb[reform((array_indices(phi, zeros))[1, *])]

          ;; Find vmax and set error = 1. / vmax
          vmax = ajs_vmax(table, vars, zmin=zmin, zmax=zmax, area=area, $
                          q0=q0, q1=q1)
          phierr[zeros] = 1. / vmax * err_multip
      ENDIF ELSE BEGIN
          ;; Consider other variables apart from absmag and sb
          IF debug GE 2 THEN message, 'Using absmag, sb and ' $
                                      + strtrim(n_other_vars, 2) $
                                      + ' other variables', /inf

          ;; Create table of galaxies at centre of each bin (zero or nonzero)
          nbins = 10
          table = replicate(table, n_elements(absmag) * n_elements(sb) $
                            * nbins ^ n_other_vars)

          ;; Assign data to table
          table.coeffs = rebin(coeffs, 5, n_elements(table))

          ;; Dimensions of n-d space of bin centres
          dims = [n_elements(absmag), n_elements(sb), $
                  intarr(n_other_vars) + nbins]
          
          ;; Bin centres for other variables
          var_arr = fltarr(n_other_vars, nbins)
          FOR i = 0, n_other_vars - 1 DO BEGIN
              var_arr[i, *] = ajs_linspace(vars[i].abs_min, vars[i].abs_max, $
                                           nbins)
          ENDFOR 

          ;; Assign values to table
          FOR j = 0, n_elements(vars) - 1 DO BEGIN 
              CASE j OF
                  0: arr_tmp = absmag
                  1: arr_tmp = sb
                  ELSE: arr_tmp = var_arr[j - 2: *]
              ENDCASE
              table.(where( $
                 tag_names(table) EQ strupcase(vars[j].abs_name))) $
                 = reform(arr_tmp[(array_indices( $
                 dims, lindgen(n_elements(table)), $
                 /dimensions))[j, *]])
          ENDFOR

          IF debug GE 2 THEN message, 'Done assigning table values', /inf
          
          ;; Find vmax and set error = 1. / vmax
          vmax = ajs_vmax(table, vars, zmin=zmin, zmax=zmax, area=area, $
                          q0=q0, q1=q1)
          vmax = reform(vmax, dims)
          FOR i = 0, n_other_vars - 1 DO BEGIN
              ;; Find maximum over one dimension
              vmax = max(vmax, dimension=3)
          ENDFOR
          phierr[zeros] = 1. / vmax[zeros] * err_multip
      ENDELSE
      IF debug GE 2 THEN print, alog10(phierr), format='(' $
                                + strtrim(n_elements(absmag), 2) + 'I2)'
  ENDELSE
END


;+
; This function finds the best-fitting Choloniewski function for
; input arrays of absolute magnitudes, surface brightness, phi
; and phi_err.
; :Returns:
;    Choloniewski parameters [mstar,alpha,phistar,sbstar,sigmasb,beta]
; :Params:
;    absmag : in, required
;       [n1] Bin centres for absolute magnitude
;    sb : in, required
;       [n2] Bin centres for surface brightness
;    phi : in, required
;       [n1, n2] Value of space density at absmag and sb
;    phierr : in, out, required
;       [n1, n2] If any values of phierr are zero, either raises an error
;       message or estimates (and changes input) values using vars etc.
;    params_init : in, optional
;       Initial guess for 
;       [mstar,alpha,phistar,sbstar,sigmasb,beta] (or default)
; :Keywords:
;    txtfile : in, optional
;       Output text file for Choloniewski parameters
;    xrange : in, optional, type=fltarr(2)
;       Range in m1 for fitting
;    yrange : in, optional, type=fltarr(2)
;       Range in m2 for fitting
;    vars : in, optional
;       For estimating errors on bins with zero phi or phierr, for Choloniewski
;       fit. Array of structures with fields {abs_name:} giving the column
;       name in table corresponding to the absolute value, {band:} giving the
;       name of the band, {type:} giving 'M'agnitudes, 'S'urface brightness
;       or 'R'adius and {app_min:, app_max: abs_min: abs_max:} giving the
;       limits. The array must be ordered so that vars[0] corresponds
;       to absmag and vars[1] corresponds to sb.
;    zmin : in, optional
;       Minimum redshift. Use with vars when performing Choloniewski
;       fit with phi = 0 and phierr = 0 in some bins
;    zmax : in, optional
;       Maximum redshift. Use with vars when performing Choloniewski
;       fit with phi = 0 and phierr = 0 in some bins
;    area : in, optional
;       Area in square degrees. Use with vars when performing Choloniewski
;       fit with phi = 0 and phierr = 0 in some bins
;    q0 : in, optional
;       Evolution parameters: for ajs_vmax
;    q1 : in, optional
;       Evolution parameters: for ajs_vmax
;    cov_mat : out, optional
;       Covariance matrix
;    _REF_EXTRA : in, optional
;       Extra keywords to pass to mpfit2dfun
; :Uses:
;    mpfit2d (idlutils)
;    ajs_choloniewski
; :History:
;    17 Jan 2008 Created, Anthony Smith
;
;    4 Apr 2008 Added xrange and yrange
;-
FUNCTION ajs_choloniewski_fit, absmag, sb, phi, phierr, params_init, $
                               txtfile=txtfile, xrange=xrange, yrange=yrange, $
                               vars=vars, zmin=zmin, $
                               zmax=zmax, area=area, q0=q0, q1=q1, $
                               cov_mat=cov_mat, _REF_EXTRA=e
  compile_opt idl2

  debug = ajs_debug()
  IF debug GE 1 THEN message, 'Choloniewski fit', /inf

  IF n_elements(params_init) EQ 0 THEN $
     params_init=[median(absmag), -1, 0.01, $
                  median(sb), 1, 1]
  params_init = double(params_init)

  ;; Bins with zero phierr: estimate phierr from Vmax
  nans = where(finite(phierr, /nan), n_nans)
  IF n_nans GT 0 THEN $
     phierr[nans] = 0
  IF min(phierr) LE 0 THEN BEGIN
      ajs_choloniewski_fit_phierr, absmag, sb, $
                                   phi, phierr, $
                                   vars=vars, zmin=zmin, zmax=zmax, $
                                   area=area, $
                                   q0=q0, q1=q1
  ENDIF

  ;; Check xrange & yrange are [min, max] not [max, min]
  IF n_elements(xrange) GT 1 THEN $
     IF xrange[0] GT xrange[1] THEN $
        xrange = reverse(xrange)
  IF n_elements(yrange) GT 1 THEN $
     IF yrange[0] GT yrange[1] THEN $
        yrange = reverse(yrange)

  ;; Restrict fit to bins within xrange and yrange
  IF n_elements(xrange) GT 1 THEN $
     fitx = where(absmag GE xrange[0] AND absmag LE xrange[1]) $
  ELSE $
     fitx = indgen(n_elements(absmag))
  IF n_elements(yrange) GT 1 THEN $
     fity = where(sb GE yrange[0] AND sb LE yrange[1]) $
  ELSE $
     fity = indgen(n_elements(sb))
  absmag_fit = absmag[fitx]
  sb_fit = sb[fity]
  phi_fit = (phi[fitx, *])[*, fity]
  phierr_fit = (phierr[fitx, *])[*, fity]
  
  ;; Make 2d arrays from 1d inputs: Array[n_elements(absmag), n_elements(sb)]
  absmag_arr = absmag_fit # (sb_fit * 0 + 1)
  sb_arr = (absmag_fit * 0 + 1) # sb_fit
  
  params=mpfit2dfun('ajs_choloniewski', $
                    absmag_arr, $
                    sb_arr, $
                    phi_fit, $
                    phierr_fit, $
                    params_init, /quiet, covar=cov_mat, perror=perror, $
                    _STRICT_EXTRA=e)

  IF keyword_set(txtfile) THEN BEGIN
      ;; Text file
      openw, unit, txtfile, /get_lun
      printf, unit, params
      free_lun, unit
  ENDIF

  RETURN, params
END


;+
; Test ajs_choloniewski_fit
;-
PRO ajs_choloniewski_fit_test
  compile_opt idl2

  params_true = [-23.5, -1, 0.01, 17, 1, 0.5]
  n_mag = 32
  n_sb = 24
  m = ajs_linspace(-20, -26, n_mag, /bincentres)
  sb = ajs_linspace(20, 14, n_sb, /bincentres)
  chol = ajs_choloniewski(m, sb, params_true)
  cholerr = chol / 10

  IF 1 THEN BEGIN 
      ;; Simulate magnitude-limited sample

      ;; Table, with one galaxy at each abs_mag
      coeffs = [1.89821e-05, 2.37252e-09, 1.97939e-06, 7.54209e-08, $
                1.93300e-08]
      table = replicate(create_struct('coeffs', coeffs, 'absmag', 0.0), n_mag)
      table.absmag = m
      
      ;; Limits
      vars = [{abs_name:'absmag', band:'K', type:'M', app_min:12, app_max:16}]
      
      ;; Calculate Vmax
      vmax = ajs_vmax(table, vars, zmin=0.01, zmax=0.3, area=500)
      vmax_bbd = rebin(vmax, n_mag, n_sb)
      
      ;; Number of galaxies at each value of absmag & sb
      n_gals = chol * vmax_bbd
      cholerr = sqrt(chol / vmax_bbd) ; = chol / sqrt(n_gals)
      min_n = 100
      IF min(n_gals) LT min_n THEN BEGIN
          chol[where(n_gals LT min_n)] = 0
          cholerr[where(n_gals LT min_n)] = 1000. / vmax_bbd[where(n_gals $
                                                                   LT min_n)]
      ENDIF   
  ENDIF
 
  ;; Fit
                                ;  params_init = params_true * 1.1
  params = ajs_choloniewski_fit(m, sb, chol, cholerr, cov_mat=cov_mat)
                                ;, params_init)
  perror = sqrt(diag_matrix(cov_mat))
  print,'True/Fitted parameters, +/- error:'
  print, transpose([[params_true], [params], [perror]])
  window, 0
  ajs_bbd_plot, m, sb, chol, phibbderr=cholerr, choloniewski=params, band='K'

  ;; Adjust errors to give 1-sigma in chi^2
  k = n_elements(cholerr)
  chi2_bin = sqrt(2. / k)
  cholerr_chi2 = sqrt(chi2_bin) * cholerr

  window, 1
  ajs_bbd_plot, m, sb, chol + cholerr_chi2, choloniewski=params, band='K'
  window, 2
  ajs_bbd_plot, m, sb, chol - cholerr_chi2, choloniewski=params, band='K'
END
