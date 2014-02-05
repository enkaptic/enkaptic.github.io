; docformat = 'rst'
;+
; Calculate Vmax given galaxy data and limits
;-

;+
; Return zminmax for magnitudes
;-
FUNCTION ajs_vmax_zminmax_mag, z, dm_k, var, table_row
  compile_opt idl2, hidden
  debug = ajs_debug()

  ;; If array not monotonic, apply smoothing and hope that works
  IF ~ ajs_is_monotonic(dm_k) THEN BEGIN
      FOR i = 2, 20 DO BEGIN
          IF ajs_is_monotonic(smooth(dm_k, i)) THEN BEGIN
              dm_k = smooth(dm_k, i)
              BREAK
          ENDIF
      ENDFOR
      IF debug GE 3 THEN BEGIN
          IF ajs_is_monotonic(dm_k) THEN $
             msg2 = 'smooth window = ' + strtrim(i, 2) $
          ELSE $
             msg2 = 'not fixed'
          message, 'Interpolation arr not monotonic: ' + msg2, /inf
      ENDIF
  ENDIF
  abs_value = table_row.(where(tag_names(table_row) $
                               EQ strupcase(var.abs_name)))
  IF var.app_min - abs_value LT min(dm_k) THEN $
     z1 = 0 $                   ; Outside interpolation range
  ELSE $
     z1 = interpol(z, dm_k, var.app_min - abs_value, /quadratic)
  IF var.app_max - abs_value GT max(dm_k) THEN $
     z2 = !values.f_infinity $  ; Outside interpolation range
  ELSE $
     z2 = interpol(z, dm_k, var.app_max - abs_value, /quadratic)
  IF debug GE 3 THEN BEGIN 
      z_interp = $
         interpol(z, dm_k, $
                  table_row.(where(tag_names(table_row) EQ $
                                   strupcase(var.app_name))) - abs_value, $
                  /quadratic)
      message, strjoin(strtrim([z_interp, z1, z2], 2), ' '), /inf
      IF var.band EQ 'K' THEN color = 0 ELSE color = 3
      oplot, [table_row.z], [z_interp], color=color, psym=3 ; black or green
  ENDIF
  IF z2 LT 0 THEN stop
  return, [z1, z2]
END


;+
; Return zminmax for surface brightness
;-
FUNCTION ajs_vmax_zminmax_sb, z, kc, var, table_row
  compile_opt idl2, hidden
  COMMON debug_block, debug     ; 0 none, 1 basics, 2 verbose, 3 ridiculous

  sb_dimming = 10 * alog10(1 + z) + kc
  IF debug GE 3 THEN IF NOT ajs_is_monotonic(sb_dimming) THEN $
     message, 'Interpolation array not monotonic', /inf

  abs_value = table_row.(where(tag_names(table_row) $
                               EQ strupcase(var.abs_name)))
  IF var.app_min - abs_value LT min(sb_dimming) THEN $
     z1 = 0 $                   ; Outside interpolation range
  ELSE $
     z1 = interpol(z, sb_dimming, var.app_min - abs_value, /quadratic)
  IF var.app_max - abs_value GT max(sb_dimming) THEN $
     z2 = !values.f_infinity $  ; Outside interpolation range
  ELSE $
     z2 = interpol(z, sb_dimming, var.app_max - abs_value, /quadratic)
  IF debug GE 3 THEN BEGIN 
      z_interp = $
         interpol(z, sb_dimming, $
                  table_row.(where(tag_names(table_row) EQ $
                                   strupcase(var.app_name))) - abs_value, $
                  /quadratic)
      message, strjoin(strtrim([z_interp, z1, z2], 2), ' '), /inf
;;       message, $
;;          strtrim(z_interp, 2) $
;;          + ' ' + strtrim(table_row.(where(tag_names(table_row) EQ $
;;                                           strupcase(var.app_name))), 2) $
;;          + ' ' + strtrim(table_row.(where(tag_names(table_row) EQ $
;;                                           strupcase(var.abs_name))), 2), /inf
      oplot, [table_row.z], [z_interp], color=1, psym=3 ; red
  ENDIF
  return, [z1, z2]
END


;+
; Return zminmax for radius
;-
FUNCTION ajs_vmax_zminmax_rad, z, d_a, var, table_row
  compile_opt idl2, hidden
  COMMON debug_block, debug     ; 0 none, 1 basics, 2 verbose, 3 ridiculous

  kpc_over_arcsec = !PI / 180 / 3600 * 1000 * d_a
  abs_value = table_row.(where(tag_names(table_row) $
                               EQ strupcase(var.abs_name)))
  IF debug GE 3 THEN IF NOT ajs_is_monotonic(kpc_over_arcsec) THEN $
     message, 'Interpolation array not monotonic', /inf
  IF abs_value / var.app_min LT min(kpc_over_arcsec) THEN $
     z1 = 0 $                   ; Outside interpolation range
  ELSE $
     z1 = interpol(z, kpc_over_arcsec, abs_value / var.app_max, /quadratic)
  IF abs_value / var.app_max GT max(kpc_over_arcsec) THEN $
     z2 = !values.f_infinity $  ; Outside interpolation range
  ELSE $
     z2 = interpol(z, kpc_over_arcsec, abs_value / var.app_min, /quadratic)
  IF debug GE 3 THEN BEGIN 
      z_interp = $
         interpol(z, kpc_over_arcsec, $
                  abs_value / table_row.(where(tag_names(table_row) EQ $
                                               strupcase(var.app_name))), $
                  /quadratic)
      message, strjoin(strtrim([z_interp, z1, z2], 2), ' '), /inf
      oplot, [table_row.z], [z_interp], color=2, psym=3 ; blue
  ENDIF
  return, [z1, z2]
END

;+
; Return the minimum and maximum redshifts at which the galaxy would
; be visible
;-
FUNCTION ajs_vmax_zminmax, z, dm_k, kc, d_a, var, table_row
  compile_opt idl2, hidden
  COMMON debug_block, debug     ; 0 none, 1 basics, 2 verbose, 3 ridiculous

  CASE var.type OF
      'M' : zminmax = ajs_vmax_zminmax_mag(z, dm_k, var, table_row)
      'S' : zminmax = ajs_vmax_zminmax_sb(z, kc, var, table_row)
      'R' : zminmax = ajs_vmax_zminmax_rad(z, d_a, var, table_row)
  ENDCASE

  return, zminmax
END


;+
; Calculate Vmax given galaxy data and limits
; :Returns:
;    vmax: in h^-3 Mpc^3
; :Params:
;    table : in, required
;       Array of structures with fields {coeffs:, z:} and columns for absolute
;       values given in vars. Redshift z is optional: raises warning
;       if outside observable limits.
;    vars : in, required
;       Array of structures with fields {abs_name:} giving the column
;       name in table corresponding to the absolute value, {band:} giving the
;       name of the band, {type:} giving 'M'agnitudes, 'S'urface brightness
;       or 'R'adius and {app_min:, app_max:} giving the limits.
; :Keywords:
;    zmin : in, optional
;       Minimum redshift. Default = min(table.z) or 0
;    zmax : in, optional
;       Maximum redshift. Default = max(table.z) or 2
;    area : in, optional
;       Area of survey in square degrees. Default whole sky.
;    z_min_limited_by : out, optional
;       Number of galaxies with minimum redshift limited by each of
;       the variables, or by zmin
;    z_max_limited_by : out, optional
;       Number of galaxies with maximum redshift limited by each of
;       the variables, or by zmax
;    v : out, optional
;       Volume, as in V/Vmax, for each galaxy
; :History:
;    7 Mar 2008 Created, Anthony Smith
;
;    25 Jul 2008 Added v keyword (V/Vmax), AJS
;
;    8 Aug 2008 Doesn't attempt V/Vmax when table.z does not exist
;-
FUNCTION ajs_vmax, table, vars, zmin=zmin, zmax=zmax, area=area, q0=q0, $
                   q1=q1, z_min_limited_by=z_min_limited_by, $
                   z_max_limited_by=z_max_limited_by, v=v
  compile_opt idl2, hidden
  COMMON debug_block, debug     ; 0 none, 1 basics, 2 verbose, 3 ridiculous
  IF n_elements(debug) EQ 0 THEN debug = 0 
  IF debug GE 1 THEN BEGIN
      message, 'Calculating Vmax for galaxies', /inf
      message, ajs_kw_string(zmin=zmin, zmax=zmax, area=area, q0=q0, q1=q1), $
               /inf
  ENDIF
  IF n_elements(zmin) EQ 0 THEN $
     IF where(tag_names(table) EQ 'Z') NE -1 THEN $
        zmin = min(table.z) $
     ELSE $
        zmin = 0
  IF n_elements(zmax) EQ 0 THEN $
     IF where(tag_names(table) EQ 'Z') NE -1 THEN $
        zmax = max(table.z) $
     ELSE $
        zmax = 2
  IF n_elements(area) EQ 0 THEN $
     area_frac = 1 $
  ELSE $
     area_frac = double(area) / (360L ^ 2 / !PI)
  
  ;; Calculate K-correction for each galaxy at each of a range of redshifts
  vmax = dblarr(n_elements(table))
  v = dblarr(n_elements(table))
  z = ajs_linspace(zmin, zmax, 201) ; or more - 2001?
  z_min_limited_by = lonarr(n_elements(vars) + 1)
  z_max_limited_by = lonarr(n_elements(vars) + 1) 
  d_a = ajs_angdidis(double(z))
  dm = ajs_distmod(double(z))
  IF n_elements(q0) GT 0 THEN BEGIN
      q0_k = q0[8]
      q0_r = q0[2]
  ENDIF
  IF n_elements(q1) GT 0 THEN BEGIN
      q1_k = q1[8]
      q1_r = q1[2]
  ENDIF

  ;; Lump all the K-corrections together, if not too many (1e7 too many)
  do_it_the_slow_way = n_elements(z) * n_elements(table) GT 8e6
  IF NOT do_it_the_slow_way THEN BEGIN
      IF debug GE 2 THEN message, 'Lumping K-corrections together', /inf
      z_all = reform(rebin(z, n_elements(z), n_elements(table)), $
                     n_elements(z) * n_elements(table))
      coeffs_all = rebin(table.coeffs, n_elements(table[0].coeffs), $
                         n_elements(table) * n_elements(z)) 
      IF debug GE 2 THEN message, 'Made z_all and coeffs_all', /inf
      dm_k_all_k = reform(ajs_dm_k(z_all, q0=q0_k, q1=q1_k, $
                                   coeffs=coeffs_all, band='K', $
                                   rmatrix=rmatrix_k, zvals=zvals_k, /nodm), $
                          n_elements(z), n_elements(table))
      dm_k_all_r = reform(ajs_dm_k(z_all, q0=q0_r, q1=q1_r, $
                                   coeffs=coeffs_all, band='r', $
                                   rmatrix=rmatrix_r, zvals=zvals_r, /nodm), $
                          n_elements(z), n_elements(table))
      z_all = 0                 ; Free memory
      coeffs_all = 0
      IF debug GE 1 THEN message, 'Finished setting up K-corrections', /inf
  ENDIF ELSE $
     IF debug GE 2 THEN message, 'K-corrections one by one', /inf

  IF debug GE 3 THEN BEGIN
      plot, [zmin, zmax], [zmin, zmax], /nodata
      oplot, [zmin, zmax], [zmin, zmax], color=5
  ENDIF
  FOR i = 0L, n_elements(table) - 1 DO BEGIN
      IF debug GE 3 THEN message, strtrim(i, 2), /inf
      z1 = dblarr(n_elements(vars)) 
      z2 = dblarr(n_elements(vars))
      IF debug GE 3 THEN message, 'Galaxy ' + strtrim(i, 2) $
                                  + ' True redshift ' $
                                  + strtrim(table[i].z, 2), /inf
      FOR j = 0, n_elements(vars) - 1 DO BEGIN
          IF vars[j].band EQ 'K' THEN BEGIN
              IF do_it_the_slow_way THEN $
                 kc = ajs_dm_k(z, q0=q0_k, q1=q1_k, coeffs=table[i].coeffs, $
                               band='K', rmatrix=rmatrix_k, zvals=zvals_k, $
                               /nodm) $
              ELSE $
                 kc = dm_k_all_k[*, i]
              dm_k = kc + dm
          ENDIF ELSE IF vars[j].band EQ 'r' THEN BEGIN
              IF do_it_the_slow_way THEN $
                 kc = ajs_dm_k(z, q0=q0_r, q1=q1_r, coeffs=table[i].coeffs, $
                               band='r', rmatrix=rmatrix_r, zvals=zvals_r, $
                               /nodm) $
              ELSE $
                 kc = dm_k_all_r[*, i]
              dm_k = kc + dm
          ENDIF

          zminmax = ajs_vmax_zminmax(z, dm_k, kc, d_a, vars[j], table[i])
          z1[j] = zminmax[0]
          z2[j] = zminmax[1]
      ENDFOR
      z1 = max([z1, zmin], z1_subs)
      z2 = min([z2, zmax], z2_subs)
      IF where(tag_names(table) EQ 'Z') NE -1 THEN $
         IF (z1 GT table[i].z * 1.01 OR z2 LT table[i].z / 1.01) THEN BEGIN
          message, 'Galaxy ' + strtrim(i, 2) + ' should not be visible!', /inf
          message, 'z1 = ' + strtrim(z1, 2) + ', z2 = ' + strtrim(z2, 2) $
                   + ', z = ' + strtrim(table[i].z, 2) $
                   + ' [' + strtrim(z1_subs, 2) $
                   + ',' + strtrim(z2_subs, 2) + ']', /inf
      ENDIF
      IF z1 GT z2 THEN BEGIN
          IF debug GE 3 THEN message, 'Vmax = 0 for galaxy ' + strtrim(i, 2), $
                                      /inf
          z1 = z2
      ENDIF ELSE BEGIN
          ;; Which limit provided zmin and zmax?
          z_min_limited_by[z1_subs]++
          z_max_limited_by[z2_subs]++
      ENDELSE
      vmax[i] = 4 * !PI / 3 * (lf_comvol(z2) - lf_comvol(z1)) * area_frac
      IF where(tag_names(table) EQ 'Z') NE -1 THEN $
         v[i] = 4 * !PI / 3 $
                * (lf_comvol(table[i].z) - lf_comvol(z1)) * area_frac
  ENDFOR 

  return, vmax
END

