; docformat = 'rst'
;+
; This procedure calculates the number counts from an input array of
; magnitudes or log(flux)
;-


;+
; This procedure calculates the number counts from an input array of
; magnitudes or log(flux)
;
; Bins must have equal width
;
; :Returns: fltarr/dblarr
;    Returns the number of galaxies per mag or log(flux), per square degree
; :Params:
;    mag : in, required
; :Keywords:
;    area : in, optional
;       Area in square degrees
;    bincentres : in, out, optional
;       Centre of each bin in magnitude
;    err_ngals : out, optional
;       Poisson errors = n / sqrt(n)
;    ninbin : out, optional
;       Number of galaxies in each bin
;    nbins : in, optional
;       Number of bins (ignored if bincentres set)
;    xrange : in, optional
;       Minimum and maximum mag (ignored if bincentres set)
; :History:
;    5 Sep 2007 Created (Anthony Smith)
;
;    6 Mar 2008 Re-written
;
;    18 Mar 2008 Added jackknife
;
;    7 Apr 2008 Jackknife estimation of Schechter function fit
;
;    15 Apr 2008 Added luminosity density
;
;    21 Apr 2008 Added double Schechter function
;-
FUNCTION ajs_number_counts, $
   mag, $
   area=area, $
   bincentres=bincentres, $
   err_ngals=err_ngals, $
   ninbin=ninbin, $
   nbins=nbins, $
   xrange=xrange
  compile_opt idl2
  debug = ajs_debug()
  IF debug GE 1 THEN message, 'Estimating number counts', /inf

  IF n_elements(bincentres) EQ 0 THEN BEGIN
      IF n_elements(nbins) EQ 0 THEN $
         nbins = 24
      IF n_elements(xrange) GT 0 THEN $
         bincentres = ajs_linspace(min(xrange), max(xrange), nbins, $
                                   /bincentres) $
      ELSE $
         bincentres = ajs_linspace(min(mag), max(mag) + 1e-6, nbins, $
                                   /bincentres)
  ENDIF
  binsize = (max(bincentres) - min(bincentres)) / (n_elements(bincentres) - 1)
  m_min = min(bincentres) - binsize / 2.

  ngals = hist1d(double(mag), min=m_min, nbins=n_elements(bincentres), $
                 binsize=binsize, obin=obin, omin=omin, omax=omax, $
                 density=ninbin) / binsize
  err_ngals = ngals / sqrt(ninbin)

  ;; Empty bins: NaN
  empty_bins = where(ninbin EQ 0, nempty)
  IF nempty GT 0 THEN BEGIN
      ngals[empty_bins] = !values.f_nan
      err_ngals[empty_bins] = !values.f_nan
  ENDIF

  IF n_elements(area) GT 0 THEN BEGIN 
      ngals /= area
      err_ngals /= area
  ENDIF

  return, ngals
END


;+
; Test ajs_number_counts
;-
PRO ajs_number_counts_test
  compile_opt idl2

  mag = ajs_linspace(10, 20, 50)
  ngals = ajs_number_counts(mag, bincentres=bincentres, $
                            area=100, nbins=30, ninbin=ninbin)
  print, ngals
  print, 'These two should be the same:'
  print, total(ninbin)
  binsize = ((max(bincentres) - min(bincentres)) / 29)
  print, total(ngals) $
         * binsize $
         * 100

END
