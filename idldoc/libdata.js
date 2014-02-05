/* Index used for searching */
/*
   Fields used:
     url, name, type, filename, authors, routine name, comments, parameters,
     categories, and attributes
*/
title = "Anthony Smith's IDL routines";
subtitle = "Galaxies etc.";
libdata = new Array();
libdataItem = 0;



libdata[libdataItem++] = new Array("misc/ajs_add_struct_tags.html", "ajs_add_struct_tags.pro", ".pro file in misc/ directory", "ajs_add_struct_tags.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_add_struct_tags.html#ajs_add_struct_tags", "ajs_add_struct_tags", "routine in ajs_add_struct_tags.pro", "ajs_add_struct_tags.pro", "", "ajs_add_struct_tags", " Add new tags to an existing structure    ", "struct       Array of structures with values  new_tag_struct       Single structure, giving the names and data types of the new        tags  ", "          -1", "      s = [{a:1, b:2}, {a:3, b:4}]       ajs_add_struct_tags, s, {c:0.0, d:''}       help, s, /structure       print, s.c  ");
  
  

libdata[libdataItem++] = new Array("gals/ajs_angdidis.html", "ajs_angdidis.pro", ".pro file in gals/ directory", "ajs_angdidis.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_angdidis.html#ajs_angdidis", "ajs_angdidis", "routine in ajs_angdidis.pro", "ajs_angdidis.pro", "", "ajs_angdidis", " Return angular diameter distance to z    See ajs_comdis for keyword parameters  ", "_REF_EXTRAz       Redshift (or array of values)  ", "          -1", "");
  
  

libdata[libdataItem++] = new Array("gals/ajs_bbd_plot.html", "ajs_bbd_plot.pro", ".pro file in gals/ directory", "ajs_bbd_plot.pro", "", "", " This procedure produces a surface/contour plot of the  bivariate brightness distribution (number of galaxies per unit  volume as a function of absolute magnitude and surface brightness)  ", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_bbd_plot.html#ajs_bbd_plot", "ajs_bbd_plot", "routine in ajs_bbd_plot.pro", "ajs_bbd_plot.pro", "", "ajs_bbd_plot", " Main procedure.    ", "phibbderr       1-sigma errors in phibbd: converted to 1-sigma errors in        corresponding chi^2, then plotted as dotted/dashed contours. Set        /phibbderr, or supply variable, with m1, m2 to calculate then        plot errors  error_contours       Set /error_contours to force plotting of errors (switches off        if Choloniewski function is plot)  plotfile       Destination for output EPS plot  ct       Color table to use (default 1)  show_plot       Open EPS plot for display  band       Name of waveband (for axis titles)  titlextitleytitlecolorbartitlemin_log_phimax_log_phicholoniewski       If a Choloniewski fit is performed, error contours will not be        displayed, unless /error_contours is set  chol_xrange       xrange for Choloniewski fit (ajs_choloniewski_fit)  chol_yrange       yrange for Choloniewski fit (ajs_choloniewski_fit)  cov_mat       Covariance matrix of Choloniewski fit  m1       Values of first magnitude for each galaxy (if phibbd empty)  m2       Values of second magnitude for each galaxy (if phibbd empty)  weights       Weight (e.g., 1/Vmax) for each galaxy (if phibbd empty). See        ajs_vmax_bbd for more details.  jackknife       Integer for each galaxy giving the jackknife region in which        the galaxy lies  xrange       Range of plot  yrange       Range of plot  nbins1       Number of bins for BBD calculation from m2 and weights  nbins2       Number of bins for BBD calculation from m2 and weights  vars       See ajs_choloniewski_fit  zmin       See ajs_choloniewski_fit  zmax       See ajs_choloniewski_fit  area       See ajs_choloniewski_fit  q0       See ajs_choloniewski_fit  q1       See ajs_choloniewski_fit  bincentres1       Centres of absolute magnitude bins (x-axis). Optional        if calculating BBD from m1, m2 and weights  bincentres2       Centres of surface brightness bins (y-axis). Optional        if calculating BBD from m1, m2 and weights  phibbd       Value of phi at bincentres1[i], bincentres2[j]. Not required        if calculating BBD from m1, m2 and weights.    ", "          -1", "      11 Jan 2008 Created, Anthony Smith       18 Mar 2008 Added jackknife       25 Mar 2008 Added keywords for Choloniewski fit       7 Apr 2008 Choloniewski included in ajs_vmax_bbd call        ajs_bbd_plot, bincentres1, bincentres2, phibbd       ajs_bbd_plot, m1=abs_mag, m2=surface_brightness, weights=weights  ");
  
  libdata[libdataItem++] = new Array("gals/ajs_bbd_plot.html#ajs_bbd_plot_test", "ajs_bbd_plot_test", "routine in ajs_bbd_plot.pro", "ajs_bbd_plot.pro", "", "ajs_bbd_plot_test", " Test ajs_bbd_plot  ", "_REF_EXTRAshow_plot", "          -1", "");
  
  

libdata[libdataItem++] = new Array("gals/ajs_bbd_read.html", "ajs_bbd_read.pro", ".pro file in gals/ directory", "ajs_bbd_read.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_bbd_read.html#ajs_bbd_read", "ajs_bbd_read", "routine in ajs_bbd_read.pro", "ajs_bbd_read.pro", "", "ajs_bbd_read", "	This procedure reads BBD (bivariate brightness distribution) 	data from a binary file created using ajs_bbd_write.                    ", "datfilebincentres1bincentres2phibbdphibbderr", "          -1", "	Written by: Anthony Smith 	July 1994 Describe modifications 	11 Jan 2008 Created, Anthony Smith       ");
  
  

libdata[libdataItem++] = new Array("gals/ajs_bbd_write.html", "ajs_bbd_write.pro", ".pro file in gals/ directory", "ajs_bbd_write.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_bbd_write.html#ajs_bbd_write", "ajs_bbd_write", "routine in ajs_bbd_write.pro", "ajs_bbd_write.pro", "", "ajs_bbd_write", " This procedure takes arrays describing the BBD (bivariate brightness  distribution) and writes the data either to a binary file or to a  text file.    ", "datfile       Name of binary file for output (read using ajs_bbd_read)  txtfile       Name of text file for output  bincentres1bincentres2phibbdphibbderr", "          -1", "      11 Jan 2008 Created, Anthony Smith  ");
  
  

libdata[libdataItem++] = new Array("gals/ajs_choloniewski.html", "ajs_choloniewski.pro", ".pro file in gals/ directory", "ajs_choloniewski.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_choloniewski.html#ajs_choloniewski", "ajs_choloniewski", "routine in ajs_choloniewski.pro", "ajs_choloniewski.pro", "", "ajs_choloniewski", "	This function generates a Choloniewski function along an input 	array of absolute magnitudes and surface brightnesses.            	Compatible with mpfit2d.  	Magnitudes only       ", "mstaralphaphistarsbstarsigmasbbetaabsmagsbparams[mstar,alpha,phistar,sbstar,sigmasb,beta] OR input separately:   ", "          -1", " 	This function returns an rray of phi for corresponding absmag, sbr    ");
  
  

libdata[libdataItem++] = new Array("gals/ajs_choloniewski_fit.html", "ajs_choloniewski_fit.pro", ".pro file in gals/ directory", "ajs_choloniewski_fit.pro", "", "", " This function finds the best-fitting Choloniewski function for  input arrays of absolute magnitudes, surface brightness, phi  and phi_err.  ", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_choloniewski_fit.html#ajs_choloniewski_fit_phierr", "ajs_choloniewski_fit_phierr", "routine in ajs_choloniewski_fit.pro", "ajs_choloniewski_fit.pro", "", "ajs_choloniewski_fit_phierr", " Estimate phierr for bins with zero phierr  ", "varszminzmaxareaq0q1absmagsbphiphierr", "          -1", "");
  
  libdata[libdataItem++] = new Array("gals/ajs_choloniewski_fit.html#ajs_choloniewski_fit", "ajs_choloniewski_fit", "routine in ajs_choloniewski_fit.pro", "ajs_choloniewski_fit.pro", "", "ajs_choloniewski_fit", " This function finds the best-fitting Choloniewski function for  input arrays of absolute magnitudes, surface brightness, phi  and phi_err.  ", "txtfile       Output text file for Choloniewski parameters  xrange       Range in m1 for fitting  yrange       Range in m2 for fitting  vars       For estimating errors on bins with zero phi or phierr, for Choloniewski        fit. Array of structures with fields {abs_name:} giving the column        name in table corresponding to the absolute value, {band:} giving the        name of the band, {type:} giving 'M'agnitudes, 'S'urface brightness        or 'R'adius and {app_min:, app_max: abs_min: abs_max:} giving the        limits. The array must be ordered so that vars[0] corresponds        to absmag and vars[1] corresponds to sb.  zmin       Minimum redshift. Use with vars when performing Choloniewski        fit with phi = 0 and phierr = 0 in some bins  zmax       Maximum redshift. Use with vars when performing Choloniewski        fit with phi = 0 and phierr = 0 in some bins  area       Area in square degrees. Use with vars when performing Choloniewski        fit with phi = 0 and phierr = 0 in some bins  q0       Evolution parameters: for ajs_vmax  q1       Evolution parameters: for ajs_vmax  cov_mat       Covariance matrix  _REF_EXTRA       Extra keywords to pass to mpfit2dfun  absmag       [n1] Bin centres for absolute magnitude  sb       [n2] Bin centres for surface brightness  phi       [n1, n2] Value of space density at absmag and sb  phierr       [n1, n2] If any values of phierr are zero, either raises an error        message or estimates (and changes input) values using vars etc.  params_init       Initial guess for        [mstar,alpha,phistar,sbstar,sigmasb,beta] (or default)  ", "          -1", "      17 Jan 2008 Created, Anthony Smith       4 Apr 2008 Added xrange and yrange        mpfit2d (idlutils)     ajs_choloniewski        Choloniewski parameters [mstar,alpha,phistar,sbstar,sigmasb,beta]  ");
  
  libdata[libdataItem++] = new Array("gals/ajs_choloniewski_fit.html#ajs_choloniewski_fit_test", "ajs_choloniewski_fit_test", "routine in ajs_choloniewski_fit.pro", "ajs_choloniewski_fit.pro", "", "ajs_choloniewski_fit_test", " Test ajs_choloniewski_fit  ", "", "          -1", "");
  
  

libdata[libdataItem++] = new Array("gals/ajs_comdis.html", "ajs_comdis.pro", ".pro file in gals/ directory", "ajs_comdis.pro", "", "", " Return comoving distance to z, in Mpc  ", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_comdis.html#ajs_comdis_er", "ajs_comdis_er", "routine in ajs_comdis.pro", "ajs_comdis.pro", "", "ajs_comdis_er", " Integrand for comoving distance integral  ", "z", "          -1", "");
  
  libdata[libdataItem++] = new Array("gals/ajs_comdis.html#ajs_lookback_er", "ajs_lookback_er", "routine in ajs_comdis.pro", "ajs_comdis.pro", "", "ajs_lookback_er", " Integrand for lookback time integral  ", "z", "          -1", "");
  
  libdata[libdataItem++] = new Array("gals/ajs_comdis.html#ajs_comdis", "ajs_comdis", "routine in ajs_comdis.pro", "ajs_comdis.pro", "", "ajs_comdis", " Return comoving distance to z, in Mpc  ", "omega_m       Matter density (default 0.3)  omega_l       Vacuum energy density (default 0.7)  h0       Hubble constant (default 100)  c       Speed of light (default 299792.458 km/s)  slow       Set /slow to force calculation for each value (defaults to        interpolation if input array has more then 10,000 elements)  lookback_time       Set /lookback_time to return lookback time in seconds.  See        ajs_lookback_time and ajs_zage.  z       Redshift (or array of values)  ", "          -1", "      13 Mar 2008 Written, Anthony Smith       10 Oct 2008 Added lookback_time, AJS       16 Jan 2009 Corrected bug, which was affecting large input arrays  ");
  
  

libdata[libdataItem++] = new Array("gals/ajs_comvol.html", "ajs_comvol.pro", ".pro file in gals/ directory", "ajs_comvol.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_comvol.html#ajs_comvol", "ajs_comvol", "routine in ajs_comvol.pro", "ajs_comvol.pro", "", "ajs_comvol", " Return comoving volume out to redshift z (flat universe only)    See ajs_comdis for other keyword parameters  ", "area       Area in square degrees. Defaults to whole sky.  _REF_EXTRAz       Redshift (or array of values)  ", "          -1", "");
  
  

libdata[libdataItem++] = new Array("misc/ajs_corr_mat.html", "ajs_corr_mat.pro", ".pro file in misc/ directory", "ajs_corr_mat.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_corr_mat.html#ajs_corr_mat", "ajs_corr_mat", "routine in ajs_corr_mat.pro", "ajs_corr_mat.pro", "", "ajs_corr_mat", " Convert covariance matrix to correlation matrix  ", "cov_mat       Covariance matrix  ", "          -1", "      7 Apr 2008 Written, Anthony Smith        Correlation matrix  ");
  
  

libdata[libdataItem++] = new Array("misc/ajs_cov_mat_plot.html", "ajs_cov_mat_plot.pro", ".pro file in misc/ directory", "ajs_cov_mat_plot.pro", "", "", " Plot 1- and 2-sigma contours for covariance matrix  ", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_cov_mat_plot.html#ajs_cov_mat_plot", "ajs_cov_mat_plot", "routine in ajs_cov_mat_plot.pro", "ajs_cov_mat_plot.pro", "", "ajs_cov_mat_plot", " Plot 1- and 2-sigma contours for covariance matrix  ", "means       [n] mean of each parameter  labels       [n] label for plot axes  reverse       [n] Reverse parameters: 1 for yes, 0 for no  true_means       [n]  True  values to show on plots  plotfile       See ajs_plot_start  show_plot       See ajs_plot_start  cov_mat       [n, n] covariance matrix (symmetric)  ", "          -1", "      31 Mar 2008 Written, Anthony Smith       14 May 2008 Added true_means  ");
  
  libdata[libdataItem++] = new Array("misc/ajs_cov_mat_plot.html#ajs_cov_mat_plot_test", "ajs_cov_mat_plot_test", "routine in ajs_cov_mat_plot.pro", "ajs_cov_mat_plot.pro", "", "ajs_cov_mat_plot_test", " Test ajs_cov_mat_plot  ", "", "          -1", "");
  
  

libdata[libdataItem++] = new Array("misc/ajs_debug.html", "ajs_debug.pro", ".pro file in misc/ directory", "ajs_debug.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_debug.html#ajs_debug", "ajs_debug", "routine in ajs_debug.pro", "ajs_debug.pro", "", "ajs_debug", " Return an integer specifying how verbose the debugging messages  should be.    0 none, 1 basics, 2 verbose, 3 ridiculous  ", "set_global       New global debug value  debug       If debug is specified, this overrides (but does not change)        the shared value  ", "          -1", "      3 Apr 2008 Written, Anthony Smith        The returned value is stored in a common block (unless debug is specified)  ");
  
  

libdata[libdataItem++] = new Array("misc/ajs_density_plot.html", "ajs_density_plot.pro", ".pro file in misc/ directory", "ajs_density_plot.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_density_plot.html#ajs_density_plot", "ajs_density_plot", "routine in ajs_density_plot.pro", "ajs_density_plot.pro", "Anthony Smith  ", "ajs_density_plot", " This procedure takes a two-dimensional array and creates a plot displaying  the value of the array at each point of the array.    Would be nice to plot shaded rectangles, but I can't figure  that out. (See ajs_pixel_plot.)    ", "min_n       Minimum value to display text for n  overplot       Set this keyword to overplot on current axes  _REF_EXTRA       Additional keywords for plot procedure    x       Bin centres, evenly spaced  y       Bin centres, evenly spaced  n       Value in each bin  ", "          -1", "      max_n : maximum value for n (if different from n itself)       label_max : display integer value of n if n < label_max    ");
  
  

libdata[libdataItem++] = new Array("misc/ajs_digamma.html", "ajs_digamma.pro", ".pro file in misc/ directory", "ajs_digamma.pro", "", "", " Return digamma function (derivative of gamma function divided by  gamma function)  ", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_digamma.html#ajs_digamma", "ajs_digamma", "routine in ajs_digamma.pro", "ajs_digamma.pro", "", "ajs_digamma", " Return digamma function (derivative of gamma function divided by  gamma function)  ", "x       Argument for digamma function  eps       Precision for result  ", "          -1", "      15 Apr 2008 Written, Anthony Smith  ");
  
  libdata[libdataItem++] = new Array("misc/ajs_digamma.html#ajs_digamma_test", "ajs_digamma_test", "routine in ajs_digamma.pro", "ajs_digamma.pro", "", "ajs_digamma_test", " Test ajs_digamma  ", "", "          -1", "");
  
  

libdata[libdataItem++] = new Array("gals/ajs_distmod.html", "ajs_distmod.pro", ".pro file in gals/ directory", "ajs_distmod.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_distmod.html#ajs_distmod", "ajs_distmod", "routine in ajs_distmod.pro", "ajs_distmod.pro", "", "ajs_distmod", " Return distance modulus to z    See ajs_comdis for keyword parameters  ", "_REF_EXTRAz       Redshift (or array of values)  ", "          -1", "");
  
  

libdata[libdataItem++] = new Array("gals/ajs_dm_k.html", "ajs_dm_k.pro", ".pro file in gals/ directory", "ajs_dm_k.pro", "", "", " This function returns some combination of the distance modulus,  K-correction and evolution correction for a galaxy at a redshift z.  ", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_dm_k.html#ajs_dm_k", "ajs_dm_k", "routine in ajs_dm_k.pro", "ajs_dm_k.pro", "Anthony Smith  ", "ajs_dm_k", " Includes K-correction if band keyword is set, and evolution correction if q0  and/or q1 is set. K(z) and/or E(z) without the distance modulus if  /nodm is set.    Uses Michael Blanton's K-correct code, with coefficients of a  typical galaxy for the K-corrections.    ", "band       Band for K-corrections: choose from ugrizYJHK  nodm       Set /nodm to return K(z) - E(z) only (no distance modulus)  q0q1       q0 and q1 define the evolution correction, E(z) = q0 * (1 + q1        * z) * z  coeffs       kcorrect coefficients. Default to typical coefficients for a        galaxy. Array of 5 elements to give same coefficients to all        galaxies, or array of [5, n_elements(z)] to specify        coefficients for each galaxy  kc       k (+ e) corrections, without distance modulus  vmatrix       supply a variable name here to speed up multiple executions  lambda       supply a variable name here to speed up multiple executions  rmatrix       for repeated execution with same z, same band, but with        different coeffs, supply a variable name here  zvals       for repeated execution with same z, same band, but with        different coeffs, supply a variable name here  z", "          -1", "25 Jan 2008: Added evolution corrections       5 Mar 2008: Added coeffs, kc, vmatrix, lambda, rmatrix, zvals keywords  DM(z) + K(z) - E(z), where E(z) = q0 * (1 + q1 * z) * z    ");
  
  libdata[libdataItem++] = new Array("gals/ajs_dm_k.html#ajs_dm_k_test", "ajs_dm_k_test", "routine in ajs_dm_k.pro", "ajs_dm_k.pro", "", "ajs_dm_k_test", " Test ajs_dm_k  ", "", "          -1", "");
  
  

libdata[libdataItem++] = new Array("gals/ajs_double_choloniewski.html", "ajs_double_choloniewski.pro", ".pro file in gals/ directory", "ajs_double_choloniewski.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_double_choloniewski.html#ajs_double_choloniewski", "ajs_double_choloniewski", "routine in ajs_double_choloniewski.pro", "ajs_double_choloniewski.pro", "", "ajs_double_choloniewski", "	This function generates a double Choloniewski function along an input 	array of absolute magnitudes and surface brightnesses.            	Compatible with mpfit2d.  	Magnitudes only       ", "absmagsbparams[mstar,alpha,phistar,sbstar,sigmasb,beta] twice   ", "          -1", " 	This function returns an rray of phi for corresponding absmag, sbr    ");
  
  

libdata[libdataItem++] = new Array("gals/ajs_double_choloniewski_fit.html", "ajs_double_choloniewski_fit.pro", ".pro file in gals/ directory", "ajs_double_choloniewski_fit.pro", "", "", "       (NB: unlikely to be useful!) This function finds the        best-fitting double Choloniewski function for input arrays of        absolute magnitudes, surface brightness, phi and phi_err.  ", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_double_choloniewski_fit.html#ajs_double_choloniewski_fit", "ajs_double_choloniewski_fit", "routine in ajs_double_choloniewski_fit.pro", "ajs_double_choloniewski_fit.pro", "", "ajs_double_choloniewski_fit", "	This function finds the best-fitting double Choloniewski function for 	input arrays of absolute magnitudes, surface brightness, phi 	and phi_err.            	Uses mpfit2d (idlutils) 	ajs_choloniewski       ", "txtfileabsmagsbphiphierrparams_init", "          -1", " 	Result	[mstar,alpha,phistar,sbstar,sigmasb,beta]    ");
  
  libdata[libdataItem++] = new Array("gals/ajs_double_choloniewski_fit.html#ajs_choloniewski_fit_test", "ajs_choloniewski_fit_test", "routine in ajs_double_choloniewski_fit.pro", "ajs_double_choloniewski_fit.pro", "", "ajs_choloniewski_fit_test", "", "", "          -1", "");
  
  

libdata[libdataItem++] = new Array("misc/ajs_double_gauss.html", "ajs_double_gauss.pro", ".pro file in misc/ directory", "ajs_double_gauss.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_double_gauss.html#ajs_double_gauss", "ajs_double_gauss", "routine in ajs_double_gauss.pro", "ajs_double_gauss.pro", "", "ajs_double_gauss", " This function returns the sum of two Gaussians, using gauss1  (mpfit). Suitable for curve fitting.    ", "_REF_EXTRA       Extra parameters to be passed to gauss1  x       Array of X values  p       [mean1, sigma1, area1, mean2, sigma2, area2] parameters for        the two Gaussian components  ", "          -1", "11 Jan 2008	Created, Anthony Smith  ");
  
  

libdata[libdataItem++] = new Array("gals/ajs_double_schechter.html", "ajs_double_schechter.pro", ".pro file in gals/ directory", "ajs_double_schechter.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_double_schechter.html#ajs_double_schechter", "ajs_double_schechter", "routine in ajs_double_schechter.pro", "ajs_double_schechter.pro", "", "ajs_double_schechter", " This function generates a double Schechter function (single value of  M-star, but two values for alpha and phi-star) along an input array of  absolute magnitudes.    Compatible with mpfit    ", "mstaralpha_1phistar_1alpha_2phistar_2loglum       Set to /loglum to input log10 luminosity instead of absolute magnitudes  absmag       Array of absolute magnitudes  params       [mstar, alpha_1, phistar_1, alpha_2, phistar_2] OR input separately:  ", "          -1", "      21 April 2008 Written, Anthony Smith        Array of phi for corresponding absmag  ");
  
  

libdata[libdataItem++] = new Array("gals/ajs_double_schechter_fit.html", "ajs_double_schechter_fit.pro", ".pro file in gals/ directory", "ajs_double_schechter_fit.pro", "", "", " This function finds the best-fitting double Schechter function for a  given array of absolute magnitudes, phi and phi_err (alpha_2 > alpha_1)  ", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_double_schechter_fit.html#ajs_double_schechter_fit", "ajs_double_schechter_fit", "routine in ajs_double_schechter_fit.pro", "ajs_double_schechter_fit.pro", "", "ajs_double_schechter_fit", " This function finds the best-fitting double Schechter function for a  given array of absolute magnitudes, phi and phi_err (alpha_2 > alpha_1)  ", "txtfile       Output text file for Schechter parameters  loglum       Set to input log10 luminosity instead of absolute magnitudes  cov_mat       Covariance matrix  range       [min, max] absmag for fitting  lum_dens       Luminosity density.  Return        value is j = ... x 10^(0.4 M_sun) h L_sun Mpc^-3  err_lum_dens       Error in luminosity density  _REF_EXTRA       Extra keywords for mpfitfun  absmagphiphierrparams_init       Initial guess for [mstar, alpha_1, phistar_1, alpha_2, phistar_2]        (or default)  ", "          -1", "      5 Sep 2007 Created, Anthony Smith       17 Jan 2008 Added loglum keyword and text file output, AJS       4 Apr 2008 Added range keyword       15 Apr 2008 Added lum_dens       17 Jul 2008 lum_dens now gives correct result for loglum        mpfit (idlutils)       ajs_schechter        [mstar, alpha_1, phistar_1, alpha_2, phistar_2]  ");
  
  libdata[libdataItem++] = new Array("gals/ajs_double_schechter_fit.html#ajs_double_schechter_fit_test", "ajs_double_schechter_fit_test", "routine in ajs_double_schechter_fit.pro", "ajs_double_schechter_fit.pro", "", "ajs_double_schechter_fit_test", " Test ajs_schechter_fit  ", "", "          -1", "");
  
  

libdata[libdataItem++] = new Array("misc/ajs_ellipse_plot.html", "ajs_ellipse_plot.pro", ".pro file in misc/ directory", "ajs_ellipse_plot.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_ellipse_plot.html#ajs_ellipse_plot", "ajs_ellipse_plot", "routine in ajs_ellipse_plot.pro", "ajs_ellipse_plot.pro", "", "ajs_ellipse_plot", " Plot an ellipse    Wrapper on plot or oplot  ", "ang       Positional angle (degrees anticlockwise)  rho       Correlation between x and y, with a and b now interpreted as        standard deviation in x and y, to give <n_sigma>-sigma contours          Equivalently, plot an ellipse within a rectangle bounded by xc        +/- a and yc +/- b, where the ellipse touches the rectangle at        +/- (xc + rho * a, yc + b) and +/- (xc + a, yc + rho * b)  overplot       Set /overplot to overplot (uses oplot rather than plot)  n_sigma       With rho, plot <n_sigma>-sigma contour(s) (equivalent to        multiplier on a and b). Can be an array of values  plotfile       See ajs_plot_start  show_plot       See ajs_plot_start  open       See ajs_plot_start  close       See ajs_plot_start  _REF_EXTRA       Keywords for plot/oplot  a       Semimajor/minor x-axis  b       Semimajor/minor y-axis  xc       x-co-ordinate of centre of ellipse  yc       y-co-ordinate of centre of ellipse  ", "          -1", "      28 Mar 2008 Written, Anthony Smith        ajs_ellipse_plot, 1, 1       ajs_ellipse_plot, 2, 3, 10, 20       ajs_ellipse_plot, 4, 1, ang=30       ajs_ellipse_plot, 1, 1, rho=0.8, /overplot       ajs_ellipse_plot, 3, 1, rho=0.5, n_sigma=2       ajs_ellipse_plot, 3, 1, rho=-0.5, n_sigma=[1,3,2]  ");
  
  

libdata[libdataItem++] = new Array("misc/ajs_gaussnd.html", "ajs_gaussnd.pro", ".pro file in misc/ directory", "ajs_gaussnd.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_gaussnd.html#ajs_gaussnd", "ajs_gaussnd", "routine in ajs_gaussnd.pro", "ajs_gaussnd.pro", "", "ajs_gaussnd", " This function returns the probability density at a given point  _x_ for a multivariate normal distribution with specified  mean, covariance matrix and amplitude.    Compatible with mpfit  ", "amp       amplitude (optional; default=1)  mean       [n] Mean of each dimension  cov       [n,n] Covariance matrix  x       [n,N] Array of n-dimensional input values  params       [1 + n + (n*(n+1))/2] = [amp,mean,cov]   	        WHERE	amp	amplitude (scalar)   		mean	[n] Mean of each dimension   		cov	[(n*(n+1))/2] Unique components of symmetric 			covariance matrix, e.g., 00, 01, 02, 11, 12, 22 for 3-d  ", "          -1", "      7 Sep 2007: Created, Anthony Smith        Identical result:       y=ajs_gaussnd([1,2,3],mean=[0,0,0],cov=[[1,0,0],[0,1,0],[0,0,1]])       y=ajs_gaussnd([1,2,3],amp=1,mean=[0,0,0],cov=[[1,0,0],[0,1,0],[0,0,1]])       y=ajs_gaussnd([1,2,3],[1,0,0,0,1,0,0,1,0,1])        [N] probability density of specified normal distribution     at _x_  ");
  
  

libdata[libdataItem++] = new Array("misc/ajs_idldoc.html", "ajs_idldoc.pro", ".pro file in misc/ directory", "ajs_idldoc.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_idldoc.html#ajs_idldoc", "ajs_idldoc", "routine in ajs_idldoc.pro", "ajs_idldoc.pro", "Anthony Smith  ", "ajs_idldoc", " This procedure generates HTML documentation for all routines  in Anthony Smith's IDL library, using IDLdoc    Creates .tar.gz of all .pro files contained within subdirectories only.    ; docformat = 'rst' (for example) should be specified at the  top of each file, or uses idldoc default.    ", "", "          -1", "28 Jan 2008 Created    ajs_idldoc    ");
  
  

libdata[libdataItem++] = new Array("misc/ajs_is_monotonic.html", "ajs_is_monotonic.pro", ".pro file in misc/ directory", "ajs_is_monotonic.pro", "", "", " Return 1 if the array is monotonic, 0 if not  ", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_is_monotonic.html#ajs_is_monotonic", "ajs_is_monotonic", "routine in ajs_is_monotonic.pro", "ajs_is_monotonic.pro", "", "ajs_is_monotonic", " Return 1 if the array is monotonic, 0 if not    ", "array", "          -1", "      12 Mar 2008 Created, Anthony Smith        print, ajs_is_monotonic([1,2,3])          1       print, ajs_is_monotonic([1,1,3])          0       print, ajs_is_monotonic([3,2,1])          1  ");
  
  libdata[libdataItem++] = new Array("misc/ajs_is_monotonic.html#ajs_is_monotonic_test", "ajs_is_monotonic_test", "routine in ajs_is_monotonic.pro", "ajs_is_monotonic.pro", "", "ajs_is_monotonic_test", " Test ajs_is_monotonic  ", "", "          -1", "");
  
  

libdata[libdataItem++] = new Array("misc/ajs_jackknife.html", "ajs_jackknife.pro", ".pro file in misc/ directory", "ajs_jackknife.pro", "", "", " This function returns the jackknife method-estimated standard  deviation (and optionally the jackknife mean) of n resamplings  of a quantity, where each of the n estimates of the quantity  is generated by removing (1/n) of the original sample data points.  ", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_jackknife.html#ajs_jackknife", "ajs_jackknife", "routine in ajs_jackknife.pro", "ajs_jackknife.pro", "", "ajs_jackknife", " This function returns the jackknife method-estimated standard  deviation (and optionally the jackknife mean) of n resamplings  of a quantity, where each of the n estimates of the quantity  is generated by removing (1/n) of the original sample data points.  ", "original_estimate       Supply this value to get bias correction (see below)  bias_correction       Bias-corrected estimate is x_c = x + b where x is the original        estimate and b is the bias_correction  jackknife_mean       Mean of the jackknife samples (NB different from        bias-corrected estimate)  jackknife_estimates       [n] or [m, n] array of n estimates of a quantity, or n        estimates of m quantities  ", "          -1", "      24 Jan 2008 Created, Anthony Smith       7 Apr 2008 Added covariance output for [m, n] input        Return value is the jackknife standard deviation or, if     jackknife_estimates is [m, n], returns jackknife covariance matrix [m, m]  ");
  
  libdata[libdataItem++] = new Array("misc/ajs_jackknife.html#ajs_jackknife_test", "ajs_jackknife_test", "routine in ajs_jackknife.pro", "ajs_jackknife.pro", "", "ajs_jackknife_test", " Test ajs_jackknife  ", "", "          -1", "");
  
  

libdata[libdataItem++] = new Array("misc/ajs_kw_string.html", "ajs_kw_string.pro", ".pro file in misc/ directory", "ajs_kw_string.pro", "", "", " This function returns a string specifying values of keywords  ", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_kw_string.html#ajs_kw_string_values", "ajs_kw_string_values", "routine in ajs_kw_string.pro", "ajs_kw_string.pro", "", "ajs_kw_string_values", " Return the string  ", "_EXTRAempty_stringnames", "          -1", "");
  
  libdata[libdataItem++] = new Array("misc/ajs_kw_string.html#ajs_kw_string", "ajs_kw_string", "routine in ajs_kw_string.pro", "ajs_kw_string.pro", "", "ajs_kw_string", " This function returns a string specifying values of keywords  ", "_REF_EXTRAempty_string       String to return for undefined keywords (default  empty )  ", "          -1", "      12 Mar 2008 Written, Anthony Smith        print, ajs_kw_string(a=2, b=b)       A=2 B=empty  ");
  
  

libdata[libdataItem++] = new Array("gals/ajs_lf_plot.html", "ajs_lf_plot.pro", ".pro file in gals/ directory", "ajs_lf_plot.pro", "", "", " This procedure produces a plot of the luminosity function.  ", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_lf_plot.html#ajs_lf_plot_xrange", "ajs_lf_plot_xrange", "routine in ajs_lf_plot.pro", "ajs_lf_plot.pro", "", "ajs_lf_plot_xrange", " Choose and set xrange  ", "bincentresloglumschechterdouble_schechter", "          -1", "");
  
  libdata[libdataItem++] = new Array("gals/ajs_lf_plot.html#ajs_lf_plot_yrange", "ajs_lf_plot_yrange", "routine in ajs_lf_plot.pro", "ajs_lf_plot.pro", "", "ajs_lf_plot_yrange", " Choose and set yrange  ", "phiinput_logphiplot_logphischechterdouble_schechter", "          -1", "");
  
  libdata[libdataItem++] = new Array("gals/ajs_lf_plot.html#ajs_lf_plot_legend", "ajs_lf_plot_legend", "routine in ajs_lf_plot.pro", "ajs_lf_plot.pro", "", "ajs_lf_plot_legend", " Add legend text to the plot  ", "legend_textlegend_posarg_legend_posschechterdouble_schechtercolorpsym_REF_EXTRA", "          -1", "");
  
  libdata[libdataItem++] = new Array("gals/ajs_lf_plot.html#ajs_lf_plot", "ajs_lf_plot", "routine in ajs_lf_plot.pro", "ajs_lf_plot.pro", "", "ajs_lf_plot", " This procedure produces a plot of the luminosity function.    Data may be entered in four forms:    (1) Values at particular absolute magnitudes (bincentres, phi[, err_phi])    (2) (double) Schechter function parameters (always plotted if supplied)    (3) Array of galaxy absolute magnitudes and weights (e.g., 1/Vmax) -  ignored if (1) set    (4) From a file, to be read using ajs_lf_read - ignored if (1) or (3) set    In addition, a (double) Schechter function fit can be performed.    ", "bincentresphierr_phineg_err_phi       Negative error, for asymmetric  pos_err_phi       Positive error, for asymmetric  lower_limit_where       Indices of values of phi where the upper error bar is a lower        limit (i.e., plot arrowhead and empty circle, if psym=16 set)  schechter       Plot Schechter function. Set to undefined variable to        calculate best-fitting Schechter function and return the        parameters. Or set /schechter to calculate and print the        parameters  sch_range       Range of absolute magnitudes for (double) Schechter fit        (ajs_schechter_fit)  double_schechter       Plot double Schechter function. Set to undefined variable to        calculate best-fitting double Schechter function and return the        parameters. Or set /double_schechter to calculate and print the        parameters  double_cov_mat       Covariance matrix of double Schechter function fit  abs_mag       Array of absolute magnitudes (use with weights)  weights       Array of weights (e.g., 1/Vmax) for estimating luminosity        function  jackknife       Used with absmag, this is an integer for each galaxy        giving the jackknife region in which the galaxy lies  nbins       Number of bins for abs_mag & weights (use with xrange or bincentres)  datfile       Input data: binary format file created by ajs_lf_write  overplot       Set /overplot to overplot  plotfile       EPS file name  show_plot       Set /show_plot to open the EPS file for viewing  open       Set /open to open the EPS file and leave it for overplotting  close       Set /close to close the EPS plotfile  xrangeyrangeinput_logphi       Set /input_logphi to indicate data values are log(phi) rather        than phi  plot_logphi       Set /plot_logphi to plot log(phi) rather than phi (redundant        if overplotting)  loglum       Set /loglum for log luminosities rather than magnitudes  ylogpsym       Uses D. Fanning's symcat (e.g., 16 = filled circle) for        non-standard psym values  symcat_thick       Thick keyword for symcat  band_name       For the x-axis: M_<band_name> - 5 log h  colorxstyle       Default /xstyle  ystyle       Default /ystyle  title       Title for plot  xtitle       Title for x-axis  ytitle       Title for y-axis  legend_pos       [xpos, ypos] for legend text. Output is position for next        legend label, for overplotting. varies from [0, 0] = [left,        bottom] to [1, 1] = [right, top].  legend_text       Legend text  cov_mat       Covariance matrix of Schechter function fit  _REF_EXTRA       Extra keywords for plots    ", "          -1", "      13 Feb 2008 Created, Anthony Smith       5 Mar 2008 err_phi no longer required       6 Mar 2008 Added abs_mag, weights and nbins       18 Mar 2008 Added jackknife       25 Mar 2008 Added legend_text and legend_pos       4 Apr 2008 Added sch_range (ajs_schechter_fit range)       7 Apr 2008 Schechter function info included in call to ajs_vmax_lf       21 Apr 2008 double Schechter function input       23 Apr 2008 psym now uses D. Fanning's symcat routine;     added lower_limit_where keyword       1 May 2008 Uses ajs_plot_pos for legend position        bincentres = [-20, -22, -24]       phi = [1e-2, 1e-2, 1e-3]       err_phi = [1e-2, 1e-3, 1e-4]       ajs_lf_plot, bincentres=bincentres, phi=phi       ajs_lf_plot, bincentres=bincentres, phi=phi, err_phi=err_phi       ajs_lf_plot, bincentres=bincentres, phi=phi, err_phi=err_phi, psym=16       ajs_lf_plot, schechter=[-24,-1,1e-2]       ajs_lf_plot, schechter=[-24,-1,1e-2], color=2, bincentres=bincentres, phi=phi, err_phi=err_phi, xrange=[-16,-26], psym=16, plotfile='~/lf.eps', /show_plot    ");
  
  

libdata[libdataItem++] = new Array("gals/ajs_lf_read.html", "ajs_lf_read.pro", ".pro file in gals/ directory", "ajs_lf_read.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_lf_read.html#ajs_lf_read", "ajs_lf_read", "routine in ajs_lf_read.pro", "ajs_lf_read.pro", "", "ajs_lf_read", " This procedure reads LF (luminosity function) data from a binary  file creaed using ajs_lf_write.    ", "datfilebincentresphilfphilferr", "          -1", "11 Jan 2008 Created, Anthony Smith  ");
  
  

libdata[libdataItem++] = new Array("gals/ajs_lf_write.html", "ajs_lf_write.pro", ".pro file in gals/ directory", "ajs_lf_write.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_lf_write.html#ajs_lf_write", "ajs_lf_write", "routine in ajs_lf_write.pro", "ajs_lf_write.pro", "", "ajs_lf_write", " This procedure takes arrays describing the LF (luminosity function)  and writes the data either to a binary file or to a text file.    ", "datfiletxtfilebincentresphilfphilferr", "          -1", "11 Jan 2008 Created, Anthony Smith  ");
  
  

libdata[libdataItem++] = new Array("misc/ajs_linspace.html", "ajs_linspace.pro", ".pro file in misc/ directory", "ajs_linspace.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_linspace.html#ajs_linspace", "ajs_linspace", "routine in ajs_linspace.pro", "ajs_linspace.pro", "", "ajs_linspace", "	This function returns an array of points equally (linearly) 	spaced between two extremes (as Python's numpy.linspace or Matlab).                    ", "stepbincentresstartstopnum", "          -1", "	Result = ajs_linspace(10, 20) 	Returns array of floating point numbers, or double-precision if 	either start or stop is double precision    ");
  
  

libdata[libdataItem++] = new Array("misc/ajs_logspace.html", "ajs_logspace.pro", ".pro file in misc/ directory", "ajs_logspace.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_logspace.html#ajs_logspace", "ajs_logspace", "routine in ajs_logspace.pro", "ajs_logspace.pro", "", "ajs_logspace", "	This function returns an array of points equally (logarithmically) 	spaced between two extremes (as Python's numpy.logspace or Matlab).                    ", "basestepstartstopnum", "          -1", "	Result = ajs_logspace(-10, 15) 	Returns array of floating point numbers, or double-precision if 	either start or stop is double precision    ");
  
  

libdata[libdataItem++] = new Array("gals/ajs_lookback_time.html", "ajs_lookback_time.pro", ".pro file in gals/ directory", "ajs_lookback_time.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_lookback_time.html#ajs_lookback_time", "ajs_lookback_time", "routine in ajs_lookback_time.pro", "ajs_lookback_time.pro", "", "ajs_lookback_time", " Return lookback time to z, in years    See ajs_comdis for keyword parameters  ", "h0       Hubble constant (default 70)  _REF_EXTRAz       Redshift (or array of values)  ", "          -1", "      10 Oct 2008 Written, Anthony Smith  ");
  
  

libdata[libdataItem++] = new Array("gals/ajs_lumdis.html", "ajs_lumdis.pro", ".pro file in gals/ directory", "ajs_lumdis.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_lumdis.html#ajs_lumdis", "ajs_lumdis", "routine in ajs_lumdis.pro", "ajs_lumdis.pro", "", "ajs_lumdis", " Return Luminosity distance to z    See ajs_comdis for keyword parameters  ", "_REF_EXTRAz       Redshift (or array of values)  ", "          -1", "");
  
  

libdata[libdataItem++] = new Array("gals/ajs_number_counts.html", "ajs_number_counts.pro", ".pro file in gals/ directory", "ajs_number_counts.pro", "", "", " This procedure calculates the number counts from an input array of  magnitudes or log(flux)  ", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_number_counts.html#ajs_number_counts", "ajs_number_counts", "routine in ajs_number_counts.pro", "ajs_number_counts.pro", "", "ajs_number_counts", " This procedure calculates the number counts from an input array of  magnitudes or log(flux)    Bins must have equal width    ", "area       Area in square degrees  bincentres       Centre of each bin in magnitude  err_ngals       Poisson errors = n / sqrt(n)  ninbin       Number of galaxies in each bin  nbins       Number of bins (ignored if bincentres set)  xrange       Minimum and maximum mag (ignored if bincentres set)  mag", "          -1", "      5 Sep 2007 Created (Anthony Smith)       6 Mar 2008 Re-written       18 Mar 2008 Added jackknife       7 Apr 2008 Jackknife estimation of Schechter function fit       15 Apr 2008 Added luminosity density       21 Apr 2008 Added double Schechter function  fltarr/dblarr     Returns the number of galaxies per mag or log(flux), per square degree  ");
  
  libdata[libdataItem++] = new Array("gals/ajs_number_counts.html#ajs_number_counts_test", "ajs_number_counts_test", "routine in ajs_number_counts.pro", "ajs_number_counts.pro", "", "ajs_number_counts_test", " Test ajs_number_counts  ", "", "          -1", "");
  
  

libdata[libdataItem++] = new Array("gals/ajs_number_counts_plot.html", "ajs_number_counts_plot.pro", ".pro file in gals/ directory", "ajs_number_counts_plot.pro", "", "", " This procedure produces a plot of the number counts.  ", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_number_counts_plot.html#ajs_number_counts_plot_xrange", "ajs_number_counts_plot_xrange", "routine in ajs_number_counts_plot.pro", "ajs_number_counts_plot.pro", "", "ajs_number_counts_plot_xrange", " Choose and set xrange  ", "bincentreslogflux", "          -1", "");
  
  libdata[libdataItem++] = new Array("gals/ajs_number_counts_plot.html#ajs_number_counts_plot_yrange", "ajs_number_counts_plot_yrange", "routine in ajs_number_counts_plot.pro", "ajs_number_counts_plot.pro", "", "ajs_number_counts_plot_yrange", " Choose and set yrange  ", "ngalsinput_logngalsplot_logngalserr_ngalspos_err_ngalsneg_err_ngals", "          -1", "");
  
  libdata[libdataItem++] = new Array("gals/ajs_number_counts_plot.html#ajs_number_counts_plot_legend", "ajs_number_counts_plot_legend", "routine in ajs_number_counts_plot.pro", "ajs_number_counts_plot.pro", "", "ajs_number_counts_plot_legend", " Add legend text to the plot  ", "legend_textlegend_posarg_legend_poscolorpsym_REF_EXTRA", "          -1", "");
  
  libdata[libdataItem++] = new Array("gals/ajs_number_counts_plot.html#ajs_number_counts_plot", "ajs_number_counts_plot", "routine in ajs_number_counts_plot.pro", "ajs_number_counts_plot.pro", "", "ajs_number_counts_plot", " This procedure produces a plot of the number counts.    ", "bincentres       x-position (e.g., magnitude, flux, log-flux)  ngals       Number density, per square degree, per unit (x_unit)  err_ngals       Uncertainty in ngals  neg_err_ngals       Negative error in ngals, for asymmetric  pos_err_ngals       Positive error in ngals, for asymmetric  mag       Array of magnitudes (or log(flux) with /logflux)  area       Area in square degrees (for use with mag)  nbins       Number of bins for mag (use with xrange or bincentres)  overplot       Set /overplot to overplot  plotfile       EPS file name  show_plot       Set /show_plot to open the EPS file for viewing  open       Set /open to open the EPS file and leave it for overplotting  close       Set /close to close the EPS plotfile  xrangeyrangeinput_logngals       Set /input_logngals to indicate data values are log(n) rather        than n  plot_logngals       Set /plot_logngals to plot log(n) rather than n (redundant        if overplotting)  logflux       Set /logflux for log flux rather than magnitudes  ylogpsym       Uses D. Fanning's symcat (e.g., 16 = filled circle) for        non-standard psym values  symcat_thick       Thick keyword for symcat  colorxstyle       Default xstyle=1  ystyle       Default ystyle=1  title       Title for plot  xtitle       Title for x-axis  ytitle       Title for y-axis  legend_pos       [xpos, ypos] for legend text. Output is position for next        legend label, for overplotting. varies from [0, 0] = [left,        bottom] to [1, 1] = [right, top].  legend_text       Legend text  euclid       Set /euclid to subtract Euclidean slope, or euclid=x to then        divide by x (specify axis ranges and titles explicitly)  _REF_EXTRA       Extra keywords for plots    ", "          -1", "      29 Oct 2008 Created, Anthony Smith       16 Jan 2009 Added euclid keyword, AJS        bincentres = [10, 11, 12]       ngals = [10, 100, 1000]       err_ngals = [1, 10, 100]       ajs_number_counts_plot, bincentres=bincentres, ngals=ngals,     err_ngals=err_ngals    ");
  
  

libdata[libdataItem++] = new Array("misc/ajs_open_file.html", "ajs_open_file.pro", ".pro file in misc/ directory", "ajs_open_file.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_open_file.html#ajs_open_file", "ajs_open_file", "routine in ajs_open_file.pro", "ajs_open_file.pro", "", "ajs_open_file", " Open a file using the default system application.    Works on Macs only!    ", "file       Path to the file to be opened.  ", "          -1", "");
  
  

libdata[libdataItem++] = new Array("misc/ajs_pause.html", "ajs_pause.pro", ".pro file in misc/ directory", "ajs_pause.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_pause.html#ajs_pause", "ajs_pause", "routine in ajs_pause.pro", "ajs_pause.pro", "", "ajs_pause", " This procedure invites the user to press return before continuing.  ", "", "          -1", "");
  
  

libdata[libdataItem++] = new Array("misc/ajs_pixel_plot.html", "ajs_pixel_plot.pro", ".pro file in misc/ directory", "ajs_pixel_plot.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_pixel_plot.html#ajs_pixel_plot", "ajs_pixel_plot", "routine in ajs_pixel_plot.pro", "ajs_pixel_plot.pro", "", "ajs_pixel_plot", " Plot (large number of) points as pixels (rectangles), with colour  (white->black) representing density    ", "binsize1       Size of bin in x-direction  binsize2       Size of bin in y-direction  nbins1       Number of bins in x-direction (ignored if xrange and binsize1        set; default 50)  nbins2       Number of bins in y-direction (ignored if yrange and binsize2        set; default 50)  xrange       Range for plot (if not overplot) and range for histogram  yrange       Range for plot (if not overplot) and range for histogram  overplot       Set /overplot to overplot  ct       Colour table: -1 = use current, -2 = fixed colour (color keyword),        -3 (or not specified) = default (white -> black)        If color keyword set, and ct != -2, generate colour table from color  color       Colour to use instead of grayscale/black  _REF_EXTRA       Extra keywords for plot  x       x-coordinate of points  y       y-coordinate of points  ", "          -1", "      18 Apr 2008 Written, Anthony Smith        ajs_plot_radec, ra, dec, /nodata, xtitle='RA / degrees',     ytitle='dec / degrees'       ajs_pixel_plot, ra, sin(dec / !RADEG), nbins1=200, nbins2=200,     /overplot, ct=40       Sample output (click to enlarge):       <a href= http://astronomy.sussex.ac.uk/~anthonys/images/sdss_radec.png ><img src= http://astronomy.sussex.ac.uk/~anthonys/images/sdss_radec.png  width= 200 ></a>  ");
  
  

libdata[libdataItem++] = new Array("misc/ajs_plot_pos.html", "ajs_plot_pos.pro", ".pro file in misc/ directory", "ajs_plot_pos.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_plot_pos.html#ajs_plot_pos", "ajs_plot_pos", "routine in ajs_plot_pos.pro", "ajs_plot_pos.pro", "", "ajs_plot_pos", " Return the location on the plot in data co-ordinates given input  value(s) between 0 and 1.    Like using xyouts or plots with /normal but between the axis limits  rather than the whole extent of the plot area  ", "x_in       x-position (as x parameter; ignored if x present)  y_in       y-position (as y parameter; ignored if y present)  xrange       Overrides !x.crange  yrange       Overrides !y.crange  xlog       Use with xrange to specify logarithmic axis. If xrange is not        specified, reads from !x.type  ylog       Use with yrange to specify logarithmic axis. If yrange is not        specified, reads from !y.type  x       x-position, between 0 (left-hand side) and 1 (right-hand        side). Or [N, 2] values for [x, y]  y       y-position, between 0 (bottom) and 1 (top). Ignored if x is        [N, 2]  ", "          -1", "      30 Apr 2008 Written, Anthony Smith        plot, [0], xrange=[-10, 20], yrange=[10, 15]       print, ajs_plot_pos(0.5)              5.0000000       print, ajs_plot_pos(0.5, 0.8)              5.0000000       14.000000       print, ajs_plot_pos(y_in=0.8)              14.000000       print, ajs_plot_pos([0.3, 0.5, 0.7], [0.2, 0.4, 0.6])            -0.99999964       5.0000000       11.000000              11.000000       12.000000       13.000000       print, ajs_plot_pos(y_in=0.8, yrange=[0.1, 0.2])            0.180000       print, ajs_plot_pos(y_in=0.5, yrange=[1, 100], /ylog)             10.0000        Position in data units corresponding to input x, y, [x, y], [x],     [y] or [[x], [y]]  ");
  
  

libdata[libdataItem++] = new Array("gals/ajs_plot_radec.html", "ajs_plot_radec.pro", ".pro file in gals/ directory", "ajs_plot_radec.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_plot_radec.html#ajs_plot_radec", "ajs_plot_radec", "routine in ajs_plot_radec.pro", "ajs_plot_radec.pro", "", "ajs_plot_radec", " Plot ra and dec on equal-area plot    Wrapper on plot/oplot, (1) plots sin(y) instead of y, but labelling  the axes to show y, (2) chooses appropriate ranges in multiples of  15s, 30s, 45s or 60s (3) reverses the x-axis, unless xrange is  specified    /overplot is equivalent to oplot, x, sin(y / !RADEG)    ", "overplot       Set /overplot to overplot  xrange       As for plot (reversed if not specified)  yrange       As for plot  xstyle       As for plot  ystyle       As for plot  _REF_EXTRA       Extra keywords for [o]plot  x       ra, degrees  y       dec, degrees  ", "          -1", "      17 Apr 2008 Written, Anthony Smith        See ajs_pixel_plot for an example  ");
  
  

libdata[libdataItem++] = new Array("misc/ajs_plot_start.html", "ajs_plot_start.pro", ".pro file in misc/ directory", "ajs_plot_start.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_plot_start.html#ajs_plot_start", "ajs_plot_start", "routine in ajs_plot_start.pro", "ajs_plot_start.pro", "", "ajs_plot_start", " Start plot (device/EPS settings, if required)    Double thickness for all lines  ", "plotfile       Name of EPS file. If EPS file required but not specified,        set /eps and it defaults to '/tmp/idl_plot.eps'  eps       Set /eps for EPS file  _REF_EXTRA       Keywords for device  ", "          -1", "      1 Apr 2008 Written, Anthony Smith       17 Apr 2008 Replaced open and show_plot keywords with eps  ");
  
  

libdata[libdataItem++] = new Array("misc/ajs_plot_stop.html", "ajs_plot_stop.pro", ".pro file in misc/ directory", "ajs_plot_stop.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_plot_stop.html#ajs_plot_stop", "ajs_plot_stop", "routine in ajs_plot_stop.pro", "ajs_plot_stop.pro", "", "ajs_plot_stop", " Stop plot (device/EPS settings, if required)    ", "plotfile       Name of EPS file. If EPS file required but not specified,        defaults to '/tmp/idl_plot.eps'  show_plot       Close EPS file (for writing) and display the file (for viewing)  open       Leave the plot open, unless /close is set  close       Set /close to close an already open EPS file.  ", "          -1", "");
  
  

libdata[libdataItem++] = new Array("misc/ajs_random_sort.html", "ajs_random_sort.pro", ".pro file in misc/ directory", "ajs_random_sort.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_random_sort.html#ajs_random_sort", "ajs_random_sort", "routine in ajs_random_sort.pro", "ajs_random_sort.pro", "", "ajs_random_sort", " Return random indices, or sort an array into random order    ajs_random_sort(n) returns a random permutation of the integers 0 to  (n - 1)    ajs_random_sort(n, min=min, max=max) returns n (different) random  integers between min and (max - 1)    ajs_random_sort(array) returns the elements of array in a random order    ", "min       Used with a single integer input, this is the minimum array        index that may be returned.  max       Used with a single integer input, this number ** - 1 ** is the        maximum array index that may be returned.  seed       Seed for randomu  input       Either (1) A single integer giving the number of array indices        to return or (2) an input array to be sorted into random order  ", "          -1", "      22 Apr 2008 Written, Anthony Smith        print, ajs_random_sort([1, 2, 3])              1           2           3       print, ajs_random_sort(4)              3           1           0           2       print, ajs_random_sort(4, min=10, max=15)             12          10          13          14  ");
  
  

libdata[libdataItem++] = new Array("gals/ajs_schechter.html", "ajs_schechter.pro", ".pro file in gals/ directory", "ajs_schechter.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_schechter.html#ajs_schechter", "ajs_schechter", "routine in ajs_schechter.pro", "ajs_schechter.pro", "", "ajs_schechter", " This function generates a Schechter function along an input array of  absolute magnitudes.    Compatible with mpfit    ", "mstaralphaphistarloglum       Set to /loglum to input log10 luminosity instead of absolute magnitudes  absmag       Array of absolute magnitudes  params       [mstar,alpha,phistar] OR input separately:  ", "          -1", "      5 Sep 2007 Created, Anthony Smith       17 Jan 2008 Added loglum keyword, AJS        Array of phi for corresponding absmag  ");
  
  

libdata[libdataItem++] = new Array("gals/ajs_schechter_fit.html", "ajs_schechter_fit.pro", ".pro file in gals/ directory", "ajs_schechter_fit.pro", "", "", " This function finds the best-fitting Schechter function for a  given array of absolute magnitudes, phi and phi_err  ", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_schechter_fit.html#ajs_schechter_fit", "ajs_schechter_fit", "routine in ajs_schechter_fit.pro", "ajs_schechter_fit.pro", "", "ajs_schechter_fit", " This function finds the best-fitting Schechter function for a  given array of absolute magnitudes, phi and phi_err  ", "txtfile       Output text file for Schechter parameters  loglum       Set to input log10 luminosity instead of absolute magnitudes  cov_mat       Covariance matrix  range       [min, max] absmag for fitting  lum_dens       Luminosity density.  Return        value is j = ... x 10^(0.4 M_sun) h L_sun Mpc^-3  err_lum_dens       Error in luminosity density  _REF_EXTRA       Extra keywords for mpfitfun  absmagphiphierrparams_init       Initial guess for [mstar,alpha,phistar] (or default)  ", "          -1", "      5 Sep 2007 Created, Anthony Smith       17 Jan 2008 Added loglum keyword and text file output, AJS       4 Apr 2008 Added range keyword       15 Apr 2008 Added lum_dens        mpfit (idlutils)       ajs_schechter        [mstar, alpha, phistar]  ");
  
  libdata[libdataItem++] = new Array("gals/ajs_schechter_fit.html#ajs_schechter_fit_test", "ajs_schechter_fit_test", "routine in ajs_schechter_fit.pro", "ajs_schechter_fit.pro", "", "ajs_schechter_fit_test", " Test ajs_schechter_fit  ", "", "          -1", "");
  
  

libdata[libdataItem++] = new Array("gals/ajs_vmax.html", "ajs_vmax.pro", ".pro file in gals/ directory", "ajs_vmax.pro", "", "", " Calculate Vmax given galaxy data and limits  ", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_vmax.html#ajs_vmax_zminmax_mag", "ajs_vmax_zminmax_mag", "routine in ajs_vmax.pro", "ajs_vmax.pro", "", "ajs_vmax_zminmax_mag", " Return zminmax for magnitudes  ", "zdm_kvartable_row", "          -1", "");
  
  libdata[libdataItem++] = new Array("gals/ajs_vmax.html#ajs_vmax_zminmax_sb", "ajs_vmax_zminmax_sb", "routine in ajs_vmax.pro", "ajs_vmax.pro", "", "ajs_vmax_zminmax_sb", " Return zminmax for surface brightness  ", "zkcvartable_row", "          -1", "");
  
  libdata[libdataItem++] = new Array("gals/ajs_vmax.html#ajs_vmax_zminmax_rad", "ajs_vmax_zminmax_rad", "routine in ajs_vmax.pro", "ajs_vmax.pro", "", "ajs_vmax_zminmax_rad", " Return zminmax for radius  ", "zd_avartable_row", "          -1", "");
  
  libdata[libdataItem++] = new Array("gals/ajs_vmax.html#ajs_vmax_zminmax", "ajs_vmax_zminmax", "routine in ajs_vmax.pro", "ajs_vmax.pro", "", "ajs_vmax_zminmax", " Return the minimum and maximum redshifts at which the galaxy would  be visible  ", "zdm_kkcd_avartable_row", "          -1", "");
  
  libdata[libdataItem++] = new Array("gals/ajs_vmax.html#ajs_vmax", "ajs_vmax", "routine in ajs_vmax.pro", "ajs_vmax.pro", "", "ajs_vmax", " Calculate Vmax given galaxy data and limits  ", "zmin       Minimum redshift. Default = min(table.z) or 0  zmax       Maximum redshift. Default = max(table.z) or 2  area       Area of survey in square degrees. Default whole sky.  q0q1z_min_limited_by       Number of galaxies with minimum redshift limited by each of        the variables, or by zmin  z_max_limited_by       Number of galaxies with maximum redshift limited by each of        the variables, or by zmax  v       Volume, as in V/Vmax, for each galaxy  table       Array of structures with fields {coeffs:, z:} and columns for absolute        values given in vars. Redshift z is optional: raises warning        if outside observable limits.  vars       Array of structures with fields {abs_name:} giving the column        name in table corresponding to the absolute value, {band:} giving the        name of the band, {type:} giving 'M'agnitudes, 'S'urface brightness        or 'R'adius and {app_min:, app_max:} giving the limits.  ", "          -1", "      7 Mar 2008 Created, Anthony Smith       25 Jul 2008 Added v keyword (V/Vmax), AJS       8 Aug 2008 Doesn't attempt V/Vmax when table.z does not exist        vmax: in h^-3 Mpc^3  ");
  
  

libdata[libdataItem++] = new Array("gals/ajs_vmax_bbd.html", "ajs_vmax_bbd.pro", ".pro file in gals/ directory", "ajs_vmax_bbd.pro", "", "", " This procedure calculates the 1/Vmax bivariate brightness  distribution from input arrays of two quantities (e.g., absolute  magnitude and surface brightness) and weights (e.g., 1/Vmax)  ", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_vmax_bbd.html#ajs_vmax_bbd", "ajs_vmax_bbd", "routine in ajs_vmax_bbd.pro", "ajs_vmax_bbd.pro", "", "ajs_vmax_bbd", " This procedure calculates the 1/Vmax bivariate brightness  distribution from input arrays of two quantities (e.g., absolute  magnitude and surface brightness) and weights (e.g., 1/Vmax)    Bins in each dimension must have equal width    ", "bincentres1       Centre of each bin in m1  bincentres2       Centre of each bin in m2  err_phi       Poisson errors = phi / sqrt(n)  ninbin       Number of galaxies in each bin  nbins1       Number of bins in m1 (ignored if bincentres set)  nbins2       Number of bins in m2 (ignored if bincentres set)  xrange       Minimum and maximum m1 (ignored if bincentres1 set)  yrange       Minimum and maximum m2 (ignored if bincentres2 set)  jackknife       Integer for each galaxy giving the jackknife region in which        the galaxy lies  choloniewski       Choloniewski fit  chol_xrange       xrange for Choloniewski fit (ajs_choloniewski_fit)  chol_yrange       yrange for Choloniewski fit (ajs_choloniewski_fit)  cov_mat       Covariance matrix of Choloniewski fit  vars       See ajs_choloniewski_fit  zmin       See ajs_choloniewski_fit  zmax       See ajs_choloniewski_fit  area       See ajs_choloniewski_fit  q0       See ajs_choloniewski_fit  q1       See ajs_choloniewski_fit  m1m2weights       E.g., 1./Vmax. Set to single value to give all galaxies the        same weight, or leave empty to give all galaxies a weight of 1.  ", "          -1", "      6 Mar 2008 Created (Anthony Smith)       18 Mar 2008 Added jackknife       7 Apr 2008 Added jackknife Choloniewski fit  fltarr/dblarr     Returns phi, value of space density  ");
  
  libdata[libdataItem++] = new Array("gals/ajs_vmax_bbd.html#ajs_vmax_bbd_test", "ajs_vmax_bbd_test", "routine in ajs_vmax_bbd.pro", "ajs_vmax_bbd.pro", "", "ajs_vmax_bbd_test", " Test ajs_vmax_lf  ", "", "          -1", "");
  
  

libdata[libdataItem++] = new Array("gals/ajs_vmax_lf.html", "ajs_vmax_lf.pro", ".pro file in gals/ directory", "ajs_vmax_lf.pro", "", "", " This procedure calculates the 1/Vmax luminosity function from input  arrays of absolute magnitude and weights (e.g., 1/Vmax)  ", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_vmax_lf.html#ajs_vmax_lf", "ajs_vmax_lf", "routine in ajs_vmax_lf.pro", "ajs_vmax_lf.pro", "", "ajs_vmax_lf", " This procedure calculates the 1/Vmax luminosity function from input  arrays of absolute magnitude and weights (e.g., 1/Vmax)    Bins must have equal width    ", "bincentres       Centre of each bin in absolute magnitude  err_phi       Poisson errors = phi / sqrt(n)  ninbin       Number of galaxies in each bin  nbins       Number of bins (ignored if bincentres set)  xrange       Minimum and maximum absmag (ignored if bincentres set)  jackknife       Integer for each galaxy giving the jackknife region in which        the galaxy lies  lum_dens       Luminosity density.  Return value is sum(lum * weight),        j = ... x 10^(0.4 M_sun) h L_sun Mpc^-3  err_lum_dens       Error on luminosity density, calculated from errors on phi, or        from jackknife resampling  loglum       For Schechter fits: /loglum indicates log(luminosity) inputs  schechter       Returns best-fitting Schechter function (using jackknife if present)  cov_mat       Returns covariance matrix of Schechter function fit (using        jackknife if present)  sch_range       Range of absolute magnitudes for (double) Schechter function fit  sch_lum_dens       As lum_dens, but for the Schechter function fit  sch_err_lum_dens       Error on Schechter luminosity density  double_schechter       Returns best-fitting double Schechter function (using        jackknife if present)  double_cov_mat       Returns covariance matrix of double Schechter function fit (using        jackknife if present)  double_sch_lum_dens       As lum_dens, but for the double Schechter function fit  double_sch_err_lum_dens       Error on double Schechter luminosity density  absmagweights       E.g., 1./Vmax. Set to a single value to give all galaxies that        weight, or leave blank to give all galaxies a weight of 1.  ", "          -1", "      5 Sep 2007 Created (Anthony Smith)       6 Mar 2008 Re-written       18 Mar 2008 Added jackknife       7 Apr 2008 Jackknife estimation of Schechter function fit       15 Apr 2008 Added luminosity density       21 Apr 2008 Added double Schechter function  fltarr/dblarr     Returns phi, value of luminosity function  ");
  
  libdata[libdataItem++] = new Array("gals/ajs_vmax_lf.html#ajs_vmax_lf_test", "ajs_vmax_lf_test", "routine in ajs_vmax_lf.pro", "ajs_vmax_lf.pro", "", "ajs_vmax_lf_test", " Test ajs_vmax_lf  ", "", "          -1", "");
  
  

libdata[libdataItem++] = new Array("gals/ajs_vmax_ngals_multiplot_4d.html", "ajs_vmax_ngals_multiplot_4d.pro", ".pro file in gals/ directory", "ajs_vmax_ngals_multiplot_4d.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_vmax_ngals_multiplot_4d.html#ajs_vmax_ngals_multiplot_4d", "ajs_vmax_ngals_multiplot_4d", "routine in ajs_vmax_ngals_multiplot_4d.pro", "ajs_vmax_ngals_multiplot_4d.pro", "", "ajs_vmax_ngals_multiplot_4d", "	This procedure provides a visualisation of a four-dimensional 	space density (number of galaxies per unit volume per unit A, B, C, 	D - absolute magnitudes, etc.). For a four-dimensional parameter 	space (luminosities, SBs, 	etc), fixing the value of two of these parameters, plot 	contours of constant Vmax along with the number of objects in 	each bin, for the other two parameters.  Do this for all 	combinations of parameters and all bins (six plots in all).   	Resets multiplot settings (multiplot,/default)   	Parameters 1, 2 & 4 reversed to give faint -> bright 	Parameter 3 not reversed (radius)   ", "phi_arrngalsallvmax_arrnbinsvmaxbincentresvmax1bincentresvmax2bincentresvmax3bincentresvmax4show_plotsplotfile_baseaxis_labelsvarminvarmaxnbinsbincentres1bincentres2bincentres3bincentres4", "          -1", "	Use calling sequence verbatim 	Six plots (EPS files)   ");
  
  

libdata[libdataItem++] = new Array("misc/ajs_widget_choice.html", "ajs_widget_choice.pro", ".pro file in misc/ directory", "ajs_widget_choice.pro", "", "", "	This function displays a widget asking for one or more options 	to be selected.  ", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("misc/ajs_widget_choice.html#ajs_widget_choice_event", "ajs_widget_choice_event", "routine in ajs_widget_choice.pro", "ajs_widget_choice.pro", "", "ajs_widget_choice_event", "", "ev", "          -1", "");
  
  libdata[libdataItem++] = new Array("misc/ajs_widget_choice.html#ajs_widget_choice", "ajs_widget_choice", "routine in ajs_widget_choice.pro", "ajs_widget_choice.pro", "", "ajs_widget_choice", "	This function displays a widget asking for one or more options 	to be selected.   ", "valuesArray of strings   ", "Widgets", "	result = AJS_WIDGET_CHOICE(['Hello', 'Boo']) 	This function returns the selected items from Values (array of 	strings).   ");
  
  

libdata[libdataItem++] = new Array("gals/ajs_zage.html", "ajs_zage.pro", ".pro file in gals/ directory", "ajs_zage.pro", "", "", "", "", "          -1", "");
  
  
  libdata[libdataItem++] = new Array("gals/ajs_zage.html#ajs_zage", "ajs_zage", "routine in ajs_zage.pro", "ajs_zage.pro", "", "ajs_zage", " Convert redshift to age of Universe, in years    See ajs_comdis for keyword parameters  ", "omega_momega_lh0       Hubble constant (default 70)  _REF_EXTRAz       Redshift (or array of values)  ", "          -1", "      10 Oct 2008 Written, Anthony Smith  ");
  
  

