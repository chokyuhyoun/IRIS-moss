function si_iv_fit, wave, spec, init=init, dn2phot=dn2phot, plot=plot
  if n_elements(init) eq 0 then init = [30., si_cen, 1d-2, 0.]
  if n_elements(dn2phot) eq 0 then dn2phot = 4.  ; Please check using iris_get_response(time)
  
  w_th_si = si_cen/3d8*sqrt(8.*alog(2.)*1.38d-23*10d0^(4.9)/(28.0855*1.6605d-27))  ; in angstrom
  w_inst = 0.026  ; in angstrom
  si_cen = 1402.77d
  si_dr = si_cen + 1.*[-1, 1]
  
  cenp = where((wave ge si_dr[0]) and (wave le si_dr[1]))
  wavep = wave[cenp]
  specp = spec[cenp]
  
  si_err = sqrt(specp*dn2phot) > 0)/dn2phot + 1.
  si_err[si_err lt 1.] = 1e5
  fit_res = gaussfit(wavep, specp, coeff, nterms=4, measure_error=si_err, estimates=init, yerror=error)
  v_d = 3d5*(coeff[1]-si_cen)/si_cen
  w_fwhm = 2.*sqrt(2.*alog(2.))*coeff[2]          ;https://iris.lmsal.com/itn38/analysis_lines_iris.html#line-fitting
  w_nth = sqrt(w_fwhm^2. - w_th_si^2. - w_inst^2.)
  chisq = total(((fit_res-spec)/si_err)^2.)
  res = {coeff:coeff, w_nth:w_nth, v_d:v_d, wavep:wavep, fit_res=fit_res, chisq:chisq}
  
  if keyword_set(plot) then begin
    p01 = plot(wave, spec)
    p02 = plot(wavep, fit_res, '-r', transp=50, over=p01)
  endif
  return, res
end