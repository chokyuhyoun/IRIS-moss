function si_iv_fit, wave, spec, init=init, dn2phot=dn2phot, plot=plot
  
  si_cen = 1402.77d
  si_dr = si_cen + 2.*[-1, 1]

  w_th_si = si_cen/3d8*sqrt(8.*alog(2.)*1.38d-23*10d0^(4.9)/(28.0855*1.6605d-27))  ; in angstrom
  w_inst = 0.026  ; in angstrom

  if n_elements(init) eq 0 then init = [30., si_cen, 5d-2, 0.]
  if n_elements(dn2phot) eq 0 then dn2phot = 4.  ; Please check using iris_get_response(time)
  
  cenp = where((wave ge si_dr[0]) and (wave le si_dr[1]))
  wavep = wave[cenp]
  specp = spec[cenp] > 0.
  
  si_err = (sqrt(specp*dn2phot))/dn2phot + 1.
  si_err[where(si_err lt 1., /null)] = 1e5
  fit_res = gaussfit(wavep, specp, coeff, nterms=4, measure_error=si_err, estimates=init, yerror=error)
  v_d = 3d5*(coeff[1]-si_cen)/si_cen
  w_fwhm = 2.*sqrt(2.*alog(2.))*coeff[2]          ;https://iris.lmsal.com/itn38/analysis_lines_iris.html#line-fitting
  w_nth = sqrt(w_fwhm^2. - w_th_si^2. - w_inst^2.)
  v_nth = w_nth*3d5/si_cen
  chisq = total(((fit_res-spec)/si_err)^2.)
  res = {coeff:coeff, w_nth:w_nth, v_nth:v_nth, v_d:v_d, wr:si_dr, chisq:chisq}
  
;  if keyword_set(plot) then begin
;    p01 = plot(wave, spec)
;    p02 = plot(wavep, fit_res, '-2r', /over, transp=50)
;  endif
;stop  
  return, res
end