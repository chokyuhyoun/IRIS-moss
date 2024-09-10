function si_iv_fit, wave, spec, init=init, dn2phot=dn2phot, plot=plot, den_cal=den_cal, spec_bin=spec_bin
  ;; positive Doppler velocity --> downward motion (redshift)
  ;; input spectrum should be normalized by exposure time, spectral pixel binning, spatial pixel binning
  ;; to calculate the Si IV total intensity, we need to restore the spectral binning. 
  
  if n_elements(den_cal) eq 0 then den_cal = 1
  if n_elements(spec_bin) eq 0 then spec_bin = 2
  res = {coeff:replicate(!values.d_nan, 4), w_nth:!values.d_nan, v_nth:!values.d_nan, v_d:!values.d_nan, $
         wr:replicate(!values.d_nan, 2), chisq:!values.d_nan, i_tot:!values.d_nan, log_den:!values.d_nan}
  if total(where(finite(spec))) eq -1 then return, res
  if n_elements(dn2phot) eq 0 then dn2phot = 4.  ; Please check using iris_get_response(time)
  
  si_cen = 1402.77d
  si_dr = si_cen + 2.*[-1, 1]
  si_err = (sqrt(spec*dn2phot))/dn2phot + 1.
  w_th_si = si_cen/3d8*sqrt(8.*alog(2.)*1.38d-23*10d0^(4.9)/(28.0855*1.6605d-27))  ; in angstrom
  ;; ~ 0.053 angstrom (https://iris.lmsal.com/itn38/diagnostics.html --> 0.05)
  w_inst = 0.026 ; in angstrom   
  
  ooe_fac = 1./(2.*sqrt(alog(2))) ; 1/e factor for w_{1/e} 
  w_th_si = w_th_si*ooe_fac  ; ~ 0.032 A
  w_inst = w_inst*ooe_fac    ; ~ 0.016 A

  cenp = where((wave ge si_dr[0]) and (wave le si_dr[1]))
  wavep = wave[cenp]
  specp = spec[cenp] > 0.
  dum = min(abs(wavep-si_cen), cen_ind)
  
  si_errp = (sqrt(specp*dn2phot))/dn2phot + 1.
  si_errp[where(si_errp lt 1., /null)] = 1e5
;  if n_elements(init) eq 0 then init = [30., si_cen, 5d-2, 0.] 
;  if n_elements(init) eq 0 then init = [specp[cen_ind], si_cen, 5d-2, 0.] ; for observational data
  if n_elements(init) eq 0 then init = [max(specp, ind), wavep[ind], 5d-2, 0.] ; for simulation data
  fit_res = gaussfit(wavep, specp, coeff, nterms=4, measure_error=si_errp, estimates=init, yerror=error)
  v_d = 3d5*(coeff[1]-si_cen)/si_cen
  w_1e = sqrt(2)*coeff[2]          ;https://iris.lmsal.com/itn38/diagnostics.html w_{1/e}
  w_nth = sqrt(w_1e^2. - w_th_si^2. - w_inst^2.)
  v_nth = w_nth*3d5/si_cen
  chisq = total((specp-fit_res)^2./fit_res, /nan)
  i_tot = total(specp, /nan) * spec_bin  ; Because input spec is normalized by spatial pixel and spectal pixel
;  stop
  if den_cal eq 0 then den = !values.d_nan else begin
    o4_1 = 1399.77
    o4_2 = 1401.16
    range = 0.5
  
    bg_dr = where(wave ge 1404. and wave le 1404.5, /null)
    sigma0 = stddev(spec[bg_dr], /nan)
    sigma1 = sqrt(sigma0^2.*n_elements(o4_1_dr))
    sigma2 = sqrt(sigma0^2.*n_elements(o4_2_dr))

    o4_1_dr = where(abs(wave-o4_1) le 0.5*range)
    o4_2_dr = where(abs(wave-o4_2) le 0.5*range)
    o4_1_int = total(spec[o4_1_dr], /nan)
    o4_2_int = total(spec[o4_2_dr], /nan)

    if (o4_1_int lt 3.*sigma1) or (o4_2_int lt 3.*sigma2) then begin
      den = !values.d_nan
;      print, 'a' 
    endif else begin
      ratio = o4_1_int/o4_2_int                  
      rat0=[0.16710843, 0.1670735, 0.16709631, 0.16718450, 0.16736760, 0.16770638, 0.16831132, 0.16937367, 0.17121214,$
            0.17433385, 0.17949909, 0.18777250, 0.20053728, 0.2193597, 0.24536467, 0.2779118, 0.3135369, 0.34697879,$
            0.3740375, 0.39339120, 0.40603911, 0.41382271, 0.4184366, 0.42111109, 0.42264118] ; from iris_ne_oiv.pro
      den0= [1.00000e+07, 1.77828e+07, 3.16228e+07, 5.62341e+07, 1.00000e+08, 1.77828e+08, 3.16228e+08, 5.62341e+08, 1.00000e+09,$
            1.77828e+09, 3.16228e+09, 5.62341e+09, 1.00000e+10, 1.77828e+10, 3.16228e+10, 5.62341e+10, 1.00000e+11, 1.77828e+11,$
            3.16228e+11, 5.62341e+11, 1.00000e+12, 1.77828e+12, 3.16228e+12, 5.62341e+12, 1.00000e+13]
      den = (ratio ge min(rat0) and ratio le max(rat0)) ? interpol(alog10(den0), rat0, ratio, /nan) : !values.d_nan 
    endelse
  endelse
  
  res.coeff = coeff
  res.w_nth = w_nth
  res.v_nth = v_nth
  res.v_d = v_d
  res.wr = si_dr
  res.chisq = chisq
  res.i_tot = i_tot
  res.log_den = den
;  stop
;  if keyword_set(plot) then begin
;    p01 = plot(wave, spec)
;    p02 = plot(wavep, fit_res, '-2r', /over, transp=50)
;  endif

  return, res
end