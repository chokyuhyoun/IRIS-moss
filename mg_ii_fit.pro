function mg_ii_fit, wave, spec, init=init, dn2phot=dn2phot, plot=plot
; please check iris_mgfit.pro
;PURPOSE-generates a vector of a double gaussian + wing line profile
;        described in schmit et al 2015 ApJ
;INPUTS: X-vector of the abscissa
;        P- 9 element vector that contains the gaussian parameters
;P[0]=Background intensity
;P[1]=wing slope
;P[2]=centroid of linear wing profile
;P[3]=amplitude of positive gaussian
;P[4]=centroid of positive gaussian
;P[5]=width^2 of the positive gaussian
;P[6]=fractional amplitude of negative gaussian
;P[7]=centroid of negative gaussian;
;P[8]=width^2 of negative gaussian
;
;OUTPUTS: Vector of same size as X with the line profile

  if n_elements(init) eq 0 then init = [0, 1d0, mean(mg_dr), 1d3, mean(mg_dr), 0.04, 0.95, mean(mg_dr), 0.01]
  if n_elements(dn2phot) eq 0 then dn2phot = 18.  ; Please check using iris_get_response(time)

  mg_h_cen = 2803.5310d0
  mg_k_cen = 2796.3501d0
  mg_triplet = 2798.823d0
  mg_h_dr = mg_h_cen + 2.*[-1, 1]
  mg_k_dr = mg_k_cen + 2.*[-1, 1]

  lims = {value:0., fixed:0, limited:[0, 0], limits:[0., 0.]}
  lims = replicate(lims, 9)
  lims[0].limited[0] = 1 & lims[0].limits[0] = 0d                      ; background
  lims[1].limited[*] = 1 & lims[1].limits = [1e-5, 1d1]                ; linear slope
  lims[2].limited[*] = 1 & lims[2].limits = mean(wave) + 0.5*[-1, 1]   ; centeroid of linear slope
  lims[3].limited[0] = 1 & lims[3].limits[0] = 150d0                   ; positive Gaussian amplitude
  lims[4].limited[*] = 1 & lims[4].limits = mean(wave) + 0.6*[-1, 1]   ; centeroid of positive Gaussian
  lims[6].limited[*] = 1 & lims[6].limits = [0, 20]                    ; fractional of negative Gaussian amplitude
  lims[7].limited[*] = 1 & lims[7].limits = mean(wave) + 0.4*[-1, 1]   ; centeroid of negative Gaussian
  lims[8].limited[1] = 1 & lims[8].limits[1] = 0.1
  
  res = {name:'Mg II k', coeff:0, v_d_3:0, wavep:0, fit_res:0, chisq:0}
  res = [res, res]
  res[1].name = 'Mg II h'
  for i=0, 1 do begin
    mg_cen = ([mg_k_cen, mg_h_cen])[i]
    mg_dr = ([[mg_k_dr], [mg_h_dr]])[*, i]
    cenp = where((wave ge mg_dr[0]) and (mg_wave le mg_dr[1]))
    wavep = mg_wave[cenp]
    wavepn = dindgen(1d4)*(mg_wavep[-1]-mg_wavep[0])*1d-4 + mg_wavep[0]
    specp = spec[cenp]

    mg_err = sqrt(specp*dn2phot)/dn2phot + 1.
    mg_err[mg_err lt 1.] = 1e5
    p = mpfitfun('iris_mgfit', wavep, specp, mg_err, init, $
                  parinfo=lims, perror=g, /quiet, $
                  maxiter=400, status=st, ftol=1d-9)
    fit_res = iris_mgfit(wavepn, p)              
    mg_rvpos = iris_postfit_gen(wavepn, fit_res)
    v_d_3 = (mg_rvpos[8]-mg_cen)*3d5/mg_cen
    chisq = total((specp-iris_mgfit(wavep, p)/mg_err)^2.)
    res[i].coeff = p
    res[i].v_d_3 = v_d_3
    res[i].wavep = wavepn
    res[i].fit_res = fit_res
    res[i].chisq = chisq
  endfor
  
  if keyword_set(plot) then begin
    p01 = plot(wave, spec)
    for i=0, 1 do p02 = plot(res[i].wavep, res[i].fit_res, '-r', transp=50, over=p01)
  endif
  return, res
end