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

  if n_elements(dn2phot) eq 0 then dn2phot = 18.  ; Please check using iris_get_response(time)

  mg_h_cen = 2803.5310d0
  mg_k_cen = 2796.3501d0
  mg_triplet = 2798.823d0

  mg_h_dr = mg_h_cen + 2.*[-1, 1]
  mg_k_dr = mg_k_cen + 2.*[-1, 1]
  
  res = {name:'', coeff:fltarr(9), v_d_3:0., wr:[0., 0.], chisq:0., cen:0., emiss:1}
  res = replicate(res, 3)
  res[2].name = 'Mg II triplet'
  for i=0, 1 do begin
    mg_cen = ([mg_k_cen, mg_h_cen])[i]
    mg_dr = ([[mg_k_dr], [mg_h_dr]])[*, i]
    cenp = where((wave ge mg_dr[0]) and (wave le mg_dr[1]))
    wavep = wave[cenp]
    wavepn = dindgen(1d4)*(wavep[-1]-wavep[0])*1d-4 + wavep[0]
    specp = spec[cenp] > 0.

    mg_err = sqrt(specp*dn2phot)/dn2phot + 1.
    mg_err[where(mg_err lt 1., /null)] = 1e5

    if n_elements(init) eq 0 then init0 = [0, 1d0, mean(mg_dr), 1d3, mean(mg_dr), 0.04, 0.95, mean(mg_dr), 0.01] $
                             else init0 = init
    lims = {value:0., fixed:0, limited:[0, 0], limits:[0., 0.]}
    lims = replicate(lims, 9)
    lims[0].limited[0] = 1 & lims[0].limits[0] = 0d                      ; background
    lims[1].limited[*] = 1 & lims[1].limits = [1e-5, 1d1]                ; linear slope
    lims[2].limited[*] = 1 & lims[2].limits = mean(mg_dr) + 0.5*[-1, 1]   ; centeroid of linear slope
    lims[3].limited[0] = 1 & lims[3].limits[0] = 150d0                   ; positive Gaussian amplitude
    lims[4].limited[*] = 1 & lims[4].limits = mean(mg_dr) + 0.6*[-1, 1]   ; centeroid of positive Gaussian
    lims[6].limited[*] = 1 & lims[6].limits = [0, 20]                    ; fractional of negative Gaussian amplitude
    lims[7].limited[*] = 1 & lims[7].limits = mean(mg_dr) + 0.4*[-1, 1]   ; centeroid of negative Gaussian
    lims[8].limited[1] = 1 & lims[8].limits[1] = 0.1
    p = mpfitfun('iris_mgfit', wavep, specp, mg_err, init0, $
                  parinfo=lims, perror=g, quiet=1, $
                  maxiter=400, status=st, ftol=1d-9)
    fit_res = iris_mgfit(wavepn, p)              
    mg_rvpos = iris_postfit_gen(wavepn, fit_res)
    v_d_3 = (mg_rvpos[8]-mg_cen)*3d5/mg_cen
    chisq = total((specp-iris_mgfit(wavep, p)/mg_err)^2.)
    res[i].coeff = p
    res[i].v_d_3 = (mg_rvpos[10] eq 1) ? v_d_3 : !values.f_nan
    res[i].wr = mg_dr
    res[i].chisq = chisq
    res[i].name = (['Mg II k', 'Mg II h'])[i]
    res[i].cen = mg_cen
;    stop 
  endfor
  
  mg_trip_dr = [2797.5, 2802.5]
  mg_trip_part = where(wave ge mg_trip_dr[0] and wave le mg_trip_dr[1])
  coeff = poly_fit(wave[mg_trip_part], spec[mg_trip_part], 2)
  partp = where((wave ge mg_triplet-0.25) and (wave le mg_triplet+0.25))
  res[2].wr = mg_trip_dr
  res[2].emiss = total((spec - poly(wave, coeff))[partp]) 
  res[2].cen = mg_triplet
;  stop
;  if keyword_set(plot) then begin
;    p01 = plot(wave, spec)
;    for i=0, 1 do begin
;      dum = [res[i].wr[0]:res[i].wr[1]:1d-4]
;      p02 = plot(dum, iris_mgfit(dum, res[i].coeff), '-r', transp=50, over=p01)
;    endfor
;  endif
;stop
  return, res
end