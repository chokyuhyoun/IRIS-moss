dir = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, dir
obj_dir = file_search(dir, 'ar?????', /test_dir)
sol_x = !null
sol_y = !null
times = !null
si_v = !null
si_nth = !null
si_peak = !null
cur_fwhm = !null
mg_h3v = !null
mg_k3v = !null
e_den = !null
sg_ind = !null
sg_dy = !null
ar_no = !null
mg_trip = !null
for i=0, n_elements(obj_dir)-1 do begin
  ar_num = fix(strmid((strsplit(obj_dir[i], '/', /extract))[-1], 2, 5))
  sav_files = file_search(obj_dir[i], '/*moss_event*.sav')
  
  for j=0, n_elements(sav_files)-1 do begin
    print, j, sav_files[j]
    restore, sav_files[j]
    ar_no = [ar_no, replicate(ar_num, eout.moss_num)] 
    sol_x = [sol_x, reform(eout.sg_phy[0, *])]
    sol_y = [sol_y, reform(eout.sg_phy[1, *])]
    times = [times, reform(eout.sg_phy[2, *])]
    si_peak = [si_peak, reform(eout.si_iv_fit_res.coeff[0, *])]
    si_v = [si_v, eout.si_iv_fit_res.v_d]
    si_nth = [si_nth, eout.si_iv_fit_res.v_nth]
    cur_fwhm = [cur_fwhm, eout.curve_fwhm]
    mg_k3v = [mg_k3v, reform(eout.mg_ii_fit_res[0, *].v_d_3)]
    mg_h3v = [mg_h3v, reform(eout.mg_ii_fit_res[1, *].v_d_3)]    
    mg_trip = [mg_trip, reform(eout.mg_ii_fit_res[2, *].emiss)]
    dum = eout.sg_ind
    dum[2, *] = j
    sg_ind = [[sg_ind], [dum]]   ; [y_pix, time in a dataset, dataset]
    sg_dy = [[sg_dy], replicate(eout.sg_dy, eout.moss_num)]

    restore, file_dirname(sav_files[j])+'/emcube.sav'
    dem_t = anytim(demstr.time)
    dem_t1 = rebin(dem_t, n_elements(dem_t), eout.moss_num)
    moss_t1 = rebin(eout.aia_phy[2, *], n_elements(dem_t), eout.moss_num)
    t_diff = abs(dem_t1 - moss_t1)
    t_diff_min = min(t_diff, dim=1, t_min_ind0)
    t_prev_ind = reform((array_indices(t_diff, t_min_ind0))[0, *]) - 2 > 0
    for k=0, eout.moss_num-1 do begin
      e_den = [e_den, total(demstr[t_prev_ind[k]].dem[eout.aia_ind[0, k], eout.aia_ind[1, k], *])]
    endfor
    eout = 0    
  endfor
endfor
event_no = intarr(n_elements(times))
no = 0
for ii=1, n_elements(event_no)-1 do begin
  same_event = where((abs(sg_ind[0, 0:ii-1] - sg_ind[0, ii]) le 2./sg_dy[ii]) and $ ; along the slit, within 2"
                     (abs(sg_ind[1, 0:ii-1] - sg_ind[1, ii]) le 1.) and $ ; successive sit-and-stare scan
                     (sg_ind[2, 0:ii-1] eq sg_ind[2, ii]), count) ; in same observation program
  if count eq 0 then no += 1
  event_no[ii] = no
endfor
pars = [[ar_no], [sol_x], [sol_y], [times], [si_v], [si_nth], [si_peak], [cur_fwhm], $
        [mg_h3v], [mg_k3v], [mg_trip], [e_den]]
par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'si_v', 'si_nth', 'si_peak', 'cur_fwhm', $
             'mg_h3v', 'mg_k3v', 'mg_trip', 'e_den']
par_titles = ['AR no', 'Solar X (arcsec)', 'Solar Y (arcsec)', 'Time', 'v$_{D, Si IV}$ (km s$^{-1}$)', $
              'v$_{nth, Si IV}$ (km s$^{-1}$)', 'I$_{Peak, Si IV}$ (DN)', 'Event Duration (s)', $
              'v$_{D, Mg II h3}$ (km $s^{-1}$)', 'v$_{D, Mg II k3}$ (km $s^{-1}$)', 'EW$_{Mg II triplet} (\AA$)', $
              'Emission Measure ($10^{26} cm^{-5}$)']
sz = [n_elements(event_no), max(event_no)+1, n_elements(par_names)]
cast = make_array(sz[0], sz[1], value = !values.f_nan)
cast[indgen(sz[0]), event_no] = 1
cast_all = rebin(cast, sz[0], sz[1], sz[2])
pars_all = rebin(reform(pars, sz[0], 1, sz[2]), sz[0], sz[1], sz[2]) * cast_all
pars_ev_mean = mean(pars_all, dim=1, /nan)
pars_ev_std = stddev(pars_all, dim=1, /nan)
pars_ev_std[where(~finite(pars_ev_std))] = 0.

com0 = !null
com1 = !null
for i=0, n_elements(par_names)-1 do begin
  com0 = [com0, par_names[i]+':pars_ev_mean[*, '+string(i, f='(i0)')+']']
  com1 = [com1, par_names[i]+':pars_ev_std[*, '+string(i, f='(i0)')+']']
endfor
dum = execute('pars_event_mean = {'+strjoin(com0, ', ')+'}')
dum = execute('pars_event_std = {'+strjoin(com1, ', ')+'}')

save, pars, par_names, pars_event_mean, pars_event_std, event_no, par_titles, $
      filename=dir + '/moss_param_event_total.sav'



end

