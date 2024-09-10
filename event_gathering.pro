dir = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, dir
obj_dir = file_search(dir, 'ar?????', /test_dir)
sol_x = !null
sol_y = !null
times = !null
si_v = !null
si_nth = !null
si_amp = !null
si_i_tot = !null
log_den = !null
cur_fwhm = !null
mg_h3v = !null
mg_k3v = !null
;em = !null
sg_ind = !null
sg_phy = !null
sg_dy = !null
ar_no = !null
mg_trip = !null
si_err = !null
mg_k_err = !null
mg_h_err = !null
mg_k_centr = !null
mg_h_centr = !null
mg_k_tot = !null
mg_h_tot = !null

pre_si_v = !null
pre_si_nth = !null
pre_si_amp = !null
pre_si_i_tot = !null
pre_log_den = !null
pre_mg_h3v = !null
pre_mg_k3v = !null
pre_mg_k_centr = !null
pre_mg_h_centr = !null
pre_mg_k_tot = !null
pre_mg_h_tot = !null
;pre_em = !null
pre_mg_trip = !null



aft_si_v = !null
aft_si_nth = !null
aft_si_amp = !null
aft_si_i_tot = !null
aft_log_den = !null
aft_mg_h3v = !null
aft_mg_k3v = !null
aft_mg_k_centr = !null
aft_mg_h_centr = !null
aft_mg_k_tot = !null
aft_mg_h_tot = !null


;pre_em = !null
pre_mg_trip = !null
aft_mg_trip = !null

fe18_avg = !null

aaa1 = !null
aaa2 = !null
aaa3 = !null
si_dlam = !null

temp_array = 5.5+findgen(21)*0.1
tr_t_ind = where(temp_array le 6.0)
if 0 then begin
  for i=0, n_elements(obj_dir)-1 do begin
    if strmatch(obj_dir[i], '*12965*') then continue
    ar_num = fix(strmid((strsplit(obj_dir[i], '/', /extract))[-1], 2, 5))
    sav_files = file_search(obj_dir[i], '/*moss_event*.sav')
    
    for j=0, n_elements(sav_files)-1 do begin
      print, j, sav_files[j]
      restore, sav_files[j]
      ar_no = [ar_no, replicate(ar_num, eout.moss_num)] 
      sol_x = [sol_x, reform(eout.sg_phy[0, *])]
      sol_y = [sol_y, reform(eout.sg_phy[1, *])]
      times = [times, reform(eout.sg_phy[2, *])]
      si_amp = [si_amp, reform(eout.si_iv_fit_res.coeff[0, *])]
      si_v = [si_v, eout.si_iv_fit_res.v_d]
      si_nth = [si_nth, eout.si_iv_fit_res.v_nth]
      si_i_tot = [si_i_tot, eout.si_iv_fit_res.i_tot]
      
      si_iv_ind = (where(strmatch(eout.line_id, '*Si IV*'), /null))[-1]
      si_wv = eout.wave_list[si_iv_ind]      
      si_dlam = [si_dlam, replicate(mean(si_wv[1:*]-si_wv[0:-2], /nan), eout.moss_num)]
      
;      print, n_elements(si_i_tot)
      log_den = [log_den, eout.si_iv_fit_res.log_den]
      cur_fwhm = [cur_fwhm, eout.curve_fwhm]
      mg_k3v = [mg_k3v, reform(eout.mg_ii_fit_res[0, *].v_d_3)]
      mg_h3v = [mg_h3v, reform(eout.mg_ii_fit_res[1, *].v_d_3)]    
      mg_trip = [mg_trip, reform(eout.mg_ii_fit_res[2, *].emiss)]
      si_err = [si_err, reform(eout.si_iv_fit_res.chisq)]
      mg_k_err = [mg_k_err, reform(eout.mg_ii_fit_res[0, *].chisq)]
      mg_h_err = [mg_h_err, reform(eout.mg_ii_fit_res[1, *].chisq)]
      mg_k_centr = [mg_k_centr, reform(eout.mg_ii_fit_res[0, *].centr)]
      mg_h_centr = [mg_h_centr, reform(eout.mg_ii_fit_res[1, *].centr)]
      mg_k_tot = [mg_k_tot, reform(eout.mg_ii_fit_res[0, *].emiss)]
      mg_h_tot = [mg_h_tot, reform(eout.mg_ii_fit_res[1, *].emiss)]

      pre_si_v = [pre_si_v, eout.pre_si_iv_fit_res.v_d]
      pre_si_nth = [pre_si_nth, eout.pre_si_iv_fit_res.v_nth]
      pre_si_amp = [pre_si_amp, reform(eout.pre_si_iv_fit_res.coeff[0, *])]
      pre_si_i_tot = [pre_si_i_tot, eout.pre_si_iv_fit_res.i_tot]
      pre_log_den = [pre_log_den, eout.pre_si_iv_fit_res.log_den]
      pre_mg_k3v = [pre_mg_k3v, reform(eout.pre_mg_ii_fit_res[0, *].v_d_3)]
      pre_mg_h3v = [pre_mg_h3v, reform(eout.pre_mg_ii_fit_res[1, *].v_d_3)]    
      pre_mg_trip = [pre_mg_trip, reform(eout.pre_mg_ii_fit_res[2, *].emiss)]
      pre_mg_k_centr = [pre_mg_k_centr, reform(eout.pre_mg_ii_fit_res[0, *].centr)]
      pre_mg_h_centr = [pre_mg_h_centr, reform(eout.pre_mg_ii_fit_res[1, *].centr)]
      pre_mg_k_tot = [pre_mg_k_tot, reform(eout.pre_mg_ii_fit_res[0, *].emiss)]
      pre_mg_h_tot = [pre_mg_h_tot, reform(eout.pre_mg_ii_fit_res[1, *].emiss)]

      aft_si_v = [aft_si_v, eout.aft_si_iv_fit_res.v_d]
      aft_si_nth = [aft_si_nth, eout.aft_si_iv_fit_res.v_nth]
      aft_si_amp = [aft_si_amp, reform(eout.aft_si_iv_fit_res.coeff[0, *])]
      aft_si_i_tot = [aft_si_i_tot, eout.aft_si_iv_fit_res.i_tot]
      aft_log_den = [aft_log_den, eout.aft_si_iv_fit_res.log_den]
      aft_mg_k3v = [aft_mg_k3v, reform(eout.aft_mg_ii_fit_res[0, *].v_d_3)]
      aft_mg_h3v = [aft_mg_h3v, reform(eout.aft_mg_ii_fit_res[1, *].v_d_3)]
      aft_mg_trip = [aft_mg_trip, reform(eout.aft_mg_ii_fit_res[2, *].emiss)]
      aft_mg_k_centr = [aft_mg_k_centr, reform(eout.aft_mg_ii_fit_res[0, *].centr)]
      aft_mg_h_centr = [aft_mg_h_centr, reform(eout.aft_mg_ii_fit_res[1, *].centr)]
      aft_mg_k_tot = [aft_mg_k_tot, reform(eout.aft_mg_ii_fit_res[0, *].emiss)]
      aft_mg_h_tot = [aft_mg_h_tot, reform(eout.aft_mg_ii_fit_res[1, *].emiss)]

      aaa1 = [aaa1, reform(eout.mg_ii_fit_res[2, *].coeff[2])]
      aaa2 = [aaa2, reform(eout.pre_mg_ii_fit_res[2, *].coeff[2])]
      aaa3 = [aaa3, reform(eout.aft_mg_ii_fit_res[2, *].coeff[2])]
      
      dum = eout.sg_ind
      dum[2, *] = j
      sg_ind = [[sg_ind], [dum]]   ; [y_pix, time in a dataset, dataset]
      sg_phy = [[sg_phy], [eout.sg_phy]]
      sg_dy = [[sg_dy], replicate(eout.sg_dy, eout.moss_num)] ; to define the event!

;      restore, file_dirname(sav_files[j])+'/emcube.sav'
      restore, file_dirname(sav_files[j])+'/Fe_XVIII_cube.sav'
      aia_sz = size(fe18)
      for k=0, eout.moss_num-1 do begin
;        t_ind = eout.aia_ind[2, k]
;        em = [em, total(demstr[eout.aia_ind[2, k]].dem[eout.aia_ind[0, k], eout.aia_ind[1, k], tr_t_ind])]  
;        pre_em = [pre_em, total(demstr[eout.pre_aia_ind[2, k]].dem[eout.pre_aia_ind[0, k], eout.pre_aia_ind[1, k], tr_t_ind])]
;        f18_ind0 = (eout.aia_ind[2, k]-10) > 0
;        f18_ind1 = (eout.aia_ind[2, k]+10) < aia_sz[3]
;        fe18_part = fe18[*, *, f18_ind0:f18_ind1] 
;        fe18_avg = [fe18_avg, mean(fe18_part, /nan)/(f18_ind1-f18_ind0+1.)] ; [DN/pix/s]  ; [fe18] = DN/s
         fe18_avg = [fe18_avg, mean(fe18[*, *, eout.aia_ind[2, k]])]
;        print, eout.pre_sg_ind[1, k], eout.sg_ind[1, k], eout.aft_sg_ind[1, k]
      endfor
      eout = 0    
    endfor
  endfor
;  stop
  thres = 0.99
  si_err[where(~finite(si_err))] = max(si_err, /nan)
  si_pdf = histogram(si_err, location=si_xbin, binsize=1d3)
  si_cdf = total(si_pdf, /cum)/n_elements(si_err)
  si_thres = (si_xbin[where(min(abs(si_cdf-thres)) eq abs(si_cdf-thres))])[0]
   
  mg_k_err[where(~finite(mg_k_err))] = max(mg_k_err, /nan)
  mg_k_pdf = histogram(mg_k_err, location=mg_k_xbin, binsize=1d3)
  mg_k_cdf = total(mg_k_pdf, /cum)/n_elements(mg_k_err)
  mg_k_thres = (mg_k_xbin[where(min(abs(mg_k_cdf-thres)) eq abs(mg_k_cdf-thres))])[0]

  mg_h_err[where(~finite(mg_h_err))] = max(mg_h_err, /nan)
  mg_h_pdf = histogram(mg_h_err, location=mg_h_xbin, binsize=1d3)
  mg_h_cdf = total(mg_h_pdf, /cum)/n_elements(mg_h_err)
  mg_h_thres = (mg_h_xbin[where(min(abs(mg_h_cdf-thres)) eq abs(mg_h_cdf-thres))])[0]
  
  print, n_elements(si_err) ; total # = 5741
  real = where((si_err lt si_thres) and (mg_k_err lt mg_k_thres) and (mg_h_err lt mg_h_thres)) 
  sg_ind = sg_ind[*, real]
  sg_dy = sg_dy[real]
  si_err = si_err[real]
  mg_k_err = mg_k_err[real]
  mg_h_err = mg_h_err[real]
 
  event_no = intarr(n_elements(si_err))
  no = 0
  for ii=1, n_elements(event_no)-1 do begin
    same_event = where((abs(sg_phy[0, 0:ii-1] - sg_phy[0, ii]) le 1./sg_dy[ii]) and $ ; along the slit, within 1"
                      (abs(sg_ind[1, 0:ii-1] - sg_ind[1, ii]) le 1.) and $ ; successive sit-and-stare scan
                      (sg_ind[2, 0:ii-1] eq sg_ind[2, ii]), count) ; in same observation program
    if count eq 0 then no += 1
    event_no[ii] = no
  endfor
  
  pixel_no = intarr(n_elements(si_err))
  no = 0
  for ii=1, n_elements(pixel_no)-1 do begin
    same_pixel = where((sg_ind[0, 0:ii-1] eq sg_ind[0, ii]) and $  ;; same slit pos
                       abs((sg_ind[1, 0:ii-1] - sg_ind[1, ii]) le 1.) and $  ;; successive sit-and-stare scan
                       (sg_ind[2, 0:ii-1] eq sg_ind[2, ii]), count) ;; in same observational program
    if count eq 0 then no += 1
    pixel_no[ii] = no                   
  endfor
  
  par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'cur_fwhm', $
               'si_v', 'pre_si_v', 'aft_si_v', $ 
               'si_nth', 'pre_si_nth', 'aft_si_nth', $ 
               'si_amp', 'pre_si_amp', 'aft_si_amp', $
               'si_i_tot', 'pre_si_i_tot', 'aft_si_i_tot', $ 
               'log_den', 'pre_log_den', 'aft_log_den', $
               'mg_h3v', 'pre_mg_h3v', 'aft_mg_h3v', $ 
               'mg_k3v', 'pre_mg_k3v', 'aft_mg_k3v', $
               'mg_h_centr', 'pre_mg_h_centr', 'aft_mg_h_centr', $
               'mg_k_centr', 'pre_mg_k_centr', 'aft_mg_k_centr', $
               'mg_k_tot', 'pre_mg_k_tot', 'aft_mg_k_tot', $
               'mg_h_tot', 'pre_mg_h_tot', 'aft_mg_h_tot', $
               'mg_trip', 'pre_mg_trip', 'aft_mg_trip', $
;               'em', 'pre_em', $
               'fe18_avg']
  par_list = !null
  for ii=0, n_elements(par_names)-1 do par_list = [par_list, '['+par_names[ii]+']']  
  dum = execute('pars = ['+strjoin(par_list, ', ')+']') 

  par_dr = [[0, 0], [-1d3, 1d3], [-1d3, 1d3], [0, 0], [0, 60], $
            [-50, 50], [-50, 50], [-50, 50], $ 
            [0, 100], [0, 100], [0, 100], $ 
            [0, 100], [0, 100], [0, 100], $
            [0, 4d3], [0, 4d3], [0, 4d3], $
            [9, 13], [9, 13], [9, 13], $
            [-15, 15], [-15, 15], [-15, 15], $
            [-15, 15], [-15, 15], [-15, 15], $
            [-15, 15], [-15, 15], [-15, 15], $
            [-15, 15], [-15, 15], [-15, 15], $
            [0, 1d4], [0, 1d4], [0, 1d4], $
            [0, 1d4], [0, 1d4], [0, 1d4], $
            [-0.25, 0.35], [-0.25, 0.35], [-0.25, 0.35], $
;            [0, 500], [0, 500], $
            [-0.5, 15]]   

  par_titles = ['AR no', 'Solar X (arcsec)', 'Solar Y (arcsec)', 'Time (UT)', 'Event Duration (s)', $
                ['', 'Prev. ', 'After ']+'v$_{D, Si IV}$ (km s$^{-1}$)', $
                ['', 'Prev. ', 'After ']+'v$_{nth, Si IV}$ (km s$^{-1}$)', $
                ['', 'Prev. ', 'After ']+'I$_{Amp, Si IV}$ (DN pix$_x^{-1}$ pix$_\lambda^{-1}$ s$^{-1}$)', $
                ['', 'Prev. ', 'After ']+'I$_{tot, Si IV}$ (DN pix$_x^{-1}$ s$^{-1}$)', $
                ['', 'Prev. ', 'After ']+'Log $n_e$ (cm$^{-3}$)', $
                ['', 'Prev. ', 'After ']+'v$_{D, Mg II h3}$ (km $s^{-1}$)', $
                ['', 'Prev. ', 'After ']+'v$_{D, Mg II k3}$ (km $s^{-1}$)', $
                ['', 'Prev. ', 'After ']+'v$_{centr, Mg II h}$ (km $s^{-1}$)', $
                ['', 'Prev. ', 'After ']+'v$_{centr, Mg II k}$ (km $s^{-1}$)', $
                ['', 'Prev. ', 'After ']+'I$_{tot, Mg II k}$ (DN pix$_x^{-1}$ s$^{-1}$)', $
                ['', 'Prev. ', 'After ']+'I$_{tot, Mg II h}$ (DN pix$_x^{-1}$ s$^{-1}$)', $
                ['', 'Prev. ', 'After ']+'EW$_{Mg II triplet} (\AA$)', $
;                'EM(logT=5.5-6.0) ($10^{26} cm^{-5}$)', 'Prev. EM(logT=5.5-6.0) ($10^{26} cm^{-5}$)', $
                'I$_{Fe XVIII}$ (DN pix$^{-1} s^{-1}$)']
  pars = pars[real, *]
endif else restore, dir + '/moss_param_event_total.sav'

sz = [n_elements(event_no), max(event_no)+1, n_elements(par_names)]   ;[total_spectra_no, event_no, par_no]
cast = make_array(sz[0], sz[1], value = !values.f_nan)
cast[indgen(sz[0]), event_no] = 1
cast_all = rebin(cast, sz[0], sz[1], sz[2])
pars_all = rebin(reform(pars, sz[0], 1, sz[2]), sz[0], sz[1], sz[2]) * cast_all

; average of each event
pars_ev_mean = mean(pars_all, dim=1, /nan)
pars_ev_std = stddev(pars_all, dim=1, /nan)
pars_ev_std[where(~finite(pars_ev_std))] = 0.


; maximum of Si IV I tot in each event
ref_ind = !null
si_i_tot = pars[*, where(par_names eq 'si_i_tot')]
for ii=0, max(event_no) do begin
  event_ind = where(event_no eq ii, count, /null)
  ind0 = (count lt 2 ) ? event_ind[0] : where(si_i_tot eq max(si_i_tot[event_ind]))
  ref_ind = [ref_ind, ind0]
endfor
pars_ev_peak = pars[ref_ind, *]

;; maximum in each event
dum = max(abs(pars_all), /nan, dim=1, index)
pars_ev_max = dum*sgn(pars_all[index])

stop
;; Si IV I tot peak in each pixel
ref_ind = !null
si_i_tot = pars[*, where(par_names eq 'si_i_tot')]
ref_ind_mg = !null
mg_k_tot = pars[*, where(par_names eq 'mg_k_tot')]
for ii=0, max(pixel_no) do begin
  pixel_ind = where(pixel_no eq ii, count , /null)
  ind0 = (count lt 2 ) ? pixel_ind[0] : where(si_i_tot eq max(si_i_tot[pixel_ind]))
  if ind0 lt min(pixel_ind) or ind0 gt max(pixel_ind) then stop
  ref_ind = [ref_ind, ind0]

  ind1 = (count lt 2 ) ? pixel_ind[0] : where(mg_k_tot eq max(mg_k_tot[pixel_ind]))
  if ind0 lt min(pixel_ind) or ind0 gt max(pixel_ind) then stop
  ref_ind_mg = [ref_ind_mg, ind1]

endfor
pars_pix_peak = pars[ref_ind, *]

if 1 then begin ;; mg analysis using Mg k peak moment
  mg_params = where(strmatch(par_names, '*mg_*'))
  pars_pix_peak[*, mg_params] = (pars[ref_ind_mg, *])[*, mg_params]
endif
;stop

var_names = ['pars', 'pars_ev_mean', 'pars_ev_std', 'pars_ev_peak', $
             'pars_ev_max', 'pars_pix_peak']

for j=0, n_elements(var_names)-1 do begin
  com0 = !null
  for i=0, n_elements(par_names)-1 do begin
    com0 = [com0, par_names[i]+':'+var_names[j]+'[*, '+string(i, f='(i0)')+']']
  endfor
  dum = execute(var_names[j]+' = {'+strjoin(com0, ', ')+'}')
endfor

save, pars, par_names, var_names, pars_ev_mean, pars_ev_std, pars_ev_peak, pars_ev_max, $
      pars_pix_peak, event_no, pixel_no, par_titles, si_err, mg_k_err, mg_h_err, par_dr, si_dlam, $
      filename=dir + '/moss_param_event_total_mg.sav'
end

