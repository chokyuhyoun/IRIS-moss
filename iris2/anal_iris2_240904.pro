function mask_ltau, im, ltau1, ltau2, color=color
  if n_elements(ltau1) eq 0 then ltau1 = -3
  if n_elements(ltau2) eq 0 then ltau2 = -6
  if n_elements(color) eq 0 then color='gray'
  im.getdata, image, x, y
  nx_zeros = fltarr(n_elements(x))
  mask = []
  for i=0, 1 do begin
    for j=-1, 1, 2 do begin
      ltau = ([ltau1, ltau2])[i]
      mask1 = polygon([x, reverse(x)], $
                      [im.yr[i] + nx_zeros, $
                      ltau + nx_zeros], /data, over=im, $
                      color=color, fill_color=color, $
                      pattern_orient=45*j, pattern_spacing=20)
      mask = [mask, mask1]
    endfor  
  endfor
;  stop                    
  return, mask                    
end


path = '/Users/khcho/Desktop/IRIS-moss-main'
iris2_path = path+'/iris2'
cd, iris2_path


restore, '/Users/khcho/Desktop/IRIS-moss-main/moss_param_event_total_mg2.sav'

;par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'cur_fwhm', $
;  'si_v', 'pre_si_v', 'aft_si_v', $
;  'si_nth', 'pre_si_nth', 'aft_si_nth', $
;  'si_amp', 'pre_si_amp', 'aft_si_amp', $
;  'si_i_tot', 'pre_si_i_tot', 'aft_si_i_tot', $
;  'log_den', 'pre_log_den', 'aft_log_den', $
;  'mg_h3v', 'pre_mg_h3v', 'aft_mg_h3v', $
;  'mg_k3v', 'pre_mg_k3v', 'aft_mg_k3v', $
;  'mg_h_centr', 'pre_mg_h_centr', 'aft_mg_h_centr', $
;  'mg_k_centr', 'pre_mg_k_centr', 'aft_mg_k_centr', $
;  'mg_k_tot', 'pre_mg_k_tot', 'aft_mg_k_tot', $
;  'mg_h_tot', 'pre_mg_h_tot', 'aft_mg_h_tot', $
;  'mg_trip', 'pre_mg_trip', 'aft_mg_trip', $
;  ;               'em', 'pre_em', $
;  'fe18_avg']
;
;save, pars, pars_struct, par_names, var_names, pars_ev_mean, pars_ev_std, pars_ev_peak, pars_ev_max, $
;      pars_pix_peak, event_no, pixel_no, par_titles, si_err, mg_k_err, mg_h_err, par_dr, si_dlam, $
;      pars_ev_peak_ind, pars_pix_peak_ind, sg_ind, pre_sg_ind, sg_file, curve_dt, $
;      filename=dir + '/moss_param_event_total_mg2.sav'
;

;iris2_moss0 = {input_filename:iris2model.input_filename, $
;              output_filename:iris2model.output_filename, $
;              ltau:iris2model.ltau, $
;              model:model, $
;              model_dimension_info:iris2model.model_dimension_info, $
;              physical_variable:iris2model.physical_variable, $
;              uncertainty:uncert, $
;              chi2:chi2, $
;              mu:iris2model.mu, $
;              mean_model:mean_model, $
;              peak_ind:fix(0.5*(nt-1)), $
;              dt:dt, $
;              time_from_peak:t_arr}

sg_files = sg_file[pars_pix_peak_ind]
sg_ind = sg_ind[*, pars_pix_peak_ind] ;[y_pix, time in a dataset, dataset]
pre_sg_ind = pre_sg_ind[*, pars_pix_peak_ind] ;[y_pix, time in a dataset, dataset]


restore, iris2_path+'/iris2_moss.sav'

; iris2_moss.model = [time, ltau, parameters, events]
; parameters : T [k],  v_los [cm/s],  v_turb [cm/s],  n_e [cm^-3]
ltau = iris2_moss[0].ltau
t_arr = iris2_moss[0].time_from_peak
;mean_temp = iris2_moss[0].mean_model[*, 0]
before_temp = reform(mean(iris2_moss.model[0:arr_eq(t_arr, -30), *, 0, *], dim=1, /nan))
event_temp = reform(iris2_moss.model[*, *, 0, *])
sz = size(event_temp)

before_temp_3d = rebin(reform(before_temp, 1, sz[2], sz[3]), sz[1], sz[2], sz[3])
temp_diff = event_temp - before_temp_3d
mean_temp_diff = mean(temp_diff, dim=3, /nan)

im01 = image_(mean_temp_diff, t_arr, ltau, $
              rgb_table=3, title='Temperature Variation', $
              xtitle='Time from Si IV peak (s)', $
              ytitle='Log($\tau$)', yr=reverse(minmax(ltau)), cb=cb01) 
cb01.title='T - T$_{0}$  (K)'
mask01 = mask_ltau(im01)

mean_temp_diff_ratio = mean(temp_diff/before_temp_3d, dim=3, /nan)
im02 = image_(mean_temp_diff_ratio*1d2, t_arr, ltau, $
              rgb_table=34, title='Percentage of Temp. enhancement', $
              xtitle='Time from Si IV peak (s)', $
              ytitle='Log($\tau$)', yr=reverse(minmax(ltau)), cb=cb02)
cb02.title='(T - T$_{0}$) / T$_{0}$  (%)'
mask02 = mask_ltau(im02)

mean_v_los = mean(iris2_moss.model[*, *, 1, *], dim=3, /nan)*1d-5
im03 = image_(mean_v_los, t_arr, ltau, $
              rgb_table=62, $
              min=0, max=2, $
              title='Line of Sight Velocity', $
              xtitle='Time from Si IV peak (s)', $
              ytitle='Log($\tau$)', yr=reverse(minmax(ltau)), cb=cb03)
cb03.title='V$_{los}$  (km $s^{-1}$)'
mask03 = mask_ltau(im03)

mean_v_turb = mean(iris2_moss.model[*, *, 2, *], dim=3, /nan)*1d-5
im04 = image_(mean_v_los, t_arr, ltau, $
              rgb_table=4, min=0, max=2, $
              title='Microturbulent Velocity', $
              xtitle='Time from Si IV peak (s)', $
              ytitle='Log($\tau$)', yr=reverse(minmax(ltau)), cb=cb04)
cb04.title='V$_{turb}$  (km $s^{-1}$)'
mask04 = mask_ltau(im04)

before_ne = reform(mean(iris2_moss.model[0:arr_eq(t_arr, -30), *, 3, *], dim=1, /nan))
event_ne = reform(iris2_moss.model[*, *, 3, *])
before_ne_3d = rebin(reform(before_ne, 1, sz[2], sz[3]), sz[1], sz[2], sz[3])
ne_diff = event_ne - before_ne_3d
mean_ne_diff = mean(ne_diff, dim=3, /nan)
mean_ne_diff_ratio = mean(ne_diff/before_ne_3d, dim=3, /nan)
im05 = image_(mean_ne_diff_ratio*1d2, t_arr, ltau, $
              rgb_table=34, title='Percentage of Electron Density enhancement', $
              xtitle='Time from Si IV peak (s)', $
              ytitle='Log($\tau$)', yr=reverse(minmax(ltau)), cb=cb05, $
              min=0, max=100)
cb05.title='(n$_e$ - n$_{e, 0}$) / n$_{e, 0}$  (%)'
mask05 = mask_ltau(im05)

 
end