;; Run IRIS^2
;

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

sg_files = sg_file[pars_pix_peak_ind]
sg_ind = sg_ind[*, pars_pix_peak_ind] ;[y_pix, time in a dataset, dataset]
pre_sg_ind = pre_sg_ind[*, pars_pix_peak_ind] ;[y_pix, time in a dataset, dataset]
curve_dt = curve_dt[pars_pix_peak_ind]

end_ind = uniq(sg_files)
start_ind = [0, end_ind[0:-2]+1]
if 0 then begin ; Run IRIS^2. Automatically save the results in "iris2model_*.sav". 
  for i=0, n_elements(end_ind)-1 do begin
    sg_files0 = sg_files[end_ind[i]]
    iris2model = iris2(sg_files0)  ; automatically save the data
  endfor
endif

iris2_files = file_search(iris2_path, 'iris2model_*.sav')

dur = 120. ; in sec
dt = 0.1 ; in sec
nt = dur/dt + 1
t_arr = (findgen(nt) - 0.5*(nt-1))*dt

iris2_moss = []
for j=0, n_elements(iris2_files)-1 do begin
  restore, iris2_files[j]  ; iris2model
  
  ; iris2model.model = [time, y, ltau, prameters]
  ; parameters : T [k],  v_los [cm/s],  v_turb [cm/s],  n_e [cm^-3]
  original_model = iris2model.model
  ; find if model has 0 temperature at some locations. = failed  
  min_temp = min(original_model[*, *, *, 0], dim=3)
  loc = array_indices(min_temp, where(min_temp eq 0))
  failed_ind_t = reform(loc[0, *])
  failed_ind_y = reform(loc[1, *])
  
  modified_model = original_model
  for ii=0, n_elements(failed_ind_t)-1 do begin
    modified_model[failed_ind_t[ii], failed_ind_y[ii], *, *] = !values.f_nan
  endfor
  modified_uncertainty = iris2model.uncertainty
  modified_uncertainty[where(~finite(modified_model))] = !values.f_nan
  mean_model = mean(mean(modified_model, dim=1, /nan), dim=1, /nan)
  for i=start_ind[j], end_ind[j] do begin 
    print, i;, iris2_files[j]  
    aux_data = readfits(iris2model.input_filename, hdr, ext=iris2model.header.nwin+1, /sil)
    time = reform(aux_data[fxpar(hdr, 'time'), *])
    if sg_files[i] ne iris2model.input_filename then stop
    ref_tind = sg_ind[1, i]
    ref_time = time[ref_tind]
    rel_time = time - ref_time
    nobs = n_elements(time)
    
    rel_time2 = rebin(rel_time, nobs, nt)
    t_arr2 = rebin(transpose(t_arr), nobs, nt)
    dum = min(abs(t_arr2 - rel_time2), t_index0, dim=1) 
    t_index = reform((array_indices(t_arr2, t_index0))[0, *])
    ind_edge = uniq(t_index)
    av_ind_cad = median(ind_edge[1:*] - ind_edge[0:-2])
    t_start_ind = (ind_edge[0] le av_ind_cad) ? 0 : (ind_edge[0]-av_ind_cad)
    t_end_ind = ((nt - ind_edge[-2]) le av_ind_cad) ? nt-1 : ind_edge[-2]+av_ind_cad
    
    model = reform(modified_model[t_index, sg_ind[0, i], *, *])
    model[0:t_start_ind, *, *] = !values.f_nan  
    model[t_end_ind:-1, *, *] = !values.f_nan
  
  
    uncert =  reform(modified_uncertainty[t_index, sg_ind[0, i], *, *])
    uncert[0:t_start_ind, *, *] = !values.f_nan
    uncert[t_end_ind:-1, *, *] = !values.f_nan
   
    chi2 = iris2model.chi2[t_index, sg_ind[0, i]]
    chi2[0:t_start_ind, *, *] = !values.f_nan
    chi2[t_end_ind:-1, *, *] = !values.f_nan
;stop  

; for test
if 0 then begin
p01 = plot(rel_time, modified_model[*, sg_ind[0, i], 20, 0], '-+', xr=[-100, 100])
p02 = plot(t_arr, model[*, 20, 0], 'r', over=p01)
endif

    if i eq 0 then begin  
      iris2_moss0 = {input_filename:iris2model.input_filename, $
                    output_filename:iris2model.output_filename, $
                    ltau:iris2model.ltau, $
                    model:model, $
                    model_dimension_info:iris2model.model_dimension_info, $
                    physical_variable:iris2model.physical_variable, $
                    uncertainty:uncert, $
                    chi2:chi2, $
                    mu:iris2model.mu, $
  ;                  header:iris2model.header, $
                    mean_model:mean_model, $
                    peak_ind:fix(0.5*(nt-1)), $
                    dt:dt, $
                    time_from_peak:t_arr}
      iris2_moss = replicate(iris2_moss0, n_elements(sg_files))
    endif else begin
      iris2_moss[i].input_filename = iris2model.input_filename
      iris2_moss[i].output_filename = iris2model.output_filename
      iris2_moss[i].ltau = iris2model.ltau
      iris2_moss[i].model = model
      iris2_moss[i].model_dimension_info = iris2model.model_dimension_info
      iris2_moss[i].physical_variable = iris2model.physical_variable
      iris2_moss[i].uncertainty = uncert
      iris2_moss[i].chi2 = chi2
      iris2_moss[i].mu = iris2model.mu
      iris2_moss[i].mean_model = mean_model
  ;    iris2_moss[i].header = iris2model.header
    endelse
  endfor
endfor

save, iris2_moss, filename=iris2_path+'/iris2_moss.sav'
 
end