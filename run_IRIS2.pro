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

end_ind = uniq(sg_files)
start_ind = [0, end_ind[0:-2]+1]
if 0 then begin
  for i=0, n_elements(end_ind)-1 do begin
    sg_files0 = sg_files[end_ind[i]]
    iris2model = iris2(sg_files0)
  endfor
endif

temp = !null
v_los = !null
v_turb = !null
n_e = !null

u_temp = !null
u_v_los = !null
u_v_turb = !null
u_n_e = !null
chi2 = !null

pre_temp = !null
pre_v_los = !null
pre_v_turb = !null
pre_n_e = !null

pre_u_temp = !null
pre_u_v_los = !null
pre_u_v_turb = !null
pre_u_n_e = !null
pre_chi2 = !null

aft_temp = !null
aft_v_los = !null
aft_v_turb = !null
aft_n_e = !null

aft_u_temp = !null
aft_u_v_los = !null
aft_u_v_turb = !null
aft_u_n_e = !null
aft_chi2 = !null


iris2_files = file_search(iris2_path, 'iris2model_*.sav')
for i=0, n_elements(end_ind)-1 do begin
  restore, iris2_files[i] ; iris2model
  for j=start_ind[i], end_ind[i] do begin
    if sg_files[j] ne iris2model.input_filename then stop    
    temp = [[temp], [reform(iris2model.model[sg_ind[1, j], sg_ind[0, j], *, 0])]]
    v_los = [[v_los], [reform(iris2model.model[sg_ind[1, j], sg_ind[0, j], *, 1])]]
    v_turb = [[v_turb], [reform(iris2model.model[sg_ind[1, j], sg_ind[0, j], *, 2])]]
    n_e = [[n_e], [reform(iris2model.model[sg_ind[1, j], sg_ind[0, j], *, 3])]]

    u_temp = [[u_temp], [reform(iris2model.uncertainty[sg_ind[1, j], sg_ind[0, j], *, 0])]]
    u_v_los = [[u_v_los], [reform(iris2model.uncertainty[sg_ind[1, j], sg_ind[0, j], *, 1])]]
    u_v_turb = [[u_v_turb], [reform(iris2model.uncertainty[sg_ind[1, j], sg_ind[0, j], *, 2])]]
    u_n_e = [[u_n_e], [reform(iris2model.uncertainty[sg_ind[1, j], sg_ind[0, j], *, 3])]]

    chi2 = [chi2, reform(iris2model.chi2[sg_ind[1, j], sg_ind[0, j]])]

;    if total(iris2model.chi2[pre_sg_ind[1, j], pre_sg_ind[0, j]]) eq 0 then $
;      pre_sg_ind[1, j] -= 1
    pre_temp = [[pre_temp], [reform(iris2model.model[sg_ind[1, j]-1, sg_ind[0, j], *, 0])]]
    pre_v_los = [[pre_v_los], [reform(iris2model.model[sg_ind[1, j]-1, sg_ind[0, j], *, 1])]]
    pre_v_turb = [[pre_v_turb], [reform(iris2model.model[sg_ind[1, j]-1, sg_ind[0, j], *, 2])]]
    pre_n_e = [[pre_n_e], [reform(iris2model.model[sg_ind[1, j]-1, sg_ind[0, j], *, 3])]]

    pre_u_temp = [[pre_u_temp], [reform(iris2model.uncertainty[sg_ind[1, j]-1, sg_ind[0, j], *, 0])]]
    pre_u_v_los = [[pre_u_v_los], [reform(iris2model.uncertainty[sg_ind[1, j]-1, sg_ind[0, j], *, 1])]]
    pre_u_v_turb = [[pre_u_v_turb], [reform(iris2model.uncertainty[sg_ind[1, j]-1, sg_ind[0, j], *, 2])]]
    pre_u_n_e = [[pre_u_n_e], [reform(iris2model.uncertainty[sg_ind[1, j]-1, sg_ind[0, j], *, 3])]]
    pre_chi2 = [pre_chi2, reform(iris2model.chi2[sg_ind[1, j]-1, sg_ind[0, j]])]

    aft_temp = [[aft_temp], [reform(iris2model.model[sg_ind[1, j]+1, sg_ind[0, j], *, 0])]]
    aft_v_los = [[aft_v_los], [reform(iris2model.model[sg_ind[1, j]+1, sg_ind[0, j], *, 1])]]
    aft_v_turb = [[aft_v_turb], [reform(iris2model.model[sg_ind[1, j]+1, sg_ind[0, j], *, 2])]]
    aft_n_e = [[aft_n_e], [reform(iris2model.model[sg_ind[1, j]+1, sg_ind[0, j], *, 3])]]

    aft_u_temp = [[aft_u_temp], [reform(iris2model.uncertainty[sg_ind[1, j]+1, sg_ind[0, j], *, 0])]]
    aft_u_v_los = [[aft_u_v_los], [reform(iris2model.uncertainty[sg_ind[1, j]+1, sg_ind[0, j], *, 1])]]
    aft_u_v_turb = [[aft_u_v_turb], [reform(iris2model.uncertainty[sg_ind[1, j]+1, sg_ind[0, j], *, 2])]]
    aft_u_n_e = [[aft_u_n_e], [reform(iris2model.uncertainty[sg_ind[1, j]+1, sg_ind[0, j], *, 3])]]
    aft_chi2 = [aft_chi2, reform(iris2model.chi2[sg_ind[1, j]+1, sg_ind[0, j]])]
  endfor
endfor
ltau = iris2model.ltau

save, temp, v_los, v_turb, n_e, $
      u_temp, u_v_los, u_v_turb, u_n_e, chi2, ltau, $
      pre_temp, pre_v_los, pre_v_turb, pre_n_e, $
      pre_u_temp, pre_u_v_los, pre_u_v_turb, pre_u_n_e, pre_chi2, $
      aft_temp, aft_v_los, aft_v_turb, aft_n_e, $
      aft_u_temp, aft_u_v_los, aft_u_v_turb, aft_u_n_e, aft_chi2, $
      pars_pix_peak, par_names, var_names, $
      filename=iris2_path+'/collection.sav'
 
end