cd, '/Users/khcho/Desktop/IRIS-moss-main'

model_dir = file_search('/Users/khcho/Desktop/IRIS-moss-main/simulation/rh_output_kh/*', /test_dir)
model_name = ['E2', 'E2+', 'E3', 'E3+', 'E1', 'E1+', 'E6', 'H2', 'E4', 'H1', 'E5', 'C1', 'C1+']

for i=0, n_elements(model_dir)-1 do begin
  si_file = file_search(model_dir[i]+'/si_iv_fit_res.sav')
  mg_file = file_search(model_dir[i]+'/mg_ii_fit_res.sav')
  if si_file eq '' or mg_file eq '' then stop
  restore, si_file
  restore, mg_file
  dum = max(sim_si4_res.i_tot, si4_max)
  dum = max(sim_mg2_res[0, *].emiss, mg2_max)
  dir_name = (strsplit(model_dir[i], '/', /extract))[-1]
;  if model_name[i] eq 'E5' then stop
  if i eq 0 then begin
    sim_res0 = {mg_res:sim_mg2_res[*, mg2_max], $
                si_res:sim_si4_res[si4_max], $
                si_max_ind:si4_max, $
                model_name:model_name[i], $
                dir_name:dir_name}
    sim_res = replicate(sim_res0, n_elements(model_name))
  endif else begin
    sim_res[i].mg_res = sim_mg2_res[*, mg2_max]
    sim_res[i].si_res = sim_si4_res[si4_max]
    sim_res[i].si_max_ind = si4_max
    sim_res[i].model_name = model_name[i]
    sim_res[i].dir_name = dir_name
  endelse

endfor
save, sim_res, filename='/Users/khcho/Desktop/IRIS-moss-main/simulation/sim_res.sav'

end