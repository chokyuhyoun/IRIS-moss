pro iris_fit_maintain, sub_dir, show=show 

dir = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, dir
sav_files0 = file_search(dir, '*moss_event*.sav', /fully)
sav_files = sav_files0[where(strmatch(sav_files0, '*'+sub_dir+'*'), /null)]

si_cen = 1402.77d
w_th_si = si_cen/3d8*sqrt(8.*alog(2.)*1.38d-23*10d0^(4.9)/(28.0855*1.6605d-27))  ; in angstrom
w_inst = 0.026 ; in angstrom
ooe_fac = 1./(2.*sqrt(alog(2))) ; 1/e factor

;stop
for i=0, n_elements(sav_files)-1 do begin
  print, sav_files[i]
  restore, sav_files[i], /relax
  si_iv_ind = (where(strmatch(eout.line_id, '*Si IV*'), /null))[-1]
  mg_ii_ind = (where(strmatch(eout.line_id, '*Mg II k 2796*'), exist_mg_ii))[0]
  
  si_array = !null
  mg_array = !null
  for j=0, eout.moss_num-1 do begin
    print, j
;    j=27
;stop
;    si_iv_fit_res0 = si_iv_fit(eout.wave_list[Si_IV_ind], $
;                              (eout.data_list[si_iv_ind])[*, j]/eout.spec_bin/eout.spat_bin, $
;                               spec_bin = eout.spec_bin) 
;    eout.si_iv_fit_res[j] = si_iv_fit_res0

    mg_ii_fit_res0 = mg_ii_fit(eout.wave_list[mg_ii_ind], (eout.data_list[mg_ii_ind])[*, j])
    eout.mg_ii_fit_res[*, j] = mg_ii_fit_res0
    mg_ii_fit_res1 = mg_ii_fit(eout.wave_list[mg_ii_ind], (eout.pre_data_list[mg_ii_ind])[*, j])
    eout.pre_mg_ii_fit_res[*, j] = mg_ii_fit_res1
    mg_ii_fit_res2 = mg_ii_fit(eout.wave_list[mg_ii_ind], (eout.aft_data_list[mg_ii_ind])[*, j])
    eout.aft_mg_ii_fit_res[*, j] = mg_ii_fit_res2
;;
  endfor
;  eout0 = rem_tag(eout, 'si_iv_fit_res')
;  eout = add_tag(eout0, si_array, 'si_iv_fit_res')
  save, eout, filename = sav_files[i]
endfor
eout = 0
eout0 = 0

;make_aia_sji_movie, sub_dir, show=show, /moss_only
end

iris_fit_maintain, 'ar'
end