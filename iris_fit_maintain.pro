pro iris_fit_maintain, sub_dir, show=show 

dir = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, dir
sav_files0 = file_search(dir, '*moss_event*.sav', /fully)
sav_files = sav_files0[where(strmatch(sav_files0, '*'+sub_dir+'*'), /null)]
;stop
for i=0, n_elements(sav_files)-1 do begin
  print, sav_files[i]
  restore, sav_files[i], /relax
  si_iv_ind = (where(strmatch(eout.line_id, '*Si IV*'), /null))[-1]
  mg_ii_ind = (where(strmatch(eout.line_id, '*Mg II k 2796*'), exist_mg_ii))[0]
  
  si_array = !null
  for j=0, eout.moss_num-1 do begin
;    j=27
;stop
    si_iv_fit_res0 = si_iv_fit(eout.wave_list[Si_IV_ind], (eout.data_list[si_iv_ind])[*, j])
    si_array = [si_array, si_iv_fit_res0]

    mg_ii_fit_res0 = mg_ii_fit(eout.wave_list[mg_ii_ind], (eout.data_list[mg_ii_ind])[*, j])
    eout.mg_ii_fit_res[*, j] = mg_ii_fit_res0

  endfor
  eout0 = rem_tag(eout, 'si_iv_fit_res')
  eout = add_tag(eout0, si_array, 'si_iv_fit_res')
  save, eout, filename = sav_files[i]
endfor
eout = 0
eout0 = 0

;make_aia_sji_movie, sub_dir, show=show, /moss_only
end

