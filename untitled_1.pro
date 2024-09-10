dir = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, dir
sav_files0 = file_search(dir, '*moss_event*.sav', /fully)
sav_files = sav_files0[where(strmatch(sav_files0, '*ar*'), /null)]
;stop
for i=0, n_elements(sav_files)-1 do begin
  print, sav_files[i]
  restore, sav_files[i], /relax
  si_iv_ind = (where(strmatch(eout.line_id, '*Si IV*'), /null))[-1]
  si_wave = eout.wave_list[Si_IV_ind]
  print, mean(si_wave[1:*] - si_wave[0:-2], /nan)
  eout = 0
endfor
end