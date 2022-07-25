out_dir_base = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, out_dir_base
out_dir = file_search(out_dir_base, 'ar12524/*', /test_dir, /fully)
for i=0, n_elements(out_dir)-1 do begin
;  i = 4
  print, out_dir[i]+' processing.....'
  find_iris_spectra_moss, out_dir[i]
;  stop
endfor

;sav_files = file_search(dir, '*2016*/*event*.sav', /fully)
;make_aia_sji_movie, 'ar12524'
end

