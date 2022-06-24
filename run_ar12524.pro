out_dir_base = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, out_dir_base
out_dir = file_search(out_dir_base, '2016*', /test_dir, /fully)
for i=0, n_elements(out_dir)-1 do begin
;  i = 1
  print, out_dir[i]+' processing.....'
  find_iris_spectra_moss, out_dir[i]
endfor

end

