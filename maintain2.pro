out_dir_base = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, out_dir_base
sav_files = file_search(out_dir_base, 'ar12529/*/*moss_event*.sav', /fully)
for i=0, n_elements(sav_files)-1 do begin
  restore, sav_files[i]
  dd = iris_obj(eout.sg_files)
  sg_dy = dd->getresy()
  eout = add_tag(temporary(eout), sg_dy, 'sg_dy')
;  stop
  save, eout, filename=sav_files[i]
  obj_destroy, dd 
endfor

end