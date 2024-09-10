pro run_ar_num, ar_num, movie=movie
  if ~keyword_set(movie) then movie = 0
  str = 'https://www.lmsal.com/hek/hcr?cmd=search-events3&outputformat=json&'+$
        'startTime=2013-07-20T00:00&stopTime=2022-07-15T00:00&target=AR&maxrasterStepsize=0&'+$
        'maxcadMeanAsrun=15&specWindows=Mg+II+k+2796,Mg+II+h+2803,Si+IV+1403&'+$
        'hasData=true&hideMostLimbScans=true&obsDesc='$
        +string(ar_num, f='(i0)')    
  hcr=ssw_hcr_query(str)
  
  out_dir_base = '/Users/khcho/Desktop/IRIS-moss-main/'
  cd, out_dir_base
  out_dir0 = out_dir_base + 'ar'+string(ar_num, f='(i0)')
  file_mkdir, out_dir0
  iris_dir = !null
  
  for i=0, n_elements(hcr)-1 do begin
  ;  i = 6
    umodes = (strsplit(hcr[i].umodes, ',', /extract))[0]
    strput, umodes, '.', strlen(umodes)-5
    iris_dir = [iris_dir, file_dirname(umodes)]
  endfor  
  order = sort(iris_dir)
  iris_dir = iris_dir[order]
  nhcr = hcr[order]
;  stop
  for i=0, n_elements(iris_dir)-1 do begin
  ;for i=2, 2 do begin
    aia_dir = iris_dir[i]+'/aia'  
    print, iris_dir[i]+' is processing.....'
    out_dir = out_dir0+path_sep()+strmid(iris_dir[i], 30, 15)
    file_mkdir, out_dir
    find_iris_spectra_moss, out_dir, iris_dir=iris_dir[i], aia_dir=aia_dir
  endfor
  obj_destroy, nhcr
  obj_destroy, hcr
  
  if movie then make_aia_sji_movie, 'ar'+string(ar_num, f='(i0)'), /moss_only
end

run_ar_num, 12415
end