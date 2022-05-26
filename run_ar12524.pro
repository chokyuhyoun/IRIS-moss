ar_no = '11890' ;'12524'
hcr=ssw_hcr_query('https://www.lmsal.com/hek/hcr?cmd=search-events3&outputformat=json&startTime=2013-07-20T00:00&stopTime=2022-05-11T00:00&target=AR&hasData=true&obsDesc='+ar_no+'&limit=200')
sit_list = !null
rep = !null
for i=0, n_elements(hcr)-1 do begin
;  if strmatch(hcr[i].goal, '*sit-and-stare*') eq 1 then sit_list = [sit_list, i]
  rep = [rep, hcr[i].iris_repeats]
endfor
nhcr = hcr[where(rep gt 10)]
stop  

dir = '/Users/khcho/Desktop/IRIS-moss-main'
nhcr = file_search(dir, '201603*', /test_dir, /fully)
cd, dir
out_dir = dir + '/output' + ar_no
file_mkdir, out_dir
for i=0, n_elements(nhcr)-1 do begin
;  i = 0
;  umodes = (strsplit(nhcr[i].umodes, ',', /extract))[0]
;  strput, umodes, '.', strlen(umodes)-5
;  iris_dir = file_dirname(umodes)
;  aia_dir = iris_dir + '/aia'
;  print, nhcr[i].starttime+' processing.....'
  
  iris_dir = nhcr[i]
  aia_dir = iris_dir
  out_dir = iris_dir
  
  find_iris_spectra_moss, out_dir, iris_dir=iris_dir, aia_dir=aia_dir
endfor

end

