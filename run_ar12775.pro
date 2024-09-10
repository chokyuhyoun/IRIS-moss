hcr=ssw_hcr_query('https://www.lmsal.com/hek/hcr?cmd=search-events3&outputformat=json&startTime=2020-10-09T00:00&stopTime=2020-10-17T00:00&target=AR&minexpMax=4&maxrasterCadMeanAsrun=40&hasData=true&hideMostLimbScans=true&limit=200');sit_list = !null
sit_list = !null
for i=0, n_elements(hcr)-1 do begin
  if strmatch(hcr[i].goal, '*sit-and-stare*') eq 1 then sit_list = [sit_list, i]
endfor
hcr = hcr[sit_list]
;nhcr = hcr[sit_list]
;nhcr = hcr

out_dir_base = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, out_dir_base
out_dir0 = out_dir_base + 'ar12775'
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
;stop
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
end

