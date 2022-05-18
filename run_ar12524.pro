;hcr=ssw_hcr_query('https://www.lmsal.com/hek/hcr?cmd=search-events3&outputformat=json&startTime=2013-07-20T00:00&stopTime=2022-05-11T00:00&target=AR&hasData=true&obsDesc=12524&limit=200')
;sit_list = !null
;for i=0, n_elements(hcr)-1 do begin
;  if strmatch(hcr[i].goal, '*sit-and-stare*') eq 1 then sit_list = [sit_list, i]
;endfor
;nhcr = hcr[sit_list]

out_dir_base = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, out_dir_base
out_dir = file_search(out_dir_base, '2016*', /test_dir, /fully)
for i=0, n_elements(out_dir)-1 do begin
;  i = 0
;  umodes = (strsplit(nhcr[i].umodes, ',', /extract))[0]
;  strput, umodes, '.', strlen(umodes)-5
;  dir0 = file_dirname(umodes)
;  dum = strsplit(dir0, '/', /ext)
;  dum[2] = 'level2_compressed'
;  dir = '/'+strjoin(dum, '/')
;  files = file_search(dir, '*.gz', /fully)
;  out_dir_sub = strmid(dum[-1], 0, 15) 
;  out_dir = out_dir_base + out_dir_sub + '/'
;  file_mkdir, out_dir
;  nfiles = file_search(out_dir, '*.fits', /fully, count=n)
;  if n eq 0 then begin
;    for j=0, n_elements(files)-1 do begin
;      if (strmatch(files[j], '*tar*') eq 1) then begin
;        spawn, 'tar -xvzf ' + files[j] + ' -C ' + out_dir
;      endif else begin
;;        stop
;        name = file_basename(files[j])
;        spawn, 'gzip -dc ' + files[j] + ' >> ' $
;                + out_dir + strmid(name, 0, strlen(name)-3)
;      endelse
;    endfor
;  endif
;  print, out_dir_sub+' processing.....'
  find_iris_spectra_moss, out_dir[i]
endfor

end

