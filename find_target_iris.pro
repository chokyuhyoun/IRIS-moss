ar_num = [11800:13050]
tot_hours = !null
if 1 then begin
  for i=0, n_elements(ar_num)-1 do begin
    str = 'https://www.lmsal.com/hek/hcr?cmd=search-events3&outputformat=json&'+$
          'startTime=2013-07-20T00:00&stopTime=2022-07-15T00:00&target=AR&maxrasterStepsize=0&'+$
          'maxcadMeanAsrun=15&specWindows=Mg+II+k+2796,Mg+II+h+2803,Si+IV+1403&'+$
          'hasData=true&hideMostLimbScans=true&obsDesc='$
          +string(ar_num[i], f='(i0)')
    hcr=ssw_hcr_query(str)
    hours = 0.
    for j=0, n_elements(hcr)-1 do begin & hours += hcr[j].duration_hours & print, j, hcr[j].duration_hours & endfor
    tot_hours = [tot_hours, hours]
    print, i, n_elements(ar_num), ar_num[i], hours
  ;  stop
  endfor
  save, str, ar_num, tot_hours, filename='IRIS_ar_obs_info.sav'
endif else restore, 'IRIS_ar_obs_info.sav'
hcr = 0
order = reverse(sort(tot_hours))
for i=0, 20 do print, ar_num[order[i]], tot_hours[order[i]]
  
end