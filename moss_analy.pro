dir = '/Users/khcho/Desktop/IRIS-moss-main'
sav_files = file_search(dir, '*moss_event*.sav')
sp_num = 0

;0   1335.71   C II 1336
;1   1349.43   Fe XII 1349
;2   1351.17   1351
;3   1355.60   O I 1356
;4   1393.78   Si IV 1394
;5   1402.77   Si IV 1403
;6   2786.52   2786
;7   2796.20   Mg II k 2796
;8   2831.33   2831

id_list = ['C II 1336', 'Fe XII 1349', '1351', 'O I 1356', 'Si IV 1394', $
           'Si IV 1403', '2786', 'Mg II k 2796', '2831']


for j=0, n_elements(id_list)-1 do begin
  j = 7
  spec = !null
  real_ind  = !null
  for i=0, n_elements(sav_files)-1 do begin
;    i = 9
    restore, sav_files[i]
    sji_1400_ind = where(strmatch(eout.sji_files, '*1400*') eq 1)
    read_iris_l2, eout.sji_files[sji_1400_ind], sji_index, sji_data, /sil
    real_1400 = sji_data[where(sji_data ge 0)]
    sji_thres = 3.*stddev(real_1400) + mean(real_1400)
    moss_sel = eout.sji_value gt sji_thres 
  
    ind = (where(strmatch(eout.line_id, id_list[j]) eq 1, count))[0]
    if count eq 0 then continue
    if i eq 0 then begin
      wave = eout.wave_list[ind]
    endif
    data = eout.data_list[ind]
    if (size(data))[1] eq n_elements(wave) then begin
      dum = data
    endif else begin
      dum = congrid(data, n_elements(wave), (size(data))[2])
    endelse
    spec = [[spec], [dum]]
    real_ind = [real_ind, moss_sel]
;    stop
  endfor
  if n_elements(spec) eq 0 then continue
  img1 = spec
  img3 = rebin(transpose(real_ind), 10, n_elements(real_ind))*1d3
  img1[0:9, *] = img3
    
  w1 = window(dim=[8d2, 6d2] )
  im01 = image_kh(img1, wave, findgen(n_elements(real_ind)), $
                  aspect_ratio=0, /current, min=0, max=stddev(spec), $
                  xtitle='Wavelength ($\AA$)', ytitle='Pixel #', $
                  title='Moss Spectra at '+id_list[j]+' Window')
  im01.xtickdir = 1
  im01.ytickdir = 1
  stop
  w1.save, id_list[j]+'.png', resol=200
  w1.close

endfor
end