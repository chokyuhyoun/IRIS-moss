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
;  j = 7
  spec = !null
  for i=0, n_elements(sav_files)-1 do begin
    restore, sav_files[i]
    read_iris_l2, eout.sg_files[0], index, /sil
    ind = (where(eout.line_id eq id_list[j], count))[0]
    if count eq 0 then continue
    if i eq 0 then begin
      wave = eout.wave_list[ind]
      cdelt2 = index[0].cdelt2
    endif
    data = eout.data_list[ind]
    if (size(data))[1] eq n_elements(wave) then begin
      dum = data
    endif else begin
      dum = congrid(data, n_elements(wave), (size(data))[2])
    endelse
    spec = [[spec], [dum/index[0].cdelt2*cdelt2]]

    print, index[0].obs_desc, index[0].cdelt1, index[0].cdelt2
;    stop
  endfor
  if n_elements(spec) eq 0 then continue
  w1 = window(dim=[8d2, 6d2] )
  im01 = image_kh(spec, wave, $
                  aspect_ratio=0, /current, min=0, max=3.*stddev(spec), $
                  xtitle='Wavelength ($\AA$)', ytitle='Pixel #', $
                  title='Moss Spectra at '+id_list[j]+' Window')
  im01.xtickdir = 1
  im01.ytickdir = 1

;  stop
  w1.save, dir+path_sep()+id_list[j]+'.png', resol=200
  w1.close

endfor
end