dir = '/Users/khcho/Desktop/IRIS-moss-main'
sav_files = file_search(dir, '*moss_event*.sav')
sp_num = 0



for i=0, n_elements(sav_files)-1 do begin
;  i = 9
  restore, sav_files[i]
  sji_1400_ind = where(strmatch(eout.sji_files, '*1400*') eq 1)
  read_iris_l2, eout.sji_files[sji_1400_ind], sji_index, sji_data, /sil
  real_1400 = sji_data[where(sji_data ge 0)]
  sji_thres = 3.*stddev(real_1400) + mean(real_1400)
  moss_sel = eout.sji_value gt sji_thres 

  if total(moss_sel) eq 0 then continue
  aia193_file = eout.aia_files[where(strmatch(eout.aia_files, '*193*'))]
  read_iris_l2, aia193_file, aia193_index, aia193_data, /sil
  
  aia_lct, rr0, gg0, bb0, wavelnth=193
  iris_lct, sji_index[0], rr2, gg2, bb2
  
  aia_indt = eout.aia_ind[2, *]
  indt = uniq(aia_indt)
  aia_num = aia_indt[indt]
  sji_num = eout.sji_ind[uniq(aia_indt)]
  
  for j=0, n_elements(aia_num)-1 do begin
    part = where(aia_indt eq aia_indt[indt[j]])
    if total(moss_sel[part]) eq 0 then continue
    part0 = where((aia_indt eq aia_indt[indt[j]]) and (moss_sel eq 0))
    part1 = where((aia_indt eq aia_indt[indt[j]]) and (moss_sel eq 1))

    selx0 = eout.sg_phy[0, part0]
    sely0 = eout.sg_phy[1, part0]
    selx1 = eout.sg_phy[0, part1]
    sely1 = eout.sg_phy[1, part1]
    
    w1 = window(dim=[8d2, 4d2])
    xr = eout.sg_phy[0, indt[j]] + 10.*[-1, 1]
    yr = eout.sg_phy[1, indt[j]] + 10.*[-1, 1]
    get_xp_yp, aia193_index[aia_indt[j]], aia_xp, aia_yp
    aia_xp -= eout.aia_shift[0, indt[j]]
    aia_yp -= eout.aia_shift[1, indt[j]]
    im01 = image(aia193_data[*, *, aia_indt[j]], aia_xp, aia_yp, $
                 /current, axis=2, pos=[0.1, 0.2, 0.5, 0.9], $
                 xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)', $
                 title='AIA 193 $\AA$ '+aia193_index[aia_indt[j]].date_obs, $
                 xr=xr, yr=yr, $
                 rgb_table=[[rr0], [gg0], [bb0]], min=0, max=7e3)
    get_xp_yp, sji_index[eout.sji_ind[indt[j]]], sji_xp, sji_yp
    im02 = image(sji_data[*, *, eout.sji_ind[indt[j]]], sji_xp, sji_yp, /current, $
                 axis=2, pos=[0.5, 0.2, 0.9, 0.9], $
                 xtitle='Solar X (arcsec)', $ 
                 title='SJI 1400 $\AA$ '+sji_index[eout.sji_ind[indt[j]]].date_obs, $
                 xr=xr, yr=yr, rgb_table=[[rr2], [gg2], [bb2]], min=0, max=5e2) 
    im02.axes[1].showtext = 0  
    zcube0 = float(eout.zcube[*, *, aia_indt[j]])
    zcube0[where(zcube0 eq 0)] = !values.f_nan           
    im011 = image(zcube0, aia_xp, aia_yp, over=im01, $
                  rgb_table=30, transp=30)
    im021 = image(zcube0, aia_xp, aia_yp, over=im02, $
                  rgb_table=30, min=0, max=254, transp=30)

    
    p01 = plot(selx0, sely0, 'o', color='k', transp=80, over=im02)
    p02 = plot(selx1, sely1, 'o', color='k', transp=30, over=im02)               
    stop                 
    w1.save, string(i, f='(i0)')+string(j, f='(i0)')+'.png', resol=200
    w1.close
  endfor
  
  
  
endfor

end