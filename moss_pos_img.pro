dir = '/Users/khcho/Desktop/IRIS-moss-main'
cd, dir
sav_files = file_search(dir, '*moss_event*.sav')
sp_num = 0
w2 = window()
p02 = plot(/test, /nodata, /current)

for i=0, n_elements(sav_files)-1 do begin
;  i = 2
  restore, sav_files[i]
  sp_num += n_elements(eout.sji_ind)
  sji_1400_ind = where(strmatch(eout.sji_files, '*1400*') eq 1)
  read_iris_l2, eout.sji_files[sji_1400_ind], sji_index, sji_data, /sil

  aia193_file = eout.aia_files[where(strmatch(eout.aia_files, '*193*'))]
  read_iris_l2, aia193_file, aia193_index, aia193_data, /sil
  
  aia_lct, rr0, gg0, bb0, wavelnth=193
  iris_lct, sji_index[0], rr2, gg2, bb2
  
  aia_indt = eout.aia_ind[2, *]
  indt = uniq(aia_indt)
  aia_num = aia_indt[indt]
  sji_num = eout.sji_ind[indt]
  
  for j=0, n_elements(aia_num)-1 do begin
;    j=2
    part = where(aia_indt eq aia_num[j])

    selx0 = eout.sg_phy[0, part]
    sely0 = eout.sg_phy[1, part]
    
    w1 = window(dim=[8d2, 4d2])
    xr = eout.sg_phy[0, indt[j]] + 10.*[-1, 1]
    yr = eout.sg_phy[1, indt[j]] + 10.*[-1, 1]
    get_xp_yp, aia193_index[aia_num[j]], aia_xp, aia_yp
    aia_xp -= eout.aia_shift[0, indt[j]]
    aia_yp -= eout.aia_shift[1, indt[j]]
    im01 = image_kh(aia193_data[*, *, aia_num[j]], aia_xp, aia_yp, $
                 /current, axis=2, pos=[0.12, 0.15, 0.47, 0.90], $
                 xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)', $
                 title='AIA 193 $\AA$ '+strmid(aia193_index[aia_num[j]].date_obs, 0, 19), $
                 xr=xr, yr=yr, xtickdir=1, ytickdir=1, $
                 rgb_table=[[rr0], [gg0], [bb0]], min=1e3, max=7e3)
    get_xp_yp, sji_index[eout.sji_ind[indt[j]]], sji_xp, sji_yp
    sji_img = sji_data[*, *, eout.sji_ind[indt[j]]]
    sji_img = iris_intscale(temporary(sji_img), sji_index[0])
    im02 = image_kh(sji_img, sji_xp, sji_yp, /current, $
                 axis=2, pos=[0.62, 0.15, 0.97, 0.9], $
                 xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)', $ 
                 title='SJI 1400 $\AA$ '+strmid(sji_index[eout.sji_ind[indt[j]]].date_obs, 0, 19), $
                 xr=xr, yr=yr, xtickdir=1, ytickdir=1, $
                 rgb_table=[[rr2], [gg2], [bb2]])  
    zcube0 = float(eout.zcube[*, *, aia_num[j]])
    zcube0[where(zcube0 eq 0)] = !values.f_nan           
    im011 = image_kh(zcube0, aia_xp, aia_yp, over=im01, $
                  rgb_table=30, transp=30)
    im021 = image_kh(zcube0, aia_xp, aia_yp, over=im02, $
                  rgb_table=30, transp=30)

    p01 = plot(selx0, sely0, 'o', color='y', sym_thick=2, transp=50, over=im02)
    t01 = text(im01.pos[0]+0.03, im01.pos[3]-0.09, 'Moss Positions', $
               color=im011.rgb_table[*, 0], transp=im011.transp, $
               font_name='margun gothic', font_size=15, font_style=1)
    t02 = text(im02.pos[0]+0.03, im02.pos[3]-0.09, 'Spectral Positions', $
               color=p01.color, transp=p01.transp, $
               font_name='margun gothic', font_size=15, font_style=1)           
;    stop                 
    w1.save, string(i, f='(i0)')+string(j, f='(i0)')+'.png', resol=200
    w1.close
  endfor
  for k=0, n_elements(eout.sji_ind)-1 do p03 = plot(eout.sji_curve[*, k], over=p02)
  
endfor

end