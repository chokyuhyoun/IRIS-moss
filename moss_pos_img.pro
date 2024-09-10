dir = '/Users/khcho/Desktop/IRIS-moss-main'
cd, dir
sav_files = file_search(dir, '*moss_event*.sav')
sp_num = 0
;w2 = window()
;p02 = plot(/test, /nodata, /current)

for i=0, n_elements(sav_files)-1 do begin
;    i = 1
  restore, sav_files[i], /relax
  sp_num += n_elements(eout.sji_fwhm)
  sji_1400_ind = where(strmatch(eout.sji_files, '*1400*') eq 1)
  read_iris_l2, eout.sji_files[sji_1400_ind], sji_index, sji_data, /sil

  aia193_file = eout.aia_files[where(strmatch(eout.aia_files, '*193*'))]
  read_iris_l2, aia193_file, aia193_index, aia193_data, /sil

  aia_lct, rr0, gg0, bb0, wavelnth=193
  iris_lct, sji_index[0], rr2, gg2, bb2

  aia_indt = eout.aia_ind[2, *]
  indt = uniq(aia_indt)
  aia_num = aia_indt[indt]
  sji_num = eout.sji_ind[2, indt]
  
  wv1_ind = (where(strmatch(eout.line_id, '*1336*')))[0]
  wv2_ind = (where(strmatch(eout.line_id, '*1403*')))[0]
  wv3_ind = (where(strmatch(eout.line_id, '*2796*')))[0]

  for j=0, n_elements(aia_num)-1 do begin
;        j=0
    aia_num_j = indt[j]
    part = where(aia_indt eq aia_num[j])

    selx0 = eout.sg_phy[0, part]
    sely0 = eout.sg_phy[1, part]

    w1 = window(dim=[12d2, 8d2], buffer=1)
    xr = eout.sg_phy[0, aia_num_j] + 10.*[-1, 1]
    yr = eout.sg_phy[1, aia_num_j] + 10.*[-1, 1]
    get_xp_yp, aia193_index[aia_num[j]], aia_xp, aia_yp
    aia_xp -= eout.aia_shift[0, aia_num_j]
    aia_yp -= eout.aia_shift[1, aia_num_j]
    im01 = image_kh(aia193_data[*, *, aia_num[j]], aia_xp, aia_yp, $
                    /current, pos=[0.08, 0.56, 0.31, 0.92], $
                    xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)', $
                    title='AIA 193 $\AA$ '+strmid(aia193_index[aia_num[j]].date_obs, 0, 19), $
                    xr=xr, yr=yr, xtickdir=1, ytickdir=1, $
                    rgb_table=[[rr0], [gg0], [bb0]], min=1e3, max=7e3)
    get_xp_yp, sji_index[eout.sji_ind[2, aia_num_j]], sji_xp, sji_yp
    sji_img = sji_data[*, *, eout.sji_ind[2, aia_num_j]]
    sji_img = iris_intscale(temporary(sji_img), sji_index[0])
    im02 = image_kh(sji_img, sji_xp, sji_yp, /current, $
                    pos=im01.pos+[0.33, 0, 0.33, 0], $
                    xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)', $
                    title='SJI 1400 $\AA$ '+strmid(sji_index[eout.sji_ind[2, aia_num_j]].date_obs, 0, 19), $
                    xr=xr, yr=yr, xtickdir=1, ytickdir=1, $
                    rgb_table=[[rr2], [gg2], [bb2]])
    zcube0 = float(eout.zcube[*, *, aia_num[j]])
    zcube0[where(zcube0 eq 0)] = !values.f_nan
    im011 = image_kh(zcube0, aia_xp, aia_yp, over=im01, $
                    rgb_table=30, transp=30)
    im021 = image_kh(zcube0, aia_xp, aia_yp, over=im02, $
                    rgb_table=30, transp=30)

    p02 = symbol(selx0, sely0, 'o', sym_color='y', sym_size=0.7, sym_thick=2, $
                 sym_transp=50, target=im02, /data)
    t01 = text(im01.pos[0]+0.01, im01.pos[3]-0.03, 'Moss Positions', $
                color=im011.rgb_table[*, 0], transp=im011.transp, $
                font_name='margun gothic', font_size=13, font_style=1)
    t02 = text(im02.pos[0]+0.01, im02.pos[3]-0.03, $
              'Spectral Positions (n='+string(n_elements(part), f='(i2)')+')', $
                color=p02.sym_color, transp=p02.sym_transp, $
                font_name='margun gothic', font_size=13, font_style=1)
    p03 = plot(/test, /nodata, pos=im01.pos+[0.66, 0, 0.66, 0], /current, $
               xr=[-70, 70], xtitle='Time from moss detection (s)', $
               title='SJI 1400 $\AA$ Light Curve', $
               font_name='margun gothic', font_size=13, font_style=1)
    p04 = plot(/test, /nodata, pos=im01.pos+[0, -0.5, 0, -0.5], /current, $
               xr=minmax(eout.wave_list[wv1_ind]), xtitle='Wavelength ($\AA$)', $
               title=eout.line_id[wv1_ind], $
                font_name='margun gothic', font_size=13, font_style=1)
    p05 = plot(/test, /nodata, pos=im02.pos+[0, -0.5, 0, -0.5], /current, $
               xr=minmax(eout.wave_list[wv2_ind]), xtitle='Wavelength ($\AA$)', $
               title=eout.line_id[wv2_ind],  $
                font_name='margun gothic', font_size=13, font_style=1)
    p06 = plot(/test, /nodata, pos=p03.pos+[0, -0.5, 0, -0.5], /current, $
               xr=minmax(eout.wave_list[wv3_ind]), xtitle='Wavelength ($\AA$)', $
               title=eout.line_id[wv3_ind], $
                font_name='margun gothic', font_size=13, font_style=1)
    sg_obs = !null
    for k=0, n_elements(part)-1 do begin
      part_k = part[k]
      p031 = plot((eout.sji_curve_ind[*, part_k] - eout.sji_ind[2, part_k])*eout.sji_dt $
                  +eout.sji_phy[part_k]-eout.aia_phy[2, part_k], $
                  eout.sji_curve[*, part_k], over=p03)
      sg_obs0 = eout.sg_phy[2, part_k] - eout.aia_phy[2, part_k]
      sg_obs = [sg_obs, sg_obs0]         
                  
      p041 = plot(eout.wave_list[wv1_ind], (eout.data_list[wv1_ind])[*, part_k], over=p04)
      p051 = plot(eout.wave_list[wv2_ind], (eout.data_list[wv2_ind])[*, part_k], over=p05)
      p061 = plot(eout.wave_list[wv3_ind], (eout.data_list[wv3_ind])[*, part_k], over=p06)
    endfor
    p04.yr = [-10, p04.yr[1]]
    p05.yr = [-10, p05.yr[1]]
    p06.yr = [-10, p06.yr[1]]
    p032 = plot([0, 0], p03.yr, $
                '--', over=p03, color=(im011.rgb_table)[*, 0], transp=im011.transp)
    p033 = plot(sg_obs, (p03.yr[1]*0.9)*(fltarr(n_elements(sg_obs))+1), 'r|', over=p03)
    p060 = plot(2798.8*[1, 1], p06.yr, '--', over=p06)
;    stop
    w1.save, string(i, f='(i02)')+string(j, f='(i02)')+'.png', resol=200
    w1.close
  endfor

endfor

end