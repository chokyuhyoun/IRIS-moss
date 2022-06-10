function data_range, data, sig=sig
  if n_elements(sig) eq 0 then sig=3
  real = data[where(data ge 0)]
  min = (mean(real) - sig*stddev(real)) > 0.
  max = (mean(real) + sig*stddev(real)) < max(real)
  return, [min, max]
end

dir = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, dir
sav_files = file_search(dir, '*_*/*event*.sav', /fully)

hfov = 50
side = 0.41
gap = 0.03
imx0 = 0.12
imy0 = 0.53
sig=5
;for i=1, n_elements(sav_files)-1 do begin
  i = 0
  print, sav_files[i]
  movie_name = strmid(sav_files[i], 0, strlen(sav_files[i])-3)+'mp45'
  t0 = systime(/sec)
  sav_dir = file_dirname(sav_files[i])+'/imgs'
  file_mkdir, sav_dir

  restore, sav_files[i], /relax
  aia94_file = eout.aia_files[where(strmatch(eout.aia_files, '*94*'))]
  aia193_file = eout.aia_files[where(strmatch(eout.aia_files, '*193*'))]
  aia1700_file = eout.aia_files[where(strmatch(eout.aia_files, '*1700*'))]
  aia171_file = eout.aia_files[where(strmatch(eout.aia_files, '*171*'))]
  aia211_file = eout.aia_files[where(strmatch(eout.aia_files, '*211*'))]
  aia1600_file = eout.aia_files[where(strmatch(eout.aia_files, '*1600*'))]
  
  read_iris_l2, aia94_file, aia94_index, aia94_data, /sil
  read_iris_l2, aia193_file, aia193_index, aia193_data, /sil
  read_iris_l2, aia1700_file, aia1700_index, aia1700_data, /sil
  read_iris_l2, aia171_file, aia171_index, aia171_data, /sil
  read_iris_l2, aia211_file, aia211_index, aia211_data, /sil
  read_iris_l2, aia1600_file, aia1600_index, /sil
  read_iris_l2, eout.sji_files[0], sji_index, sji_data, /sil
  fe_xviii_data = aia94_data - aia171_data/450. - aia211_data/120.
  
;  fe_xviii_dr = data_range(fe_xviii_data, sig=sig)
;  aia193_dr = data_range(aia193_data, sig=sig)
;  aia1700_dr = data_range(aia1700_data, sig=sig) 
   
  aia94_time = anytim(aia94_index.date_obs)
  aia1700_time = anytim(aia1700_index.date_obs)
  aia193_time = anytim(aia193_index.date_obs)
  aia1600_time = anytim(aia1600_index.date_obs)
  sji_time = anytim(sji_index.date_obs)
 
  aia_lct, r94, g94, b94, wavelnth=94, /load
  aia_lct, r193, g193, b193, wavelnth=193, /load
  aia_lct, r1700, g1700, b1700, wavelnth=1700, /load
  iris_lct, sji_index[0], rsji, gsji, bsji
  
  w1 = window(dim=[8d2, 8d2], buffer=1)
  im01 = image_kh(hanning(2, 2), /current, $
                  pos=[imx0, imy0, imx0+side, imy0+side], $
                  ytitle='Solar Y (arcsec)', $
                  xtickdir=1, ytickdir=1, xshowtext=0, $
                  rgb_table=[[r94], [g94], [b94]])                  
  im02 = image_kh(hanning(2, 2), /current, $
                  pos=im01.pos + [side+gap, 0, side+gap, 0], $
                  xtickdir=1, ytickdir=1, xshowtext=0, yshowtext=0, $
                  rgb_table=[[r193], [g193], [b193]])
  im03 = image_kh(hanning(2, 2), /current, $
                  pos=im01.pos - [0, side+gap, 0, side+gap], $
                  xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)', $
                  xtickdir=1, ytickdir=1, $
                  rgb_table=[[r1700], [g1700], [b1700]])
  im04 = image_kh(hanning(2, 2), /current, $
                  pos=im01.pos + [side+gap, -side-gap, side+gap, -side-gap], $
                  xtitle='Solar X (arcsec)', $
                  xtickdir=1, ytickdir=1, yshowtext=0, $
                  rgb_table=[[rsji], [gsji], [bsji]])
  im011 = image_kh(hanning(2, 2), over=im01, rgb_table=34, transp=0)
  im021 = image_kh(hanning(2, 2), over=im02, rgb_table=34, transp=0)
  im031 = image_kh(hanning(2, 2), over=im03, rgb_table=34, transp=0)
  im041 = image_kh(hanning(2, 2), over=im04, rgb_table=34, transp=0)
  p04 = plot([0, 1], [0, 1], over=im04, thick=2)
  t01 = text(im01.pos[0]+0.02, im01.pos[3]-0.02, '', /current, $
             font_size=12, font_name='malgun gothic', font_style=1, $
             vertical_align=1, font_color='white')
  t02 = text(im02.pos[0]+0.02, im02.pos[3]-0.02, '', /current, $
             font_size=12, font_name='malgun gothic', font_style=1, $
             vertical_align=1, font_color='white')
  t03 = text(im03.pos[0]+0.02, im03.pos[3]-0.02, '', /current, $
             font_size=12, font_name='malgun gothic', font_style=1, $
             vertical_align=1, font_color='white')
  t04 = text(im04.pos[0]+0.02, im04.pos[3]-0.02, '', /current, $
             font_size=12, font_name='malgun gothic', font_style=1, $
             vertical_align=1, font_color='white')             
  dum = min(abs(aia193_time-sji_time[0]), i_ind)
  dum = min(abs(aia193_time-sji_time[-1]), f_ind)

  video = idlffvideowrite(movie_name)
  framerate = 20
  wdims = w1.dimensions
  stream = video.addvideostream(wdims[0], wdims[1], framerate)
  percent = !null
  for j=i_ind, f_ind do begin
    dum = floor(10.*(j-i_ind)/(f_ind-i_ind))*10
    if dum ne percent then begin
      print, string(dum, f='(i3)')+' %'
      percent = dum
    endif
    dum = min(abs(aia94_time - aia193_time[j]), match94)
    dum = min(abs(sji_time - aia193_time[j]), matchsji)
    dum = min(abs(aia1700_time - aia193_time[j]), match1700)
    
    get_xp_yp, sji_index[matchsji], sji_xp, sji_yp
    get_xp_yp, aia94_index[match94], aia94_xp, aia94_yp
    get_xp_yp, aia193_index[j], aia193_xp, aia193_yp
    get_xp_yp, aia1700_index[match1700], aia1700_xp, aia1700_yp
    
;    aia94_xp -= poly(aia94_time[match94]-aia1600_time[0], eout.aia_shift[*, 0])
;    aia94_yp -= poly(aia94_time[match94]-aia1600_time[0], eout.aia_shift[*, 1])
;    aia193_xp -= poly(aia193_time[j]-aia1600_time[0], eout.aia_shift[*, 0])
;    aia193_yp -= poly(aia193_time[j]-aia1600_time[0], eout.aia_shift[*, 1])
;    aia1700_xp -= poly(aia1700_time[match1700]-aia1600_time[0], eout.aia_shift[*, 0])
;    aia1700_yp -= poly(aia1700_time[match1700]-aia1600_time[0], eout.aia_shift[*, 1])

    im01.xr = sji_index[matchsji].crval1 + hfov*[-1, 1]
    im01.yr = sji_index[matchsji].crval2 + hfov*[-1, 1]
    im02.xr = im01.xr
    im02.yr = im01.yr
    im03.xr = im01.xr
    im03.yr = im01.yr
    im04.xr = im01.xr
    im04.yr = im01.yr
    
    sji_img = iris_intscale(sji_data[*, *, matchsji], sji_index[matchsji])

;    im01.setdata, fe_xviii_data[*, *, match94], aia94_xp, aia94_yp
;    im02.setdata, aia193_data[*, *, j], aia193_xp, aia193_yp
;    im03.setdata, aia1700_data[*, *, match1700], aia1700_xp, aia1700_yp
;    im04.setdata, sji_img, sji_xp, sji_yp
    setdata_hi_res, im01, fe_xviii_data[*, *, match94], aia94_xp, aia94_yp
    setdata_hi_res, im02, aia193_data[*, *, j], aia193_xp, aia193_yp
    setdata_hi_res, im03, aia1700_data[*, *, match1700], aia1700_xp, aia1700_yp
    setdata_hi_res, im04, sji_img, sji_xp, sji_yp

    moss_img = float(eout.zcube[*, *, j])
    moss_img[where(moss_img eq 0)] = !values.f_nan
;    im011.setdata, moss_img, aia193_xp, aia193_yp
;    im021.setdata, moss_img, aia193_xp, aia193_yp
;    im031.setdata, moss_img, aia193_xp, aia193_yp
;    im041.setdata, moss_img, aia193_xp, aia193_yp
    setdata_hi_res, im011, moss_img, aia193_xp, aia193_yp
    setdata_hi_res, im021,moss_img, aia193_xp, aia193_yp
    setdata_hi_res, im031,moss_img, aia193_xp, aia193_yp
    setdata_hi_res, im041,moss_img, aia193_xp, aia193_yp

    p04.setdata, sji_xp[sji_index[matchsji].sltpx1ix*[1, 1]], im04.yr
        
    im01.min = 0
    im01.max = 70
    im02.min = 1e2
    im02.max = 7e3
    im03.min = 1e2
    im03.max = 3e3
    im011.min = 0
    im021.min = 0
    im031.min = 0
    im041.min = 0
    im011.max = 1
    im021.max = 1
    im031.max = 1
    im041.max = 1

    t01.string = 'Fe XVIII '+strmid(aia94_index[match94].date_obs, 0, 19)
    t02.string = 'AIA 193$\AA$ '+strmid(aia193_index[j].date_obs, 0, 19)
    t03.string = 'AIA 1700$\AA$ '+strmid(aia1700_index[match1700].date_obs, 0, 19)
    t04.string = 'SJI '+string(eout.sji_wave, f='(i4)')+'$\AA$ '$
                 + strmid(sji_index[matchsji].date_obs, 0, 19)
;    w1.save, sav_dir+path_sep()+string(j, f='(i04)')+'.png', resol=100
    timestamp = video.put(stream, w1.copywindow())
;    stop
  endfor
;  img_files = file_search(sav_dir, '*.png', /fully)
;  ffmpeg, img_files, 20, output=movie_name
  video.cleanup
  print, 'It took '+string((systime(/sec)-t0)/60., f='(f6.2)')+' min'
;endfor


end