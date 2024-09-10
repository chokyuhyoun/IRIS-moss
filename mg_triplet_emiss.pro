dir = '/Users/khcho/Desktop/IRIS-moss-main'
cd, dir

sav_files = file_search(dir, '*moss_event*.sav')
restore, sav_files[3], /relax

ny = 9
nt = 5
;sg_indy = 213 + indgen(ny) - (ny-1)/2
sg_indy = 215 + indgen(ny) - (ny-1)/2
sg_indx = 1338 + indgen(nt) - (nt-1)/2

col_vec = (colortable(39))[254.*findgen(ny)/(ny-1), *]

dd = iris_obj(eout.sg_files[0])
sg_xp = dd->getxpos()
sg_yp = dd->getypos()
sg_time = anytim(tai2utc(dd->ti2tai(), /stime))
iris_resp = iris_get_response((dd->ti2utc())[0]) 
sg_spec_all = dd->getvar(8, /load)
sg_spec = sg_spec_all[*, sg_indy, sg_indx]

read_iris_l2, eout.sji_files[0], sji_index, sji_data
sji_time = anytim(sji_index.date_obs)
iris_lct, sji_index[0], rr2, gg2, bb2
sji_match = !null
for i=0, nt-1 do begin
  dum = min(abs(sji_time-sg_time[sg_indx[i]]), sji_match0)
  sji_match = [sji_match, sji_match0]
endfor
sji_match = sji_match[uniq(sji_match)]


w01 = window(dim=[16d2, 8d2])
p01 = !null
p02 = !null
for i=0, n_elements(sji_match)-1 do begin
  sji_img = iris_intscale(sji_data[*, *, sji_match[i]], sji_index[sji_match[i]])
  get_xp_yp, sji_index[sji_match[i]], sji_xp, sji_yp
  p01n = image_kh(sji_img, sji_xp, sji_yp, /current, $
                  pos=[0.01, 0.55, 0.17, 0.95] + [0.19, 0, 0.19, 0]*i + [0.05, 0, 0.05, 0], $
                  title=sji_index[sji_match[i]].date_obs, $
                  xtitle='Solar X (arcsec)', ytitle=(i eq 0) ? 'Solar Y (arcsec)' : '', $
                  rgb_table=[[rr2], [gg2], [bb2]], xtickdir=1, ytickdir=1, $
                  xr=sg_xp[sg_indx[2]]+5*[-1, 1], yr=sg_yp[sg_indy[2]]+5*[-1, 1])               
  if i ne 0 then p01n.yshowtext=0
  for j=0, ny-1 do begin
    p01m = symbol([sg_xp[sg_indx[0]]], [sg_yp[sg_indy[j]]], 'o', target=p01n, /data, $
                sym_color=reform(col_vec[j, *]), sym_size=0.7, sym_thick=2, sym_transp=0)
  endfor
  p01 = [p01, p01n]
  if sg_indx[i] eq 1338 then begin
    xind = interpol(findgen(n_elements(sji_xp)), sji_xp, sg_xp[sg_indx[i]])
    yind = interpol(findgen(n_elements(sji_yp)), sji_yp, sg_yp[sg_indy])
    light_curve_ind = [sji_match[i]-10:sji_match[i]+10]
    light_curve = interpolate(sji_data, xind, yind, light_curve_ind, /grid)
  endif
endfor


for k=0, nt-1 do begin
  p02n = plot(eout.wave_list[8], sg_spec[*, 0, 0], /current, /nodata, $
    pos=p01[0].pos + [0.19, 0, 0.19, 0]*k - [0, 0.5, 0, 0.5], $
    xtitle='Wavelength ($\AA$)', ytitle=(k eq 0) ? 'Intensity (DN)' : '', $
    xr=[2798.3, 2799.7], yr=[0, 100], xtickinterval=0.5, $
    title=anytim(sg_time[sg_indx[k]], /ccsds), $
    font_name='malgun gothic', font_size=13, font_style=1)
  if k ne 0 then p02n.yshowtext=0
  for j=0, ny-1 do begin
    p02m = plot(eout.wave_list[8], sg_spec[*, j, k], over=p02n, $
      color=reform(col_vec[j, *]), transp=0)
    p02m1 = plot(2798.75*[1, 1], p02n.yr, over=p02n, '--', color='gray', transp=50)
    p02m2 = plot(2798.82*[1, 1], p02n.yr, over=p02n, '--', color='gray', transp=50)
  endfor
  p02n.yr = [-10, p02n.yr[1]]
  p02 = [p02, p02n]
endfor

light_curve = reform(light_curve)
p03 = plot(sji_time[light_curve_ind], light_curve[0, *], /nodata, /current, $
           pos=[p02[3].pos[0]+0.05, p01[1].pos[1], p02[-1].pos[2], p01[1].pos[3]], $
           title = '2016-03-20 10:00 (UT)', $ 
           xtitle='Time (mm:ss)', ytitle='SJI 1400 $\AA$ Intensity (DN)', $ 
           font_name='malgun gothic', font_size=13, font_style=1)
           
for j=0, ny-1 do begin
  p03n = plot(sji_time[light_curve_ind], light_curve[j, *], over=p03, '-', $
              color=reform(col_vec[j, *]))
endfor

p03.xtickname = strmid(anytim(p03.xtickval, /ccsds), 14, 5)
p03.xtickinterval = 20

for k=0, nt-1 do begin
  p03n = plot(sg_time[sg_indx[k]]*[1, 1], p03.yr, over=p03, $
              '-2', color='gray', transp=50)
endfor
w01.save, 'Mg_II_Emiss.png', resol=200
end