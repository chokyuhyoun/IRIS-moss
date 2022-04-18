;; test aia - iris co-alignment 
;


dir = '/Users/khcho/Desktop/moss_test2'
sav_file = file_search(dir, '*.sav')
restore, sav_file
aia_map = eout.event_map
sji_map = eout.sji_map

;aia_tai = anytim(eout.aia_time, /tai)

aia_file = eout.aia_files(where(strmatch(eout.aia_files, '*1600.fits')))
sji_file = eout.sji_files(where(strmatch(eout.sji_files, '*1400.fits')))
read_iris_l2, aia_file, aia_index, aia_data
read_iris_l2, sji_file, sji_index, sji_data
aia_sz = size(aia_data)
sji_sz = size(sji_data)

sg_time = anytim(eout.sg_time)
aia_time = anytim(aia_index.date_obs)
sji_time = anytim(sji_index.date_obs)

aia_stack = !null
sji_stack = !null

;for i=0, n_elements(eout.sg_time)-1 do begin
  i = 0
  dum = min(abs(sg_time[i] - aia_time), match_aia)
  dum = min(abs(sg_time[i] - sji_time), match_sji)
  
  aia_img = aia_data[*, *, match_aia]
  aia_xpdata = (findgen(aia_sz[1])-(aia_sz[1]-1)*0.5)*aia_index[match_aia].cdelt1 $
                + aia_index[match_aia].crval1 - eout.xshift_iris2aia[match_aia]
  aia_ypdata = (findgen(aia_sz[2])-(aia_sz[2]-1)*0.5)*aia_index[match_aia].cdelt2 $
                + aia_index[match_aia].crval2 - eout.yshift_iris2aia[match_aia]
  sji_xpdata = (findgen(sji_sz[1])-(sji_sz[1]-1)*0.5)*sji_index[match_sji].cdelt1 $
                + sji_index[match_sji].crval1
  sji_ypdata = (findgen(sji_sz[2])-(sji_sz[2]-1)*0.5)*sji_index[match_sji].cdelt2 $
                + sji_index[match_sji].crval2
                  
  sji_int_xpix = interpol(findgen(sji_sz[1]), sji_xpdata, aia_xpdata)
  sji_int_ypix = interpol(findgen(sji_sz[2]), sji_ypdata, aia_ypdata)
  
  sji_corr_aia = interpolate(sji_data[*, *, match_sji], $
                             sji_int_xpix, sji_int_ypix, /grid)
  sji_stack = [sji_stack, sji_corr_aia[*]]
  aia_stack = [aia_stack, aia_img[*]]
;endfor

del = alignoffset(aia_img, sji_corr_aia)

aia_title = 'AIA 1600 $\AA$ '+aia_index[match_aia].date_obs
sji_title = 'SJI 1400 $\AA$ '+sji_index[match_sji].date_obs
w1 = window(dimension=[8e2, 8e2])
im02 = image(aia_img, aia_xpdata, aia_ypdata, axis=2, /current, $
             xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)', $
             title=aia_title, max=1e3, position=[0.15, 0.15, 0.9, 0.9], $
             font_size=15)
im03 = image(sji_data[*, *, match_sji], sji_xpdata, sji_ypdata, $
             min=0, max=3e2, /over)
t02 = text(0.2, 0.85, $
           'Misalignment : ('+ string(del[0]*aia_index[0].cdelt1, f='(f5.2)')+'", ' $
           + string(del[1]*aia_index[0].cdelt2, f='(f5.2)')+'")', $
           font_size=15, color='white', current=w1) 
             
video = idlffvideowrite('align4.mp4')
wdims = w1.dimensions
stream = video.addvideostream(wdims[0], wdims[1], 5)
for i=0, 50 do begin
  if (i mod 10) lt 5 then begin
    im03.hide = 1
    im02.title = aia_title
  endif else begin
    im03.hide = 0
    im03.title = sji_title
  endelse
  timestemp = video.put(stream, im02.copywindow())
endfor
video.cleanup


real = where(sji_stack gt 0 and aia_stack gt 0)
aia_stack = aia_stack[real]
sji_stack = sji_stack[real]

h2d = hist_2d(aia_stack, sji_stack)
cc = correlate(aia_stack, sji_stack)
coeff = poly_fit(aia_stack, sji_stack, 1)
im33 = image(h2d, axis=2, position=[0.15, 0.15, 0.9, 0.9], max=2e1, $
             xrange=[0, 6e2], yrange=[0, 5e2], rgb_table=39, aspect_ratio=0, $
             xtitle='AIA 1600 $\AA$ DN', ytitle='SJI 1400 $\AA$ DN', $
             title='AIA 1600 $\AA$ Vs SJI 1400 $\AA$', font_size=15)
p33 = plot([0, 1d4], poly([0, 1d4], coeff), '--w2', /over)             
text1 = text(0.18, 0.82, 'cc = '+string(cc, f='(f5.2)'), color='white', font_size=15, $
             /over, /nor)
text2 = text(0.18, 0.77, 'y = '+string(coeff[1], f='(f5.2)')+' x + '+string(coeff[0], f='(f5.2)'), $
             color='white', font_size=15, /over, /nor)
im33.save, 'align4.png', resol=200
end