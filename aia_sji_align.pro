pro get_xp_yp, index, xp, yp
  xp = (findgen(index.naxis1) - 0.5*(index.naxis1-1))*index.cdelt1 + index.crval1
  yp = (findgen(index.naxis2) - 0.5*(index.naxis2-1))*index.cdelt2 + index.crval2
end


pro aia_sji_align, aia_img, aia_index, sji_img, sji_index, del, $
                   aia_part=aia_part, sji_corr_part=sji_corr_part, cor0=cor0

aia_sz = size(aia_img)
sji_sz = size(sji_img)
;aia_wcs0 = fitshead2wcs(aia_index)
;aia_wcs = aia_wcs0
;sji_wcs = fitshead2wcs(sji_index)

get_xp_yp, aia_index, aia_xp0, aia_yp0
get_xp_yp, sji_index, sji_xp0, sji_yp0 

del = [0., 0.]
del0 = [0., 0.]
if keyword_set(cor0) ne 0 then run=2 else run=1
cor0 = !null
for ii=0, run do begin
  del = del + del0
;  aia_wcs.crval = aia_wcs0.crval - del
;  aia_coord = wcs_get_coord(aia_wcs)
;  aia_pixel = wcs_get_pixel(sji_wcs, aia_coord)
;  sji_corr_aia = reform(interpolate(sji_img, aia_pixel[0, *, *], aia_pixel[1, *, *], $
;                                    missing=-200))
;  xmin = min(where(sji_corr_aia[*, 0.5*aia_sz[2]] gt 0), max=xmax)
;  ymin = min(where(sji_corr_aia[0.5*aia_sz[1], *] gt 0), max=ymax)
;  aia_part = aia_img[xmin:xmax-1, ymin:ymax-1]
;  sji_corr_part = sji_corr_aia[xmin:xmax-1, ymin:ymax-1]
;  del0 = alignoffset(aia_part, sji_corr_part, cor)

  sji_xp1 = interpol(findgen(aia_sz[1]), aia_xp0-del[0], sji_xp0)
  sji_yp1 = interpol(findgen(aia_sz[2]), aia_yp0-del[1], sji_yp0)
  sji_corr_aia = interpolate(aia_img, sji_xp1, sji_yp1, /grid)
  del0 = alignoffset(sji_corr_aia, sji_img, cor)
  
  cor0 = [cor0, cor]
;  print, del, cor
;  stop
endfor
del = del*[aia_index.cdelt1, aia_index.cdelt2]
;stop
end


dir = '/Users/khcho/Desktop/moss_test2'
sji_files = file_search(dir, 'iris_l2_*SJI*.fits')
aia_files = file_search(dir, 'aia_l2*.fits')
dum = where(strmatch(sji_files, '*1400*'), n)
if n then begin
  read_iris_l2, sji_files[dum], sji_index, sji_data, /sil
  print, 'Use IRIS SJI_1400A for co-alignment'
endif else begin
  read_iris_l2, sji_files[0], sji_index, sji_data, /sil
  print, 'No IRIS SJI_1400A. Use IRIS '+sji_index.tdesc1+' for co-alignment'
endelse
dum = where(strmatch(aia_files, '*1600*'))
read_iris_l2, aia_files[dum], aia1600_index, aia1600_data, /sil
aia_time = anytim(aia1600_index.date_obs)
sji_time = anytim(sji_index.date_obs)
dels = !null
for i=0, n_elements(sji_index)-1 do begin
  dum = min(abs(sji_time[i]-aia_time), match)
  aia_img = aia1600_data[*, *, match]
  sji_img = sji_data[*, *, i]
  aia_sji_align, aia_img, aia1600_index[match], sji_img, sji_index[match], del
  dels = [[dels], [del]]
endfor
end


