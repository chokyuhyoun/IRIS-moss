;pro find_iris_spectra_moss, dir, eout, outdir=outdir

dir = '/Users/khcho/Desktop/moss_test2'
ti = systime(/sec)

; DEFAULT SETTING
cd, current=current_dir
this_file = file_which('find_iris_moss_events3.pro')
if this_file eq '' then this_file = file_search('~/desktop', 'find_iris_moss_events3.pro', /fully)
cd, file_dirname(this_file)

if n_elements(outdir) eq 0 then outdir = dir
sg_files = file_search(dir, 'iris_l2_*raster*.fits')
sji_files = file_search(dir, 'iris_l2_*SJI*.fits')
aia_files = file_search(dir, 'aia_l2*.fits')
aia_dir = file_dirname(aia_files[0])
save_dir = file_dirname(sji_files[0])

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
dum = where(strmatch(aia_files, '*193*'))
read_iris_l2, aia_files[dum], aia193_index, aia193_data, /sil

; MAKE TIME ARRAY
sji_time = anytim(sji_index.date_obs)
aia1600_time = anytim(aia1600_index.date_obs)
aia193_time = anytim(aia193_index.date_obs)
aia_sz = size(aia193_data)
sji_sz = size(sji_data)

;; pointing correction
;aia_shift = !null  ; in arcsec, relative to sji 
;for i=0, n_elements(sji_index)-1 do begin
;  dum = min(abs(sji_time[i]-aia1600_time), match)
;  aia_img = aia1600_data[*, *, match]
;  sji_img = sji_data[*, *, i]
;  aia_sji_align, aia_img, aia1600_index[match], sji_img, sji_index[match], del, /cor0
;  aia_shift = [[aia_shift], [del]]
;endfor
;stop

if n_elements(sg_files) eq 0 then begin
  message, 'No IRIS data found. Please check the directory'
endif

; RUN MOSSVAR ON AIA CUBE TO DETECT EVENTS
mossvar, dir, '',mash,4,4,/read,/cube_data
if mash['status'] ne 'data ok' then begin
  message, 'Failed to find moss structure from AIA data.'
  stop
endif
mossvar,dir,'',mash,4,4,/moss,/cnet,/fexviii,/loop,/filter1700
if keyword_set(movie) then mossvar,dir,'',mash,4,4,/variability,/cleanflare,/movie
if ~keyword_set(movie) then mossvar,dir,'',mash,4,4,/variability,/cleanflare
zcube = mash['zmatch no loop']

; VARIABLE TO SAVE
sji_ind = !null ; t index for IRIS SJI cube
aia1600_ind = !null 
aia_ind = !null ; [x, y, t] index for IRIS / AIA cube
aia_phy = !null  ; [x arcsec, y arcsec, sec from the beginning of AIA observation]
aia_shift = !null  ; in arcsec, compare to sji
sg_ind = !null ; 
sg_spec = !null ; [lmabda, moss index]
sg_phy = !null

; TOLERANCE (moss from AIA 193A and SG, spatio-temporally) 
tol_dist = 1.5*aia193_index[0].cdelt1
tol_t_sec = 0.5*aia193_index[0].cdelt3


for i=0, n_elements(sg_files)-1 do begin ; raster file no
  read_iris_l2, sg_files[i], sg_index, sg_data, /sil
  sg_time = anytim(sg_index.date_obs)
  dd = iris_obj(sg_files[i])
  sg_xp = dd -> getxpos()
  sg_yp = dd -> getypos()
  
  for j=0, n_elements(sg_xp)-1 do begin ; raster scan no
; time check
    sg_aia193_tdiff = min(abs(sg_time[j]-aia193_time), match_aia193)
    sg_aia1600_tdiff = min(abs(sg_time[j]-aia1600_time), match_aia1600)
    sg_sji_tdiff = min(abs(sg_time[j]-sji_time), match_sji)
    if sg_aia193_tdiff gt tol_t_sec then continue

; align correction
    aia_img = aia1600_data[*, *, match_aia1600]
    sji_img = sji_data[*, *, match_sji]
    aia_sji_align, aia_img, aia1600_index[match_aia1600], $
                   sji_img, sji_index[match_sji], del;, /cor0
    
;    aia_wcs = fitshead2wcs(aia193_index[match_aia193])
;    aia_wcs.crval = aia_wcs.crval - del
;    aia_coord = wcs_get_coord(aia_wcs)
;    aia193_xxp = reform(aia_coord[0, *, *])
;    aia193_yyp = reform(aia_coord[1, *, *])
    get_xp_yp, aia193_index[match_aia193], aia_xp, aia_yp
    aia_xp -= del[0]
    aia_yp -= del[1]
    aia193_xxp = rebin(aia_xp, n_elements(aia_xp), n_elements(aia_yp))
    aia193_yyp = rebin(transpose(aia_yp), n_elements(aia_xp), n_elements(aia_yp))

;FIND IT!    
    moss_img = zcube[*, *, match_aia193]
    moss_ind = where(moss_img eq 1)        
    moss_xp = aia193_xxp[moss_ind]
    moss_yp = aia193_yyp[moss_ind]
    
    nx = n_elements(moss_xp)
    ny = n_elements(sg_yp)
    moss_xxp = rebin(moss_xp, nx, ny)
    moss_yyp = rebin(moss_yp, nx, ny)
    sg_xxp = rebin(transpose(replicate(sg_xp[j], ny)), nx, ny)
    sg_yyp = rebin(transpose(sg_yp), nx, ny)
    
    dist = sqrt((sg_xxp - moss_xxp)^2. + (sg_yyp - moss_yyp)^2.)
    min_dist = min(dist, dim=1, closest_ind)
    ind0 = where(min_dist le tol_dist, count)
;    stop
    if count ne 0 then begin
      real_ind1 = closest_ind[ind0]
      real_ind = array_indices(dist, real_ind1)

      spec_yind = reform(real_ind[1, *])
      moss_ind_spec = array_indices(moss_img, moss_ind[real_ind[0, *]])

      for k=0, n_elements(real_ind1)-1 do begin
        sg_spec = [[sg_spec], [sg_data[*, spec_yind[k], j]]]
        sg_ind = [[sg_ind], [spec_yind[k], j, i]]
        sg_phy = [[sg_phy], [sg_xp[j], sg_yp[spec_yind[k]], sg_time[j]]]
        sji_ind = [sji_ind, match_sji]
        aia1600_ind = [aia1600_ind, match_aia1600]
        aia_shift = [[aia_shift], [del]]
        aia_ind = [[aia_ind], [moss_ind_spec[*, k], match_aia193]]
        aia_phy = [[aia_phy], $
          [moss_xp[real_ind[0, k]], moss_yp[real_ind[0, k]], aia193_time[match_aia193]]]
      endfor ; every sg pixels
;      stop          
    endif                  
  endfor ; raster scan
endfor ; sg files

eout = {sg_spec:sg_spec, sg_ind:sg_ind, sg_phy:sg_phy, sji_ind:sji_ind, $
        aia_ind:aia_ind, aia_phy:aia_phy, aia1600_ind:aia1600_ind, zcube:zcube, $ 
        sg_files:sg_files, sji_files:sji_files, aia_files:aia_files, $
        aia_shift:aia_shift}

save_filename = save_dir+'/moss_event_'+strmid(file_basename(sji_files[0]), 8, 15)+'.sav'
save, eout, filename=save_filename

;; verify the moss position --> save images
hfov = 25
w1 = window(dim=[12d2, 5d2])
aia_lct, rr0, gg0, bb0, wavelnth=193
aia_lct, rr1, gg1, bb1, wavelnth=1600
iris_lct, sji_index[0], rr2, gg2, bb2
im01 = image_kh(aia193_data[*, *, 0], /current, position=[0.08, 0.15, 0.4, 0.85], $
                xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)', title='AIA 193', $
                rgb_table=[[rr0], [gg0], [bb0]], min=0, max=5e3, $
                xr=[-hfov, hfov], yr=[-hfov, hfov], aspect_ratio=1)
im02 = image_kh(aia193_data[*, *, 0], /current, $
                position=[im01.pos[2], 0.15, im01.pos[2]+0.3, 0.85], $
                xtitle='Solar X (arcsec)', title='AIA 1600', $
                rgb_table=[[rr1], [gg1], [bb1]], min=0, max=5e3, $
                xr=[-hfov, hfov], yr=[-hfov, hfov], aspect_ratio=1)
im02.axes[1].showtext = 0
im03 = image_kh(sji_data[*, *, 0], /current, $
                position=[im02.pos[2], 0.15, im02.pos[2]+0.3, 0.85], $
                xtitle='Solar X (arcsec)', title='SJI', $
                rgb_table=[[rr2], [gg2], [bb2]], min=0, max=1e3,$ 
                xr=[-hfov, hfov], yr=[-hfov, hfov], aspect_ratio=1)
im03.axes[1].showtext = 0

con02 = contour(aia193_data[*, *, 0], color='k', c_value=[4d3], c_label_show=0, over=im02)
con03 = contour(aia1600_data[*, *, 0], color='k', c_value=[5e2], c_label_show=0, over=im03)

im011 = image(zcube[*, *, 0], over=im01, rgb_table=30, min=0, max=1, transp=30)
im021 = image(zcube[*, *, 0], over=im02, rgb_table=30, min=0, max=1, transp=30)
im031 = image(zcube[*, *, 0], over=im03, rgb_table=30, min=0, max=1, transp=30)
p02 = plot([sg_phy[0, 0]], [sg_phy[1, 0]], '+', sym_size=3, sym_thick=2, over=im03)                
t0 = text(0.02, 0.05, 'Moss Position', color=(im011.rgb_table)[*, 1], transp=im011.transp, $
          font_size=15, font_name='malgun_gothic')
t02 = text(im02.pos[0]+0.02, im02.pos[1]+0.05, 'AIA 193 $\AA$', font_size=13)
t03 = text(im03.pos[0]+0.02, im03.pos[1]+0.05, 'AIA 1600 $\AA$', font_size=13)
for i=0, (size(sg_spec))[2]-1 do begin
;  i = 0
  sjiz = sji_ind[i]
  get_xp_yp, sji_index[sjiz], sji_xp, sji_yp
  im03.setdata, sji_data[*, *, sjiz], sji_xp, sji_yp
  im03.title = 'SJI '+string(sji_index[sjiz].twave1, f='(i0)')+' $\AA$ ' $
                + strmid(sji_index[sjiz].date_obs, 0, 19)
  im03.max = 1e3
  
  aiaz = aia_ind[2, i]
  get_xp_yp, aia193_index[aiaz], aia_xp, aia_yp
  aia_xp -= aia_shift[0, i]
  aia_yp -= aia_shift[1, i]
  im01.setdata, aia193_data[*, *, aiaz], aia_xp, aia_yp
  im01.title = 'AIA 193 $\AA$ '+strmid(aia193_index[aiaz].date_obs, 0, 19)
  im01.max = 5e3
  
  im02.setdata, aia1600_data[*, *, aia1600_ind[i]], aia_xp, aia_yp
  im02.title = 'AIA 1600 $\AA$ '+strmid(aia1600_index[aia1600_ind[i]].date_obs, 0, 19)
  im02.max = 5e2
  
  con02.setdata, aia193_data[*, *, aiaz], aia_xp+0.6, aia_yp+0.6
  con03.setdata, aia1600_data[*, *, aia1600_ind[i]], aia_xp+0.6, aia_yp+0.6
   
  dum = float(zcube[*, *, aiaz])
  dum[where(dum eq 0)] = !values.f_nan
  im011.setdata, dum, aia_xp, aia_yp
  im021.setdata, dum, aia_xp, aia_yp
  im031.setdata, dum, aia_xp, aia_yp
  
  xr = sg_phy[0, i] + [-1, 1]*hfov
  yr = sg_phy[1, i] + [-1, 1]*hfov
  im01.xr = xr
  im01.yr = yr
  im02.xr = xr
  im02.yr = yr
  im03.xr = xr
  im03.yr = yr
  
  p02.setdata, [eout.sg_phy[0, i]], [eout.sg_phy[1, i]]  
;stop
  file_mkdir, save_dir+'/imgs'
  w1.save, save_dir+'/imgs/'+string(i, f='(i03)')+'.png', resol=200
endfor

print, 'It takes '+string((systime(/sec)-ti)/6d1, f='(f4.1)')+' s'
end

