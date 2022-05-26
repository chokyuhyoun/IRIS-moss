pro find_iris_spectra_moss, out_dir, eout, iris_dir=iris_dir, aia_dir=aia_dir

;dir = '/Users/khcho/Desktop/IRIS-moss-main/20160318_175434'
ti = systime(/sec)

; DEFAULT SETTING
cd, current=current_dir
this_file = file_which('find_iris_spectra_moss.pro')
if this_file eq '' then this_file = file_search('~/desktop', 'find_iris_moss_events3.pro', /fully)
cd, file_dirname(this_file)

if n_elements(aia_dir) eq 0 then aia_dir = out_dir
if n_elements(iris_dir) eq 0 then iris_dir = out_dir
sg_files = file_search(iris_dir, 'iris_l2_*raster*.fits')
sji_files = file_search(iris_dir, 'iris_l2_*SJI*.fits')
aia_files = file_search(aia_dir, 'aia_l2*.fits')
save_filename = out_dir+'/moss_event_'+strmid(file_basename(sji_files[0]), 8, 15)+'.sav'

dum = where(strmatch(sji_files, '*1400*'), n)
if n then begin
  read_iris_l2, sji_files[dum], sji_index, sji_data, /sil
  print, 'Use IRIS SJI_1400A for co-alignment'
endif else begin
  read_iris_l2, sji_files[0], sji_index, sji_data, /sil
  print, 'No IRIS SJI_1400A. Use IRIS '+sji_index.tdesc1+' for co-alignment'
endelse

; Check SJI data & fill up the missed time
dum = median(median(sji_data, dim=2), dim=1)
miss = where(dum lt 0, miss_n)
for i=0, miss_n-1 do begin
  sji_data[*, *, miss[i]] = 0.5*(sji_data[*, *, miss[i]-1] + sji_data[*, *, miss[i]+1])
  print, 'Filled '+string(miss[i], f='(i0)')+'th SJI data'
endfor

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
t_lim = 60.
t_lim_sji = floor(t_lim*2 / sji_index[0].cdelt3)

; TOLERANCE (moss from AIA 193A and SG, spatio-temporally)
;tol_dist = 1.5*aia193_index[0].cdelt1
tol_dist = 1 ; in arcsec
tol_t_sec = 0.5*aia193_index[0].cdelt3

if n_elements(sg_files) eq 0 then begin
  message, 'No IRIS data found. Please check the directory'
endif

;; RUN MOSSVAR ON AIA CUBE TO DETECT EVENTS
;if ~file_exist(save_filename) then begin
;  mossvar, aia_dir, '', mash, 4, 4,/read, /cube_data
;  if mash['status'] ne 'data ok' then begin
;    message, 'Failed to find moss structure from AIA data.'
;    stop
;  endif
;  mossvar, aia_dir, '', mash, 4, 4, /moss, /cnet, /fexviii, /loop, /filter1700
;  mossvar, aia_dir, '', mash, 4, 4, /variability, /cleanflare
;  zcube = mash['zmatch no loop']
;  
;;  stop
;endif else begin
if file_exist(save_filename) then begin
  restore, save_filename
  zcube = eout.zcube

  ; VARIABLES TO SAVE
  sg_ind = !null ;; [ypix, scan #, file #]
  sg_phy = !null
  sg_expt = !null
  sji_ind = !null ; t index for IRIS SJI cube
  sji_value = !null
  sji_peak_pos = !null
  sji_fwhm = !null
  sji_i_0 = !null
  sji_curve = !null
  sji_curve_ind = !null
  sji_peak_val = !null
  sji_phy = !null
  aia1600_ind = !null
  aia_ind = !null ; [x, y, t] index for IRIS / AIA cube
  aia_phy = !null  ; [x arcsec, y arcsec, sec from the beginning of AIA observation]
  aia_shift = !null  ; in arcsec, compare to sji
  aia193_value = !null
  
  for i=0, n_elements(sg_files)-1 do begin ; raster file no
    dd = iris_obj(sg_files[i])
    sg_xp = dd->getxpos()
    sg_yp = dd->getypos()
    sg_time = anytim(tai2utc(dd->ti2tai(0), /stime))
    iris_resp = iris_get_response((dd->ti2utc())[0])
    if i eq 0 then begin
      nwin = dd->getnwin()
      line_id = dd->getline_id()
      data_list = list(len=nwin)
      wave_list = list(len=nwin)
      for ii=0, nwin-1 do begin
        wave_list[ii] = dd->getlam(ii)
        sg_expt = [sg_expt, (dd->getexp(ii))[0]]
      endfor
    endif
  
    for j=0, n_elements(sg_xp)-1 do begin ; raster scan no
  ; time check
      sg_aia193_tdiff = min(abs(sg_time[j]-aia193_time), match_aia193)
      sg_aia1600_tdiff = min(abs(sg_time[j]-aia1600_time), match_aia1600)
      sg_sji_tdiff = min(abs(sg_time[j]-sji_time), match_sji)
      if sg_aia193_tdiff gt tol_t_sec then continue ; No matching time
      moss_img = zcube[*, *, match_aia193]
      moss_ind = where(moss_img eq 1, count)
      if count eq 0 then continue        ; No moss in AIA 193 A
  
  ; align correction
      aia_img = aia1600_data[*, *, match_aia1600]
      sji_img = sji_data[*, *, match_sji]
      aia_sji_align, aia_img, aia1600_index[match_aia1600], $
                     sji_img, sji_index[match_sji], del;, /cor0
  
      get_xp_yp, aia193_index[match_aia193], aia_xp, aia_yp
      aia_xp -= del[0]
      aia_yp -= del[1]
      aia193_xxp = rebin(aia_xp, n_elements(aia_xp), n_elements(aia_yp))
      aia193_yyp = rebin(transpose(aia_yp), n_elements(aia_xp), n_elements(aia_yp))
  
  ; check distances among moss positions and sg pixels
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
  
      if count eq 0 then continue   ;; No pixels within 1 arcsec
      real_ind1 = closest_ind[ind0]
      real_ind = array_indices(dist, real_ind1)
  
      spec_yind = reform(real_ind[1, *])
      moss_ind_spec = array_indices(moss_img, moss_ind[real_ind[0, *]])
  
      for k=0, n_elements(real_ind1)-1 do begin ; sg pixels within tolerance
        xpos = sg_xp[j]
        ypos = sg_yp[spec_yind[k]]
        get_xp_yp, sji_index[match_sji], sji_xp, sji_yp
        sji_xpix = interpol(findgen(n_elements(sji_xp)), sji_xp, xpos)
        sji_ypix = interpol(findgen(n_elements(sji_yp)), sji_yp, ypos)
  ; check light curve in SJI 1400 : variavility < 60s, I-I_0 > 3 sigma
        sji_curve_ind0 = [match_sji-t_lim_sji:match_sji+t_lim_sji]
        sji_curve0 = interpolate(sji_data, sji_xpix, sji_ypix, sji_curve_ind0, /grid)
        sji_curve0 = reform(sji_curve0) / sji_index[0].exptime * iris_resp.dn2phot_sji[1]
        sji_fwhm0 = !null
        var_chk, sji_curve0, sji_index, sji_i_00, sji_peak_val0, sji_peak_pos0, sji_fwhm0
        if n_elements(sji_fwhm0) eq 0 then continue
        sji_peak_pos0 += match_sji-t_lim_sji
  
  ; save variables
        for l=0, nwin-1 do begin
          dum = (dd->getvar(l))[*, spec_yind[k], j]
          data_list[l] = [[data_list[l]], [(dd->descale_array(dum))/ sg_expt[l]]]
        endfor
        sg_ind = [[sg_ind], [spec_yind[k], j, i]]  ; [ypix, scan #, file #]
        sg_phy = [[sg_phy], [xpos, ypos, sg_time[j]]]
        sji_ind = [[sji_ind], [sji_xpix, sji_ypix, match_sji]]
        sji_phy = [sji_phy, sji_time[match_sji]]
        sji_curve = [[sji_curve], [sji_curve0]]
        sji_curve_ind = [[sji_curve_ind], [sji_curve_ind0]]
        sji_peak_pos = [sji_peak_pos, sji_peak_pos0]
        sji_peak_val = [sji_peak_val, sji_peak_val0]
        sji_i_0 = [sji_i_0, sji_i_00]
        sji_fwhm = [sji_fwhm, sji_fwhm0]
        sji_value = [sji_value, sji_curve[n_elements(sji_curve0)/2]]
        aia1600_ind = [aia1600_ind, match_aia1600]
        aia_shift = [[aia_shift], [del]]
        aia_ind = [[aia_ind], [moss_ind_spec[*, k], match_aia193]]
        aia_phy = [[aia_phy], $
          [moss_xp[real_ind[0, k]], moss_yp[real_ind[0, k]], aia193_time[match_aia193]]]
  
  
        aia_xpix = interpol(findgen(n_elements(aia_xp)), aia_xp, xpos)
        aia_ypix = interpol(findgen(n_elements(aia_yp)), aia_yp, ypos)
        aia193_value = [aia193_value, interpolate(aia193_data[*, *, match_aia193], aia_xpix, aia_ypix)]
      endfor ;sg pixels
  
    endfor ; raster scan
    obj_destroy, dd
  endfor ; sg files

  if n_elements(data_list[0]) ne 0 then begin
    eout = {data_list:data_list, wave_list:wave_list, line_id:line_id, $
            sg_ind:sg_ind, sg_phy:sg_phy, $
            sji_ind:sji_ind, sji_phy:sji_phy, sji_value:sji_value, sji_peak_pos:sji_peak_pos, $
            sji_fwhm:sji_fwhm, sji_i_0:sji_i_0, sji_curve:sji_curve, $
            sji_curve_ind:sji_curve_ind, sji_dt:sji_index[0].cdelt3, $
            aia1600_ind:aia1600_ind, aia_ind:aia_ind, aia_phy:aia_phy, aia_shift:aia_shift, $
            aia193_value:aia193_value, zcube:zcube, sji_wave:sji_index[0].twave1, $
            sg_files:sg_files, sji_files:sji_files, aia_files:aia_files}
    save, eout, filename=save_filename
  endif else begin
    print, 'No results'
    save, sg_spec, filename=out_dir+'/No_event_'+strmid(file_basename(sji_files[0]), 8, 15)+'.sav'
  endelse
endif


print, 'It takes '+string((systime(/sec)-ti)/6d1, f='(f4.1)')+' min'
end
