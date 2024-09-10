pro find_iris_spectra_moss, out_dir, eout, iris_dir=iris_dir, aia_dir=aia_dir

;out_dir = '/Users/khcho/Desktop/IRIS-moss-main/20160320_110924'
print, systime()
ti = systime(/sec)
obj_pix_num = 0
moss_num = 0

; DEFAULT SETTING
cd, current=current_dir
this_file = file_which('find_iris_spectra_moss.pro')
if this_file eq '' then this_file = file_search('~/desktop', 'find_iris_moss_events3.pro', /fully)
cd, file_dirname(this_file)

if n_elements(aia_dir) eq 0 then aia_dir = out_dir
if n_elements(iris_dir) eq 0 then iris_dir = out_dir
sg_files = file_search(iris_dir, 'iris_l2_*raster*.fits')
sji_files = file_search(iris_dir, 'iris_l2_*SJI*.fits')
aia_files = file_search(aia_dir, 'aia_l2*.fits', count=n)
if n eq 0 then begin
  print, 'No avaliable AIA files'
  save, filename=out_dir+'No_AIA_files.sav'
  return
endif

save_filename = out_dir+'/moss_event_'+strmid(file_basename(sji_files[0]), 8, 15)+'.sav'

dum = where(strmatch(sji_files, '*1400*') eq 1, n)
if n gt 0 then begin
  read_iris_l2, sji_files[dum], sji_index, sji_data, /sil
  print, 'Use IRIS SJI_1400A for co-alignment'
  sji_time = anytim(sji_index.date_obs)

  ; Check SJI data & fill up the missed time
  dum = median(median(sji_data, dim=2), dim=1)
  miss = where(dum lt 0, miss_n)
  for i=0, miss_n-1 do begin
    sji_data[*, *, miss[i]] = 0.5*(sji_data[*, *, miss[i]-1] + sji_data[*, *, miss[i]+1])
    print, 'Filled '+string(miss[i], f='(i0)')+'th SJI data'
  endfor
  
  dum = where(strmatch(aia_files, '*1600.fits'))
  read_iris_l2, aia_files[dum], aia1600_index, aia1600_data, /sil
  
  aia1600_time = anytim(aia1600_index.date_obs)

  ; co-alignment between AIA 1600 A and SJI 1400 A - 
  del1 = !null
  cor1 = !null
  dum = min(abs(aia1600_time - sji_time[0]), mn)
  dum = min(abs(aia1600_time - sji_time[-1]), mx)
  for i=mn, mx do begin
    dum = min(abs(aia1600_time[i]-sji_time), match)
    aia_img = aia1600_data[*, *, i]
    sji_img = sji_data[*, *, match]
    aia_sji_align, aia_img, aia1600_index[i], $
                   sji_img, sji_index[match], del0, cor0=cor0
    del1 = [[del1], [del0]]  
    cor1 = [cor1, cor0]
  ;  stop
  endfor
  tt = aia1600_time[mn:mx]
  unvalid = where(cor1 lt 0.4)  ;  empirical value
  if n_elements(unvalid)*1./n_elements(cor1) gt 0.3 then begin
    del1[*] = 0.
  endif else begin
    del1[*, unvalid] = !values.f_nan
  endelse
  
  aia_shift = [[reform(del1[0, *])], [reform(del1[1, *])], [tt]] ; [delx, dely, time]
endif else begin
  print, 'No IRIS SJI_1400A. Use IRIS header information for co-alignment'
  read_iris_l2, sji_files[0], sji_index, sji_data, /sil
  sji_time = anytim(sji_index.date_obs)
  aia_shift = [[fltarr(n_elements(sji_time))], [fltarr(n_elements(sji_time))], [sji_time]]
endelse

dum = where(strmatch(aia_files, '*193.fits'))
read_iris_l2, aia_files[dum], aia193_index, aia193_data, /sil

; MAKE TIME ARRAY
aia193_time = anytim(aia193_index.date_obs)

aia_sz = size(aia193_data)
sji_sz = size(sji_data)
t_lim = 150.
t_lim_sji = floor(t_lim / sji_index[0].cdelt3)

; TOLERANCE (moss from AIA 193A and SG, spatio-temporally)
;tol_dist = 1.5*aia193_index[0].cdelt1
tol_dist = 1 ; in arcsec
tol_t_sec = 0.5*aia193_index[0].cdelt3

if n_elements(sg_files) eq 0 then begin
  print, 'No IRIS data found. Please check the directory'
  stop
endif

;; RUN MOSSVAR ON AIA CUBE TO DETECT EVENTS
  mossvar, aia_dir, '', mash, 4, 4,/read, /cube_data
  if mash['status'] ne 'data ok' then begin
    print, 'Failed to find moss structure from AIA data.'
    stop
  endif
  mossvar, aia_dir, '', mash, 4, 4, /moss, /cnet, /fexviii, /loop, /filter1700, $
          /demfilter, /rundem, dem_save_path=out_dir
  mossvar, aia_dir, '', mash, 4, 4, /variability, /cleanflare

  zcube = mash['zmatch no loop']
  err_t = where(total(total(zcube, 1), 1) gt 100, /null)
  zcube[*, *, err_t] = 0
  
  ; VARIABLES TO SAVE
  sg_ind = !null ;; [ypix, scan #, file #]
  sg_phy = !null
  sg_expt = !null
  sji_ind = !null ; [xpix, ypix, t] index for IRIS SJI cube
  sji_phy = !null ; sji time
  curve_fwhm = !null
  curve_peak_ind = !null
  si_iv_curve = !null
  si_iv_time = !null
  pre_sg_ind = !null
  aft_sg_ind = !null
  
  aia_ind = !null ; [x, y, t] index for IRIS / AIA cube
  aia_phy = !null  ; [x arcsec, y arcsec, sec from the beginning of AIA observation]
  aia193_value = !null
  pre_aia_ind = !null
  
  si_iv_fit_res = !null
  mg_ii_fit_res = !null
  pre_si_iv_fit_res = !null
  pre_mg_ii_fit_res = !null
  aft_si_iv_fit_res = !null
  aft_mg_ii_fit_res = !null

  for i=0, n_elements(sg_files)-1 do begin ; raster file no
    dd = iris_obj(sg_files[i])
    sg_xp = dd->getxpos()
    sg_yp = dd->getypos()
    sg_time = anytim(tai2utc(dd->ti2tai(), /stime))
    iris_resp = iris_get_response((dd->ti2utc())[0])
    sg_dy = dd->getresy()
    if i eq 0 then begin
      nwin = dd->getnwin()
      line_id = dd->getline_id()
      data_list = list(len=nwin)
      pre_data_list = list(len=nwin)
      aft_data_list = list(len=nwin)
      wave_list = list(len=nwin)
      for ii=0, nwin-1 do begin
        wave_list[ii] = dd->getlam(ii)
        sg_expt = [[sg_expt], [(dd->getexp(iwin=ii))]]
      endfor
    endif
    Si_IV_ind = (where(strmatch(line_id, '*Si IV*'), /null))[-1]
    mg_ii_ind = (where(strmatch(line_id, '*2796*'), /null))[0]
    sit_and_stare = dd->getsit_and_stare()
    sg_dt = mean(sg_time[1:*] - sg_time[0:-2])
    t_lim_sg = floor(t_lim*2 / sg_dt)
    if (Si_IV_ind ne -1) and sit_and_stare then begin
      Si_IV_int_ind = where(wave_list[Si_IV_ind] gt 1402.77-2. and wave_list[Si_IV_ind] lt 1402.77+2.)
    endif else continue ;; only for SI IV obs. and sit-and-stare
    
    spat_bin =  (dd->binning_region('FUV'))[0]
    spec_bin = (dd->binning_spectral(si_iv_ind))[0]   
;    stop
    for j=0, n_elements(sg_xp)-1 do begin ; raster scan no
  ; time check
      sg_aia193_tdiff = min(abs(sg_time[j]-aia193_time), match_aia193)
;      sg_aia1600_tdiff = min(abs(sg_time[j]-aia1600_time), match_aia1600)
      sg_sji_tdiff = min(abs(sg_time[j]-sji_time), match_sji)
      if sg_aia193_tdiff gt tol_t_sec then continue ; No matching time
      moss_img = zcube[*, *, match_aia193]
      moss_ind = where(moss_img eq 1, count)
      if count eq 0 then continue        ; No moss in AIA data
  
  ; align correction  
      get_xp_yp, aia193_index[match_aia193], aia_xp, aia_yp
      aia_xp -= interpol(aia_shift[*, 0], aia_shift[*, 2], aia193_time[match_aia193], /nan)
      aia_yp -= interpol(aia_shift[*, 1], aia_shift[*, 2], aia193_time[match_aia193], /nan)
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
      get_xp_yp, sji_index[match_sji], sji_xp, sji_yp
      
      for k=0, n_elements(real_ind1)-1 do begin ; sg pixels within tolerance
        obj_pix_num += 1
        xpos = sg_xp[j]
        ypos = sg_yp[spec_yind[k]]
        sji_xpix = interpol(findgen(n_elements(sji_xp)), sji_xp, xpos)
        sji_ypix = interpol(findgen(n_elements(sji_yp)), sji_yp, ypos)
    ; check light curve in SJI 1400 : variavility < 60s, I-I_0 > 3 sigma
        if (SI_IV_ind ne -1) and sit_and_stare then begin
          si_iv_curve0 = !null
          ind00 = 0 > (j-t_lim_sg)
          ind01 = (n_elements(sg_time)-1) < (j+t_lim_sg)
          si_iv_time0 = sg_time[ind00:ind01]
          curve_dt = sg_dt
          for l=ind00, ind01 do begin
            dum = reform((dd->getvar(SI_IV_ind))[*, spec_yind[k], l])
            dum = (dd->descale_array(dum))/ sg_expt[l, si_iv_ind]
            si_iv_curve0 = [si_iv_curve0,  total(dum[Si_IV_int_ind])] 
          endfor
          t_lim_ind = t_lim_sg
        endif else begin             
          ind00 = 0 > match_sji-t_lim_sji
          ind01 = (n_elements(sji_time)-1) < (match_sji+t_lim_sji)
          sji_curve_ind0 = [ind00:ind01]
          sji_curve0 = total(total(sji_data[sji_xpix-1:sji_xpix+1, sji_ypix-1:sji_ypix+1, sji_curve_ind0], 1), 1)
          si_iv_curve0 = reform(sji_curve0) / sji_index[0].exptime * iris_resp.dn2phot_sji[1]
          si_iv_time0 = sji_time[sji_curve_ind0]
          curve_dt = sji_index[0].cdelt3
          t_lim_ind = t_lim_sji
        endelse
        if n_elements(si_iv_curve0) lt 2*t_lim_ind + 1 then begin ;; center of si_iv_curve0 = j
          deficit = 2*t_lim_ind+1-n_elements(si_iv_curve0)
          if match_sji lt 0.5*n_elements(sji_index) then begin
            si_iv_curve0 = [fltarr(deficit)*!values.f_nan, si_iv_curve0]
            si_iv_time0 = [fltarr(deficit)*!values.f_nan, si_iv_time0]
          endif else begin
            si_iv_curve0 = [si_iv_curve0, fltarr(deficit)*!values.f_nan]
            si_iv_time0 = [si_iv_time0, fltarr(deficit)*!values.f_nan]
          endelse
        endif
        curve_fwhm0 = !null
;        if total(finite(si_iv_curve0)) lt 0.5*n_elements(si_iv_curve0) then stop
        var_chk, si_iv_curve0, curve_dt, curve_fwhm0, i_0, curve_peak_ind0, before_ind, after_ind
;        stop
        if n_elements(curve_fwhm0) eq 0 then continue
        moss_num += 1
        j_new = j - floor(300./curve_dt)+1 + curve_peak_ind0
        j0 = (j_new+before_ind) > 0  ;; before_ind : negative
        j1 = (j_new+after_ind) < (n_elements(sg_xp)-1)
        
        if si_iv_ind ge 0 then begin
          dum = (dd->getvar(si_iv_ind))[*, spec_yind[k], j_new]
          si_spec = (dd->descale_array(dum)) *1. / sg_expt[j_new, si_iv_ind] / spec_bin / spat_bin  
          si_iv_fit_res0 = si_iv_fit(wave_list[Si_IV_ind], si_spec, spec_bin=spec_bin)
          
          dum = (dd->getvar(si_iv_ind))[*, spec_yind[k], j0]
          pre_si_spec = (dd->descale_array(dum)) *1. / sg_expt[j0, si_iv_ind] / spec_bin / spat_bin   
          pre_si_iv_fit_res0 = si_iv_fit(wave_list[Si_IV_ind], pre_si_spec, spec_bin=spec_bin)
          
          dum = (dd->getvar(si_iv_ind))[*, spec_yind[k], j1]
          aft_si_spec = (dd->descale_array(dum)) *1. / sg_expt[j1, si_iv_ind] / spec_bin / spat_bin 
          aft_si_iv_fit_res0 = si_iv_fit(wave_list[Si_IV_ind], aft_si_spec, spec_bin=spec_bin)
        endif else si_iv_fit_res0 = 0
        
        if mg_ii_ind ge 0 then begin
          dum = (dd->getvar(mg_ii_ind))[*, spec_yind[k], j_new]
          mg_spec = (dd->descale_array(dum)) / sg_expt[j_new, mg_ii_ind]
          mg_ii_fit_res0 = mg_ii_fit(wave_list[mg_ii_ind], mg_spec)

          dum = (dd->getvar(mg_ii_ind))[*, spec_yind[k], j0]
          pre_mg_spec = (dd->descale_array(dum)) / sg_expt[j0, mg_ii_ind]
          pre_mg_ii_fit_res0 = mg_ii_fit(wave_list[mg_ii_ind], pre_mg_spec)

          dum = (dd->getvar(mg_ii_ind))[*, spec_yind[k], j1]
          aft_mg_spec = (dd->descale_array(dum)) / sg_expt[j1, mg_ii_ind]
          aft_mg_ii_fit_res0 = mg_ii_fit(wave_list[mg_ii_ind], aft_mg_spec)
        endif else mg_ii_fit_res0 = [0, 0]
;        stop
  ; save variables
        for l=0, nwin-1 do begin  ; for spectral windows
          dum = (dd->getvar(l))[*, spec_yind[k], j_new]
          data_list[l] = [[data_list[l]], [(dd->descale_array(dum))/ sg_expt[j_new, l]]]
          dum = (dd->getvar(l))[*, spec_yind[k], j0]
          pre_data_list[l] = [[pre_data_list[l]], [(dd->descale_array(dum))/ sg_expt[j0, l]]]
          dum = (dd->getvar(l))[*, spec_yind[k], j1]
          aft_data_list[l] = [[aft_data_list[l]], [(dd->descale_array(dum))/ sg_expt[j1, l]]]

        endfor
        sg_ind = [[sg_ind], [spec_yind[k], j_new, i]]  ; [ypix, scan #, file #]
        pre_sg_ind = [[pre_sg_ind], [spec_yind[k], j0, i]]  
        aft_sg_ind = [[aft_sg_ind], [spec_yind[k], j1, i]]
        sg_phy = [[sg_phy], [xpos, ypos, sg_time[j_new]]]
        sji_ind = [[sji_ind], [sji_xpix, sji_ypix, match_sji]]
        sji_phy = [sji_phy, sji_time[match_sji]]
        si_iv_curve = [[si_iv_curve], [si_iv_curve0]]
        si_iv_time = [[si_iv_time], [si_iv_time0]]
        curve_fwhm = [curve_fwhm, curve_fwhm0]
        curve_peak_ind = [curve_peak_ind, curve_peak_ind0]
        aia_ind = [[aia_ind], [moss_ind_spec[*, k], match_aia193]]
        pre_aia_ind = [[pre_aia_ind], [moss_ind_spec[*, k], (match_aia193-curve_fwhm0/12.)>0]]
        aia_phy = [[aia_phy], $
          [moss_xp[real_ind[0, k]], moss_yp[real_ind[0, k]], aia193_time[match_aia193]]]
  
        aia_xpix = interpol(findgen(n_elements(aia_xp)), aia_xp, xpos)
        aia_ypix = interpol(findgen(n_elements(aia_yp)), aia_yp, ypos)
        aia193_value = [aia193_value, interpolate(aia193_data[*, *, match_aia193], aia_xpix, aia_ypix)]
        si_iv_fit_res = [si_iv_fit_res, si_iv_fit_res0]
        mg_ii_fit_res = [[mg_ii_fit_res], [mg_ii_fit_res0]]
        pre_si_iv_fit_res = [pre_si_iv_fit_res, pre_si_iv_fit_res0]
        pre_mg_ii_fit_res = [[pre_mg_ii_fit_res], [pre_mg_ii_fit_res0]]
        aft_si_iv_fit_res = [aft_si_iv_fit_res, aft_si_iv_fit_res0]
        aft_mg_ii_fit_res = [[aft_mg_ii_fit_res], [aft_mg_ii_fit_res0]]
;      stop
      endfor ;sg pixels
;    if strmatch(aia193_index[match_aia193].date_obs, '*09:42:29*') eq 1 then stop
    endfor ; raster scan
    obj_destroy, dd
  endfor ; sg files

;  stop
  if n_elements(data_list[0]) ne 0 then begin
    dum = file_search(out_dir+'/No_event_'+strmid(file_basename(sji_files[0]), 8, 15)+'.sav', $
                      count=prev_exist)
    if prev_exist ge 1 then file_delete, dum
    eout = {data_list:data_list, wave_list:wave_list, line_id:line_id, $
            sg_ind:sg_ind, sg_phy:sg_phy, sg_dy:sg_dy, pre_sg_ind:pre_sg_ind, aft_sg_ind:aft_sg_ind, $
            sji_ind:sji_ind, sji_phy:sji_phy, $
            curve_fwhm:curve_fwhm, si_iv_curve:si_iv_curve, si_iv_time:si_iv_time, $
            curve_dt:curve_dt, curve_peak_ind:curve_peak_ind, $
            aia_ind:aia_ind, aia_phy:aia_phy, aia_shift:aia_shift, $
            aia193_value:aia193_value, zcube:zcube, sji_wave:sji_index[0].twave1, $
            sg_files:sg_files, sji_files:sji_files, aia_files:aia_files, $
            si_iv_fit_res:si_iv_fit_res, mg_ii_fit_res:mg_ii_fit_res, $
            obj_pix_num:obj_pix_num, moss_num:moss_num, $
            pre_si_iv_fit_res:pre_si_iv_fit_res, pre_mg_ii_fit_res:pre_mg_ii_fit_res, $
            pre_data_list:pre_data_list, pre_aia_ind:pre_aia_ind, $
            aft_si_iv_fit_res:aft_si_iv_fit_res, aft_mg_ii_fit_res:aft_mg_ii_fit_res, $
            aft_data_list:aft_data_list, spat_bin:spat_bin, spec_bin:spec_bin}

    save, eout, filename=save_filename
  endif else begin
    print, 'No results'
    dum = file_search(save_filename, count=prev_exist)
    if prev_exist ge 1 then file_delete, dum
    eout = {wave_list:wave_list, line_id:line_id, $
            aia_shift:aia_shift, $
            zcube:zcube, sji_wave:sji_index[0].twave1, $
            sg_files:sg_files, sji_files:sji_files, aia_files:aia_files, $
            obj_pix_num:obj_pix_num, moss_num:moss_num}
    save, eout, filename=out_dir+'/No_event_'+strmid(file_basename(sji_files[0]), 8, 15)+'.sav'
  endelse
  mash = 0
  data_list = 0
  pre_data_list = 0
  aft_data_list = 0
  wave_list = 0
;endif
print, out_dir + ' finished. It takes '+string((systime(/sec)-ti)/6d1, f='(f4.1)')+' min'

end