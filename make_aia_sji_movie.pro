pro make_aia_sji_movie, sav_files

dir = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, dir
;sav_files = file_search(dir, '*2020*/*event*.sav', /fully)

win_dim = [16d2, 10d2]
side = 330
gap = 30
imx0 = 100
imy0 = 640
sig = 5
nplot = 5
py0 = imy0 - gap - side
py1 = imy0 + side 
px0 = imx0 + 2.*side + gap + 50
pxside1 = 110.
pxside2 = pxside1*5.
px1 = win_dim[0]-30
pyside = (py1-py0)/nplot

wave_list = ['94', '193', '1700', '171', '211', '1600', 'sji']
wv_inds = [0, 1, 2, 6]
wv_drange = [[0, 30], [1e2, 5e3], [1e2, 3e3], [0, 255]]
si_cen = 1402.77d
mg_h_cen = 2803.5310d0
mg_k_cen = 2796.3501d0
mg_triplet = 2798.823d0
si_dr = si_cen + 1.*[-1, 1]
mg_h_dr = mg_h_cen + 2.*[-1, 1]
mg_k_dr = mg_k_cen + 2.*[-1, 1]
w_th_si = si_cen/3d8*sqrt(8.*alog(2.)*1.38d-23*10d0^(4.9)/(28.0855*1.6605d-27))  ; in angstrom
w_inst = 0.026  ; in angstrom

for i=0, n_elements(sav_files)-1 do begin
;for i=12, 12 do begin
  print, sav_files[i]
  movie_name = strmid(sav_files[i], 0, strlen(sav_files[i])-3)+'mp4'
  t0 = systime(/sec)
  sav_dir = file_dirname(sav_files[i])+'/imgs'
  file_mkdir, sav_dir

  restore, sav_files[i], /relax
  Si_win_ind = (where(eout.line_id eq 'Si IV 1403', si_win))[0]
  Mg_win_ind = (where(eout.line_id eq 'Mg II k 2796', mg_win))[0]
  data = list(len=n_elements(wave_list))
  indices = list(len=n_elements(wave_list))
  times = list(len=n_elements(wave_list))
  for j = 0, n_elements(wave_list)-1 do begin
    if wave_list[j] eq 'sji' then begin
      file = eout.sji_files[0]
    endif else begin
      file = eout.aia_files[where(strmatch(eout.aia_files, '*_'+wave_list[j]+'.fits'))]
    endelse
    read_iris_l2, file, index0, data0, /sil
    sz = size(data0)
    data[j] = data0 / rebin(reform(index0.exptime, 1, 1, sz[3]), sz[1], sz[2], sz[3])
    indices[j] = index0
    times[j] = anytim(index0.date_obs)
  endfor
  fe_xviii_data = data[0] - data[3]/450. - data[4]/120.
  fe_xviii_curve = total(total(fe_xviii_data, 1), 1)/(size(data[1]))[1]/(size(data[1]))[2]
  moss_num_curve = total(total(eout.zcube, 1), 1)
  fovx = (indices[-1])[0].fovx
  fovy = (indices[-1])[0].fovy
  hfov = 0.5*(fovx > fovy)
  titles = ['Fe XVIII ', 'AIA 193$\AA$ ', 'AIA 1700$\AA$ ', $
            'SJI '+string(eout.sji_wave, f='(i4)')+'$\AA$ ']

;  stop
  w1 = window(dim=win_dim, buffer=1)
  whole_win = plot(indgen(2), /current, /nodata, axis_style=0, pos=[0, 0, 1, 1], $
                   xr=[0, 1], yr=[0, 1])
  im01 = !null
  im011 = !null
  t01 = !null
  for jj=0, 1 do begin
    for ii=0, 1 do begin
      kk = ii + jj*2
      if kk eq 3 then begin
        iris_lct, (indices[-1])[0], rr, gg, bb 
      endif else begin
        aia_lct, rr, gg, bb, wavelnth=wave_list[wv_inds[kk]], /load
      endelse
      dum0 = image_kh(hanning(2, 2), /current, /dev, $
                      pos=[imx0+ii*(gap+side),      imy0-jj*(gap+side), $
                           imx0+ii*(gap+side)+side, imy0-jj*(gap+side)+side], $
                      xtitle = (jj eq 1) ? 'Solar X (arcsec)' : '', $
                      ytitle = (ii eq 0) ? 'Solar Y (arcsec)' : '', $
                      xtickdir=1, ytickdir=1, xticklen=0.03, yticklen=0.03, $
                      rgb_table=[[rr], [gg], [bb]], background_color='black')
      if ii eq 1 then dum0.yshowtext = 0
      if jj eq 0 then dum0.xshowtext = 0
      if (ii eq 1) and (jj eq 1) then p01 = plot([0, 1], [0, 1], over=dum0, thick=2) 
      dum1 = image_kh(hanning(2, 2), over=dum0, rgb_table=34, transp=0)
      dum2 = text(dum0.pos[0]+0.01, dum0.pos[3]-0.02, '', /current, $
                     font_size=12, font_name='malgun gothic', font_style=1, $
                     vertical_align=1, font_color='white')
      im01 = [im01, dum0]
      im011 = [im011, dum1]
      t01 = [t01, dum2]                               
    endfor
  endfor
    
  p03 = objarr(nplot)
  p031 = objarr(nplot)
  t031 = objarr(nplot)
  
  p04 = objarr(nplot)
  p041 = objarr(nplot, 2)
  t041 = objarr(nplot, 2)
  t040 = objarr(nplot)
  
  for l=0, nplot-1 do begin
    p03[l] = plot(indgen(2), /hide, /current, /dev, $
                   pos=[px0, py1-pyside*(l+1), px0+pxside1, py1-pyside*l], $
                   title=(l eq 0) ? 'Si IV 1403' : '', $
                   font_size=12, font_name='malgun gothic', font_style=1, $
                   xr=si_cen+1.*[-1, 1], xticklen=0.08, yticklen=0, $
                   xtitle='Wavelength ($\AA$)', xtickinterval=1., xminor=1, yshowtext=0, $
                   yminor=1, xstyle=1)
    if l ne nplot-1 then p03[l].xshowtext = 0
    p031[l] = plot(indgen(2), /hide, '-2r', transp=40, over=p03[l])
    t031[l] = text(p03[l].pos[0]+0.005, p03[l].pos[3]-0.007, '  ', $
                   font_color=p031[l].color, transp=p031[l].transp, $
                   font_size=9, vertical_align=1)

    p04[l] = plot(indgen(2), /hide, /current, /dev, $
                   pos=[px1-pxside2, py1-pyside*(l+1), px1, py1-pyside*l], $
                   title=(l eq 0) ? 'Mg II k 2796' : '', $
                   font_size=12, font_name='malgun gothic', font_style=1, $
                   xr=[2795d, 2805d]+[-0.2, 0.2], xtitle='Wavelength ($\AA$)', $
                   yr=[0, 1], xticklen=0.08, yticklen=0, yshowtext=0, yminor=1, $
                   xstyle=1)
    if l ne nplot-1 then p04[l].xshowtext = 0           
    for m=0, 1 do begin
      p041[l, m] = plot(indgen(2), /hide, '-2r', transp=40, over=p04[l])
      par_pos = ([1133., 1345.]/win_dim[0])[m]
      t041[l, m] = text(par_pos, p04[l].pos[3]-0.07, '  ', $
                        font_color=p041[l, m].color, transp=p041[l, m].transp, $
                        font_size=10, vertical_align=1)
    endfor
    t040[l] = text(mean(p04[l].pos[[0, 2]]), p04[l].pos[3]-0.05, '  ', $
                  font_size=12, font_name='malgun gothic', font_style=1, $
                  align=0.5)    
  endfor
  
  iris_resp = iris_get_response((times[0])[0])
  no_obs_t = !null
  if (~si_win) then begin
    dum0 = text(mean((p03[0].pos)[[0, 2]]), mean((p03[0].pos)[[1, 3]]), $
                'No obs.', /normal, align=0.5, vertical_align=0.5, $
                font_size=12, font_name='malgun gothic', font_style=1)
    no_obs_t = [no_obs_t, dum0]            
  endif else begin
    si_wave = eout.wave_list[si_win_ind]
    si_coeff = !null
    si_dn2phot = iris_resp.dn2phot_sg[0]
  endelse
  si_cen_pos = (p03[0].convertcoord(si_cen, 0, /data, /to_normal))[0]
  p_si_cen = plot(si_cen_pos*[1, 1], [p03[-1].pos[1], p03[0].pos[3]], over=whole_win, $
                  '--k', transp=50)

  if (~mg_win) then begin
    dum0 = text(mean((p03[0].pos)[[0, 2]]), mean((p03[0].pos)[[1, 3]]), $
                'No obs.', /normal, align=0.5, vertical_align=0.5, yr=[0, 1], $
                font_size=12, font_name='malgun gothic', font_style=1)
    no_obs_t = [no_obs_t, dum0]            
  endif else begin
    mg_wave = eout.wave_list[mg_win_ind]
    mg_coeff = !null
    mg_dn2phot = iris_resp.dn2phot_sg[1] 
  endelse
  mg_trip_pos = (p04[0].convertcoord(mg_triplet, 0, /data, /to_normal))[0]
  p_mg_trip = plot(mg_trip_pos*[1, 1], [p04[-1].pos[1], p04[0].pos[3]], over=whole_win, $
                   '--k', transp=50)
  mg_h_pos = (p04[0].convertcoord(mg_h_cen, 0, /data, /to_normal))[0]
  p_mg_h = plot(mg_h_pos*[1, 1], [p04[-1].pos[1], p04[0].pos[3]], over=whole_win, $
                   '--k', transp=50)
  mg_k_pos = (p04[0].convertcoord(mg_k_cen, 0, /data, /to_normal))[0]
  p_mg_k = plot(mg_k_pos*[1, 1], [p04[-1].pos[1], p04[0].pos[3]], over=whole_win, $
                   '--k', transp=50)
  
  i_time = (times[-1])[0]
  f_time = (times[-1])[-1]
;  if (indices[-1])[0].cdelt3 gt (indices[0])[0].cdelt3 then begin
    dum = min(abs(times[1]-i_time), ind1)
    dum = min(abs(times[1]-f_time), ind2)
    t_arr = (times[1])[ind1:ind2]
;  endif else t_arr = times[-1]

  dum = floor((times[0]-0.1)/600d0)*600d0
  xtickv = (dum[uniq(dum)])[1:*]
  p05 = plot(times[0], fe_xviii_curve, /current, /dev, $
             pos=[imx0, 60, win_dim[0]-80, 200], $
             xstyle=1, xtitle=strmid(anytim((times[0])[0], /ccsds), 0, 10)+' (UT)', $
             xtickv=xtickv, xtickname=strmid(anytim(xtickv, /ccsds), 11, 5), $
             ytitle='Fe XVIII Light Curve', $
             font_size=12, font_name='malgun gothic', font_style=1)
  p05.axes[3].hide = 1
  p05yr = p05.yr
  p051 = fillplot((times[0])[0:ind1], [[fe_xviii_curve[0:ind1]], [fltarr(ind1+1)+p05.yr[0]]], $
                  over=p05, fill_color='gray', color='gray')
  p052 = plot(times[0], moss_num_curve, /current, pos=p05.pos, $
              color='red', /histogram, ytitle='# of Moss', yminor=0, $
              font_size=12, font_name='malgun gothic', font_style=1, xstyle=1)
  for ii=0, 2 do p052.axes[ii].hide=1
  p052.axes[3].showtext = 1
  p052.axes[3].color='red'
  p05.xstyle=1
  p05.yr = p05yr

  if strmatch(sav_files[i], '*moss_event*') then begin
    sg_ind = intarr(n_elements(eout.curve_fwhm))
    for ii = 0, n_elements(eout.curve_fwhm)-1 do begin
      dum = min(abs(t_arr - eout.sg_phy[2, ii]), dum1)
      sg_ind[ii] = dum1 
      p053 = plot(eout.sg_phy[2, ii]*[1, 1], [(p05yr[1]-p05yr[0])*[0.6, 0.8]+p05yr[0]], over=p05, $
                  '-2r')
    endfor
  endif else sg_ind = !null
;  stop
;================================================================
  percent = !null
  match0 = intarr(4) - 1
  match1 = intarr(4)
  for j=0, n_elements(t_arr)-1 do begin
;  for j=163, 108 do begin
    dum = floor(10.*j/(n_elements(t_arr)-1))*10
    if dum ne percent then begin
      print, string(dum, f='(i3)')+' %'
      percent = dum
    endif

    moss = where(sg_ind eq j, count)
;    if count eq 0 then continue

    dum = min(abs(times[-1] - t_arr[j]), match_sji)
    dum = min(abs(times[1] - t_arr[j]), match_aia)
    dum = min(abs(times[2] - t_arr[j]), match_1600)
    if match_aia gt (size(fe_xviii_data))[3]-1 then continue
    match1 = [match_aia, match_aia, match_1600, match_sji]
    cenx = (indices[-1])[match_sji].crval1
    ceny = (indices[-1])[match_sji].crval2
    xr = cenx + hfov*[-1, 1]
    yr = ceny + hfov*[-1, 1]
    for k=0, 3 do begin   ; AIA & SJI
      wv_ind = wv_inds[k]
      if match0[k] eq match1[k] then continue 
      match0[k] = match1[k]
      match = match1[k]
      cur_data = (k eq 0) ? fe_xviii_data[*, *, match] : (data[wv_ind])[*, *, match] 
      get_xp_yp, (indices[wv_ind])[match], xp, yp
      if k ne 3 then begin
        xp -= interpol(eout.aia_shift[*, 0], eout.aia_shift[*, 2], (times[wv_ind])[match], /nan)
        yp -= interpol(eout.aia_shift[*, 1], eout.aia_shift[*, 2], (times[wv_ind])[match], /nan)        
      endif else cur_data = iris_intscale(temporary(cur_data), (indices[wv_ind])[match])
      im01[k].xr = xr
      im01[k].yr = yr
      setdata_hi_res, im01[k], cur_data, xp, yp
      im01[k].min = wv_drange[0, k]
      im01[k].max = wv_drange[1, k]      
      if k eq 1 then begin
        moss_img = float(eout.zcube[*, *, match])
        moss_img[where(moss_img eq 0)] = !values.f_nan
        aia193_xp = xp
        aia193_yp = yp
      endif
      if k eq 3 then begin
        p01.setdata, xp[(indices[-1])[match_sji].sltpx1ix*[1, 1]], yr
        for l=0, 3 do begin
          setdata_hi_res, im011[l], moss_img, aia193_xp, aia193_yp
          im011[l].min = 0
          im011[l].max = 1
        endfor
      endif
      t01[k].string = titles[k] + strmid((indices[wv_ind])[match].date_obs, 0, 19)
    endfor

;-------------------------------    

    page = 0

plot_spectra :
    rep = (count-1)/nplot
    png_name = (count gt 0) ? sav_dir+path_sep()+string(j, f='(i05)')+'_moss'+string(page, f='(i01)')+'.png' $
      : sav_dir+path_sep()+string(j, f='(i05)')+'.png'

    for l = 0, count-1 do begin  ; Spectrograph data
      if l ge 5 then continue
      if si_win then begin
        temp_data = (eout.data_list[si_win_ind])[*, moss[l]]
        cenp = where((si_wave ge si_dr[0]) and (si_wave le si_dr[1]))
        si_wavep = si_wave[cenp]
        temp_datap = temp_data[cenp]
        p03[l].setdata, si_wavep, temp_datap
        p03[l].ystyle =2 
        p03[l].yr = [-3, p03[l].yr[1]*1.2]
        p03[l].axes[1].showtext=1
        p03[l].yticklen=0.07
  ;      p03[l].ytickval = p03[l].ytickval[0:-2]
        si_err = sqrt((temp_datap*si_dn2phot) > 0)/si_dn2phot + 1.
        si_err[si_err lt 1.] = 1e5
        fit_res = gaussfit(si_wavep, temp_datap, coeff, nterms=4, measure_error=si_err, estimates=[30., si_cen, 1d-2, 0.])
        p031[l].setdata, si_wavep, fit_res
        w_fwhm = 2.*sqrt(2.*alog(2.))*coeff[2]
        w_nth = sqrt(w_fwhm^2. - w_th_si^2. - w_inst^2.)
        t031[l].string = 'A = '+string(coeff[0], f='(f0.1)')+' DN!c'+ $
                         '$v_D$ = '+string(3d5*(coeff[1]-si_cen)/si_cen, f='(f0.1)')+' km s$^{-1}$!c'+ $
                         '$v_{nth}$ = '+string(w_nth*3d5/si_cen, f='(f0.1)')+' km s$^{-1}$'
        p03[l].hide = 0
        p031[l].hide = 0
        t031[l].hide = 0
        si_coeff = [[si_coeff], [coeff]]
      endif ; si_win

      if mg_win then begin 
        temp_data = (eout.data_list[mg_win_ind])[*, moss[l]]
        dum = where((mg_wave ge p04[l].xr[0]) and (mg_wave le p04[l].xr[1]))
        p04[l].setdata, mg_wave[dum], temp_data[dum]
        p04[l].ystyle = 2
        p04[l].yr = [-10., p04[l].yr[1]]
        p04[l].axes[1].showtext = 1
        p04[l].yticklen = 0.03
  ;      p04[l].ytickval = p04[l].ytickval[0:-2]
        t040[l].string = strmid(anytim(eout.sg_phy[2, moss[l]], /ccsds), 0, 19)+'!c'+$
                        '(' +string(eout.sg_phy[0, moss[l]], f='(f7.2)')+$
                        ', '+string(eout.sg_phy[1, moss[l]], f='(f7.2)')+')'
        t040[l].hide = 0
        p04[l].hide = 0        
        for m=0, 1 do begin
          mg_dr = ([[mg_k_dr], [mg_h_dr]])[*, m]
          cenp = where((mg_wave ge mg_dr[0]) and (mg_wave le mg_dr[1]))
          mg_wavep = mg_wave[cenp]
          mg_wavepn = dindgen(1d4)*(mg_wavep[-1]-mg_wavep[0])*1d-4 + mg_wavep[0]
          temp_datap = temp_data[cenp]
          mg_err = sqrt(temp_datap*mg_dn2phot)/mg_dn2phot + 1.
          mg_err[mg_err lt 1.] = 1e5
          lims = {value:0., fixed:0, limited:[0, 0], limits:[0., 0.]}
          lims = replicate(lims, 9)
          lims[0].limited[0] = 1 & lims[0].limits[0] = 0d    ; background
          lims[1].limited[*] = 1 & lims[1].limits = [1e-5, 1d1]   ; linear slope
          lims[2].limited[*] = 1 & lims[2].limits = mean(mg_dr) + 0.5*[-1, 1]   ; centeroid of linear slope
          lims[3].limited[0] = 1 & lims[3].limits[0] = 150d0   ; positive Gaussian amplitude
          lims[4].limited[*] = 1 & lims[4].limits = mean(mg_dr) + 0.6*[-1, 1]   ; centeroid of positive Gaussian
          lims[6].limited[*] = 1 & lims[6].limits = [0, 20] ; fractional of negative Gaussian amplitude
          lims[7].limited[*] = 1 & lims[7].limits = mean(mg_dr) + 0.4*[-1, 1]  ; centeroid of negative Gaussian
          lims[8].limited[1] = 1 & lims[8].limits[1] = 0.1
          p_init = [0, 1d0, mean(mg_dr), 1d3, mean(mg_dr), 0.04, 0.95, mean(mg_dr), 0.01]
          p = mpfitfun('iris_mgfit', mg_wavep, temp_datap, mg_err, p_init, $
                       parinfo=lims, perror=g, /quiet, $
                       maxiter=400, status=st, ftol=1d-9)             
          p041[l, m].setdata, mg_wavepn, iris_mgfit(mg_wavepn, p)
          mg_rvpos = iris_postfit_gen(mg_wavepn, iris_mgfit(mg_wavepn, p))
          par_pos = ([1133., 1343.]/win_dim[0])[m]
          mg_cen = ([mg_k_cen, mg_h_cen])[m]
          kh = (m eq 0) ? 'k' : 'h'
;          print, p
          t041[l, m].string = '$v_{'+kh+'3}$ = '+string((mg_rvpos[8]-mg_cen)*3d5/mg_cen, f='(f0.2)')+' km s$^{-1}$'
          p041[l, m].hide = 0
          t041[l, m].hide = 0 
;          stop
        endfor  
      endif               
    endfor
    p051.setdata, (times[0])[0:ind1+j], [[fe_xviii_curve[0:ind1+j]], [fltarr(ind1+j+1)+p05.yr[0]]]

    w1.save, png_name, resol=130
    
    if count ne 0 then begin
      for l=0, nplot-1 do begin
        p03[l].yticklen = 0
        p03[l].yshowtext = 0
        p03[l].hide = 1
        p031[l].hide = 1
        t031[l].hide = 1
        p04[l].yshowtext = 0
        p04[l].yticklen = 0
        p04[l].hide = 1
        t040[l].hide = 1
        for m=0, 1 do begin
          p041[l, m].hide = 1
          t041[l, m].hide = 1
        endfor
      endfor
    endif

    if rep ne 0 then begin
      moss = moss[nplot:*]
      count = n_elements(moss)
      dum = where(sg_ind eq j, dum1) 
      page += 1

      goto, plot_spectra
    endif 
  endfor ; j
  w1.close
;  video.cleanup
  png_files = file_search(sav_dir, '*.png')
  ffmpeg, png_files, output = strmid(sav_files[i], 0, strlen(sav_files[i])-3)+'mp4'
  print, 'It took '+string((systime(/sec)-t0)/60., f='(f0.2)')+' min'
endfor ; i
data = 0
indices = 0
end