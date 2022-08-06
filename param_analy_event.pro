dir = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, dir
ar_num = 'ar12529'
dum = file_search(dir, '*event*.sav')
sav_files = dum[where(strmatch(dum, '*'+ar_num+'*'))]
out_dir = '/'+strjoin((strsplit(sav_files[0], '/', /extract))[0:-3], '/')+'/'

t0 = !null
t1 = !null
tr = !null
tt = !null
ind0 = !null
ind1 = !null
si_v = !null
si_nth = !null
si_peak = !null
mg_k_v = !null
mg_h_v = !null
mg_trip = !null
moss_num = !null
err_si = !null
err_mg_k = !null
err_mg_h = !null
aia_phy = !null
phase = !null
sg_ind = !null
sg_dy = !null
dem_res = !null
if 1 then begin
  for i=0, n_elements(sav_files)-1 do begin
    restore, sav_files[i], /relax
    t0 = [t0, eout.aia_shift[0, 2]]
    t1 = [t1, eout.aia_shift[-1, 2]]
    tr = [tr, 0.5*(t0[-1]+t1[-1])]  
    if eout.moss_num gt 0 then begin
      ind0 = [ind0, n_elements(tt)]
      tt = [tt, reform(eout.sg_phy[2, *])]
      si_peak = [si_peak, reform(eout.si_iv_fit_res.coeff[0, *])]
      si_v = [si_v, eout.si_iv_fit_res.v_d]
      si_nth = [si_nth, eout.si_iv_fit_res.v_nth]
      mg_k_v = [mg_k_v, reform(eout.mg_ii_fit_res[0, *].v_d_3)]
      mg_h_v = [mg_h_v, reform(eout.mg_ii_fit_res[1, *].v_d_3)]
      ind1 = [ind1, n_elements(tt)-1]
      err_si = [err_si, eout.si_iv_fit_res.chisq]
      err_mg_h = [err_mg_h, reform(eout.mg_ii_fit_res[0, *].chisq)]
      err_mg_k = [err_mg_k, reform(eout.mg_ii_fit_Res[1, *].chisq)]
      mg_trip = [mg_trip, reform(eout.mg_ii_fit_res[2, *].emiss)]
      aia_phy = [[aia_phy], [eout.aia_phy]]
      cen = (size(eout.si_iv_curve))[1]/2
      phase = [phase, sgn(cen - 0.5 - eout.curve_peak_ind)]  ; -1 = early, 1 = late
      dum = eout.sg_ind
      dum[2, *] = i
      sg_ind = [[sg_ind], [dum]]   ; [y_pix, time in a dataset, dataset]
      sg_dy = [[sg_dy], replicate(eout.sg_dy, eout.moss_num)]
      
      dir_name = file_dirname(sav_files[i])
      restore, dir_name+'/emcube.sav'
      dem_t = anytim(demstr.time)
      dem_t1 = rebin(dem_t, n_elements(dem_t), eout.moss_num)
      moss_t1 = rebin(eout.aia_phy[2, *], n_elements(dem_t), eout.moss_num)
      t_diff = abs(dem_t1 - moss_t1)
      t_diff_min = min(t_diff, dim=1, t_min_ind0)
      t_prev_ind = reform((array_indices(t_diff, t_min_ind0))[0, *]) - 2 > 0
      for j=0, eout.moss_num-1 do begin
        dem_res = [[dem_res], [reform(demstr[t_prev_ind[j]].dem[eout.aia_ind[0, j], eout.aia_ind[1, j], *])]]
      endfor
    endif
    moss_num = [moss_num, eout.moss_num]
  endfor
  event_no = intarr(total(moss_num))
  no = 0
  for i=1, n_elements(event_no)-1 do begin
    same_event = where((abs(sg_ind[0, 0:i-1] - sg_ind[0, i]) le 2./sg_dy[i]) and $
      (abs(sg_ind[1, 0:i-1] - sg_ind[1, i]) le 1.) and $
      (sg_ind[2, 0:i-1] eq sg_ind[2, i]), count)
    if count eq 0 then no += 1
    event_no[i] = no
  endfor

  temp_arr = 5.5+findgen(21)*0.1

  titles = ['Si IV Peak Intensity (Photons pix$^{-1}$ s$^{-1}$)', $
            'Si IV Nonthermal Velocity (km s$^{-1}$)', $
            'Si IV Doppler Velocity (km s$^{-1}$)', $
            'Mg II k3 Doppler Velocity (km s$^{-1}$)', $
            'Mg II h3 Doppler Velocity (km s$^{-1}$)']
  pars = [[si_peak], [si_nth], [si_v], [mg_k_v], [mg_h_v], [total(dem_res, 1)]]
  par_names = ['si_peak', 'si_nth', 'si_v', 'mg_k_v', 'mg_h_v', 'e_den']
  
  sz = [n_elements(event_no), max(event_no)+1, n_elements(par_names)]
  cast = make_array(sz[0], sz[1], value = !values.f_nan)
  cast[indgen(sz[0]), event_no] = 1
  cast_all = rebin(cast, sz[0], sz[1], sz[2])
  pars_all = rebin(reform(pars, sz[0], 1, sz[2]), sz[0], sz[1], sz[2]) * cast_all
  pars_ev_mean = mean(pars_all, dim=1, /nan)
  pars_ev_std = stddev(pars_all, dim=1, /nan)
  pars_ev_std[where(~finite(pars_ev_std))] = 0.
  t_ev = (reform(aia_phy[2, *]))[uniq(event_no)]
  save, t0, t1, tr, tt, ind0, ind1, si_peak, si_v, si_nth, mg_k_v, mg_h_v, mg_trip, $
        err_mg_h, err_mg_k, err_si, event_no, cast, pars_ev_mean, pars_ev_std, t_ev, $
        dem_res, $
        filename=out_dir+'gethering_moss_param.sav'
  eout = 0 
endif else restore, out_dir+'gethering_moss_param.sav'

;----------------------------------





stop
;

;--------------------------------------------------------
w1 = window(dim=[8d2, 8d2])
p011 = objarr(max(event_no)+1)
p012 = objarr(max(event_no)+1, 3)
p01 = plot3d(pars_ev_mean[*, 0], pars_ev_mean[*, 1], pars_ev_mean[*, 2], $
            /current, /nodata, clip=0, axis=2, $
            ytitle='Si V$_{nth}$', ztitle='Si V$_D$', xtitle='Si peak Int.', $
            font_name='malgun gothic', font_size=13, font_style=1)
for i=0, 2 do begin
  p01.axes[i].showtext = 1
  p01.axes[i].text_orient = 90
endfor
gether_sav_file = file_search(dir, '*gether*.sav')
for j=0, n_elements(gether_sav_file)-1 do begin
  restore, gether_sav_file[j]
  col = (['red', 'blue'])[j]
  for i=0, max(event_no) do begin ; event 
;    col = reform((colortable(40))[(bytscl(t_ev))[i], *])
    p011[i] = plot3d([pars_ev_mean[i, 0]], [pars_ev_mean[i, 1]], [pars_ev_mean[i, 2]], $
                     over=p01, sym_object=orb(), sym_color=col, sym_transparency=60, $
                     sym_size = (total(cast, 1, /nan)/30 + 0.5)[i])
    p012[i, 0] = plot3d(pars_ev_mean[i, 0]+pars_ev_std[i, 0]*[1, -1], $
                        pars_ev_mean[i, 1]*[1, 1], $
                        pars_ev_mean[i, 2]*[1, 1], $
                        over=p01, color='gray', transp=80)
    p012[i, 1] = plot3d(pars_ev_mean[i, 0]*[1, 1], $
                        pars_ev_mean[i, 1]+pars_ev_std[i, 1]*[1, -1], $
                        pars_ev_mean[i, 2]*[1, 1], $
                        over=p01, color='gray', transp=80)
    p012[i, 2] = plot3d(pars_ev_mean[i, 0]*[1, 1], $
                        pars_ev_mean[i, 1]*[1, 1], $
                        pars_ev_mean[i, 2]+pars_ev_std[i, 2]*[1, -1], $
                        over=p01, color='gray', transp=80)      
  ;stop                   
  endfor
endfor          
stop
;-----------------------------------------------------------------------------------
if 0 then begin 
  for i=0, 4 do begin
    win_dim = [10d2, 6d2]
    w1 = window(dim=win_dim, buffer=1)
    
    par = pars[*, i]    
    xr0 = floor(t0[0]/86400d0)*86400d0
    xr1 = ceil(t1[-1]/86400d0)*86400d0
    xtickv = dindgen((xr1-xr0)/86400d0/2+1)*86400d0*2 + xr0
  
    p01 = plot(tt, par, '.', /current, pos=[0.1, 0.1, 0.7, 0.7], $
               xtitle='Time (2016)', ytitle=titles[i], $
               xstyle=1, ystyle=1, xr=[xr0, xr1], xtickv=xtickv, xminor=3, $
               xtickname=strmid(anytim(xtickv, /ccsds), 5, 5), $
               font_name='malgun gothic', font_size=13, font_style=1)
  
    m = !null
    sig = !null
    tmid = !null
    for j=0, n_elements(ind0)-1 do begin
      m = [m, mean(par[ind0[j]:ind1[j]], /nan)]
      sig = [sig, stddev(par[ind0[j]:ind1[j]], /nan)]
      tmid = [tmid, mean(tt[ind0[j]:ind1[j]], /nan)] 
    endfor
    p03 = errorplot(tmid, m, sig, over=p01, $
                    symbol='o', sym_filled=1, sym_color='red', sym_fill_color='red', $
                    errorbar_color='blue', errorbar_thick=2, yr=p01.yr, ystyle=1)
    tt00 = [0, dblarr(n_elements(t0)*2)]
    par00 = [0, dblarr(n_elements(t0)*2)]
    for j=0, n_elements(t0)-1 do begin
      p02 = fillplot([t0[j], 0.5*(t0+t1)[j], t1[j]], $
                    [[p01.yr[1]*[1, 1, 1]], [p01.yr[0]*[1, 1, 1]]], $
                    over=p01, fill_color='light_gray', color='light_gray')
      p02.order, /send_to_back
      tt00[2*j] = t0[j]
      tt00[2*j+1] = t1[j]
      par00[2*j] = moss_num[j]
    endfor
    p04 = plot(tt00, par00, /current, $
              pos=[p01.pos[0], p01.pos[3]+0.01, p01.pos[2], 0.9], $
              xr=p01.xr, xshowtext=0, xtickv=xtickv, ytitle='# of Moss', $
              hist=1, fill_color='blue', fill_background=1, color='k', $
              font_name='malgun gothic', font_size=13, font_style=1)
    
    hist = histogram(par, loc=xbin, nbin=20)
    p05 = barplot(xbin, hist, /hor, /current, $
                  pos=[p01.pos[2]+0.01, p01.pos[1], 0.9, p01.pos[2]], $
                  yr=p01.yr, ystyle=1, yshowtext=0, xtitle='# of Moss', $
                  fill_color='blue', color='k', $
                  font_name='malgun gothic', font_size=13, font_style=1)
    t05 = text(p05.pos[2]-0.02, p05.pos[3]-0.02, /normal, $
               'm = '+string(mean(par, /nan), f='(f0.2)')+$
               '!c$\sigma$ = '+string(stddev(par, /nan), f='(f0.2)'), $
               vertical_align=1, align=1, font_size=10)
;    stop                  
    w1.save, out_dir+par_names[i]+'.png', resol=200 
    w1.close
  endfor
endif

;-----------------------------------------------------------------

if 0 then begin
  for i=0, 4 do begin
    for j=i+1, 4 do begin
  
      param1 = pars[*, i]
      param2 = pars[*, j]
      real = where(finite(param1) and finite(param2))
      w2 = window(dim=[8d2, 8d2], buffer=0)
      if i lt 3 then err1 = err_si else err1 = (i eq 3) ? err_mg_k : err_mg_h 
      if j lt 3 then err2 = err_si else err2 = (j eq 3) ? err_mg_k : err_mg_h
      err = 0.5*(bytscl(err1) + bytscl(err2))
      
      p21 = scatterplot(param1, param2, /current, pos=[0.15, 0.15, 0.7, 0.7], $
                        rgb_table=colortable(49, /rev), magnitude=err, symbol='o', sym_filled=1, $ 
                        sym_size=0.2, transp=20, $
                        xtitle=titles[i], ytitle=titles[j], $
                        font_name='malgun gothic', font_size=13, font_style=1)
      p21.xr = p21.xr
      p21.yr = p21.yr
      p212 = plot([0, 0], p21.yr, '--2', color='gray', transp=50, over=p21)
      p213 = plot(p21.xr, [0, 0], '--2', color='gray', transp=50, over=p21)
      fit_res = linfit(param1[real], param2[real])
      corr = correlate(param1[real], param2[real])
      p211 = plot(p21.xr, poly(p21.xr, fit_res), over=p21, $
                  '--r', transp=30)

      t211 = text(p21.pos[2]-0.02, p21.pos[3]-0.02, /normal, /current, $
                  'n = ' + string(n_elements(param1), f='(i0)')+ $
                  '!c y = ' + string(fit_res[1], f='(f0.2)') + 'x ' + string(fit_res[0], f='(f+-0.2)') +$
                  '!c cc = '+ string(corr, f='(f0.2)'), $
                  vertical_align=1, align=1, font_size=10)
      hist1 = histogram(param1, loc=xbin, nbin=20)
      p22 = barplot(xbin, hist1, /current, $
                    pos=[p21.pos[0], p21.pos[3]+0.01, p21.pos[2], 0.95], $
                    xshowtext=0, xr=p21.xr, xstyle=1, ytitle='# of moss', $
                    font_name='malgun gothic', font_size=13, font_style=1)
      t22 = text(p22.pos[2]-0.02, p22.pos[3]-0.02, /normal, /current, $
                 'm = '+string(mean(param1, /nan), f='(f0.1)')+ $
                 '!c $\sigma$ = '+string(stddev(param1, /nan), f='(f0.1)'), $
                 vertical_align=1, align=1, font_size=10)
      hist2 = histogram(param2, loc=ybin, nbin=20)
      p23 = barplot(ybin, hist2, /current, /horizontal, $
                    pos=[p21.pos[2]+0.01, p21.pos[1], 0.95, p21.pos[3]], $
                    yshowtext=0, yr=p21.yr, ystyle=1, xtitle='# of moss', $
                    font_name='malgun gothic', font_size=13, font_style=1)
      t23 = text(p23.pos[2]-0.02, p23.pos[3]-0.02, /normal, /current, $
                    'm = '+string(mean(param2, /nan), f='(f0.1)')+ $
                    '!c $\sigma$ = '+string(stddev(param2, /nan), f='(f0.1)'), $
                    vertical_align=1, align=1, font_size=10)
;      stop
      w2.save, out_dir+'corr_'+par_names[i]+'_vs_'+par_names[j]+'.png', resol=200
      w2.close 
    endfor
  endfor
endif
;-----------------------------------------------------------------

if 0 then begin
  w3 = window(dim=[8d2, 8d2])
  param1 = mg_k_v - mg_h_v
  param2 = si_v - mg_k_v 
  err = (bytscl(err_si) + bytscl(err_mg_k) + bytscl(err_mg_h))/3.
  
  p31 = scatterplot(param1, param2, /current, pos=[0.15, 0.15, 0.7, 0.7], $
                    rgb_table=colortable(49, /rev), magnitude=err, symbol='o', sym_filled=1, $
                    sym_size=0.2, transp=20, $
                    xtitle='$v_{k3} - v_{h3}$', ytitle='$v_{Si IV} - v_{k3}$', $
                    font_name='malgun gothic', font_size=13, font_style=1)
  p31.xr = p31.xr
  p31.yr = p31.yr
  p312 = plot([0, 0], p31.yr, '--2', color='gray', transp=50, over=p31)  
  p313 = plot(p31.xr, [0, 0], '--2', color='gray', transp=50, over=p31)
  
  fit_res = linfit(param1[real], param2[real])
  corr = correlate(param1[real], param2[real])
  p311 = plot(p31.xr, poly(p31.xr, fit_res), over=p31, '--r', transp=30)
  t311 = text(p31.pos[2]-0.02, p31.pos[3]-0.02, /normal, /current, $
              'n = ' + string(n_elements(param1), f='(i0)')+ $
              '!c y = ' + string(fit_res[1], f='(f0.2)') + 'x ' + string(fit_res[0], f='(f+-0.2)') +$
              '!c cc = '+ string(corr, f='(f0.2)'), $
              vertical_align=1, align=1, font_size=10)
  hist1 = histogram(param1, loc=xbin, nbin=20)
  p32 = barplot(xbin, hist1, /current, $
                pos=[p31.pos[0], p31.pos[3]+0.01, p31.pos[2], 0.95], $
                xshowtext=0, xr=p31.xr, xstyle=1, ytitle='# of moss', $
                font_name='malgun gothic', font_size=13, font_style=1)
  t32 = text(p32.pos[2]-0.02, p32.pos[3]-0.02, /normal, /current, $
            'm = '+string(mean(param1, /nan), f='(f0.1)')+ $
            '!c $\sigma$ = '+string(stddev(param1, /nan), f='(f0.1)'), $
            vertical_align=1, align=1, font_size=10)
  hist2 = histogram(param2, loc=ybin, nbin=20)
  p33 = barplot(ybin, hist2, /current, /horizontal, $
                pos=[p31.pos[2]+0.01, p31.pos[1], 0.95, p31.pos[3]], $
                yshowtext=0, yr=p31.yr, ystyle=1, xtitle='# of moss', $
                font_name='malgun gothic', font_size=13, font_style=1)
  t33 = text(p33.pos[2]-0.02, p33.pos[3]-0.02, /normal, /current, $
              'm = '+string(mean(param2, /nan), f='(f0.1)')+ $
              '!c $\sigma$ = '+string(stddev(param2, /nan), f='(f0.1)'), $
              vertical_align=1, align=1, font_size=10)
  ;    stop
  w3.save, out_dir+'corr_3_velocities.png', resol=200
  w3.close
endif

;------------------------------------------------------------------------------------

if 0 then begin
  param1 = mg_k_v
  param2 = mg_h_v
  k_ge_0 = (param1 ge 0)
  h_ge_0 = (param2 ge 0)
  h_ge_k = (param2 ge param1)
  h_ge_nek = (param2 ge -param1) 
  quad11 = (h_ge_0 and ~h_ge_k)
  quad12 = (k_ge_0 and h_ge_k)
  quad21 = (~k_ge_0 and h_ge_nek)
  quad22 = (h_ge_0 and ~h_ge_nek)
  quad31 = (~h_ge_0 and h_ge_k)
  quad32 = (~k_ge_0 and ~h_ge_k)
  quad41 = (k_ge_0 and ~h_ge_nek)
  quad42 = (~h_ge_0 and h_ge_nek)
  quad_all = float([[quad11], [quad12], [quad21], [quad22], [quad31], [quad32], [quad41], [quad42]])
  nan_pos = where(~finite(param1) or ~finite(param2))
  quad_all[nan_pos, *] = !values.f_nan 
  
  w4 = window(dim=[8d2, 8d2], buffer=0)
  
  p41 = scatterplot(param1, param2, rgb_table=33, magnitude = phase+0, /current, pos=[0.15, 0.15, 0.7, 0.7], $
                    symbol='.', sym_filled=1, sym_color='r', $
                    sym_size=2, transp=70, xr=[-9, 9], yr=[-9, 9], $
                    xtitle=titles[3], ytitle=titles[4], $
                    font_name='malgun gothic', font_size=13, font_style=1)
  p41.xr = p41.xr
  p41.yr = p41.yr
  p412 = plot([0, 0], p41.yr, '--2', color='gray', transp=50, over=p41)
  p413 = plot(p41.xr, [0, 0], '--2', color='gray', transp=50, over=p41)
  p414 = plot(p41.xr, p41.yr, '--2', color='gray', transp=50, over=p41)
  p415 = plot(p41.xr, -p41.yr, '--2', color='gray', transp=50, over=p41)

  for j=0, 7 do begin
    r = 7.
    theta = !dtor*(360./8.*j+22.5)
    num = n_elements(param1[where(quad_all[*, j] eq 1, /null)])
    t4n = text(r*cos(theta), r*sin(theta), $
               string(num, f='(i0)')+' ('+string(num*100./n_elements(param1), f='(f0.1)')+'%)', $
               target=p41, /data, color='black', vertical=0.5, align=0.5)
  endfor
  
  hist41 = histogram(param1[where(phase lt 0)], loc=xbin, min=-10, max=10, binsize=1)
  hist42 = histogram(param1[where(phase ge 0)], loc=xbin, min=-10, max=10, binsize=1)
  p421 = barplot(xbin, hist41, /current, $
                pos=[p41.pos[0], p41.pos[3]+0.01, p41.pos[2], 0.95], $
                fill_color=p41.rgb_table[*, 0], transp=p41.transp, $
                xshowtext=0, xr=p41.xr, xstyle=1, ytitle='# of moss', $
                font_name='malgun gothic', font_size=13, font_style=1)
  p422 = barplot(xbin, hist41+hist42, /current, over=p421, transp=p41.transp, $
                fill_color=p41.rgb_table[*, -1], bottom_value=hist41)              
  t42 = text(p421.pos[2]-0.02, p421.pos[3]-0.02, /normal, /current, $
            'm = '+string(mean(param1, /nan), f='(f0.1)')+ $
            '!c $\sigma$ = '+string(stddev(param1, /nan), f='(f0.1)'), $
            vertical_align=1, align=1, font_size=10)
  hist43 = histogram(param2[where(phase lt 0)], loc=ybin, min=-10, max=10, binsize=1)
  hist44 = histogram(param2[where(phase ge 0)], loc=ybin, min=-10, max=10, binsize=1)
  
  p431 = barplot(ybin, hist43, /current, /horizontal, $
                  pos=[p41.pos[2]+0.01, p41.pos[1], 0.95, p41.pos[3]], $
                  fill_color=p41.rgb_table[*, 0], transp=p41.transp, $
                  yshowtext=0, yr=p41.yr, ystyle=1, xtitle='# of moss', $
                  font_name='malgun gothic', font_size=13, font_style=1)
  p432 = barplot(ybin, hist43+hist44, /current, over=p431, /horiz, transp=p41.transp, $
                 fill_color=p41.rgb_table[*, -1], bottom_value=hist43)
  t43 = text(p431.pos[2]-0.02, p431.pos[3]-0.02, /normal, /current, $
              'm = '+string(mean(param2, /nan), f='(f0.1)')+ $
              '!c $\sigma$ = '+string(stddev(param2, /nan), f='(f0.1)'), $
              vertical_align=1, align=1, font_size=10)
  ;      stop
  t441 = text(p421.pos[1]+0.03, p421.pos[3]-0.03, 'Early phase', /normal, /current, $
             color=p41.rgb_table[*, 0], transp=p41.transp)
  t442 = text(p421.pos[1]+0.03, p421.pos[3]-0.05, 'Late phase', /normal, /current, $
             color=p41.rgb_table[*, -1], transp=p41.transp)             
  w4.save, out_dir+'mg_k_h_anal.png', resol=200
;  w4.close
endif
;---------------------------------------------------------------------------------

if 0 then begin
  dist_sq = aia_phy[0, *]^2. + aia_phy[1, *]^2
  cos_theta = sqrt(1d0-(dist_sq/960.^2.))
  w5 = window(dim=[8d2, 8d2])
  p5 = scatterplot(cos_theta, si_nth, /current, pos=[0.15, 0.15, 0.9, 0.9], $
                   ytitle='Si IV nonthermal velocity (km s$^{-1}$)', $
                   xtitle='$\mu$ = Cos $\theta$', $
                   sym_filled=1, symbol='o', sym_size=0.2, sym_transp=70, transp=70, $
                   magnitude=tt, rgb_table=39, $
                   font_name='malgun gothic', font_size=13, font_style=1)
                   
  w5.save, out_dir+'Si_IV_nth_mu.png', resol=200
endif

;---------
if 1 then begin
  sym_size = 0.5
  col1 = reform((colortable(33))[0, *])
  col2 = reform((colortable(33))[-1, *])
  w6 = window(dim=[1d3, 1d3])
  p61 = scatterplot3d(si_nth, si_v, si_peak, $
              /current, xtitle='Si V$_{nth}$', ytitle='Si V$_D$', ztitle='Si peak Int.', zr=[-20, 600], $
              sym_object=orb(), sym_size=sym_size, clip=0, transp=70, $
               rgb_table=39, magnitude=findgen(n_elements(phase)))
;  p62 = plot3d(si_nth[where(phase gt 0)], si_v[where(phase gt 0)], si_peak[where(phase gt 0)], $
;              over=p61, linestyle='None', sym_object=orb(), clip=0, transp=70)
  

endif
end