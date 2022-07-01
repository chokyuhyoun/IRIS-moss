dir = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, dir
sav_files = file_search(dir, '*2020*/*event*.sav', /fully)
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
moss_num = !null
err_si = !null
err_mg_k = !null
err_mg_h = !null

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
      err_mg_h = [err_si, reform(eout.mg_ii_fit_res[0, *].chisq)]
      err_mg_k = [err_si, reform(eout.mg_ii_fit_Res[1, *].chisq)]
    endif
    moss_num = [moss_num, eout.moss_num]
  endfor
  save, t0, t1, tr, tt, ind0, ind1, si_peak, si_v, si_nth, mg_k_v, mg_h_v, filename='gethering_moss_param.sav'
  eout =0 
endif else restore, 'gethering_moss_param.sav'

titles = ['Si IV Peak Intensity (Photons pix$^{-1}$ s$^{-1}$)', $
          'Si IV Doppler Velocity (km s$^{-1}$)', $
          'Si IV Nonthermal Velocity (km s$^{-1}$)', $
          'Mg II k3 Doppler Velocity (km s$^{-1}$)', $
          'Mg II h3 Doppler Velocity (km s$^{-1}$)']
pars = [[si_peak], [si_v], [si_nth], [mg_k_v], [mg_h_v]]
par_names = ['si_peak', 'si_v', 'si_nth', 'mg_k_v', 'mg_h_v']  

;-----------------------------------------------------------------------------------
if 1 then begin 
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
  ;  stop                  
    w1.save, out_dir+par_names[i]+'.png', resol=200 
    w1.close
  endfor
endif

;-----------------------------------------------------------------

for i=0, 4 do begin
  for j=i+1, 4 do begin

    param1 = pars[*, i]
    param2 = pars[*, j]
    real = where(finite(param1) and finite(param2))
    w2 = window(dim=[8d2, 8d2], buffer=1)
    if i lt 3 then err1 = err_si else err1 = (i eq 3) ? err_mg_k : err_mg_h 
    if j lt 3 then err2 = err_si else err2 = (j eq 3) ? err_mg_k : err_mg_h
    err = 0.5*(bytscl(err1) + bytscl(err2))
    
    p01 = scatterplot(param1, param2, /current, pos=[0.15, 0.15, 0.7, 0.7], $
                      rgb_table=colortable(49, /rev), magnitude=err, symbol='o', sym_filled=1, $ 
                      sym_size=0.2, transp=20, $
                      xtitle=titles[i], ytitle=titles[j], $
                      font_name='malgun gothic', font_size=13, font_style=1)
    fit_res = linfit(param1[real], param2[real])
    corr = correlate(param1[real], param2[real])
    p011 = plot(p01.xr, poly(p01.xr, fit_res), over=p01, $
                '--r', transp=30)
    t011 = text(p01.pos[2]-0.02, p01.pos[3]-0.02, /normal, /current, $
                'n = ' + string(n_elements(param1), f='(i0)')+ $
                '!c y = ' + string(fit_res[1], f='(f0.2)') + 'x ' + string(fit_res[0], f='(f+-0.2)') +$
                '!c cc = '+ string(corr, f='(f0.2)'), $
                vertical_align=1, align=1, font_size=10)
    hist1 = histogram(param1, loc=xbin, nbin=20)
    p02 = barplot(xbin, hist1, /current, $
                  pos=[p01.pos[0], p01.pos[3]+0.01, p01.pos[2], 0.95], $
                  xshowtext=0, xr=p01.xr, xstyle=1, ytitle='# of moss', $
                  font_name='malgun gothic', font_size=13, font_style=1)
    t02 = text(p02.pos[2]-0.02, p02.pos[3]-0.02, /normal, /current, $
               'm = '+string(mean(param1, /nan), f='(f0.1)')+ $
               '!c $\sigma$ = '+string(stddev(param1, /nan), f='(f0.1)'), $
               vertical_align=1, align=1, font_size=10)
    hist2 = histogram(param2, loc=ybin, nbin=20)
    p03 = barplot(ybin, hist2, /current, /horizontal, $
                  pos=[p01.pos[2]+0.01, p01.pos[1], 0.95, p01.pos[3]], $
                  yshowtext=0, yr=p01.yr, ystyle=1, xtitle='# of moss', $
                  font_name='malgun gothic', font_size=13, font_style=1)
    t03 = text(p03.pos[2]-0.02, p03.pos[3]-0.02, /normal, /current, $
                  'm = '+string(mean(param2, /nan), f='(f0.1)')+ $
                  '!c $\sigma$ = '+string(stddev(param2, /nan), f='(f0.1)'), $
                  vertical_align=1, align=1, font_size=10)
;    stop
    w2.save, out_dir+'corr_'+par_names[i]+'_vs_'+par_names[j]+'.png', resol=200
    w2.close 
  endfor
endfor


end