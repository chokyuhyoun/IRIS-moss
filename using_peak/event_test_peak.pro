dir = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, dir

gather_files = file_search(dir, 'moss_param_event_total.sav')
restore, gather_files

save_dir = dir+'using_peak/'
file_mkdir, save_dir
how_to_select = 'peak pixel !c in each event'
overall = 1
all_scatter = 1
all_density = 1
;par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'si_v', 'si_nth', 'si_peak', 'si_i_tot', 'log_den', 'cur_fwhm', $
;             'mg_h3v', 'mg_k3v', 'mg_trip', 'e_den']
par_xr = [[-1000, 1000], [-50, 50], [0, 200], [0, 200], [0, 4d3], [9, 13], $
          [0, 60], [-30, 30], [-30, 30], [-0.3, 0.3], [0, 600]] 
params = [1, [4:n_elements(par_names)-1]]

ars = fix(pars_event_peak.ar_no[uniq(pars_event_peak.ar_no)])
col = transpose((colortable(6, ncol=n_elements(ars)+1))[1:*, *])
transp = 30

if overall then begin
  w1_sz = [10d2, 12d2]
  nx = 2.
  ny = ceil(n_elements(params)/nx)
  xmar = 100
  ymar = 100
  xgap = 100
  ygap = 75
  p01_xs = (w1_sz[0] - 2*xmar - (nx-1)*xgap)/nx
  p01_ys = (w1_sz[1] - 2*ymar - (ny-1)*ygap)/ny
  p01_x0 = xmar+indgen(nx)*(p01_xs + xgap)
  p01_y0 = w1_sz[1] - ymar - (indgen(ny)+1)*(p01_ys+ygap) + ygap
  w1 = window(dim=w1_sz)
  p00 = plot(indgen(2), /nodata, axis=0, /current, pos=[0, 0, 1, 1], xr=[0, 1], yr=[0, 1])
  p01 = objarr(n_Elements(params))
  for j=0, ny-1 do begin
    for i=0, nx-1 do begin
      k0 = j*nx + i
      if k0 ge n_elements(params) then continue
      k = params[k0]
      hist0 = 0
      p01[k0] = plot(indgen(2), /nodata, /current, /dev, $
                    pos=[p01_x0[i], p01_y0[j], p01_x0[i]+p01_xs, p01_y0[j]+p01_ys], $
                    xtitle=par_titles[k], ytitle='Number', font_size=13)
      for l=0, n_elements(ars)-1 do begin
        selected = pars_event_peak.(k)[where(pars_event_peak.ar_no eq ars[l])]
        hist = histogram(selected, nbins=50, loc=xbin, $
                         min=par_xr[0, k0], max=par_xr[1, k0])
        p011 = barplot(xbin, hist0+hist, bottom_value=hist0, over=p01[k0], $
                       fill_color=col[*, l], linestyle=' ', xstyle=1, transp=transp)
        hist0 = hist0 + hist
        if k0 ne 0 then begin
          mean_pos = p01[k0].convertcoord(mean(selected, /nan), p01[k0].yr[1], /data, /to_normal)
          p012 = plot(mean_pos[0]*[1, 1], mean_pos[1]+[0.005, 0.01], $
                      over=p00, color=col[*, l], thick=2, transp=transp)
        endif
      endfor
      p01[k0].xr = par_xr[*, k0]
      if k0 ne 0 then begin
        t01 = text(p01[k0].pos[2]-0.01, p01[k0].pos[3]-0.01, /current, $
                   'm = '+string(mean(pars_event_peak.(k), /nan), f='(f0.2)')+'!c'+ $
                   '$\sigma$ = '+string(stddev(pars_event_peak.(k), /nan), f='(f0.2)'), $
                   vertical_align=1, align=1, font_size=12, transp=transp)
        p013 = plot(mean(pars_event_peak.(k), /nan)*[1, 1], p01[k0].yr, over=p01[k0], $
                    ':k2', transp=50) 
      endif
    endfor
  endfor
  for l=0, n_elements(ars)-1 do begin
    t011 = text(p01[-2].pos[0], p01[-1].pos[3]-0.02*l, $
              'AR '+string(ars[l], f='(i0)')+' : '+$
              string(total(pars_event_peak.ar_no eq ars[l]), f='(i4)'), $
               /current, color=col[*, l], font_size=13, align=0, vertical_align=1, transp=transp)
  endfor
  t010 = text(0.01, 0.99, how_to_select, align=0, vertical_align=1, font_size=11)
  w1.save, save_dir+'total_histogram.png', resol=200
endif


if all_scatter then begin
  w2_sz = [13d2, 13d2]
  n_params = n_elements(params)
  xmar = [100, 100]
  ymar = [100, 100]
  p02_xs = (w2_sz[0]-total(xmar))/(n_params)
  p02_ys = (w2_sz[1]-total(ymar))/(n_params)
  par_titles2 = par_titles
  for ii=0, n_elements(par_titles2)-1 do par_titles2[ii] = strjoin(strsplit(par_titles[ii], '(', /extract), '!c!c(')
  w2 = window(dim=w2_sz)
  p02 = !null
  for j=0, n_params-1 do begin ; y direction
    for i=0, n_params-1 do begin ; x direction
      xpar = params[i]
      ypar = params[j]
      pos = [xmar[0]+i*p02_xs, w2_sz[1]-ymar[0]-(j+1)*p02_ys, $
        xmar[0]+(i+1)*p02_xs, w2_sz[1]-ymar[0]-j*p02_ys]
      p02n = plot(indgen(2), /nodata, /current, /dev, pos=pos, $
                  xr=par_xr[*, i], yr=par_xr[*, j], xshowtext=0, yshowtext=0, $
                  xminor=1, yminor=1, xmajor=3, ymajor=3, font_size=10, $
                  xtitle=par_titles2[params[i]], ytitle=par_titles2[params[j]])
      p02n.xtickval = p02n.xtickval[0:-2]
      p02n.ytickval = p02n.ytickval[0:-2]
      p02n.axes[0].showtext = (j eq n_params-1) ? 1 : 0
      p02n.axes[1].showtext = (i eq 0) ? 1 : 0
      p02n.axes[2].showtext = (j eq 0) ? 1 : 0
      p02n.axes[3].showtext = (i eq n_params-1) ? 1 : 0
      for l=0, n_elements(ars)-1 do begin
        xselected = pars_event_peak.(xpar)[where(pars_event_peak.ar_no eq ars[l])]
        yselected = pars_event_peak.(ypar)[where(pars_event_peak.ar_no eq ars[l])]
        p021 = plot(xselected, yselected, ' .', over=p02n, $
          color=col[*, l], transp=50)
      endfor

    endfor
  endfor
  t020 = text(0.01, 0.99, how_to_select, align=0, vertical_align=1, font_size=11)  
  w2.save, save_dir+'all_scatter.png', resol=300
endif

if all_density then begin
  w3_sz = [13d2, 13d2]
  n_params = n_elements(params)
  xmar = [100, 100]
  ymar = [100, 100]
  p03_xs = (w3_sz[0]-total(xmar))/(n_params)
  p03_ys = (w3_sz[1]-total(ymar))/(n_params)
  par_titles2 = par_titles
  for ii=0, n_elements(par_titles2)-1 do par_titles2[ii] = strjoin(strsplit(par_titles[ii], '(', /extract), '!c!c(')
  w3 = window(dim=w3_sz)
  p03= !null  
  for j=0, n_params-1 do begin ; y direction
    for i=0, n_params-1 do begin ; x direction
      xpar = params[i]
      ypar = params[j]
      xparam = pars_event_peak.(xpar)
      yparam = pars_event_peak.(ypar)
      
      nbin = 50
      xbinsize = (par_xr[1, i]-par_xr[0, i])/nbin
      ybinsize = (par_xr[1, j]-par_xr[0, j])/nbin
      h2d = hist_2d(xparam, yparam, bin1 = xbinsize, bin2 = ybinsize, $
                    min1 = par_xr[0, i], max1 = par_xr[1, i], $
                    min2 = par_xr[0, j], max2 = par_xr[1, j])
      pos = [xmar[0]+i*p03_xs, w3_sz[1]-ymar[0]-(j+1)*p03_ys, $
             xmar[0]+(i+1)*p03_xs, w3_sz[1]-ymar[0]-j*p03_ys]
      p03n = image_kh(bytscl(h2d), par_xr[0, i]+findgen(nbin+1)*xbinsize+0.5*xbinsize, $
                      par_xr[0, j]+findgen(nbin+1)*ybinsize+0.5*ybinsize, $ 
                      /current, /dev, aspect=0, hi_res = 0, $ 
                      xr=par_xr[*, i], yr=par_xr[*, j], $
                      rgb_table=22, pos=pos, xminor=1, yminor=1, xmajor=3, ymajor=3, $
                      xtitle=par_titles2[params[i]], ytitle=par_titles2[params[j]], $
                      font_style=0, font_name='Helvetica', font_size=10)
      p03n.xtickval = p03n.xtickval[0:-2]
      p03n.ytickval = p03n.ytickval[0:-2]
      p03n.axes[0].showtext = (j eq n_params-1) ? 1 : 0
      p03n.axes[1].showtext = (i eq 0) ? 1 : 0
      p03n.axes[2].showtext = (j eq 0) ? 1 : 0
      p03n.axes[3].showtext = (i eq n_params-1) ? 1 : 0
      real = where(finite(xparam) and finite(yparam))
      cc = correlate(xparam[real], yparam[real])
      t03n = text(p03n.pos[0]+0.01, p03n.pos[3]-0.01, string(cc, f='(f5.2)'), $
                  vertical_align=1, align=0)
;stop
    endfor
  endfor
  t030 = text(0.01, 0.99, how_to_select, align=0, vertical_align=1, font_size=11)
  w3.save, save_dir+'all_density.png', resol=300
endif





end
;
;
;if sin_theta_and_si_nonthermal then begin
;  w2 = window(dim=[8d2, 8d2])
;  
;  ;par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'si_v', 'si_nth', 'si_peak', 'cur_fwhm', $
;  ;             'mg_h3v', 'mg_k3v', 'mg_trip', 'e_den']
;  
;  tag_names = tag_names(pars_event_peak)
;  ind1 = 1
;  ind2 = 5
;  
;  par1 = pars_event_peak.(ind1)
;  dd = pars_event_peak.sol_x^2. + pars_event_peak.sol_y^2.
;  par1 = dd/950.^2.
;  par2 = pars_event_peak.(ind2)
;  
;  p02 = plot(indgen(2), /nodata, /current, $
;             xtitle='sin $\theta$', ytitle=par_titles[ind2], font_size=13)
;  
;  for l=0, n_elements(ars)-1 do begin
;    selected = where(pars_event_peak.ar_no eq ars[l])
;    sample1 = par1[selected]
;    sample2 = par2[selected]
;    p021 = plot(sample1, sample2, over=p02, $
;                'o', sym_size=0.5, sym_filled=1, sym_color=col[*, l], transp=90, $
;                font_size=13)
;    t021 = text(p02.pos[2]-0.2, p02.pos[3]-0.05-0.03*l, $
;                'AR '+string(ars[l], f='(i0)')+' : '+$
;                string(total(pars_event_peak.ar_no eq ars[l]), f='(i3)'), $
;                /current, color=col[*, l], font_size=13, align=0, vertical_align=1, transp=transp)
;  endfor
;  real = where(finite(par1) and finite(par2))
;  res = poly_fit(par1[real], par2[real], 1)
;  cc = correlate(par1[real], par2[real])
;  xx = [p02.xr[0]:p02.xr[1]:1d-2]
;  p023 = plot(xx, poly(xx, res), over=p02, '--3', color='gray')
;  t023 = text(p02.pos[0]+0.05, p02.pos[3]-0.05, $
;              'y = '+string(res[1], f='(f0.2)')+' x + ' $
;              +string(res[0], f='(f0.2)')+'!ccc = '+string(cc, f='(f0.2)'), $
;              font_size=13, vertical_align=1)
;              
;  w2.save, 'sin theta & si v nth.png', resol=200
;
;endif
;
;if mg_ii_comparison then begin
;  w3 = window(dim=[8d2, 8d2])
;
;  ;par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'si_v', 'si_nth', 'si_peak', 'cur_fwhm', $
;  ;             'mg_h3v', 'mg_k3v', 'mg_trip', 'e_den']
;
;  tag_names = tag_names(pars_event_peak)
;  ind1 = 8
;  ind2 = 9
;
;  par1 = pars_event_peak.(ind1)
;  par2 = pars_event_peak.(ind2)
;  real = where(finite(par1) and finite(par2))
;  par1 = par1[real]
;  par2 = par2[real]
;  
;  p03 = plot(indgen(2), /nodata, /current, $
;    xtitle=par_titles[ind1], ytitle=par_titles[ind2], font_size=13, $
;    xr=20*[-1, 1], yr=20*[-1, 1])
;  p038 = plot(p03.xr, [0, 0], ':', color='gray', transp=20, over=p03)
;  p039 = plot([0, 0], p03.yr, ':', color='gray', transp=20, over=p03)
;  
;  for l=0, n_elements(ars)-1 do begin
;    selected = where(pars_event_peak.ar_no eq ars[l])
;    sample1 = par1[selected]
;    sample2 = par2[selected]
;    p031 = plot(sample1, sample2, over=p03, $
;      'o', sym_size=0.2, sym_filled=1, sym_color=col[*, l], transp=70, $
;      font_size=13)
;    t031 = text(p03.pos[2]-0.2, p03.pos[3]-0.05-0.03*l, $
;      'AR '+string(ars[l], f='(i0)')+' : '+$
;      string(total(pars_event_peak.ar_no[real] eq ars[l]), f='(i3)'), $
;      /current, color=col[*, l], font_size=13, align=0, vertical_align=1, transp=transp)
;  endfor
;  res = poly_fit(par1, par2, 1)
;  cc = correlate(par1, par2)
;  xx = [p03.xr[0]:p03.xr[1]:1d-2]
;;  p033 = plot(xx, poly(xx, res), over=p03, '--3', color='gray')
;  t033 = text(p03.pos[0]+0.05, p03.pos[3]-0.05, $
;    'y = '+string(res[1], f='(f0.2)')+' x ' $
;    +string(res[0], f='(f0.2)')+'!ccc = '+string(cc, f='(f0.2)'), $
;    font_size=13, vertical_align=1)
;  p037 = plot(p03.xr, p03.yr, '-k', color='black', transp=0, over=p03)
;
;  w3.save, 'mg ii comparison.png', resol=200
;
;endif
;
;
;if si_dop_and_em then begin
;  w4 = window(dim=[8d2, 8d2])
;
;  ;par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'si_v', 'si_nth', 'si_peak', 'cur_fwhm', $
;  ;             'mg_h3v', 'mg_k3v', 'mg_trip', 'e_den']
;
;  tag_names = tag_names(pars_event_peak)
;  ind1 = 11
;  ind2 = 4
;
;  par1 = pars_event_peak.(ind1)
;  par2 = pars_event_peak.(ind2)
;  real = where(finite(par1) and finite(par2))
;  par1 = par1[real]
;  par2 = par2[real]
;
;  p04 = plot(indgen(2), /nodata, /current, $
;              xtitle=par_titles[ind1], ytitle=par_titles[ind2], font_size=13)
;;  p048 = plot(p04.xr, [0, 0], ':', color='gray', transp=20, over=p04)
;;  p049 = plot([0, 0], p04.yr, ':', color='gray', transp=20, over=p04)
;
;  for l=0, n_elements(ars)-1 do begin
;    selected = where(pars_event_peak.ar_no eq ars[l])
;    sample1 = par1[selected]
;    sample2 = par2[selected]
;    p041 = plot(sample1, sample2, over=p04, $
;      'o', sym_size=0.2, sym_filled=1, sym_color=col[*, l], transp=70, $
;      font_size=13)
;    t041 = text(p04.pos[2]-0.2, p04.pos[3]-0.05-0.03*l, $
;      'AR '+string(ars[l], f='(i0)')+' : '+$
;      string(total(pars_event_peak.ar_no[real] eq ars[l]), f='(i3)'), $
;      /current, color=col[*, l], font_size=13, align=0, vertical_align=1, transp=transp)
;  endfor
;  res = poly_fit(par1, par2, 1)
;  cc = correlate(par1, par2)
;  xx = [p04.xr[0]:p04.xr[1]:1d-2]
;  p043 = plot(xx, poly(xx, res), over=p04, '--3', color='gray')
;  t043 = text(p04.pos[0]+0.05, p04.pos[3]-0.05, $
;    'y = '+string(res[1], f='(f0.2)')+' x ' $
;    +string(res[0], f='(f+0.2)')+'!ccc = '+string(cc, f='(f0.2)'), $
;    font_size=13, vertical_align=1)
;
;  w4.save, 'Si iv and EM.png', resol=200
;endif
;
;if si_mg_vel then begin
;  w6 = window(dim=[8d2, 8d2])
;
;  ;par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'si_v', 'si_nth', 'si_peak', 'cur_fwhm', $
;  ;             'mg_h3v', 'mg_k3v', 'mg_trip', 'e_den']
;
;
;  tag_names = tag_names(pars_event_peak)
;  ind1 = 9
;  ind2 = 4
;
;  par1 = pars_event_peak.(ind1)
;  par2 = pars_event_peak.(ind2)
;  real = where(finite(par1) and finite(par2))
;  par1 = par1[real]
;  par2 = par2[real]
;
;
;  p06 = plot(indgen(2), /nodata, /current, $
;    xtitle=par_titles[ind1], ytitle=par_titles[ind2], font_size=13)
;  p068 = plot(p06.xr, [0, 0], ':', color='gray', transp=20, over=p06)
;  p069 = plot([0, 0], p06.yr, ':', color='gray', transp=20, over=p06)
;
;  for l=0, n_elements(ars)-1 do begin
;    selected = where(pars_event_peak.ar_no eq ars[l])
;    sample1 = par1[selected]
;    sample2 = par2[selected]
;    p061 = plot(sample1, sample2, over=p06, $
;      'o', sym_size=0.2, sym_filled=1, sym_color=col[*, l], transp=70, $
;      font_size=13)
;    t061 = text(p06.pos[2]-0.2, p06.pos[3]-0.05-0.03*l, $
;      'AR '+string(ars[l], f='(i0)')+' : '+$
;      string(total(pars_event_peak.ar_no[real] eq ars[l]), f='(i3)'), $
;      /current, color=col[*, l], font_size=13, align=0, vertical_align=1, transp=transp)
;  endfor
;  res = poly_fit(par1, par2, 1)
;  cc = correlate(par1, par2)
;  xx = [p06.xr[0]:p06.xr[1]:1d-2]
;    p063 = plot(xx, poly(xx, res), over=p06, '--3', color='gray')
;  t063 = text(p06.pos[0]+0.05, p06.pos[3]-0.05, $
;    'y = '+string(res[1], f='(f0.2)')+' x ' $
;    +string(res[0], f='(f0.2)')+'!ccc = '+string(cc, f='(f0.2)'), $
;    font_size=13, vertical_align=1)
;;  p067 = plot(p06.xr, p06.yr, '-k', color='black', transp=0, over=p06)
;  p068 = plot(p06.xr, [0, 0], ':', color='gray', transp=20, over=p06)
;  p069 = plot([0, 0], p06.yr, ':', color='gray', transp=20, over=p06)
;  w6.save, 'Si vel_ Mg II vel.png', resol=200
;
;endif
;
;
;if all then begin
;  w5 = window(dim=[8d2, 8d2])
;  params = [1, [4:n_elements(par_names)-1]]
;
;  ;par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'si_v', 'si_nth', 'si_peak', 'cur_fwhm', $
;  ;             'mg_h3v', 'mg_k3v', 'mg_trip', 'e_den']
;  for i=0, n_elements(params)-1 do begin
;    for j=i+1, n_elements(params)-1 do begin
;      ind1 = params[i]
;      ind2 = params[j]
;
;      par1 = pars_event_peak.(ind1)
;      par2 = pars_event_peak.(ind2)
;      real = where(finite(par1) and finite(par2))
;      par1 = par1[real]
;      par2 = par2[real]
;
;      p05 = plot(indgen(2), /nodata, /current, $
;        xtitle=par_titles[ind1], ytitle=par_titles[ind2], font_size=13)
;      ;  p058 = plot(p05.xr, [0, 0], ':', color='gray', transp=20, over=p05)
;      ;  p059 = plot([0, 0], p05.yr, ':', color='gray', transp=20, over=p05)
;
;      for l=0, n_elements(ars)-1 do begin
;        selected = where(pars_event_peak.ar_no eq ars[l])
;        sample1 = par1[selected]
;        sample2 = par2[selected]
;        p051 = plot(sample1, sample2, over=p05, $
;          'o', sym_size=0.2, sym_filled=1, sym_color=col[*, l], transp=70, $
;          font_size=13)
;        t051 = text(p05.pos[2]-0.2, p05.pos[3]-0.05-0.03*l, $
;          'AR '+string(ars[l], f='(i0)')+' : '+$
;          string(total(pars_event_peak.ar_no[real] eq ars[l]), f='(i3)'), $
;          /current, color=col[*, l], font_size=13, align=0, vertical_align=1, transp=transp)
;      endfor
;      res = poly_fit(par1, par2, 1)
;      cc = correlate(par1, par2)
;      xx = [p05.xr[0]:p05.xr[1]:1d-2]
;      p053 = plot(xx, poly(xx, res), over=p05, '--3', color='gray')
;      t053 = text(p05.pos[0]+0.05, p05.pos[3]-0.05, $
;        'y = '+string(res[1], f='(f0.2)')+' x ' $
;        +string(res[0], f='(f+0.2)')+'!ccc = '+string(cc, f='(f0.2)'), $
;        font_size=13, vertical_align=1)
;
;
;      w5.save, par_names[ind1]+'_'+par_names[ind2]+'.png', resol=200
;      w5.erase
;    endfor
;  endfor
;
;endif
;
;end