dir = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, dir

gather_files = file_search(dir, 'moss_param_event_total.sav')
restore, gather_files

save_dir = dir+'using_max/'
file_mkdir, save_dir
how_to_select = 'max value !cin each event'
overall = 1 
all_scatter = 1
all_density = 1
;par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'si_v', 'si_nth', 'si_max', 'si_i_tot', 'log_den', 'cur_fwhm', $
;             'mg_h3v', 'mg_k3v', 'mg_trip', 'e_den']
par_xr = [[-1000, 1000], [-50, 50], [0, 200], [0, 200], [0, 4d3], [9, 13], $
          [0, 60], [-30, 30], [-30, 30], [-0.3, 0.3], [0, 600]] 
params = [1, [4:n_elements(par_names)-1]]

ars = fix(pars_event_max.ar_no[uniq(pars_event_max.ar_no)])
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
        selected = pars_event_max.(k)[where(pars_event_max.ar_no eq ars[l])]
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
                   'm = '+string(mean(pars_event_max.(k), /nan), f='(f0.2)')+'!c'+ $
                   '$\sigma$ = '+string(stddev(pars_event_max.(k), /nan), f='(f0.2)'), $
                   vertical_align=1, align=1, font_size=12, transp=transp)
        p013 = plot(mean(pars_event_max.(k), /nan)*[1, 1], p01[k0].yr, over=p01[k0], $
                    ':k2', transp=50) 
      endif
    endfor
  endfor
  for l=0, n_elements(ars)-1 do begin
    t011 = text(p01[-2].pos[0], p01[-1].pos[3]-0.02*l, $
              'AR '+string(ars[l], f='(i0)')+' : '+$
              string(total(pars_event_max.ar_no eq ars[l]), f='(i4)'), $
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
        xselected = pars_event_max.(xpar)[where(pars_event_max.ar_no eq ars[l])]
        yselected = pars_event_max.(ypar)[where(pars_event_max.ar_no eq ars[l])]
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
      xparam = pars_event_max.(xpar)
      yparam = pars_event_max.(ypar)
      
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