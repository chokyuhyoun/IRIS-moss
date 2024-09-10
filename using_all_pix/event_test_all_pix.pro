dir = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, dir

gather_files = file_search(dir, 'moss_param_event_total.sav')
restore, gather_files

save_dir = dir+'using_all_pix/'
file_mkdir, save_dir
how_to_select = 'every pixels'
overall = 1 
all_scatter = 0 
all_density = 0
mg_ii_kh = 0
mg_si = 0
abs_siv_nth = 0 
ar_tot_peak = 0
mu_nth = 0

;one_peak = where(~finite(pars[*, 10]) or ~finite(pars[*, 11]))   
;pars = pars[one_peak, *]

;  par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'cur_fwhm', $ ; 4
;               'si_v', 'pre_si_v', 'si_nth', 'pre_si_nth', 'si_peak', 'pre_si_peak', $ ; 10
;               'si_i_tot', 'pre_si_i_tot', 'log_den', 'pre_log_den', $ ; 14
;               'mg_h3v', 'pre_mg_h3v', 'mg_k3v', 'pre_mg_k3v', $ ; 18
;               'mg_trip', 'pre_mg_trip', $ ; 20 
;               'fe18_avg'] ; 21
params = [1, 4, [5:21:2]]
par_titles2 = par_titles
for ii=0, n_elements(par_titles2)-1 do $
  par_titles2[ii] = strjoin(strsplit(par_titles[ii], '(', /extract), '!c!c(')

ars = fix((pars[*, 0])[uniq(pars[*, 0])])
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
        selected = (pars[*, k])[where(pars[*, 0] eq ars[l])]
        hist = histogram(selected, nbins=50, loc=xbin, $
                         min=par_dr[0, k], max=par_dr[1, k])
        p011 = barplot(xbin, hist0+hist, bottom_value=hist0, over=p01[k0], $
                       fill_color=col[*, l], linestyle=' ', xstyle=1, transp=transp)
        hist0 = hist0 + hist
        if k0 ne 0 then begin
          mean_pos = p01[k0].convertcoord(median(selected), p01[k0].yr[1], /data, /to_normal)
          p012 = plot(mean_pos[0]*[1, 1], mean_pos[1]+[0.005, 0.01], $
                      over=p00, color=col[*, l], thick=2, transp=transp)
        endif
      endfor
      p01[k0].xr = par_dr[*, k]
      if k0 ne 0 then begin
        t01 = text(p01[k0].pos[2]-0.01, p01[k0].pos[3]-0.01, /current, $
                   'm = '+string(mean(pars[*, k], /nan), f='(f0.2)')+'!c'+ $
                   '$\sigma$ = '+string(stddev(pars[*, k], /nan), f='(f0.2)'), $
                   vertical_align=1, align=1, font_size=12, transp=transp)
        p013 = plot(mean(pars[*, k], /nan)*[1, 1], p01[k0].yr, over=p01[k0], $
                    ':k2', transp=50) 
      endif
    endfor
  endfor
  for l=0, n_elements(ars)-1 do begin
    t011 = text(p01[-2].pos[0], p01[-1].pos[3]-0.02*l, $
              'AR '+string(ars[l], f='(i0)')+' : '+$
              string(total(pars[*, 0] eq ars[l]), f='(i4)'), $
               /current, color=col[*, l], font_size=13, align=0, vertical_align=1, transp=transp)
  endfor
  t010 = text(0.01, 0.99, how_to_select, align=0, vertical_align=1, font_size=11)
  w1.save, save_dir+'total_histogram.png', resol=200
endif

stop
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
        xselected = (pars[*, xpar])[where(pars[*, 0] eq ars[l])]
        yselected = (pars[*, ypar])[where(pars[*, 0] eq ars[l])]
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

  w3 = window(dim=w3_sz)
  p03 = !null

  for j=0, n_params-1 do begin ; y direction
    for i=0, n_params-1 do begin ; x direction
      xpar = params[i]
      ypar = params[j]
      xparam = (pars[*, xpar])
      yparam = (pars[*, ypar])

      
      
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

if mg_ii_kh then begin
  w4 = window(dim=[8d2, 8d2])
;par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'si_v', 'si_nth', 'si_peak', 'si_i_tot', 'log_den', 'cur_fwhm', $
;             'mg_h3v', 'mg_k3v', 'mg_trip', 'e_den']  
  i = 7
  j = 8
  xpar = params[i]
  ypar = params[j]
  xparam = (pars[*, xpar])
  yparam = (pars[*, ypar])
  real = where(finite(xparam) and finite(yparam))
  xparam = xparam[real]
  yparam = yparam[real]
  nbin = 200
  xbinsize = (par_xr[1, i]-par_xr[0, i])/nbin
  ybinsize = (par_xr[1, j]-par_xr[0, j])/nbin
  xp = par_xr[0, i]+findgen(nbin+1)*xbinsize+0.5*xbinsize
  yp = par_xr[0, j]+findgen(nbin+1)*ybinsize+0.5*ybinsize  
  h2d = hist_2d(xparam, yparam, bin1 = xbinsize, bin2 = ybinsize, $
                min1 = par_xr[0, i], max1 = par_xr[1, i], $
                min2 = par_xr[0, j], max2 = par_xr[1, j])
  im04 = image_kh(h2d/total(h2d), xp, yp, $
                  /current, aspect=0, pos=[0.1, 0.15, 0.85, 0.9], $ 
                  xr=[-10, 10], yr=[-10, 10], rgb_table=22, max=0.005, $
                  title='Joint Probability Distribution function of $v_{D, Mg II h3}$ and $v_{D, Mg II k3}$', $
                  xtitle=par_titles[params[i]], ytitle=par_titles[params[j]], $
                  font_style=0, font_name='Helvetica', font_size=13)
  cb04 = colorbar(target=im04, orientation=1, border=1, textpos=1, $
                  pos=im04.pos[[2, 1, 2, 3]]+[0, 0, 0.03, 0])
  p040 = plot([0, 0], im04.yr, ':', over=im04)
  p041 = plot(im04.xr, [0, 0], ':', over=im04)
  p042 = plot([-50, 50], [-50, 50], ':', over=im04)
  p0421 = plot([-50, 50], [50, -50], ':', over=im04)
  down = where(xparam lt 0 and yparam lt 0)
  res1 = reform(poly_fit(xparam[down], yparam[down], 1, yerror=err1))
;  p043 = plot([-50, 0], poly([-50, 0], res1), '-2', over=im04)
;  t043 = text(-6, -8, 'y = '+string(res1[1], f='(f4.2)')+'x'+string(res1[0], f='(f+-5.2)'), $
;              /data, font_size=13)
  t0431 = text(-9, -7, '$<v_{D, Mg II k3} - v_{D, Mg II h3}>_{med}$ = '+$
               string(median(yparam[down]-xparam[down]), f='(f-5.2)'),$
               /data, font_size=13)
  t0432 = text(-9, -8, '$<v_{D, Mg II k3}/v_{D, Mg II h3}>_{med}$ = '+$
                string(median(yparam[down]/xparam[down]), f='(f-5.2)'),$
                /data, font_size=13)
  h_rho = 200 ;  Hansteen et al. (2010) https://iopscience.iop.org/article/10.1088/0004-637X/718/2/1070/meta
  rho_ratio = exp(-40./h_rho)
  delta_e_down = 1.-(1./rho_ratio*(median(xparam[down]/yparam[down]))^2.)
  print, 'del E/E_0 (down) = '+string(delta_e_down) 

  up = where(xparam gt 0 and yparam gt 0)
  res2 = reform(poly_fit(xparam[up], yparam[up], 1, yerror=err2))
;  p044 = plot([0, 50], poly([0, 50], res2), '-2', over=im04)
;  t044 = text(6, 4.5, 'y = '+string(res2[1], f='(f4.2)')+'x'+string(res2[0], f='(f+-5.2)'), $
;              /data, font_size=13, clip=0)
  t0441 = text(0.5, 8.5, '$<v_{D, Mg II k3} - v_{D, Mg II h3}>_{med}$ = '+$
              string(median(yparam[up]-xparam[up]), f='(f-5.2)'),$
              /data, font_size=13, clip=0)
  t0442 = text(0.5, 7.5, '$<v_{D, Mg II k3}/v_{D, Mg II h3}>_{med}$ = '+$
              string(median(yparam[up]/xparam[up]), f='(f-5.2)'),$
              /data, font_size=13)
  delta_e_up = 1.-(rho_ratio*(median(yparam[up]/xparam[up]))^2.)
  print, 'del E/E_0 (up) = '+string(delta_e_up)

;  print, err1, err2
  quad11 = where(yparam gt 0. and yparam lt xparam, n11)
  quad12 = where(xparam gt 0. and yparam gt xparam, n12)
  quad21 = where(xparam lt 0. and yparam gt -xparam, n21)
  quad22 = where(yparam gt 0. and yparam lt -xparam, n22)
  quad31 = where(yparam lt 0. and yparam gt xparam, n31)
  quad32 = where(xparam lt 0. and yparam lt xparam, n32)
  quad41 = where(xparam gt 0. and yparam lt -xparam, n41)
  quad42 = where(yparam lt 0. and yparam gt -xparam, n42)
  quad_all = float([n11, n12, n21, n22, n31, n32, n41, n42])

  for ii=0, 7 do begin
    theta = !dtor*(360./8.*ii+22.5)
    t045 = text(5.*cos(theta), 5.*sin(theta), $
                string(quad_all[ii]/total(quad_all)*100., f='(f4.1)')+'%', /data, $
                align=0.5)
  endfor
  t045 = text(im04.pos[0]+0.05, im04.pos[3]-0.05, $
              'cc = '+string(correlate(xparam, yparam), f='(f4.2)'), $
              /current, /normal, font_size=13)
  w4.save, save_dir+'mg_ii_comp.png', resol=300
  w41 = window(dim=[8d2, 8d2])
  hist_down = histogram(yparam[down]-xparam[down], loc=xbin, binsize=0.1)
  p045 = plot(xbin, hist_down, 'r2', /hist, xr=[-5, 5], /current, $
              xtitle='$V_{D, k3} - V_{D, h3} (km s^{-1})$', ytitle='# of pixels', $
              font_size=13)
  down_med = median(yparam[down]-xparam[down])
  up_med = median(yparam[up]-xparam[up])
  p0451 = plot(down_med*[1, 1], p045.yr, ':2r', over=p045)            
  hist_up = histogram(yparam[up]-xparam[up], loc=xbin, binsize=0.1)
  p046 = plot(xbin, hist_up, /hist, 'b2', over=p045)
  p0461 = plot(up_med*[1, 1], p045.yr, ':2b', over=p045)
  t045 = text(p045.pos[2]-0.05, p045.pos[3]-0.1, '$Downward_{median} = $'+string(down_med, f='(f5.2)'), $
              align=1, vertical_align=1, color='red', font_size=13) 
  t046 = text(p045.pos[2]-0.05, p045.pos[3]-0.05, '$Upward_{median} = $'+string(up_med, f='(f5.2)'), $
              align=1, vertical_align=1, color='blue', font_size=13)
  p047 = plot([0, 0], p045.yr, ':k', over=p045)
  w41.save, save_dir+'diff_hist.png', resol=300
;  stop
endif 


if mg_si then begin
  w5 = window(dim=[8d2, 8d2])
  ;par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'si_v', 'si_nth', 'si_peak', 'si_i_tot', 'log_den', 'cur_fwhm', $
  ;             'mg_h3v', 'mg_k3v', 'mg_trip', 'e_den']
  i = 8
  j = 1
  xpar = params[i]
  ypar = params[j]
  xparam = (pars[*, xpar])
  yparam = (pars[*, ypar])
  real = where(finite(xparam) and finite(yparam))
  xparam = xparam[real]
  yparam = yparam[real]
  nbin = 200
  xbinsize = (par_xr[1, i]-par_xr[0, i])/nbin
  ybinsize = (par_xr[1, j]-par_xr[0, j])/nbin
  xp = par_xr[0, i]+findgen(nbin+1)*xbinsize+0.5*xbinsize
  yp = par_xr[0, j]+findgen(nbin+1)*ybinsize+0.5*ybinsize
  h2d = hist_2d(xparam, yparam, bin1 = xbinsize, bin2 = ybinsize, $
    min1 = par_xr[0, i], max1 = par_xr[1, i], $
    min2 = par_xr[0, j], max2 = par_xr[1, j])
  im05 = image_kh(h2d/total(h2d), xp, yp, $
    /current, aspect=0, pos=[0.1, 0.15, 0.85, 0.9], $
    xr=[-20, 20], yr=[-20, 20], rgb_table=22, max=0.003, $
    title='Joint Probability Distribution function of $v_{D, Mg II k3}$ and $v_{D, Si IV}$', $
    xtitle=par_titles[params[i]], ytitle=par_titles[params[j]], $
    font_style=0, font_name='Helvetica', font_size=13)
  cb04 = colorbar(target=im05, orientation=1, border=1, textpos=1, $
    pos=im05.pos[[2, 1, 2, 3]]+[0, 0, 0.03, 0])
  p050 = plot([0, 0], im05.yr, ':', over=im05)
  p051 = plot(im05.xr, [0, 0], ':', over=im05)
  p052 = plot([-50, 50], [-50, 50], ':', over=im05)
  p0521 = plot([-50, 50], [50, -50], ':', over=im05)
  down = where(xparam lt 0 and yparam lt 0)
  res1 = reform(poly_fit(xparam[down], yparam[down], 1, yerror=err1))
  ;  p043 = plot([-50, 0], poly([-50, 0], res1), '-2', over=im05)
  ;  t043 = text(-6, -8, 'y = '+string(res1[1], f='(f4.2)')+'x'+string(res1[0], f='(f+-5.2)'), $
  ;              /data, font_size=13)
  t0531 = text(-9*2, -7*2, '$<v_{D, Si IV} - v_{D, Mg II k3}>_{med}$ = '+$
    string(median(yparam[down]-xparam[down]), f='(f-5.2)'),$
    /data, font_size=13)
  t0532 = text(-9*2, -8*2, '$<v_{D, Si IV}/v_{D, Mg II k3}>_{med}$ = '+$
    string(median(yparam[down]/xparam[down]), f='(f-5.2)'),$
    /data, font_size=13)

  up = where(xparam gt 0 and yparam gt 0)
  res2 = reform(poly_fit(xparam[up], yparam[up], 1, yerror=err2))
  ;  p044 = plot([0, 50], poly([0, 50], res2), '-2', over=im05)
  ;  t044 = text(6, 4.5, 'y = '+string(res2[1], f='(f4.2)')+'x'+string(res2[0], f='(f+-5.2)'), $
  ;              /data, font_size=13, clip=0)
  t0541 = text(0.5*2, 8.5*2, '$<v_{D, Si IV} - v_{D, Mg II k3}>_{med}$ = '+$
    string(median(yparam[up]-xparam[up]), f='(f-5.2)'),$
    /data, font_size=13, clip=0)
  t0542 = text(0.5*2, 7.5*2, '$<v_{D, Si IV}/v_{D, Mg II k3}>_{med}$ = '+$
    string(median(yparam[up]/xparam[up]), f='(f-5.2)'),$
    /data, font_size=13)

  print, err1, err2
  quad11 = where(yparam gt 0. and yparam lt xparam, n11)
  quad12 = where(xparam gt 0. and yparam gt xparam, n12)
  quad21 = where(xparam lt 0. and yparam gt -xparam, n21)
  quad22 = where(yparam gt 0. and yparam lt -xparam, n22)
  quad31 = where(yparam lt 0. and yparam gt xparam, n31)
  quad32 = where(xparam lt 0. and yparam lt xparam, n32)
  quad41 = where(xparam gt 0. and yparam lt -xparam, n41)
  quad42 = where(yparam lt 0. and yparam gt -xparam, n42)
  quad_all = float([n11, n12, n21, n22, n31, n32, n41, n42])

  for ii=0, 7 do begin
    theta = !dtor*(360./8.*ii+22.5)
    t055 = text(5.*cos(theta)*2, 5.*sin(theta)*2, $
      string(quad_all[ii]/total(quad_all)*100., f='(f4.1)')+'%', /data, $
      align=0.5)
  endfor
  t055 = text(im05.pos[0]+0.05, im05.pos[3]-0.05, $
      'cc = '+string(correlate(xparam, yparam), f='(f4.2)'), $
      /current, /normal, font_size=13)
  w5.save, save_dir+'mg_si.png', resol=300
endif

if ar_tot_peak then begin
  p06_sz = 200
  xmar = 80
  ymar = 80
  w6 = window(dim=[xmar*2+p06_sz*n_elements(ars), ymar*3+p06_sz*2])
  nbin = 50
  xind = 5
  yind = 7
  xdr = [0, 50]
  ydr = [0, 50]
  xbinsize = (xdr[1]-xdr[0])/nbin
  ybinsize = (ydr[1]-ydr[0])/nbin
  for ii=0, n_elements(ars)-1 do begin
    one_ar = where(pars[*, 0] eq ars[ii])
    width = sqrt((pars[one_ar, xind]/3d5*1403d0)^2.+0.053^2.+0.026*2)
    print, minmax(pars[one_ar, xind]), minmax(width)
    h2d = hist_2d(width*pars[one_ar, xind+1], pars[one_ar, yind], $
                  bin1 = xbinsize, bin2 = ybinsize, $
                  min1 = xdr[0], max1 = xdr[1], $
                  min2 = ydr[0], max2 = ydr[1])
    pos = [xmar+p06_sz*ii, ymar*2+p06_sz, xmar+p06_sz*(ii+1), ymar*2+p06_sz*2]
    im06 = image_kh(bytscl(h2d), xdr[0]+findgen(nbin+1)*xbinsize+0.5*xbinsize, $
                    ydr[0]+findgen(nbin+1)*ybinsize+0.5*ybinsize, $
                    /current, /dev, aspect=0, hi_res = 1, pos = pos, $ 
                    xr=xdr, yr=ydr, rgb_table=22, $
                    title='AR '+string(ars[ii], f='(i0)'), $
                    xtitle=par_titles[xind], ytitle=par_titles[yind], $
                    font_style=0, font_name='Helvetica', font_size=13)
    if ii ne 0 then im06.yshowtext = 0
    if ii ne 4 then im06.xtickval = (im06.xtickval)[0:-2]         
;    stop       
  endfor

  nbin = 50
  xind = 4
  yind = 11
  xdr = par_xr[*, xind-3]
  ydr = par_xr[*, yind-3]
  xbinsize = (xdr[1]-xdr[0])/nbin
  ybinsize = (ydr[1]-ydr[0])/nbin
  for ii=0, n_elements(ars)-1 do begin
    one_ar = where(pars[*, 0] eq ars[ii])
    width = sqrt((pars[one_ar, xind]/3d5*1403d0)^2.+0.053^2.+0.026*2)
    h2d = hist_2d(pars[one_ar, xind], pars[one_ar, yind], $
                  bin1 = xbinsize, bin2 = ybinsize, $
                  min1 = xdr[0], max1 = xdr[1], $
                  min2 = ydr[0], max2 = ydr[1])
    pos = [xmar+p06_sz*ii, ymar, xmar+p06_sz*(ii+1), ymar+p06_sz]
    im06 = image_kh(bytscl(h2d), xdr[0]+findgen(nbin+1)*xbinsize+0.5*xbinsize, $
      ydr[0]+findgen(nbin+1)*ybinsize+0.5*ybinsize, $
      /current, /dev, aspect=0, hi_res = 1, pos = pos, $
      xr=xdr, yr=ydr, rgb_table=22, $
      title='AR '+string(ars[ii], f='(i0)'), $
      xtitle=par_titles[xind], ytitle=par_titles[yind], $
      font_style=0, font_name='Helvetica', font_size=13)
    if ii ne 0 then im06.yshowtext = 0
    if ii ne 4 then im06.xtickval = (im06.xtickval)[0:-2]
  endfor
endif

if mu_nth then begin
  sin_theta = sqrt(pars[*, 1]^2. + pars[*, 2]^2.)/960.
  mu = sqrt(1.-sin_theta^2.)
  real = where(finite(mu) and finite(pars[*, 5]))
  xparam = mu[real]
  yparam = pars[real, 5]
  nbin = 50 
  xdr = [0, 1d0]
  ydr = [0., 100.]
  xbinsize = (xdr[1]-xdr[0])/nbin
  ybinsize = (ydr[1]-ydr[0])/nbin
  h2d = hist_2d(xparam, yparam, bin1=xbinsize, bin2=ybinsize, $
                min1 = xdr[0], max1 = xdr[1], $
                min2 = ydr[0], max2 = ydr[1])              
  w7 = window(dim=[8d2, 8d2])
  im07 = image_kh(h2d, xdr[0]+findgen(nbin+1)*xbinsize+0.5*xbinsize, $
                  ydr[0]+findgen(nbin+1)*ybinsize+0.5*ybinsize, $
                  /current, aspect=0, hi_res = 1, $
                  xr=xdr, yr=ydr, rgb_table=22, $
                  xtitle='Cos $\theta$', ytitle=par_titles[5], $
                  font_style=0, font_name='Helvetica', font_size=13)
  w7.save, save_dir+'sin_theta_nth.png', resol=300
endif

end