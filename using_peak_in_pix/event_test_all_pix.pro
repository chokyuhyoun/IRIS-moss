dir = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, dir

gather_files = file_search(dir+'moss_param_event_total_mg2.sav')
restore, gather_files



how_to_select = 'using_peak_in_pix'
; using_max / using_mean / using_peak / using_all_pix / using_peak_in_pix
obj_arr = pars_pix_peak
; pars_ev_max / pars_ev_mean / pars_ev_peak / pars / pars_pix_peak 
save_dir = dir+how_to_select+'/'
file_mkdir, save_dir


restore, '/Users/khcho/Desktop/IRIS-moss-main/simulation/sim_res.sav'
model_order = ['C1', 'E1', 'E2', 'E3', 'C1+', 'E1+', 'E2+', 'E3+', $
               'E4', 'E5', 'E6', 'H1', 'H2']

sim_col0 = transpose(colortable(38, ncolor=10))
sim_col = fltarr(3, n_elements(model_order))
sim_col[*, 0:3] = sim_col0[*, 0:3]
sim_col[*, 4:7] = sim_col0[*, 0:3]
sim_col[*, 8:*] = sim_col0[*, 4:-2]
  
  
overall = 1
overall_fe18 = 0
peak_prev = 0
after_prev = 0
all_scatter = 0 
all_density = 0
diff_density = 0
mg_ii_kh = 1
mg_si = 1
abs_siv_nth = 1
mu_nth = 1

;one_peak = where(~finite(pars.mg_h3v) or ~finite(pars.mg_k3v))   
;pars = pars[one_peak, *]

;  par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'cur_fwhm', $ ; 4
;               'si_v', 'pre_si_v', 'si_nth', 'pre_si_nth', 'si_amp', 'pre_si_peak', $ ; 10
;               'si_i_tot', 'pre_si_i_tot', 'log_den', 'pre_log_den', $ ; 14
;               'mg_h3v', 'pre_mg_h3v', 'mg_k3v', 'pre_mg_k3v', $ ; 18
;               'mg_trip', 'pre_mg_trip', $ ; 20 
;               'fe18_avg'] ; 21

selected_par = ['sol_x', 'fe18_avg', 'log_den', 'cur_fwhm', $
                'si_i_tot', 'si_amp', 'si_nth', 'si_v', $
                'mg_k_centr', 'mg_h_centr', 'mg_trip']
params = intarr(n_elements(selected_par))
for ii=0, n_elements(selected_par)-1 do params[ii] = where(par_names eq selected_par[ii], /null)                

par_titles2 = par_titles
for ii=0, n_elements(par_titles2)-1 do $
  par_titles2[ii] = strjoin(strsplit(par_titles[ii], '(', /extract), '!c!c(')

ars = fix((obj_arr.ar_no)[uniq(obj_arr.ar_no)])
col = transpose((colortable(6, ncol=n_elements(ars)+1))[1:*, *])
transp = 30
diff_dr = [[0d0, 0], [0, 0], [-3, 3], [0, 0], $
           [-1d3, 1d3], [-40, 40], [-30, 30], [-30, 30], $
           [-10, 10], [-10, 10], [-0.1, 0.1]]

;----------------------------------------------------------------------
if overall then begin
  w1_sz = [10d2, 12d2]
  nx = 2.
  ny = ceil(n_elements(params)/nx)
  xmar = 100
  ymar = 50
  xgap = 100
  ygap = 75
  p01_xs = (w1_sz[0] - 2*xmar - (nx-1)*xgap)/nx
  p01_ys = (w1_sz[1] - 2*ymar - (ny-1)*ygap)/ny
  p01_x0 = xmar+indgen(nx)*(p01_xs + xgap)
  p01_y0 = w1_sz[1] - ymar - (indgen(ny)+1)*(p01_ys+ygap) + ygap + 10
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
                    yminor=1, $
                    xtitle=par_titles[k], ytitle='Number', font_size=13)
      for l=0, n_elements(ars)-1 do begin
        selected = (obj_arr.(k))[where(obj_arr.(0) eq ars[l])]
;        selected = (obj_arr.(k))[where(obj_arr.ar_no eq ars[l] and (~finite(obj_arr.mg_h3v) or ~finite(obj_arr.mg_k3v)))]
;        stop
        hist = histogram(selected, nbins=50, loc=xbin, $
                         min=par_dr[0, k], max=par_dr[1, k])
        p011 = barplot(xbin+(xbin[1]-xbin[0])*0.5, hist0+hist, bottom_value=hist0, over=p01[k0], $
                       fill_color=col[*, l], linestyle=' ', xstyle=1, transp=transp)
        hist0 = hist0 + hist
        if k0 ne 0 then begin
          mean_pos = p01[k0].convertcoord(median(selected), p01[k0].yr[1], /data, /to_normal)
          p012 = plot(mean_pos[0]*[1, 1], mean_pos[1]+[0.005, 0.01], $
                      over=p00, color=col[*, l], thick=1.5, transp=transp)
        endif
      endfor
      p01[k0].xr = par_dr[*, k]
      p01[k0].yr = p01[k0].yr * 1.4
      if k0 ne 0 then begin
        t01 = text(p01[k0].pos[2]-0.01, p01[k0].pos[3] - 0.01 - (k0 eq 3 ? 0 : 0.035), /current, $
                   'm = '+string(median(obj_arr.(k)), f='(f0.2)')+'!c'+ $
                   '$\sigma$ = '+string(stddev(obj_arr.(k), /nan), f='(f0.2)'), $
                   vertical_align=1, align=1, font_size=12, transp=0)
        p013 = plot(median(obj_arr.(k))*[1, 1], p01[k0].yr, over=p01[k0], $
                    ':k', thick=1.5, transp=50) 
      endif
      offset = 0
      for ii=0, n_elements(sim_res)-1 do begin
        model = where(sim_res.model_name eq model_order[ii], /null)
        xval = !null
        if k eq where(par_names eq 'si_amp') then xval = [sim_res[model].si_res.coeff[0]]
        if k eq where(par_names eq 'si_i_tot') then xval = [sim_res[model].si_res.i_tot]
        if k eq where(par_names eq 'si_nth') then xval = [sim_res[model].si_res.v_nth]
        if k eq where(par_names eq 'si_v') then xval = [sim_res[model].si_res.v_d]
        if k eq where(par_names eq 'mg_h_centr') then xval = [sim_res[model].mg_res[1].centr]
        if k eq where(par_names eq 'mg_k_centr') then xval = [sim_res[model].mg_res[0].centr]
;        if k eq where(par_names eq 'mg_h3v') then xval = [sim_res[model].mg_res[1].v_d_3]
;        if k eq where(par_names eq 'mg_k3v') then xval = [sim_res[model].mg_res[0].v_d_3]
        
        if k eq where(par_names eq 'mg_trip') then xval = [sim_res[model].mg_res[2].emiss]
        plus = strmatch(model_order[ii], '*+') ? 1 : 0 
        if n_elements(xval) eq 0 then continue
        if xval[0] ge p01[k0].xr[0] and xval[0] le p01[k0].xr[1] then begin
          p018 = plot(xval, [p01[k0].yr[1]*(0.78-plus*0.1+0.01*ii)], over=p01[k0], $
                      symbol=plus ? 1 : 24, $
                      sym_size=0.5, sym_filled=0, sym_color=sim_col[*, ii], $
                      transp=transp, sym_thick=1.5)
        endif else begin
          tem_xpos = (p01[k0].xr[1]-p01[k0].xr[0])*0.92+p01[k0].xr[0]
          tem_ypos = p01[k0].yr[1]*(0.8-offset*0.03)
          p018 = plot([tem_xpos], [tem_ypos], over=p01[k0], $
                      symbol=plus ? 1 : 24, $
                      sym_size=0.5, sym_filled=0, sym_color=sim_col[*, ii], $
                      transp=transp, sym_thick=1.5)
          p019 = arrow(tem_xpos+[0, (p01[k0].xr[1]-p01[k0].xr[0])*0.04], tem_ypos+[0, 0], target=p01[k0], $
                      /data, color=sim_col[*, ii], head_ind=1, line_thick=1.5, $
                      head_size=0.3, transp=transp)
;          p019 = text(tem_xpos, tem_ypos, '$\rightarrow$', target=p01[k0], /data, $
;                      color=sim_col[*, ii], align=0, vertical_align=0.5, transp=transp, $
;                      font_size=13)
          offset += 1
;          stop
        endelse
      endfor
      if k0 eq 4 or k0 eq 6 then t011x = p01[k0].pos[2]-0.04 else t011x = p01[k0].pos[0]+0.01
      t011 = text(t011x, p01[k0].pos[3]-0.01, '('+string(byte(97+k0))+')', $
                  font_size=13, font_style=0, align=0, vertical_align=1)
    endfor
  endfor
  for l=0, n_elements(ars)-1 do begin
    t011 = text(p01[-2].pos[0], p01[-1].pos[3]-0.02*l, $
              'AR '+string(ars[l], f='(i0)')+' : '+$
              string(total(obj_arr.ar_no eq ars[l]), f='(i4)'), $
               /current, color=col[*, l], font_size=13, align=0, vertical_align=1, transp=transp)
  endfor
  for ii=0, n_elements(model_order)-1 do begin
    plus = strmatch(model_order[ii], '*+') ? 1 : 0
    t018_xpos = p01[-2].pos[0] + 0.2 + 0.11*(ii/8) + 0.05*plus
    t018_ypos = p01[-1].pos[3] - (ii mod 4)*0.02 - (ii eq 12 ? 4 : 0)*0.02
    t018 = text(t018_xpos, t018_ypos, /current, $
                model_order[ii], font_size=13, font_color=sim_col[*, ii], $
                align=0, vertical_align=1, transp=transp)
  endfor
  stop
;  t010 = text(0.01, 0.99, how_to_select, align=0, vertical_align=1, font_size=11)
;  w1.save, save_dir+'total_histogram.pdf', resol=200
;  w1.save, save_dir+'total_histogram.pdf', page_size=w1.dimen/1d2, width=w1.dimen[0]/1d2, /bitmap
endif

;----------------------------------------------------------------------
if overall_fe18 then begin
  for nn = 0, 1 do begin
    fe18_dr = -0.5 + 4*([0., 1.]+nn)
    fe18_dr = nn eq 0 ? [-0.5, 4] : [4, 15]

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
        fe18_bin_sel = (obj_arr.(k))[where((obj_arr.fe18_avg ge fe18_dr[0]) and $
                                           (obj_arr.fe18_avg lt fe18_dr[1]), /null)]
        for l=0, n_elements(ars)-1 do begin
          selected = [(obj_arr.(k))[where((obj_arr.ar_no eq ars[l]) and $
                                         (obj_arr.fe18_avg ge fe18_dr[0]) and $
                                         (obj_arr.fe18_avg lt fe18_dr[1]), /null)]]
;          stop                               
          if n_elements(selected) ne 0 then begin
            hist = histogram(selected, nbins=50, loc=xbin, $
              min=par_dr[0, k], max=par_dr[1, k])
            p011 = barplot(xbin+(xbin[1]-xbin[0])*0.5, hist0+hist, bottom_value=hist0, over=p01[k0], $
              fill_color=col[*, l], linestyle=' ', xstyle=1, transp=transp)
            hist0 = hist0 + hist
            if k0 ne 0 then begin
              mean_pos = p01[k0].convertcoord(median(selected), p01[k0].yr[1], /data, /to_normal)
              p012 = plot(mean_pos[0]*[1, 1], mean_pos[1]+[0.005, 0.01], $
                over=p00, color=col[*, l], thick=2, transp=transp)
            endif
          endif
        endfor
        p01[k0].xr = par_dr[*, k]
        if k0 ne 0 then begin
          t01 = text(p01[k0].pos[2]-0.01, p01[k0].pos[3]-0.01, /current, $
            'm = '+string(median(fe18_bin_sel), f='(f0.2)')+'!c'+ $
            '$\sigma$ = '+string(stddev(fe18_bin_sel, /nan), f='(f0.2)'), $
            vertical_align=1, align=1, font_size=12, transp=transp)
          p013 = plot(median(fe18_bin_sel)*[1, 1], p01[k0].yr, over=p01[k0], $
            ':k2', transp=50)
        endif
      endfor
    endfor
    for l=0, n_elements(ars)-1 do begin
      t011 = text(p01[-2].pos[0], p01[-1].pos[3]-0.02*l, $
        'AR '+string(ars[l], f='(i0)')+' : '+$
        string(total(obj_arr.ar_no eq ars[l]), f='(i4)'), $
        /current, color=col[*, l], font_size=13, align=0, vertical_align=1, transp=transp)
    endfor
;    t010 = text(0.01, 0.99, how_to_select, align=0, vertical_align=1, font_size=11)
    title_string = string(nn, f='(i0)')
    title_string = nn eq 0 ? 'fe18 < 4' : 'fe18 > 4' 
;    w1.save, save_dir+'total_histogram_'+title_string+'.png', resol=200
    w1.close 
  endfor
endif

;----------------------------------------------------------------------
if peak_prev then begin
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
      if k0 le 1 or k0 eq 3 or k0 ge n_elements(params) then continue
      k = params[k0]
      hist0 = 0
      p01[k0] = plot(indgen(2), /nodata, /current, /dev, $
        pos=[p01_x0[i], p01_y0[j], p01_x0[i]+p01_xs, p01_y0[j]+p01_ys], $
        xtitle='Event - Prev. '+par_titles[k], ytitle='Number', font_size=13)
      diff = obj_arr.(k) - obj_arr.(k+1)
      for l=0, n_elements(ars)-1 do begin
        selected = diff[where(obj_arr.ar_no eq ars[l])]
        hist = histogram(selected, nbins=50, loc=xbin, $
                         min=diff_dr[0, k0], max=diff_dr[1, k0])
        p011 = barplot(xbin+(xbin[1]-xbin[0])*0.5, hist0+hist, bottom_value=hist0, over=p01[k0], $
          fill_color=col[*, l], linestyle=' ', xstyle=1, transp=transp)
        hist0 = hist0 + hist
        if k0 ne 0 then begin
          mean_pos = p01[k0].convertcoord(median(selected), p01[k0].yr[1], /data, /to_normal)
          p012 = plot(mean_pos[0]*[1, 1], mean_pos[1]+[0.005, 0.01], $
            over=p00, color=col[*, l], thick=2, transp=transp)
        endif
      endfor
      p01[k0].xr = diff_dr[*, k0]
      if k0 ne 0 then begin
        t01 = text(p01[k0].pos[2]-0.01, p01[k0].pos[3]-0.01, /current, $
          'm = '+string(median(diff), f='(f0.2)')+'!c'+ $
          '$\sigma$ = '+string(stddev(diff, /nan), f='(f0.2)'), $
          vertical_align=1, align=1, font_size=12, transp=transp)
        p013 = plot(median(diff)*[1, 1], p01[k0].yr, over=p01[k0], $
          ':k2', transp=50)
      endif
    endfor
  endfor
;  t010 = text(0.01, 0.99, how_to_select, align=0, vertical_align=1, font_size=11)
;  w1.save, save_dir+'histogram_peak_prev.png', resol=200
endif

;----------------------------------------------------------------------
if after_prev then begin
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
      if k0 le 1 or k0 eq 3 or k0 ge n_elements(params) then continue
      k = params[k0]
      hist0 = 0
      p01[k0] = plot(indgen(2), /nodata, /current, /dev, $
        pos=[p01_x0[i], p01_y0[j], p01_x0[i]+p01_xs, p01_y0[j]+p01_ys], $
        xtitle='After - Prev. '+par_titles[k], ytitle='Number', font_size=13)
      diff = obj_arr.(k+2) - obj_arr.(k+1)
      for l=0, n_elements(ars)-1 do begin
        selected = diff[where(obj_arr.ar_no eq ars[l])]
        hist = histogram(selected, nbins=50, loc=xbin, $
          min=diff_dr[0, k0], max=diff_dr[1, k0])
        p011 = barplot(xbin+(xbin[1]-xbin[0])*0.5, hist0+hist, bottom_value=hist0, over=p01[k0], $
          fill_color=col[*, l], linestyle=' ', xstyle=1, transp=transp)
        hist0 = hist0 + hist
        if k0 ne 0 then begin
          mean_pos = p01[k0].convertcoord(median(selected), p01[k0].yr[1], /data, /to_normal)
          p012 = plot(mean_pos[0]*[1, 1], mean_pos[1]+[0.005, 0.01], $
            over=p00, color=col[*, l], thick=2, transp=transp)
        endif
      endfor
      p01[k0].xr = diff_dr[*, k0]
      if k0 ne 0 then begin
        t01 = text(p01[k0].pos[2]-0.01, p01[k0].pos[3]-0.01, /current, $
          'm = '+string(median(diff), f='(f0.2)')+'!c'+ $
          '$\sigma$ = '+string(stddev(diff, /nan), f='(f0.2)'), $
          vertical_align=1, align=1, font_size=12, transp=transp)
        p013 = plot(median(diff)*[1, 1], p01[k0].yr, over=p01[k0], $
          ':k2', transp=50)
      endif
    endfor
  endfor
;  t010 = text(0.01, 0.99, how_to_select, align=0, vertical_align=1, font_size=11)
;  w1.save, save_dir+'histogram_after_prev.png', resol=200
endif

;stop
;----------------------------------------------------------------------
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
                  xr=par_dr[*, xpar], yr=par_dr[*, ypar], xshowtext=0, yshowtext=0, $
                  xminor=1, yminor=1, xmajor=3, ymajor=3, font_size=10, $
                  xtitle=par_titles2[params[i]], ytitle=par_titles2[params[j]])
      p02n.xtickval = p02n.xtickval[0:-2]
      p02n.ytickval = p02n.ytickval[0:-2]
      p02n.axes[0].showtext = (j eq n_params-1) ? 1 : 0
      p02n.axes[1].showtext = (i eq 0) ? 1 : 0
      p02n.axes[2].showtext = (j eq 0) ? 1 : 0
      p02n.axes[3].showtext = (i eq n_params-1) ? 1 : 0
      for l=0, n_elements(ars)-1 do begin
        xselected = (obj_arr.(xpar))[where(obj_arr.ar_no eq ars[l])]
        yselected = (obj_arr.(ypar))[where(obj_arr.ar_no eq ars[l])]
        p021 = plot(xselected, yselected, ' .', over=p02n, $
          color=col[*, l], transp=50)
      endfor

    endfor
  endfor
  t020 = text(0.01, 0.99, how_to_select, align=0, vertical_align=1, font_size=11)  
;  w2.save, save_dir+'all_scatter.png', resol=300
endif

;----------------------------------------------------------------------
if all_density then begin
  w3_sz = [13d2, 14d2]
  n_params = n_elements(params)
  xmar = [100, 100]
  ymar = [100, 100]
  p03_xs = (w3_sz[0]-total(xmar))/(n_params)
  p03_ys = p03_xs

  w3 = window(dim=w3_sz)
  p03 = objarr(n_params, n_params)
  p0300 = plot(indgen(2), /nodata, axis=0, position=[0, 0, 1, 1], /current, $
               xr=[0, 1], yr=[0, 1])
                
  for j=0, n_params-1 do begin ; y direction
    for i=0, n_params-1 do begin ; x direction
      xpar = params[i]
      ypar = params[j]
;      if i eq 6 then xpar = xpar + 1
;      if j eq 6 then ypar = ypar + 1
      xparam = (obj_arr.(xpar))
      yparam = (obj_arr.(ypar))
      
      nbin = 50
      xbinsize = (par_dr[1, xpar]-par_dr[0, xpar])/nbin
      ybinsize = (par_dr[1, ypar]-par_dr[0, ypar])/nbin
      h2d = hist_2d(xparam, yparam, bin1 = xbinsize, bin2 = ybinsize, $
                    min1 = par_dr[0, xpar], max1 = par_dr[1, xpar], $
                    min2 = par_dr[0, ypar], max2 = par_dr[1, ypar])
      pos = [xmar[0]+i*p03_xs, w3_sz[1]-ymar[0]-(j+1)*p03_ys, $
             xmar[0]+(i+1)*p03_xs, w3_sz[1]-ymar[0]-j*p03_ys]
      p03n = image_kh(h2d/total(h2d, /nan), $
                      par_dr[0, xpar]+findgen(nbin+1)*xbinsize+0.5*xbinsize, $
                      par_dr[0, ypar]+findgen(nbin+1)*ybinsize+0.5*ybinsize, $ 
                      /current, /dev, aspect=0, high_res = 0, max=0.01, min=0, $ 
                      xr=par_dr[*, xpar], yr=par_dr[*, ypar], $
                      rgb_table=22, pos=pos, xminor=1, yminor=1, xmajor=3, ymajor=3, $
                      xtitle=par_titles2[params[i]], ytitle=par_titles2[params[j]], $
                      font_style=0, font_name='Helvetica', font_size=10)
      p03n.xtickval = p03n.xtickval[0:-2]
      p03n.ytickval = p03n.ytickval[0:-2]
      p03n.axes[0].showtext = (j eq n_params-1) ? 1 : 0
      p03n.axes[1].showtext = (i eq 0) ? 1 : 0
      p03n.axes[2].showtext = (j eq 0) ? 1 : 0
      p03n.axes[3].showtext = (i eq n_params-1) ? 1 : 0
      for ii=0, 3 do p03n.axes[ii].tickfont_size=10
;      real = where(finite(xparam) and finite(yparam))
      real = where(xparam gt par_dr[0, xpar] and xparam lt par_dr[1, xpar] and $
                   yparam gt par_dr[0, ypar] and yparam lt par_dr[1, ypar]) 
      cc = correlate(xparam[real], yparam[real])
      t03n = text(p03n.pos[0]+0.01, p03n.pos[3]-0.01, string(cc, f='(f5.2)'), $
                  vertical_align=1, align=0)
      
      dev = 0.0
      clab = intarr(5)
      mmar = 0.005
      bthick = 2.5
      btransp = 40
      if (i eq 5 and j eq 4) or (i eq 6 and (j eq 4 or j eq 5)) or (i eq 7 and j eq 6) then begin ;(a) Si IV 
        pbord = plot(p03n.pos[[0, 0, 2, 2, 0]] + dev*[1, 1, -1, -1, 1], $
                      p03n.pos[[1, 3, 3, 1, 1]] + dev*[1, -1, -1, 1, 1], $
                      thick=bthick, color='r', over=p0300, /data, transp=btransp)
        pbord.order, /bring_to_front
        if clab[0] eq 0 then begin
          clab[0] = 1
          ctex = text(p03n.pos[2]-mmar, p03n.pos[1]+mmar, '3.2.1', $
                      color=pbord.color, align=1)
        endif
      endif
      if (i eq 9 and j eq 8) or ((i eq 8 or i eq 9) and (j eq 7)) then begin ;(b) Mg II h&k
        pbord = plot(p03n.pos[[0, 0, 2, 2, 0]] + dev*[1, 1, -1, -1, 1], $
                      p03n.pos[[1, 3, 3, 1, 1]] + dev*[1, -1, -1, 1, 1], $
                      thick=bthick, color='g', over=p0300, /data, transp=btransp)
        pbord.order, /bring_to_front
        if clab[1] eq 0 then begin
          clab[1] = 1
          ctex = text(p03n.pos[2]-mmar, p03n.pos[1]+mmar, '3.2.2', $
            color=pbord.color, align=1)
        endif
      endif
      if (i eq 6 and j eq 0) then begin ;(c) Si nth vs mu
        pbord = plot(p03n.pos[[0, 0, 2, 2, 0]] + dev*[1, 1, -1, -1, 1], $
                      p03n.pos[[1, 3, 3, 1, 1]] + dev*[1, -1, -1, 1, 1], $
                      thick=bthick, color='b', over=p0300, /data, transp=btransp)
        pbord.order, /bring_to_front
        if clab[2] eq 0 then begin
          clab[2] = 1
          ctex = text(p03n.pos[2]-mmar, p03n.pos[1]+mmar, '3.2.3', $
            color=pbord.color, align=1)
        endif
      endif
      if j eq 10 and (i eq 4 or i eq 5 or i eq 7) then begin ;(d) Mg II triplet
        pbord = plot(p03n.pos[[0, 0, 2, 2, 0]] + dev*[1, 1, -1, -1, 1], $
                      p03n.pos[[1, 3, 3, 1, 1]] + dev*[1, -1, -1, 1, 1], $
                      thick=bthick, color='orange', over=p0300, /data, transp=btransp)
        pbord.order, /bring_to_front
        if clab[3] eq 0 then begin
          clab[3] = 1
          ctex = text(p03n.pos[2]-mmar, p03n.pos[1]+mmar, '3.2.4', $
            color=pbord.color, align=1)
        endif
      endif
      if (i eq 4 and j eq 1) or (i eq 10 and j eq 1) then begin ;(e) Fe xviii
        pbord = plot(p03n.pos[[0, 0, 2, 2, 0]] + dev*[1, 1, -1, -1, 1], $
                      p03n.pos[[1, 3, 3, 1, 1]] + dev*[1, -1, -1, 1, 1], $
                      thick=bthick, color='purple', over=p0300, /data, transp=btransp)                      
        pbord.order, /bring_to_front
        if clab[4] eq 0 then begin
          clab[4] = 1
          ctex = text(p03n.pos[2]-mmar, p03n.pos[1]+mmar, '3.2.5', $
            color=pbord.color, align=1)
        endif
      endif                
      p03[i, j] = p03n
;stop 
    endfor
  endfor
  pbord.order, /bring_to_front
  dum = 200 
  cb03 = colorbar(target=p03n, orientation=0, border=1, textpos=0, major=5, /dev, $
                pos=[w3_sz[0]/2-dum, 80, w3_sz[0]/2+dum, 100], $
                title='Joint Probability Density Function')
;  t030 = text(0.01, 0.99, how_to_select, align=0, vertical_align=1, font_size=11)
;  w3.save, save_dir+'all_density.png', resol=300
;  w3.save, save_dir+'all_density.pdf', page_size=w3.dimen/1d2, width=w3.dimen[0]/1d2, /bitmap
endif

;----------------------------------------------------------------------
if diff_density then begin
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
      if (i lt 2) or (i eq 3) or (j lt 2) or (j eq 3) then continue
      xpar = params[i]
      ypar = params[j]
      xparam = (obj_arr.(xpar) - obj_arr.(xpar+1))
      yparam = (obj_arr.(ypar) - obj_arr.(ypar+1))

      nbin = 50d0
      xbinsize = (diff_dr[1, i]-diff_dr[0, i])/nbin
      ybinsize = (diff_dr[1, j]-diff_dr[0, j])/nbin
      h2d = hist_2d(xparam, yparam, bin1 = xbinsize, bin2 = ybinsize, $
                    min1 = diff_dr[0, i], max1 = diff_dr[1, i], $
                    min2 = diff_dr[0, j], max2 = diff_dr[1, j])
      pos = [xmar[0]+i*p03_xs, w3_sz[1]-ymar[0]-(j+1)*p03_ys, $
             xmar[0]+(i+1)*p03_xs, w3_sz[1]-ymar[0]-j*p03_ys]
      p03n = image_kh(bytscl(h2d), $
                      diff_dr[0, i]+findgen((size(h2d))[1])*xbinsize+0.5*xbinsize, $
                      diff_dr[0, j]+findgen((size(h2d))[2])*ybinsize+0.5*ybinsize, $
                      /current, /dev, aspect=0, high_res = 0, $
                      xr=diff_dr[*, i], yr=diff_dr[*, j], $
                      rgb_table=22, pos=pos, xminor=1, yminor=1, xmajor=3, ymajor=3, $
                      xtitle=par_titles2[params[i]], ytitle=par_titles2[params[j]], $
                      font_style=0, font_name='Helvetica', font_size=10)
      p03n.xtickval = p03n.xtickval[0:-2]
      p03n.ytickval = p03n.ytickval[0:-2]
      p03n.axes[0].showtext = (j eq 10) ? 1 : 0
      p03n.axes[1].showtext = (i eq 2) ? 1 : 0
      p03n.axes[2].showtext = (j eq 2) ? 1 : 0
      p03n.axes[3].showtext = (i eq 10) ? 1 : 0
      for ii=0, 3 do p03n.axes[ii].tickfont_size=10
      ;      real = where(finite(xparam) and finite(yparam))
      real = where(xparam gt diff_dr[0, i] and xparam lt diff_dr[1, i] and $
                    yparam gt diff_dr[0, j] and yparam lt diff_dr[1, j])
      cc = correlate(xparam[real], yparam[real])
      t03n = text(p03n.pos[0]+0.01, p03n.pos[3]-0.01, string(cc, f='(f5.2)'), $
        vertical_align=1, align=0)
      ;stop
    endfor
  endfor
;  t030 = text(0.01, 0.99, how_to_select+' + diff.', align=0, vertical_align=1, font_size=11)
;  w3.save, save_dir+'diff_density.png', resol=300
endif

if mg_ii_kh then begin
  w4 = window(dim=[8d2, 8d2])
;par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'si_v', 'si_nth', 'si_amp', 'si_i_tot', 'log_den', 'cur_fwhm', $
;             'mg_h3v', 'mg_k3v', 'mg_trip', 'e_den']  
;  xpar = (where(par_names eq 'mg_h_centr'))[0]
;  ypar = (where(par_names eq 'mg_k_centr'))[0]

  xpar = (where(par_names eq 'mg_h3v'))[0]
  ypar = (where(par_names eq 'mg_k3v'))[0]

  xparam = (obj_arr.(xpar))
  yparam = (obj_arr.(ypar))
  real = where(finite(xparam) and finite(yparam))
  xparam = xparam[real]
  yparam = yparam[real]
  nbin = 200
  dr = [-15, 15]
  xbinsize = (par_dr[1, xpar]-par_dr[0, xpar])/nbin
  ybinsize = (par_dr[1, ypar]-par_dr[0, ypar])/nbin
  xp = par_dr[0, xpar]+findgen(nbin+1)*xbinsize+0.5*xbinsize
  yp = par_dr[0, ypar]+findgen(nbin+1)*ybinsize+0.5*ybinsize  
  h2d = hist_2d(xparam, yparam, bin1 = xbinsize, bin2 = ybinsize, $
                min1 = par_dr[0, xpar], max1 = par_dr[1, xpar], $
                min2 = par_dr[0, ypar], max2 = par_dr[1, ypar])
  im04 = image_kh(h2d/total(h2d), xp, yp, $
                  /current, aspect=0, pos=[0.1, 0.15, 0.85, 0.9], $ 
                  xr=par_dr[*, xpar], yr=par_dr[*, ypar], rgb_table=22, max=0.01, $
;                  title='Joint Probability Distribution Function of $v_{centr, Mg II h}$ and $v_{centr, Mg II k}$', $
                  title='Joint Probability Distribution Function of $v_{D, Mg II h}$ and $v_{D, Mg II k}$', $
                  xtitle=par_titles[xpar], ytitle=par_titles[ypar], $
                  font_style=0, font_name='Helvetica', font_size=13)
  cb04 = colorbar(target=im04, orientation=1, border=1, textpos=1, $
                  pos=im04.pos[[2, 1, 2, 3]]+[0, 0, 0.03, 0])
  p040 = plot([0, 0], im04.yr, ':', over=im04)
  p041 = plot(im04.xr, [0, 0], ':', over=im04)
  p042 = plot([-50, 50], [-50, 50], ':', over=im04)
  p0421 = plot([-50, 50], [50, -50], ':', over=im04)
  down = where(xparam gt 0 and yparam gt 0)
  up = where(xparam lt 0 and yparam lt 0)
;  res1 = reform(poly_fit(xparam[down], yparam[down], 1, yerror=err1))
;  p043 = plot([-50, 0], poly([-50, 0], res1), '-2', over=im04)
;  t043 = text(-6, -8, 'y = '+string(res1[1], f='(f4.2)')+'x'+string(res1[0], f='(f+-5.2)'), $
;              /data, font_size=13)
;  t0431 = text(-9, -7, '$<v_{D, Mg II k3} - v_{D, Mg II h3}>_{med}$ = '+$
;               string(median(yparam[down]-xparam[down]), f='(f-5.2)'),$
;               /data, font_size=13)
;  t0432 = text(-9, -8, '$<v_{D, Mg II k3}/v_{D, Mg II h3}>_{med}$ = '+$
;                string(median(yparam[down]/xparam[down]), f='(f-5.2)'),$
;                /data, font_size=13)
;  h_rho = 200 ;  Hansteen et al. (2010) https://iopscience.iop.org/article/10.1088/0004-637X/718/2/1070/meta
;  rho_ratio = exp(-40./h_rho)
;  delta_e_down = 1.-(1./rho_ratio*(median(xparam[down]/yparam[down]))^2.)
;  print, 'del E/E_0 (down) = '+string(delta_e_down) 
;

;  res2 = reform(poly_fit(xparam[up], yparam[up], 1, yerror=err2))
;;  p044 = plot([0, 50], poly([0, 50], res2), '-2', over=im04)
;;  t044 = text(6, 4.5, 'y = '+string(res2[1], f='(f4.2)')+'x'+string(res2[0], f='(f+-5.2)'), $
;;              /data, font_size=13, clip=0)
;  t0441 = text(0.5, 8.5, '$<v_{D, Mg II k3} - v_{D, Mg II h3}>_{med}$ = '+$
;              string(median(yparam[up]-xparam[up]), f='(f-5.2)'),$
;              /data, font_size=13, clip=0)
;  t0442 = text(0.5, 7.5, '$<v_{D, Mg II k3}/v_{D, Mg II h3}>_{med}$ = '+$
;              string(median(yparam[up]/xparam[up]), f='(f-5.2)'),$
;              /data, font_size=13)
;  delta_e_up = 1.-(rho_ratio*(median(yparam[up]/xparam[up]))^2.)
;  print, 'del E/E_0 (up) = '+string(delta_e_up)

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
    t045 = text(dr[1]*0.5*cos(theta), dr[1]*0.5*sin(theta), $
                '('+string(ii+1, f='(i0)')+') '+$
                string(quad_all[ii]/total(quad_all)*100., f='(f4.1)')+'%', /data, $
                align=0.5)
  endfor
;  t045 = text(im04.pos[0]+0.05, im04.pos[3]-0.05, $
;              'cc = '+string(correlate(xparam, yparam), f='(f4.2)'), $
;              /current, /normal, font_size=13)

                       
  for ii=0, n_elements(model_order)-1 do begin
    model = where(sim_res.model_name eq model_order[ii], /null)
    if n_elements(model) eq 0 then continue
;    xval = [sim_res[model].mg_res[1].centr]
;    yval = [sim_res[model].mg_res[0].centr]
    xval = [sim_res[model].mg_res[1].v_d_3]
    yval = [sim_res[model].mg_res[0].v_d_3]

    plus = strmatch(model_order[ii], '*+') ? 1 : 0
    p48 = plot(xval, yval, over=im04, $
                symbol=plus ? 1 : 24, sym_size=1, sym_filled=0, sym_color=sim_col[*, ii], $
                transp=transp, sym_thick=2)
    t048_xpos = im04.pos[2] - 0.2 + 0.11*(ii/8) + 0.05*plus
    t048_ypos = im04.pos[1] + 0.15 - (ii mod 4)*0.02 - (ii eq 12 ? 4 : 0)*0.02
    t048 = text(t048_xpos, t048_ypos, /current, $
               model_order[ii], font_size=13, font_color=sim_col[*, ii], $
               align=0, vertical_align=1, transp=transp)

;    print, res[2, peak_time[ii]].emiss           
  endfor

;  w4.save, save_dir+'mg_ii_comp.png', resol=300


  w41 = window(dim=[8d2, 8d2])
  hist_down = histogram(abs(yparam[down])-abs(xparam[down]), loc=xbin, binsize=0.2)
  p045 = plot(xbin, hist_down, 'r2', /hist, xr=[-5, 5], /current, $
              xtitle='$|v_{cent, k3}| - |v_{cent, h3}| (km s^{-1})$', ytitle='# of pixels', $
              font_size=13)
  down_med = median(abs(yparam[down])-abs(xparam[down]))
  up_med = median(abs(yparam[up])-abs(xparam[up]))
  hist_up = histogram(abs(yparam[up])-abs(xparam[up]), loc=xbin, binsize=0.2)
  p046 = plot(xbin, hist_up, /hist, 'b2', over=p045)
  p0451 = plot(down_med*[1, 1], p045.yr, '--2r', over=p045)   
  p0461 = plot(up_med*[1, 1], p045.yr, '--2b', over=p045)
  t045 = text(p045.pos[2]-0.05, p045.pos[3]-0.1, '$Downward_{median} = $'+string(down_med, f='(f5.2)'), $
              align=1, vertical_align=1, color='red', font_size=13) 
  t046 = text(p045.pos[2]-0.05, p045.pos[3]-0.05, '$Upward_{median} = $'+string(up_med, f='(f5.2)'), $
              align=1, vertical_align=1, color='blue', font_size=13)
  p047 = plot([0, 0], p045.yr, '--2k', transp=50, over=p045)

;  w41.save, save_dir+'diff_hist.png', resol=300
;  stop
endif 


if mg_si then begin
  w5 = window(dim=[8d2, 8d2])
  ;par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'si_v', 'si_nth', 'si_amp', 'si_i_tot', 'log_den', 'cur_fwhm', $
  ;             'mg_h3v', 'mg_k3v', 'mg_trip', 'e_den']
  xpar = (where(par_names eq 'mg_k_centr'))[0]
  ypar = (where(par_names eq 'si_v'))[0]
  xparam = (obj_arr.(xpar))
  yparam = (obj_arr.(ypar))
  real = where(finite(xparam) and finite(yparam))
  xparam = xparam[real]
  yparam = yparam[real]
  nbin = 50
  xbinsize = 40./nbin
  ybinsize = 40./nbin
  xp = -20.+findgen(nbin+1)*xbinsize+0.5*xbinsize
  yp = -20.+findgen(nbin+1)*ybinsize+0.5*ybinsize
  h2d = hist_2d(xparam, yparam, bin1 = xbinsize, bin2 = ybinsize, $
                min1 = -20, max1 = 20, $
                min2 = -20, max2 = 20)
  im05 = image_kh(h2d/total(h2d), xp, yp, $
    /current, aspect=0, pos=[0.1, 0.15, 0.85, 0.9], $
    xr=[-20, 20], yr=[-20, 20], rgb_table=22, max=0.01, $
    title='Joint Probability Distribution function of $v_{D, Mg II k3}$ and $v_{D, Si IV}$', $
    xtitle=par_titles[xpar], ytitle=par_titles[ypar], $
    font_style=0, font_name='Helvetica', font_size=13)
  cb05 = colorbar(target=im05, orientation=1, border=1, textpos=1, $
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
;  w5.save, save_dir+'mg_si.png', resol=300
endif

if mu_nth then begin
  ind01 = (where(par_names eq 'sol_x'))[0]
  ind02 = (where(par_names eq 'sol_y'))[0]
  ind1 = (where(par_names eq 'si_nth'))[0]
  sin_theta = sqrt(obj_arr.(ind01)^2. + obj_arr.(ind02)^2.)/960.
  mu = sqrt(1.-sin_theta^2.)

  real = where(finite(mu) and finite(obj_arr.(ind1)))
  xparam = mu[real]
  yparam = (obj_arr.(ind1))[real]
  nbin = 50 
  xdr = [0, 1d0]
  ydr = [0., 150.]
  xbinsize = (xdr[1]-xdr[0])/nbin
  ybinsize = (ydr[1]-ydr[0])/nbin
  h2d = hist_2d(xparam, yparam, bin1=xbinsize, bin2=ybinsize, $
                min1 = xdr[0], max1 = xdr[1], $
                min2 = ydr[0], max2 = ydr[1])              
  w7 = window(dim=[5d2, 5d2])
  im07 = image_kh(h2d/total(h2d), $
                  xdr[0]+findgen((size(h2d))[1])*xbinsize+0.5*xbinsize, $
                  ydr[0]+findgen((size(h2d))[2])*ybinsize+0.5*ybinsize, $
                  /current, aspect=0, high_res = 1, $
                  pos=[80, 80, 430, 430], /dev, $
                  xr=xdr, yr=ydr, rgb_table=22, $
                  title='$\mu$ Vs. $v_{nth, Si IV}$', $
                  xtitle='$\mu$ = cos $\theta$', ytitle=par_titles[ind1], $
                  font_style=0, font_name='Helvetica', font_size=13)
  t7 = text(im07.pos[0]+0.05, im07.pos[3]-0.05, $
            'r = '+string(correlate(xparam, yparam), f='(f4.2)'), $
            font_size=13, vertical_align=1, align=0)
  cb07 = colorbar(target=im07, orientation=1, border=1, textpos=1, $
                  pos=im07.pos[[2, 1, 2, 3]]+[0, 0, 0.03, 0])                  
;  w7.save, save_dir+'mu_nth.png', resol=300
;  w7.save, save_dir+'mu_nth.pdf', page_size=w7.dimen/1d2, width=w7.dimen[0]/1d2, /bitmap
endif

end