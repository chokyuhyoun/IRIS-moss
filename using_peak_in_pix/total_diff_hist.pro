dir = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, dir

gather_files = file_search(dir, 'moss_param_event_total.sav')
restore, gather_files

how_to_select = 'using_peak_in_pix'
save_dir = dir+how_to_select+'/'
file_mkdir, save_dir


;one_peak = where(~finite(pars[*, 10]) or ~finite(pars[*, 11]))
;pars = pars[one_peak, *]

;  par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'cur_fwhm', $ ; 4
;               'si_v', 'pre_si_v', 'si_nth', 'pre_si_nth', 'si_peak', 'pre_si_peak', $ ; 10
;               'si_i_tot', 'pre_si_i_tot', 'log_den', 'pre_log_den', $ ; 14
;               'mg_h3v', 'pre_mg_h3v', 'mg_k3v', 'pre_mg_k3v', $ ; 18
;               'mg_trip', 'pre_mg_trip', $ ; 20
;               'fe18_avg'] ; 21
params = [1, 4, [5:21:2]]
diff_dr = [[0d0, 0], [0, 0], [-30, 30], [-30, 30], [-100, 100], [-10, 50], [-2., 2], [-10, 10], [-10, 10], [-0.1, 0.1]] 
par_titles2 = par_titles
for ii=0, n_elements(par_titles2)-1 do $
  par_titles2[ii] = strjoin(strsplit(par_titles[ii], '(', /extract), '!c!c(')

ars = fix((pars[*, 0])[uniq(pars[*, 0])])
col = transpose((colortable(6, ncol=n_elements(ars)+1))[1:*, *])
transp = 30

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
    if k0 le 1 then continue
    if k0 ge n_elements(params)-1 then continue
    k = params[k0]
    hist0 = 0
    p01[k0] = plot(indgen(2), /nodata, /current, /dev, $
                  pos=[p01_x0[i], p01_y0[j], p01_x0[i]+p01_xs, p01_y0[j]+p01_ys], $
                  xr = diff_dr[*, k0], $
                  xtitle=par_titles[k]+'-'+par_titles[k+1], ytitle='Number', font_size=13)
    pars_diff = pars_pix_peak[*, k] - pars_pix_peak[*, k+1]
    if k0 eq 6 then pars_diff = alog10(10.^pars_pix_peak[*, k] - 10.^pars_pix_peak[*, k+1])
    for l=0, n_elements(ars)-1 do begin
      sel_ind = where(pars_pix_peak[*, 0] eq ars[l])
      selected = pars_diff[sel_ind]
;      sel_ind = (pars[*, k])[where(pars[*, 0] eq ars[l])]
      
      hist = histogram(selected, nbins=50, loc=xbin, min=diff_dr[0, k0], max=diff_dr[1, k0])
      p011 = barplot(xbin+(xbin[1]-xbin[0])*0.5, hist0+hist, bottom_value=hist0, over=p01[k0], $
        fill_color=col[*, l], linestyle=' ', xstyle=1, transp=transp)
      hist0 = hist0 + hist
      if k0 ne 0 then begin
        mean_pos = p01[k0].convertcoord(mean(selected, /nan), p01[k0].yr[1], /data, /to_normal)
        p012 = plot(mean_pos[0]*[1, 1], mean_pos[1]+[0.005, 0.01], $
          over=p00, color=col[*, l], thick=2, transp=transp)
      endif
;      stop
    endfor
    if k0 ne 0 then begin
      t01 = text(p01[k0].pos[2]-0.01, p01[k0].pos[3]-0.01, /current, $
        'med = '+string(mean(pars_diff, /nan), f='(f0.2)')+'!c'+ $
        '$\sigma$ = '+string(stddev(pars_diff, /nan), f='(f0.2)'), $
        vertical_align=1, align=1, font_size=12, transp=transp)
      p013 = plot(mean(pars_diff, /nan)*[1, 1], p01[k0].yr, over=p01[k0], $
        ':k2', transp=50)
    endif
;    stop
  endfor
endfor
;for l=0, n_elements(ars)-1 do begin
;  t011 = text(p01[-2].pos[0], p01[-1].pos[3]-0.02*l, $
;    'AR '+string(ars[l], f='(i0)')+' : '+$
;    string(total(pars[*, 0] eq ars[l]), f='(i4)'), $
;    /current, color=col[*, l], font_size=13, align=0, vertical_align=1, transp=transp)
;endfor
t010 = text(0.01, 0.99, how_to_select, align=0, vertical_align=1, font_size=11)
w1.save, save_dir+'total_diff_histogram.png', resol=200
end