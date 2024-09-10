dir = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, dir

gather_files = file_search(dir, 'moss_param_event_total.sav')
restore, gather_files

save_dir = dir
how_to_select = 'every pixels'
;one_peak = where(~finite(pars[*, 10]) or ~finite(pars[*, 11]))
;pars = pars[one_peak, *]

;par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'si_v', 'si_nth', 'si_peak', 'si_i_tot', 'log_den', 'cur_fwhm', $
;             'mg_h3v', 'mg_k3v', 'mg_trip', 'e_den']
par_xr = [[0, 0], [-1000d0, 1000], [0, 0], [0, 0], [-50, 50], [0, 150], [0, 100], [0, 50], [9, 13], [0, 60], $
          [-15, 15], [-15, 15], [-0.3, 0.3], [0, 500]]

params = [1, 9, 4, 5, 6, 7, 10, 11, 12]

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
    if k0 ge n_elements(params) then continue
    k = params[k0]
    hist0 = 0
    p01[k0] = plot(indgen(2), /nodata, /current, /dev, $
      pos=[p01_x0[i], p01_y0[j], p01_x0[i]+p01_xs, p01_y0[j]+p01_ys], $
      xtitle=par_titles[k], ytitle='Number', font_size=13)
    for l=0, n_elements(ars)-1 do begin
      selected = (pars[*, k])[where(pars[*, 0] eq ars[l])]
      hist = histogram(selected, nbins=50, loc=xbin, $
        min=par_xr[0, k], max=par_xr[1, k])
      p011 = barplot(xbin, hist0+hist, bottom_value=hist0, over=p01[k0], $
        fill_color=col[*, l], linestyle=' ', xstyle=1, transp=transp)
      hist0 = hist0 + hist
      if k0 ne 0 then begin
        mean_pos = p01[k0].convertcoord(median(selected), p01[k0].yr[1], /data, /to_normal)
        p012 = plot(mean_pos[0]*[1, 1], mean_pos[1]+[0.005, 0.01], $
          over=p00, color=col[*, l], thick=2, transp=transp)
      endif
    endfor
    p01[k0].xr = par_xr[*, k]
    if k0 ne 0 then begin
      t01 = text(p01[k0].pos[2]-0.01, p01[k0].pos[3]-0.01, /current, $
        'm = '+string(median(pars[*, k]), f='(f0.2)')+'!c'+ $
        '$\sigma$ = '+string(stddev(pars[*, k], /nan), f='(f0.2)'), $
        vertical_align=1, align=1, font_size=12, transp=transp)
      p013 = plot(median(pars[*, k])*[1, 1], p01[k0].yr, over=p01[k0], $
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
w1.save, save_dir+'fig1.pdf', resol=200

end