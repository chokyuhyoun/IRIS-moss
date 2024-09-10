dir = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, dir

gather_files = file_search(dir+'moss_param_event_total_mg.sav')
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
  
;  par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'cur_fwhm', $ ; 4
;               'si_v', 'pre_si_v', 'si_nth', 'pre_si_nth', 'si_amp', 'pre_si_peak', $ ; 10
;               'si_i_tot', 'pre_si_i_tot', 'log_den', 'pre_log_den', $ ; 14
;               'mg_h3v', 'pre_mg_h3v', 'mg_k3v', 'pre_mg_k3v', $ ; 18
;               'mg_trip', 'pre_mg_trip', $ ; 20 
;               'fe18_avg'] ; 21

;stop

w1 = window(dim=[10d2, 9d2])
xmar = [100, 100]
ymar = [50, 80]
ymar = [50, 80]
nx = 2.
ny = 2.
xs = 320.
ys = xs
xgap = (w1.dimension[0] - total(xmar) - nx*xs)/(nx-1)
ygap = (w1.dimension[1] - total(ymar) - ny*ys)/(ny-1)

pos = fltarr(4, nx*ny)
for ii=0, nx-1 do begin
  for jj=0, ny-1 do begin
    kk = ii + jj*nx
    pos[*, kk] = [xmar[0] + ii*(xs+xgap), $
                  w1.dimension[1] - ymar[0] - (jj+1)*ys - jj*ygap, $
                  xmar[0] + ii*(xs+xgap)+xs, $
                  w1.dimension[1] - ymar[0] - jj*ys - jj*ygap]
  endfor
endfor

nbin = 100
dr = [-15d, 15d]
xbinsize = (dr[1]-dr[0])/nbin
ybinsize = (dr[1]-dr[0])/nbin
xp = dr[0]+findgen(nbin+1)*xbinsize
yp = dr[0]+findgen(nbin+1)*ybinsize


;-------- (a) h centr Vs. k centr from obs.
xpar = (where(par_names eq 'mg_h_centr'))[0]
ypar = (where(par_names eq 'mg_k_centr'))[0]
real = where(finite(obj_arr.(xpar)) and finite(obj_arr.(ypar)))
xparam = (obj_arr.(xpar))[real]
yparam = (obj_arr.(ypar))[real]
centr_hk = hist_2d(xparam, yparam, $
                   bin1 = xbinsize, bin2 = ybinsize, $
                   min1 = dr[0], max1 = dr[1], $
                   min2 = dr[0], max2 = dr[1])
im01 = image_kh(centr_hk/total(centr_hk), xp, yp, /current, /dev, $
                pos=pos[*, 0], rgb_table=22, max=0.02, $
                title='(a) v$_{centr, Mg II h}$ Vs. v$_{centr, Mg II k}$ from observations', $
                xtitle=par_titles[xpar], ytitle=par_titles[ypar], $
                font_style=0, font_name='Helvetica', font_size=13)
cb01 = colorbar(target=im01, orientation=1, border=1, textpos=1, major=5, $
                pos=im01.pos[[2, 1, 2, 3]]+[0, 0, 0.01, 0])
p010 = plot([0, 0], im01.yr, ':', over=im01)
p011 = plot(im01.xr, [0, 0], ':', over=im01)
p012 = plot([-50, 50], [-50, 50], ':', over=im01)
p0121 = plot([-50, 50], [50, -50], ':', over=im01)

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
  t015 = text(dr[1]*0.6*cos(theta), dr[1]*0.7*sin(theta), $
              '('+string(ii+1, f='(i0)')+') '+$
              string(quad_all[ii]/total(quad_all)*100., f='(f4.1)')+'%', /data, $
              align=0.5)
endfor
rere = where(xparam gt -15. and xparam lt 15. and yparam gt -15. and yparam lt 15.)
pcc = correlate(xparam[rere], yparam[rere])
t016 = text(im01.pos[0]+0.015, im01.pos[3]-0.03, 'r = '+string(pcc, f='(f4.2)'), $
            font_size=12)

;--------  (b) h_centr Vs k_centr from simulation

im02 = image_kh(centr_hk, xp, yp, /current, /nodata, /dev, $
                pos=pos[*, 1], rgb_table=22, max=0.02, $
                xtitle=par_titles[xpar], ytitle=par_titles[ypar], $
                title='(b) v$_{centr, Mg II h}$ Vs. v$_{centr, Mg II k}$ from simulations', $
                font_style=0, font_name='Helvetica', font_size=13)
p020 = plot([0, 0], im02.yr, ':', over=im02)
p021 = plot(im02.xr, [0, 0], ':', over=im02)
p022 = plot([-50, 50], [-50, 50], ':', over=im02)
p0221 = plot([-50, 50], [50, -50], ':', over=im02)

for ii=0, n_elements(model_order)-1 do begin
  model = where(sim_res.model_name eq model_order[ii], /null)
  if n_elements(model) eq 0 then continue
  xval = sim_res[model].mg_res[1].centr
  yval = sim_res[model].mg_res[0].centr

  plus = strmatch(model_order[ii], '*+') ? 1 : 0
  
  if abs(xval) gt im02.xr[1] or abs(yval) gt im02.yr[1] then begin
    how_far = sqrt(xval^2.+yval^2.)
    nv = [xval, yval]/how_far
    tt = findgen(how_far)
    xout = min(where(abs(tt*nv[0]) gt im02.xr[1]))
    yout = min(where(abs(tt*nv[1]) gt im02.yr[1]))
    out = min([xout, yout])-4
    xval1 = nv[0]*out
    yval1 = nv[1]*out
    p029 = arrow(xval1+[0, nv[0]*1.5], yval1+[0, nv[1]*1.5], target=im02, $
                 /data, color=sim_col[*, ii], head_ind=1, line_thick=1.5, $
                 head_size=0.3, transp=transp)
  endif else begin
    xval1 = xval
    yval1 = yval
  endelse
  
  p028 = plot([xval1], [yval1], over=im02, $
              symbol=plus ? 1 : 24, sym_size=0.7, sym_filled=0, sym_color=sim_col[*, ii], $
              transp=transp, sym_thick=1.5)
  t028_xpos = im02.pos[2] - 0.11 + 0.07*(ii/8) + 0.03*plus
  t028_ypos = im02.pos[1] + 0.10 - ((ii mod 4) + (ii eq 12 ? 4 : 0))*0.018
  t028 = text(t028_xpos, t028_ypos, /current, target=im02, $
              model_order[ii], font_size=11, font_color=sim_col[*, ii], $
              align=0, vertical_align=1, transp=transp)

  ;    print, res[2, peak_time[ii]].emiss
endfor

;;-------- (c) |h_centr| - |k_centr|  histogram from obs. 
;
;down = where(xparam gt 0 and yparam gt 0)
;up = where(xparam lt 0 and yparam lt 0)
;hist_down = histogram(abs(yparam[down])-abs(xparam[down]), loc=xbin, binsize=0.2)
;p031 = plot(xbin, hist_down, 'r2', /hist, xr=[-3, 3], /current, /dev, $
;           pos = pos[*, 2], title='(c) Absolute difference between two Mg II centroid', $ 
;           xtitle='$|v_{cent, Mg II k}| - |v_{cent, Mg II h}| (km s^{-1})$', ytitle='# of pixels', $
;           font_size=13)
;hist_up = histogram(abs(yparam[up])-abs(xparam[up]), loc=xbin, binsize=0.2)
;p032 = plot(xbin, hist_up, /hist, 'b2', over=p031)
;p031.yr = p031.yr*1.2
;down_med = median(abs(yparam[down])-abs(xparam[down]))
;up_med = median(abs(yparam[up])-abs(xparam[up]))
;p0311 = plot(down_med*[1, 1], p031.yr, '--2r', over=p031)
;p0321 = plot(up_med*[1, 1], p031.yr, '--2b', over=p031)
;t031 = text(p031.pos[2]-0.02, p031.pos[3]-0.01, '$Upward_{med} = $'+string(up_med, f='(f5.2)'), $
;            align=1, vertical_align=1, color='blue', font_size=12)
;t032 = text(p031.pos[2]-0.02, p031.pos[3]-0.023, '$Downward_{med} = $'+string(down_med, f='(f5.2)'), $
;            align=1, vertical_align=1, color='red', font_size=12)
;p033 = plot([0, 0], p031.yr, '--2k', transp=50, over=p031)


;--------  (d) h3 Vs k3 from obs

xpar3 = (where(par_names eq 'mg_h3v'))[0]
ypar3 = (where(par_names eq 'mg_k3v'))[0]
real = where(finite(obj_arr.(xpar3)) and finite(obj_arr.(ypar3)))
xparam3 = (obj_arr.(xpar3))[real]
yparam3 = (obj_arr.(ypar3))[real]
h2d_h3k3 = hist_2d(xparam3, yparam3, $
                  bin1 = xbinsize, bin2 = ybinsize, $
                  min1 = dr[0], max1 = dr[1], $
                  min2 = dr[0], max2 = dr[1])
im04 = image_kh(h2d_h3k3/total(h2d_h3k3), xp, yp, /current, /dev, $
                pos=pos[*, 2], rgb_table=22, max=0.02, $
                title='(c) v$_{D, Mg II h3}$ Vs. v$_{D, Mg II k3}$ from observations', $ 
                xtitle=par_titles[xpar3], ytitle=par_titles[ypar3], $
                font_style=0, font_name='Helvetica', font_size=13)
cb04 = colorbar(target=im04, orientation=1, border=1, textpos=1, major=5, $
                pos=im04.pos[[2, 1, 2, 3]]+[0, 0, 0.01, 0])
p040 = plot([0, 0], im04.yr, ':', over=im04)
p041 = plot(im04.xr, [0, 0], ':', over=im04)
p042 = plot([-50, 50], [-50, 50], ':', over=im04)
p0421 = plot([-50, 50], [50, -50], ':', over=im04)

quad11 = where(yparam3 gt 0. and yparam3 lt xparam3, n11)
quad12 = where(xparam3 gt 0. and yparam3 gt xparam3, n12)
quad21 = where(xparam3 lt 0. and yparam3 gt -xparam3, n21)
quad22 = where(yparam3 gt 0. and yparam3 lt -xparam3, n22)
quad31 = where(yparam3 lt 0. and yparam3 gt xparam3, n31)
quad32 = where(xparam3 lt 0. and yparam3 lt xparam3, n32)
quad41 = where(xparam3 gt 0. and yparam3 lt -xparam3, n41)
quad42 = where(yparam3 lt 0. and yparam3 gt -xparam3, n42)
quad_all = float([n11, n12, n21, n22, n31, n32, n41, n42])

for ii=0, 7 do begin
  theta = !dtor*(360./8.*ii+22.5)
  t035 = text(dr[1]*0.6*cos(theta), dr[1]*0.7*sin(theta), target=im04, $
              '('+string(ii+1, f='(i0)')+') '+$
              string(quad_all[ii]/total(quad_all)*100., f='(f4.1)')+'%', /data, $
              align=0.5)
endfor

rere = where(xparam3 gt -15. and xparam3 lt 15. and yparam3 gt -15. and yparam3 lt 15.)
pcc = correlate(xparam3[rere], yparam3[rere])
t036 = text(im04.pos[0]+0.015, im04.pos[3]-0.03, 'r = '+string(pcc, f='(f4.2)'), $
            font_size=12)


; ------------ (e) v_{centr, h} Vs. v_{D, h3}
xpar = (where(par_names eq 'mg_h_centr'))[0]
ypar = (where(par_names eq 'mg_h3v'))[0]
real = where(finite(obj_arr.(xpar)) and finite(obj_arr.(ypar)))
xparam = (obj_arr.(xpar))[real]
yparam = (obj_arr.(ypar))[real]
h2d_hcen_h3 = hist_2d(xparam, yparam, $
                    bin1 = xbinsize, bin2 = ybinsize, $
                    min1 = dr[0], max1 = dr[1], $
                    min2 = dr[0], max2 = dr[1])
pcc = correlate(xparam, yparam)                    
im05 = image_kh(h2d_hcen_h3/total(h2d_hcen_h3), xp, yp, /current, /dev, $
                pos=pos[*, 3], rgb_table=22, max=0.01, $
                title='(d) v$_{centr, Mg II h}$ Vs. v$_{D, Mg II h3}$ from observations', $
                xtitle=par_titles[xpar], ytitle=par_titles[xpar3], $
                font_style=0, font_name='Helvetica', font_size=13)
cb05 = colorbar(target=im05, orientation=1, border=1, textpos=1, major=5, $
                pos=im05.pos[[2, 1, 2, 3]]+[0, 0, 0.01, 0])
p050 = plot([0, 0], im05.yr, ':', over=im05)
p051 = plot(im05.xr, [0, 0], ':', over=im05)
p052 = plot([-50, 50], [-50, 50], ':', over=im05)
p0521 = plot([-50, 50], [50, -50], ':', over=im05)

rere = where(xparam gt -15. and xparam lt 15. and yparam gt -15. and yparam lt 15.)
pcc = correlate(xparam[rere], yparam[rere])
t056 = text(im05.pos[0]+0.015, im05.pos[3]-0.03, 'r = '+string(pcc, f='(f4.2)'), $
            font_size=12)



;;----------- (f) Difference between v_{centr, Mg II h} and v_{D, Mg II h3}
;
;hist_diff = histogram(yparam-xparam, loc=xbin, binsize=1.)
;p061 = plot(xbin, hist_diff, 'k2', /hist, xr=[-15, 15], /current, /dev, $
;            pos = pos[*, 5], title='(f) Difference between two Mg II h velocities', $
;            xtitle='$|v_{D, Mg II h3}| - |v_{cent, Mg II h}| (km s^{-1})$', ytitle='# of pixels', $
;            font_size=13)
;med = median(yparam-xparam)
;p061.yr = p061.yr*1.2
;p0611 = plot(med*[1, 1], p061.yr, '--2k', over=p061)
;t061 = text(p061.pos[2]-0.02, p061.pos[3]-0.01, '$median = $'+string(med, f='(f5.2)'), $
;            align=1, vertical_align=1, font_size=12)

;w1.save, save_dir+'Mg_II_anal.png', resol=300
w1.save, save_dir+'Mg_II_anal.pdf', page_size=w1.dimen/1d2, width=w1.dimen[0]/1d2, /bitmap

end