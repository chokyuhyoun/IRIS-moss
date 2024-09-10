dir = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, dir

gather_files = file_search(dir+'moss_param_event_total.sav')
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


overall = 0
overall_fe18 = 0
peak_prev = 0
after_prev = 0
all_scatter = 0
all_density = 0
diff_density = 0
mg_ii_kh = 0
mg_si = 0
abs_siv_nth = 0
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
  
  
w1_sz = [6d2, 3d2]
w1 = window(dim=[w1_sz])
k0 = 7
k = params[k0]
hist0 = 0
p01 = objarr(n_Elements(params))
p01[k0] = plot(indgen(2), /nodata, /current, /dev, $
  pos=[60, 60, 260, 260], $
  xtitle=par_titles[k], ytitle='Number', font_size=10)
for l=0, n_elements(ars)-1 do begin
  selected = (obj_arr.(k))[where(obj_arr.(0) eq ars[l])]
  ;        selected = (obj_arr.(k))[where(obj_arr.ar_no eq ars[l] and (~finite(obj_arr.mg_h3v) or ~finite(obj_arr.mg_k3v)))]
  ;        stop
  hist = histogram(selected, nbins=50, loc=xbin, $
    min=par_dr[0, k], max=par_dr[1, k])
  p011 = barplot(xbin+(xbin[1]-xbin[0])*0.5, hist0+hist, bottom_value=hist0, over=p01[k0], $
    fill_color=col[*, l], linestyle=' ', xstyle=1, transp=transp)
  hist0 = hist0 + hist
;  if k0 ne 0 then begin
;    mean_pos = p01[k0].convertcoord(median(selected), p01[k0].yr[1], /data, /to_normal)
;    p012 = plot(mean_pos[0]*[1, 1], mean_pos[1]+[0.005, 0.01], $
;      over=p00, color=col[*, l], thick=2, transp=transp)
;  endif
endfor
p01[k0].xr = par_dr[*, k]
p01[k0].yr = p01[k0].yr * 1.1
if k0 ne 0 then begin
  t01 = text(p01[k0].pos[2]-0.01, p01[k0].pos[3] - 0.2 - (k0 eq 3 ? 0 : 0.04), /current, $
    'm = '+string(median(obj_arr.(k)), f='(f0.2)')+'!c'+ $
    '$\sigma$ = '+string(stddev(obj_arr.(k), /nan), f='(f0.2)'), $
    vertical_align=1, align=1, font_size=10, transp=0)
  p013 = plot(median(obj_arr.(k))*[1, 1], p01[k0].yr, over=p01[k0], $
    ':k2', transp=50)
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
    p018 = plot(xval, [p01[k0].yr[1]*(0.82-plus*0.1+0.01*ii)], over=p01[k0], $
      symbol=plus ? 1 : 24, $
      sym_size=1.5, sym_filled=0, sym_color=sim_col[*, ii], $
      transp=transp, sym_thick=2)
  endif else begin
    tem_xpos = [(p01[k0].xr[1]-p01[k0].xr[0])*0.92+p01[k0].xr[0]]
    tem_ypos = [p01[k0].yr[1]*(0.8-offset*0.03)]
    p018 = plot(tem_xpos, tem_ypos, over=p01[k0], $
      symbol=plus ? 1 : 24, $
      sym_size=5, sym_filled=0, sym_color=sim_col[*, ii], $
      transp=transp, sym_thick=2)
    p019 = text(tem_xpos, tem_ypos, '$\rightarrow$', target=p01[k0], /data, $
      color=sim_col[*, ii], align=0, vertical_align=0.5, transp=transp, $
      font_size=10)
    offset += 1
    ;          stop
  endelse
endfor


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
ydr = [0., 100.]
xbinsize = (xdr[1]-xdr[0])/nbin
ybinsize = (ydr[1]-ydr[0])/nbin
h2d = hist_2d(xparam, yparam, bin1=xbinsize, bin2=ybinsize, $
  min1 = xdr[0], max1 = xdr[1], $
  min2 = ydr[0], max2 = ydr[1])
im07 = image_kh(h2d/total(h2d), $
  xdr[0]+findgen((size(h2d))[1])*xbinsize+0.5*xbinsize, $
  ydr[0]+findgen((size(h2d))[2])*ybinsize+0.5*ybinsize, $
  /current, aspect=0, hi_res = 1, $
  pos=[350, 60, 550, 260], /dev, $
  xr=xdr, yr=ydr, rgb_table=22, $
  title='JPDF of $\mu$ and $v_{nth, Si IV}$', $
  xtitle='$\mu$ = cos $\theta$', ytitle=par_titles[ind1], $
  font_style=0, font_name='Helvetica', font_size=10)
t7 = text(im07.pos[0]+0.02, im07.pos[3]-0.05, $
  'cc = '+string(correlate(xparam, yparam), f='(f4.2)'), $
  font_size=10, vertical_align=1, align=0)
cb07 = colorbar(target=im07, orientation=1, border=1, textpos=1, $
  pos=im07.pos[[2, 1, 2, 3]]+[0, 0, 0.01, 0])


w1.save, save_dir+'senior_review.png', resol=300
  
end    