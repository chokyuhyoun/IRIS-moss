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
  
;  par_names = ['ar_no', 'sol_x', 'sol_y', 'times', 'cur_fwhm', $ ; 4
;               'si_v', 'pre_si_v', 'si_nth', 'pre_si_nth', 'si_amp', 'pre_si_peak', $ ; 10
;               'si_i_tot', 'pre_si_i_tot', 'log_den', 'pre_log_den', $ ; 14
;               'mg_h3v', 'pre_mg_h3v', 'mg_k3v', 'pre_mg_k3v', $ ; 18
;               'mg_trip', 'pre_mg_trip', $ ; 20 
;               'fe18_avg'] ; 21

;stop

w1 = window(dim=[10d2, 4d2])
xmar = [100, 100]
ymar = [50, 80]
nx = 3.
ny = 1.
xs = 270.
ys = xs
xgap = 0
ygap = 0

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
;dr = [-15d, 15d]
xbinsize = (dr[1]-dr[0])/nbin
ybinsize = (dr[1]-dr[0])/nbin
xp = dr[0]+findgen(nbin+1)*xbinsize
yp = dr[0]+findgen(nbin+1)*ybinsize


;-------- (a) Si IV I_tot Vs EW Mg II trip 
xpar = (where(par_names eq 'si_i_tot'))[0]
ypar = (where(par_names eq 'mg_trip'))[0]
real = where(finite(obj_arr.(xpar)) and finite(obj_arr.(ypar)))
xparam = (obj_arr.(xpar))[real]
yparam = (obj_arr.(ypar))[real]

xbinsize = (par_dr[1, xpar]-par_dr[0, xpar])/nbin
ybinsize = (par_dr[1, ypar]-par_dr[0, ypar])/nbin
xp = par_dr[0, xpar]+findgen(nbin+1)*xbinsize
yp = par_dr[0, ypar]+findgen(nbin+1)*ybinsize

centr_hk = hist_2d(xparam, yparam, $
                   bin1 = xbinsize, bin2 = ybinsize, $
                   min1 = par_dr[0, xpar], max1 = par_dr[1, xpar], $
                   min2 = par_dr[0, ypar], max2 = par_dr[1, ypar])
im01 = image_kh(centr_hk/total(centr_hk), xp, yp, /current, /dev, $
                pos=pos[*, 0], rgb_table=22, $
                title='(a) I$_{tot, Si IV}$ Vs. EW$_{Mg II trip}$', $
                xtitle=par_titles[xpar], ytitle=par_titles[ypar], max=0.005, $
                font_style=0, font_name='Helvetica', font_size=13, aspect_ratio=0, $
                xtickv=findgen(4)*1d3, xticklen=0.04, yticklen=0.04, yminor=3)
;cb01 = colorbar(target=im01, orientation=1, border=1, textpos=1, major=5, $
;                pos=im01.pos[[2, 1, 2, 3]]+[0, 0, 0.01, 0])
rere = where(xparam gt par_dr[0, xpar] and xparam lt par_dr[1, xpar] and $
             yparam gt par_dr[0, ypar] and yparam lt par_dr[1, ypar])
pcc = correlate(xparam[rere], yparam[rere])
t016 = text(pos[0, 0]+20, pos[3, 0]-20, 'r = '+string(pcc, f='(f4.2)'), /dev, $
            font_size=12, vertical_align=1)


;-------- (b) Si IV I_tot Vs Si IV nth
xpar = (where(par_names eq 'si_amp'))[0]
ypar = (where(par_names eq 'mg_trip'))[0]
real = where(finite(obj_arr.(xpar)) and finite(obj_arr.(ypar)))
xparam = (obj_arr.(xpar))[real]
yparam = (obj_arr.(ypar))[real]

xbinsize = (par_dr[1, xpar]-par_dr[0, xpar])/nbin
ybinsize = (par_dr[1, ypar]-par_dr[0, ypar])/nbin
xp = par_dr[0, xpar]+findgen(nbin+1)*xbinsize
yp = par_dr[0, ypar]+findgen(nbin+1)*ybinsize

centr_hk = hist_2d(xparam, yparam, $
  bin1 = xbinsize, bin2 = ybinsize, $
  min1 = par_dr[0, xpar], max1 = par_dr[1, xpar], $
  min2 = par_dr[0, ypar], max2 = par_dr[1, ypar])
im02 = image_kh(centr_hk/total(centr_hk), xp, yp, /current, /dev, $
  pos=pos[*, 1], rgb_table=22, $
  title='(b) I$_{amp, Si IV}$ Vs. EW$_{Mg II trip}$', $
  xtitle=par_titles[xpar], ytitle=par_titles[ypar], max=im01.max, $
  font_style=0, font_name='Helvetica', font_size=13, aspect_ratio=0, yshowtext=0, $
  xtickv=findgen(5)*20, xticklen=0.04, yticklen=0.04, yminor=3)
;cb01 = colorbar(target=im01, orientation=1, border=1, textpos=1, major=5, $
;  pos=im01.pos[[2, 1, 2, 3]]+[0, 0, 0.01, 0])
rere = where(xparam gt par_dr[0, xpar] and xparam lt par_dr[1, xpar] and $
  yparam gt par_dr[0, ypar] and yparam lt par_dr[1, ypar])
pcc = correlate(xparam[rere], yparam[rere])
t016 = text(pos[0, 1]+20, pos[3, 1]-20, 'r = '+string(pcc, f='(f4.2)'), /dev, $
            font_size=12, vertical_align=1)


;-------- (c) Si IV I_amp Vs Si IV nth
xpar = (where(par_names eq 'si_v'))[0]
ypar = (where(par_names eq 'mg_trip'))[0]
real = where(finite(obj_arr.(xpar)) and finite(obj_arr.(ypar)))
xparam = (obj_arr.(xpar))[real]
yparam = (obj_arr.(ypar))[real]

xbinsize = (par_dr[1, xpar]-par_dr[0, xpar])/nbin
ybinsize = (par_dr[1, ypar]-par_dr[0, ypar])/nbin
xp = par_dr[0, xpar]+findgen(nbin+1)*xbinsize
yp = par_dr[0, ypar]+findgen(nbin+1)*ybinsize

centr_hk = hist_2d(xparam, yparam, $
  bin1 = xbinsize, bin2 = ybinsize, $
  min1 = par_dr[0, xpar], max1 = par_dr[1, xpar], $
  min2 = par_dr[0, ypar], max2 = par_dr[1, ypar])
im03 = image_kh(centr_hk/total(centr_hk), xp, yp, /current, /dev, $
  pos=pos[*, 2], rgb_table=22, $
  title='(c) V$_{D, Si IV}$ Vs. EW$_{Mg II trip}$', $
  xtitle=par_titles[xpar], ytitle=par_titles[ypar], max=im01.max, $
  font_style=0, font_name='Helvetica', font_size=13, aspect_ratio=0, yshowtext=0, $
  xticklen=0.04, yticklen=0.04, yminor=3)
cb01 = colorbar(target=im03, orientation=1, border=1, textpos=1, major=5, $
  pos=im03.pos[[2, 1, 2, 3]]+[0, 0, 0.01, 0])
rere = where(xparam gt par_dr[0, xpar] and xparam lt par_dr[1, xpar] and $
  yparam gt par_dr[0, ypar] and yparam lt par_dr[1, ypar])
pcc = correlate(xparam[rere], yparam[rere])
t016 = text(pos[0, 2]+20, pos[3, 2]-20, 'r = '+string(pcc, f='(f4.2)'), /dev, $
            font_size=12, vertical_align=1)  


;w1.save, save_dir+'Mg_trip_anal.pdf', page_size=w1.dimen/1d2, width=w1.dimen[0]/1d2, /bitmap

end