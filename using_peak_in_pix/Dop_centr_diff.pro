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

ars = fix((obj_arr.ar_no)[uniq(obj_arr.ar_no)])
col = transpose((colortable(6, ncol=n_elements(ars)+1))[1:*, *])
transp = 30

sim_col0 = transpose(colortable(38, ncolor=10))
sim_col = fltarr(3, n_elements(model_order))
sim_col[*, 0:3] = sim_col0[*, 0:3]
sim_col[*, 4:7] = sim_col0[*, 0:3]
sim_col[*, 8:*] = sim_col0[*, 4:-2]

obj_name1 = ['mg_h3v', 'mg_k3v']
obj_name2 = ['mg_h_centr', 'mg_k_centr']

w01 = window(dim=[1d3, 1d3])
for i=0, n_elements(obj_name1)-1 do begin
;  i = 0
  xpar = (where(strmatch(par_names, obj_name1[i]), /null))[0]
  ypar = (where(strmatch(par_names, obj_name2[i]), /null))[0]
  xparam = (obj_arr.(xpar))
  yparam = (obj_arr.(ypar))
  
  nbin = 50
  xbinsize = (par_dr[1, xpar]-par_dr[0, xpar])/nbin
  ybinsize = (par_dr[1, ypar]-par_dr[0, ypar])/nbin
  h2d = hist_2d(xparam, yparam, bin1 = xbinsize, bin2 = ybinsize, $
                min1 = par_dr[0, xpar], max1 = par_dr[1, xpar], $
                min2 = par_dr[0, ypar], max2 = par_dr[1, ypar])
  
  p011 = image_kh(bytscl(h2d), $
                 par_dr[0, xpar]+findgen(nbin+1)*xbinsize+0.5*xbinsize, $
                 par_dr[0, ypar]+findgen(nbin+1)*ybinsize+0.5*ybinsize, $ 
                 /current, /dev, aspect=0, hi_res = 0, $
                 pos=[100, 70+500*i, 450, 70+350+500*i], $ 
                 xr=par_dr[*, xpar], yr=par_dr[*, ypar], $
                 rgb_table=22, xminor=1, yminor=1, xmajor=3, ymajor=3, $
                 xtitle=par_titles[xpar], ytitle=par_titles[ypar], $
                 font_style=0, font_name='Helvetica', font_size=15)  
  real = where(xparam gt par_dr[0, xpar] and xparam lt par_dr[1, xpar] and $
               yparam gt par_dr[0, ypar] and yparam lt par_dr[1, ypar])
  cc = correlate(xparam[real], yparam[real])
  t03n = text(p011.pos[0]+0.01, p011.pos[3]-0.01, 'cc = '+string(cc, f='(f5.2)'), $
              vertical_align=1, align=0)
    
  hist0 = 0
  p012 = plot(indgen(2), /nodata, /current, /dev, $
              pos=[600, 150+500*i, 900, 350+500*i], $
              xtitle=par_names[ypar]+' - '+par_names[xpar]+' (km s$^{-1}$)', ytitle='Number', font_size=13)
  diff = yparam - xparam   
  for l=0, n_elements(ars)-1 do begin
    selected = diff[where(obj_arr.ar_no eq ars[l])]
    hist = histogram(selected, nbins=50, loc=xbin, $
                     min=-15, max=15)
    p014 = barplot(xbin+(xbin[1]-xbin[0])*0.5, hist0+hist, bottom_value=hist0, over=p012, $
                   fill_color=col[*, l], linestyle=' ', xstyle=1, transp=transp)
    hist0 = hist0 + hist

    mean_pos = p012.convertcoord(median(selected), p012.yr[1], /data, /to_normal)
    p013 = plot(mean_pos[0]*[1, 1], mean_pos[1]+[0.005, 0.01], $
                over=p012, color=col[*, l], thick=2, transp=transp)
  endfor
  p012.xr = [-15, 15]
  if k0 ne 0 then begin
    t01 = text(p012.pos[2]-0.01, p012.pos[3]-0.01, /current, $
              'm = '+string(median(diff), f='(f0.2)')+'!c'+ $
              '$\sigma$ = '+string(stddev(diff, /nan), f='(f0.2)'), $
              vertical_align=1, align=1, font_size=12, transp=transp)
    p013 = plot(median(diff)*[1, 1], p012.yr, over=p012, $
                ':k2', transp=50)
  endif  
endfor
w01.save, save_dir+'Dop_centroid_comp.png'
end