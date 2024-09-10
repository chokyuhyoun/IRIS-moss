path = '/Users/khcho/Desktop/IRIS-moss-main'
iris2_path = path+'/iris2'
cd, iris2_path

restore, 'collection.sav'
;save, temp, v_los, v_turb, n_e, $
;      u_temp, u_v_los, u_v_turb, u_n_e, chi2, ltau, $
;      pre_temp, pre_v_los, pre_v_turb, pre_n_e, $
;      pre_u_temp, pre_u_v_los, pre_u_v_turb, pre_u_n_e, pre_chi2, $
;      aft_temp, aft_v_los, aft_v_turb, aft_n_e, $
;      aft_u_temp, aft_u_v_los, aft_u_v_turb, aft_u_n_e, aft_chi2, $
;      pars_pix_peak, par_names, var_names, $
;      filename=iris2_path+'/collection.sav'


npix = n_elements(pars_pix_peak.ar_no)
del_t = aft_temp-temp
u_del_t = sqrt(u_temp^2. + pre_u_temp^2.)

;p01 = plot(indgen(2), /nodata)
;for i=0, npix-1 do p01n = plot(ltau, temp[*, i], over=p01, transp=99)
u_mask = abs(del_t) gt u_del_t
u_mask *= 1.
u_mask[where(u_mask eq 0)] = !values.f_nan

xmar = [100, 100]
ymar = [80, 100]
im_sz = [300, 300]
pl_sz = [im_sz[0], 100]
win_sz = [total(xmar)+im_sz[0]*3, total(ymar)+2*(im_sz[1]+pl_sz[1])]
w01 = window(dim=win_sz)

for i=0, 2 do begin
  target = (list(pre_temp, temp, aft_temp))[i]
  im01 = image_kh(target, ltau, findgen(npix), $
                  axis=2, /current, /dev, aspect_ratio=0, xstyle=1, $
                  pos=[xmar[0]+im_sz[0]*i, ymar[0], $
                       xmar[0]+im_sz[0]*(i+1), ymar[0]+im_sz[1]], $
                  rgb_table=33, min=2d3, max=2d4, xr=minmax(ltau), $
                  xtitle='Log $\tau$', ytitle='# of Pixel', $
                  font_size=13)
  if i ne 0 then im01.yshowtext=0
;  cb01 = colorbar(target=im01, pos=im01.pos[[2, 1, 2, 3]]+[0, 0, 0.01, 0], $
;                  title='$\Delta T (10^3 K)$', orient=1, textpos=1, /border)
  ;im011 = image_kh(u_mask, ltau, findgen(npix), over=im01, $
  ;                 aspect_ratio=0, rgb_table=0, transp=50)
  p01 = plot(ltau, mean(del_t, dim=2)*1d-3, '2', /current, /dev, $
              pos=[xmar[0]+im_sz[0]*i, ymar[0]+im_sz[1], $
                   xmar[0]+im_sz[0]*(i+1), ymar[0]+im_sz[1]+pl_sz[1]], $
              xr=minmax(ltau), $
              xtickformat='(a1)', ytitle='$\Delta T (10^3 K)$', font_size=13)
  p001 = plot(ltau, median(temp-pre_temp, dim=2)*1d-3, color='gray', over=p01)
  p011 = plot(p01.xr, [0, 0], ':', over=p01)
endfor
  cb01 = colorbar(target=im01, pos=im01.pos[[2, 1, 2, 3]]+[0, 0, 0.01, 0], $
                  title='$\Delta T (K)$', orient=1, textpos=1, /border)

; temp - pre_temp
;ltau_par =
;
;nbin = 50
;xbinsize = (par_xr[1, i]-par_xr[0, i])/nbin
;ybinsize = (par_xr[1, j]-par_xr[0, j])/nbin
;h2d = hist_2d(xparam, yparam, bin1 = xbinsize, bin2 = ybinsize, $
;  min1 = par_xr[0, i], max1 = par_xr[1, i], $
;  min2 = par_xr[0, j], max2 = par_xr[1, j])
;
;



end