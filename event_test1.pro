dir = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, dir

gather_files = file_search(dir, '/*/*gather*.sav')
w1 = window(dim=[10d2, 6d2])
p01 = plot(indgen(2), /nodata, /current, $
            pos=[0.15, 0.15, 0.5, 0.7], xtitle='$\rho_e$', ytitle='Si IV $v_D$', $
            xr=[0, 1d3], yr=[-50, 50], $
            font_size=12, font_name='malgun gothic', font_style=1)

col39 = colortable(39)
for i=0, n_elements(gather_files)-1 do begin
  col = reform(col39[254.*i/(n_elements(gether_files)-1), *])
  restore, gather_files[i]
;  par_names = ['sol_x', 'sol_y', 'times', 'si_v', 'si_nth', 'si_peak', 'mg_h3v', 'mg_k3v', 'e_den']
  p011 = plot(pars_ev_mean[*, -1], pars_ev_mean[*, 3], '.', over=p01, $
              color=col)
endfor



end