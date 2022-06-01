dir = '/Users/khcho/Desktop/IRIS-moss-main'
cd, dir
sav_files = file_search(dir, '*moss_event*.sav')
sp_num = 0
pix_time = list()
blue = 2797.5
red = 2802.5
MgII_trip = [2798.75, 2798.82]

div = !null
div_date = !null
fwhm = !null

for i=0, n_elements(sav_files)-1 do begin
  restore, sav_files[i], /relax
  sp_num += n_elements(eout.sji_fwhm)
  sg_indy = eout.sg_ind[0, *]
  uniq_pix = sg_indy[uniq(sg_indy, sort(sg_indy))]
  wave = eout.wave_list[8]
  data = eout.data_list[8]
  dum = min(abs(wave-blue), blue_ind)
  dum = min(abs(wave-red), red_ind)
  div = [div, n_elements(pix_time)]
  div_date = [div_date, strmid(eout.sg_files[0], 36, 15)]
  for j=0, n_elements(uniq_pix)-1 do begin
    ind = where(sg_indy eq uniq_pix[j], count)
    in_pix = intarr(count)+1
    for k=0, count-1 do begin
      wave_part = wave[blue_ind:red_ind] 
      data_part = data[blue_ind:red_ind, ind[k]]
      coeff = poly_fit(wave_part, data_part, 2, $
                       measure_error=sqrt(data_part), sigma=sigma, yerror=yerror)
      MgII_data = interpol(data_part, wave_part, MgII_trip)
      MgII_base = poly(MgII_trip, coeff)
      if total(MgII_data - MgII_base) gt yerror then begin
        in_pix[k] = 2
;        w1 = window()
;        p01 = plot(wave, data[*, ind[k]], /current)
;        p01.yr = [-10, p01.yr[1]]
;        p011 = plot(wave_part, poly(wave_part, coeff), '--', over=p01)
;        for dum=0, 1 do p012 = plot(MgII_trip[dum]*[1, 1], p01.yr, '--', over=p01, color='gray')
;        stop 
      endif
      if k eq 0 then fwhm = [fwhm, eout.sji_fwhm[ind[k]]]
    endfor
    pix_time.add, in_pix
  endfor
endfor

length = 0
for i=0, n_elements(pix_time)-1 do length = length > n_elements(pix_time[i])
table = fltarr(n_elements(pix_time), length)
table[*] = !values.f_nan
for i=0, n_elements(pix_time)-1 do begin
  table[i, 0:n_elements(pix_time[i])-1] = pix_time[i]
endfor
tot_num = n_elements(where(finite(table)))
mg_tri_emiss_num = n_elements(where(table eq 2))
w2 = window(dim=[8d2, 6d2])
im02 = image(table, axis=2, pos=[0.1, 0.4, 0.98, 0.8], /current, $
             rgb_table=54, min=0, max=2, yr=[0, 20], ytickinterval=5, $
             xtitle='Position #', ytitle='Time (pix)', $
             font_size=13, font_name='malgon gothic', font_style=1)
t021 = text(0.05, 0.9, /current, $
            'Total # of Moss Spectra : '+string(tot_num, f='(i4)'), $
            color=reform((colortable(im02.rgb_table))[128, *]) ,$
            font_size=13, font_name='malgon gothic', font_style=1)
t022 = text(0.05, 0.85, /current, $
            '# of Mg II Triplet emiss. : '+string(mg_tri_emiss_num, f='(i3)') $
            + ' ('+string(mg_tri_emiss_num*1d2/tot_num, f='(f5.2)')+' %)', $
            color=reform((colortable(im02.rgb_table))[254, *]) ,$
            font_size=13, font_name='malgon gothic', font_style=1)
for i=1, im02.yr[1]-1 do p02 = plot(im02.xr, i*[1, 1], '-w', over=im02) 
for i=1, im02.xr[1]-1 do p02 = plot(i*[1, 1], im02.yr, '-w', over=im02)

for i=0, n_elements(div)-1 do begin
  t023 = text(div[i]+0.5, -3, div_date[i], target=im02, clip=0, /data, $
              align=1, vertical_align=0.5, orientation=90)
endfor

w2.save, 'moss_stat.png', resol=200
end