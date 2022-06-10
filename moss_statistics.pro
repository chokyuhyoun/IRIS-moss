dir = '/Users/khcho/Desktop/IRIS-moss-main'
cd, dir
sav_files = file_search(dir, '*_*/*event*.sav', /fully)
sp_num = 0
pix_time = list()
blue = 2797.5
red = 2802.5
MgII_trip = [2798.75, 2798.82]


div = !null
div_date = !null
fwhm = !null

sg_moss_num = !null
sg_moss_time = !null
aia_moss_num = list(len=n_elements(sav_files))
aia_moss_time = list(len=n_elements(sav_files))
aia_waves = ['*94*', '*171*', '*211*']
fe_xviii_curve = list(len=n_elements(sav_files))
aia_pixel_area = (0.6*725)^2.  ; km^2

if 0 then begin
  for i=0, n_elements(sav_files)-1 do begin
    restore, sav_files[i], /relax
    if strmatch(sav_files[i], '*moss_event*') then begin
      sp_num += n_elements(eout.curve_fwhm)
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
          if k eq 0 then fwhm = [fwhm, eout.curve_fwhm[ind[k]]]
        endfor
        pix_time.add, in_pix
      endfor
      sg_time0 = reform(eout.sg_phy[2, *])
      sg_time1 = sg_time0[uniq(sg_time0, sort(sg_time0))]
      sg_num = intarr(n_elements(sg_time1))
      for ii=0, n_elements(sg_time1)-1 do sg_num[ii] = total(sg_time0 eq sg_time1[ii])
      sg_moss_time = [sg_moss_time, sg_time1]
      sg_moss_num = [sg_moss_num, sg_num]
    endif
    
    aia_data = list(len=n_elements(aia_waves))
    for j=0, n_elements(aia_waves)-1 do begin
      aia_file = eout.aia_files[where(strmatch(eout.aia_files, aia_waves[j]))]
      read_iris_l2, aia_file, index, data, /sil
      sz = size(data)
      aia_data[j] = data / rebin(reform(index.exptime, 1, 1, sz[3]), sz[1], sz[2], sz[3])
    endfor
    fe_xviii = aia_data[0] - aia_data[1]/450. - aia_data[2]/120.
    curve0 = !null
    num0 = !null
    for k=0, sz[3]-1 do begin
      get_xp_yp, index[k], xp, yp
      xxp = rebin(xp, sz[1], sz[2])
      yyp = rebin(transpose(yp), sz[1], sz[2])
      dist = sqrt(xxp^2. + yyp^2.)*725d
      cos_theta = dist/7d5
      curve0 = [curve0, total(fe_xviii[*, *, k]/cos_theta)/aia_pixel_area/sz[1]/sz[2]]
      num0 = [num0, total(eout.zcube[*, *, k]/cos_theta)/aia_pixel_area/sz[1]/sz[2]]
    endfor
    fe_xviii_curve[i] = curve0
    aia_moss_num[i] = num0
    aia_moss_time[i] = anytim(index.date_obs)
  ;  stop
  endfor
  save, fe_xviii_curve, aia_moss_num, aia_moss_time, sg_moss_num, sg_moss_time, $
        pix_time, div, div_date, filename='moss_stat.sav'
endif else restore, 'moss_stat.sav', /relax

length = 0
for i=0, n_elements(pix_time)-1 do length = length > n_elements(pix_time[i])
table = fltarr(n_elements(pix_time), length)
table[*] = !values.f_nan
for i=0, n_elements(pix_time)-1 do begin
  table[i, 0:n_elements(pix_time[i])-1] = pix_time[i]
endfor
tot_num = n_elements(where(finite(table)))
mg_tri_emiss_num = n_elements(where(table eq 2))

if 0 then begin
  w2 = window(dim=[8d2, 6d2])
  im02 = image(table, axis=2, pos=[0.1, 0.3, 0.98, 0.8], /current, $
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
endif 


aia_time = !null
fe_curve = !null
moss_num = !null
for i=0, n_elements(aia_moss_time)-1 do begin
;  p31 = plot(aia_moss_time[i], fe_xviii_curve[i], '.k', over=p30)
;  p32 = plot(aia_moss_time[i], aia_moss_num[i], '.r', over=p30)
  aia_time = [aia_time, aia_moss_time[i]]
  fe_curve = [fe_curve, fe_xviii_curve[i]]
  moss_num = [moss_num, aia_moss_num[i]]
endfor
w3 = window(dim=[10d2, 6d2])
p30 = plot(aia_time, fe_curve*1d6, '.', /current, $
            pos=[0.1, 0.15, 0.75, 0.9], xstyle=1, $
            yr=[-2, 8], xminor=3, $
            xtitle='Time (2016)', ytitle='Fe XVIII Light curve (DN s$^{-1}$ km$^2$)', $
            font_size=13, font_name='malgun gothic', font_style=1)
nxticks = fix((p30.xr[1]-p30.xr[0])/86400d0)
xtickv = (dindgen(nxticks) + ceil(p30.xr[0]/86400d0))*86400d0
xtickstr = strmid(anytim(xtickv, /ccsds), 5, 5)
p30.xtickval = xtickv
p30.xtickname = xtickstr
p30.axes[3].hide = 1            
p31 = plot(aia_time, moss_num*1d6, '.r', /current, axis_style=0, $
            pos=p30.pos, xr=p30.xr, yr=[0, 0.003], $
            font_size=13, font_name='malgun gothic', font_style=1)
ax31 = axis('Y', location='right', target=p31, $
            title='# of Moss (s$^{-1}$ Mm$^2$)' , color=p31.color, $
            tickfont_size=13, tickfont_name='malgun gothic', tickfont_style=1)
p32 = plot(sg_moss_time, sg_moss_num, 'o', color='blue', /current, axis_style=0, $
           pos=p30.pos, xr=p30.xr, yr=[0, 7], sym_filled=1, sym_size=0.5 )
ax32 = axis('Y', location=p30.xr[1]+(p30.xr[1]-p30.xr[0])*0.25, target=p32, $
            title='# of Moss on IRIS Slit', textpos=1, tickdir=1, color=p32.color, $
            tickfont_size=13, tickfont_name='malgun gothic', tickfont_style=1)
w3.save, 'light_curve_moss.png', resol=200            
end