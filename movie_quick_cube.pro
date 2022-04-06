PRO movie_quick_cube, mash, outfile, x0,x1,y0,y1,d1,d2, comp=comp, onehour=onehour


  ;datain - AIA data from newar.pro structures
  ;structure - output struct from mossfinder
  ;outfile - string for filename
  ;filter - threshold level of variability counts

  sd = size(mash['zmatch'])

  a193 = mash['aia_193'].data
  a171 = mash['aia_171'].data
  a94 = mash['aia_94'].data

  zmatch = mash['zmatch no loop']
  zmatch_d = mash['zmatch']
  fe18 = mash['fe18']

  date94 = mash['aia_94'].index.date_obs
  date193 = mash['aia_193'].index.date_obs

  if x1 eq 0 then x1 = sd[1]-1
  if y1 eq 0 then y1 = sd[2]-1


   video_file = outfile
   video = idlffvideowrite(video_file)
   framerate = 25
   framedims = [x1*2.2,y1+100]
   if d1 gt 0 then framedims = [d1,d2]

   stream = video.addvideostream(framedims[0], framedims[1], framerate)

   ;loadct, 1
   set_plot, 'z', /copy
   device, set_resolution=framedims, set_pixel_depth=24, decomposed=0

   nframes = sd[3]
   if keyword_set(onehour) then nframes = 300


   !p.multi =[0,2,1,0,0]

   A = FINDGEN(9) * (!PI*2/8.)
   usersym,cos(a)/3.,sin(a)/3.,/fill


for i=0, nframes-1 do begin
    ;print,i
    aia_lct,/load,wave='193'
    gamma= 0.7

    plot_image, a193[x0:x1,y0:y1,i]^gamma,position=[0.05,0.05,0.5,0.95],min=100^gamma,max=8000^gamma

    loadct,0
    for j=0,sd[2]-1 do begin
      wrow = where(zmatch[*,j,i] eq 1,wnum)
      if wnum gt 0 then oplot,[[wrow]],[[intarr(wnum)+j]],psym=8,color=254
    endfor

  if keyword_set(comp) then begin
    loadct,1
    for j=0,sd[2]-1 do begin
      wrow = where(zmatch_d[*,j,i] eq 1,wnum)
      if wnum gt 0 then oplot,[[wrow]],[[intarr(wnum)+j]],psym=8,color=160
    endfor
    endif


    tm = min(abs(anytim(date94)-anytim(date193[i])),t94)

    aia_lct,/load,wave='94'
    gamma = 0.3
    plot_image, fe18[x0:x1,y0:y1,i]^gamma,position=[0.5,0.05,0.95,0.95],ytickformat="(A1)",min=1,max=300^gamma  ;110

    loadct,0
    for j=0,sd[2]-1 do begin
      wrow = where(zmatch[*,j,i] eq 1,wnum)
      if wnum gt 0 then oplot,[[wrow]],[[intarr(wnum)+j]],psym=8,color=254
    endfor

  if keyword_set(comp) then begin
    loadct,7
    for j=0,sd[2]-1 do begin
      wrow = where(zmatch_d[*,j,i] eq 1,wnum)
      if wnum gt 0 then oplot,[[wrow]],[[intarr(wnum)+j]],psym=8,color=80
    endfor
  endif


    xyouts,0.1,0.9,date193[i],/normal,charsize=1.2,color=254,charthick=2
    xyouts,0.1,0.8,trim(i),/normal,charsize=1.2,color=254,charthick=2

	  timestamp = video.put(stream, tvrd(true=1))

   endfor

!p.multi=0
video.cleanup
device, /close
set_plot,'x'

end
