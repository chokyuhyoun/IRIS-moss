;pro find_iris_spectra_moss, dir, eout, outdir=outdir

dir = '/Users/khcho/Desktop/moss_test2'


; DEFAULT SETTING
t0 = systime(/sec)
cd, current=current_dir
this_file = file_which('find_iris_moss_events3.pro')
if this_file eq '' then this_file = file_search('~/desktop', 'find_iris_moss_events3.pro', /fully)
cd, file_dirname(this_file)

if n_elements(outdir) eq 0 then outdir = dir
sg_files = file_search(dir, 'iris_l2_*raster*.fits')
sji_files = file_search(dir, 'iris_l2_*SJI*.fits')
aia_files = file_search(dir, 'aia_l2*.fits')
aia_dir = file_dirname(aia_files[0])

read_iris_l2, sg_files, sg_index, sg_data, /sil
dum = where(strmatch(sji_files, '*1400*'), n)
if n then begin
  read_iris_l2, sji_files[dum], sji_index, sji_data, /sil
  print, 'Use IRIS SJI_1400A for co-alignment'
endif else begin
  read_iris_l2, sji_files[0], sji_index, sji_data, /sil
  print, 'No IRIS SJI_1400A. Use IRIS '+sji_index.tdesc1+' for co-alignment'
endelse
dum = where(strmatch(aia_files, '*1600*'))
read_iris_l2, aia_files[dum], aia1600_index, aia1600_data, /sil
dum = where(strmatch(aia_files, '*193*'))
read_iris_l2, aia_files[dum], aia193_index, aia193_data, /sil

; MAKE TIME ARRAY
sg_time = anytim(sg_index.date_obs)
sji_time = anytim(sji_index.date_obs)
aia1600_time = anytim(aia1600_index.date_obs)
aia193_time = anytim(aia193_index.date_obs)
aia_sz = size(aia1600_data)
sji_sz = size(sji_data)

;del = !null
;for i=0, n_elements(sji_time)-1 do begin
;  dum = min(abs(sji_time[i] - aia1600_time), match)
;  aia_xpdata = (findgen(aia_sz[1])-(aia_sz[1]-1)*0.5)*aia1600_index[match].cdelt1 $ 
;                + aia1600_index[match].crval1
;  aia_ypdata = (findgen(aia_sz[2])-(aia_sz[2]-1)*0.5)*aia1600_index[match].cdelt2 $ 
;                + aia1600_index[match].crval2
;  sji_xpdata = (findgen(sji_sz[1])-(sji_sz[1]-1)*0.5)*sji_index[i].cdelt1 $
;                + sji_index[i].crval1
;  sji_ypdata = (findgen(sji_sz[2])-(sji_sz[2]-1)*0.5)*sji_index[i].cdelt2 $
;                + sji_index[i].crval2
;
;  sji_int_xpix = interpol(findgen(sji_sz[1]), sji_xpdata, aia_xpdata)
;  sji_int_ypix = interpol(findgen(sji_sz[2]), sji_ypdata, aia_ypdata)
;
;  sji_corr_aia = interpolate(sji_data[*, *, i], $
;                             sji_int_xpix, sji_int_ypix, /grid)
;  
;  delxy = alignoffset(sji_corr_aia, aia1600_data[*, *, match])
;  del = [[del], [delxy]] 
;endfor

if n_elements(sg_files) eq 0 then begin
  message, 'No IRIS data found. Please check the directory'
endif

; RUN MOSSVAR ON AIA CUBE TO DETECT EVENTS
mossvar, dir, '',mash,4,4,/read,/cube_data
if mash['status'] ne 'data ok' then begin
  message, 'Failed to find moss structure from AIA data.'
  stop
endif
mossvar,dir,'',mash,4,4,/moss,/cnet,/fexviii,/loop,/filter1700
if keyword_set(movie) then mossvar,dir,'',mash,4,4,/variability,/cleanflare,/movie
if ~keyword_set(movie) then mossvar,dir,'',mash,4,4,/variability,/cleanflare
zcube = mash['zmatch no loop']

; VARIABLE TO SAVE INFOMATION
sji_ind = !null ; [x, y, t] index for IRIS SJI cube
sji_xyt = !null  ; [x arcsec, y arcsec, sec from the beginning of SJI observation] 
aia_ind = !null ; [x, y, t] index for IRIS / AIA cube
aia_xyt = !null  ; [x arcsec, y arcsec, sec from the beginning of AIA observation]
diff = !null ;  [x diff in arcsec, y diff in arcsec, time diff in sec]
aia_shift = !null ; [x shift, y shift in arcsec]
sg_spec = !null ; [lmabda, moss index]

; TOLERANCE (moss from AIA 193A and SG, spatio-temporally) 
tol_x_arc = 0.5*aia193_index[0].cdelt1
tol_y_arc = 0.5*aia193_index[0].cdelt2
tol_t_sec = 0.5*aia193_index[0].cdelt3

;FIND IT!
moss_num = total(total(zcube, 1), 1)
for i=0, n_elements(moss_num)-1 do begin ; AIA time
  if moss_num[i] eq 0 then continue

; FIND NEAREST SJI AND SPECTROGRAPH
  dum = min(abs(aia193_time[i]-sg_time), match_sg)
  if dum gt tol_t_sec then continue
  
  stop


    
endfor ; AIA time










end


