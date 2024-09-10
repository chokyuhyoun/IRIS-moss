PRO dem4moss_cube, mash, bin, outdir, makedem=makedem, save=save, noreturn=noreturn


;need to compile these to calculate dem
;ar_dem_timestamp.pro
;dem_blos_imgplot.pro
;imgplot.pro
;Feb 20th - Added switch for rough coalignment of the wavelength channels

;makedem - run the dem inversion
;save - save the output strucutre, warning, it may be huge
;noreturn - skip returning the dem structure to the main program hash structure



ss = size(mash['aia_193'].data)

tarray = mash['temp_array']
nt = n_elements(tarray)

tlen = ss[3]

demchannels = ['94','131','171','193','211','335']

n_demchannels = n_elements(demchannels)
dkeys = mash.Keys()
demkeys = dkeys[1:n_demchannels]

tclosest = mash['tclosest']

a94 = mash['aia_94'].data
a131 = mash['aia_131'].data
a171 = mash['aia_171'].data
a193 = mash['aia_193'].data
a211 = mash['aia_211'].data
a335 = mash['aia_335'].data

dsize = fltarr(n_demchannels)

for c=0,n_demchannels-1 do begin
	dsize[c] = n_elements(mash[demkeys[c]].index.date_obs)
endfor

print,dsize

nbins = (tlen/bin)+1

o94 = fltarr(ss[1],ss[2],nbins)
o131 = fltarr(ss[1],ss[2],nbins)
o171 = fltarr(ss[1],ss[2],nbins)
o193 = fltarr(ss[1],ss[2],nbins)
o211 = fltarr(ss[1],ss[2],nbins)
o335 = fltarr(ss[1],ss[2],nbins)

dlow = fltarr(ss[1],ss[2],nbins)
dmed = fltarr(ss[1],ss[2],nbins)
dhigh = fltarr(ss[1],ss[2],nbins)

lowpass = fltarr(ss[1],ss[2],nbins)
medpass = fltarr(ss[1],ss[2],nbins)
highpass = fltarr(ss[1],ss[2],nbins)



for t = 0,tlen-1,bin do begin

	print,t

	;nearest times to 193
	;best to use correlated data using the wavecorr switch in mossvar_fromcube

	t94 = tclosest[t,0]
	t131 = tclosest[t,1]
	t171 = tclosest[t,2]
	t193 = t
	t211 = tclosest[t,3]
	t335 = tclosest[t,4]

	img = fltarr(ss[1],ss[2],6)

	tpos = [t94,t131,t171,t193,t211,t335]
	dmax = dsize - 1

	bshoot =  dmax-tpos

	ashoot = bshoot+1


	;part to bin data in time
	if bin eq 1 then begin
		img[*,*,0] = float(a94[*,*,t94])
		img[*,*,1] = float(a131[*,*,t131])
		img[*,*,2] = float(a171[*,*,t171])
		img[*,*,3] = float(a193[*,*,t193])
		img[*,*,4] = float(a211[*,*,t211])
		img[*,*,5] = float(a335[*,*,t335])
	endif

	if bin gt 1 then begin
		print,'binned dem'
		if bshoot[0] gt bin then img[*,*,0] = total(a94[*,*,t94:t94+bin-1],3)/bin else img[*,*,0] = total(reform(a94[*,*,t94:t94+bshoot[0]],ss[1],ss[2],ashoot[0]) ,3)/ashoot[0]
		if bshoot[1] gt bin then img[*,*,1] = total(a131[*,*,t131:t131+bin-1],3)/bin else img[*,*,1] = total(reform(a131[*,*,t131:t131+bshoot[1]],ss[1],ss[2],ashoot[1]) ,3)/ashoot[1]
		if bshoot[2] gt bin then img[*,*,2] = total(a171[*,*,t171:t171+bin-1],3)/bin else img[*,*,2] = total(reform(a171[*,*,t171:t171+bshoot[2]],ss[1],ss[2],ashoot[2]) ,3)/ashoot[2]
		if bshoot[3] gt bin then img[*,*,3] = total(a193[*,*,t193:t193+bin-1],3)/bin else img[*,*,3] = total(reform(a193[*,*,t193:t193+bshoot[3]],ss[1],ss[2],ashoot[3]) ,3)/ashoot[3]
		if bshoot[4] gt bin then img[*,*,4] = total(a211[*,*,t211:t211+bin-1],3)/bin else img[*,*,4] = total(reform(a211[*,*,t211:t211+bshoot[4]],ss[1],ss[2],ashoot[4]) ,3)/ashoot[4]
		if bshoot[5] gt bin then img[*,*,5] = total(a335[*,*,t335:t335+bin-1],3)/bin else img[*,*,5] = total(reform(a335[*,*,t335:t335+bshoot[5]],ss[1],ss[2],ashoot[5]) ,3)/ashoot[5]
	endif


	;binned data
	o94[*,*,t/bin] =  img[*,*,0]
	o131[*,*,t/bin] = img[*,*,1]
	o171[*,*,t/bin] = img[*,*,2]
	o193[*,*,t/bin] = img[*,*,3]
	o211[*,*,t/bin] = img[*,*,4]
	o335[*,*,t/bin] = img[*,*,5]


	demtime = mash['aia_193'].index[t].date_obs

	if keyword_set(makedem) then begin
		temp = 0

		aia_sparse_em_solve, img, oem=oem, status=status, tolfac=1.4, tolfunc=tolfunc,symmbuff=1.0

		;plot_image,oem[*,*,10]

		singledem = create_struct('dem',oem,'time',demtime,'status',status)


		if t eq 0 then demstr = [singledem]
		if t gt 1 then demstr = [demstr,singledem]
	endif

endfor

	if keyword_set(save) then save,demstr,filename=outdir+'/emcube.sav'
	if ~keyword_set(noreturn) then mash['dem']=demstr

end