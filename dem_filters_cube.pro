PRO dem_filters_cube, mash, bin, rfilter,lowpass,medpass,highpass,dlow,dmed,dhigh,hbinning=hbinning

;dir - directory of DEMs
;nstart nend - start and end indices
;cutx - x range for background region
;cuty - y range for background region

;hbinning - to bin the high temp filter more
;bin = number of binned steps


ss = size(mash['dem'].dem)
ndems = ss[4]


dlow = fltarr(ss[1],ss[2],ndems)
dmed = fltarr(ss[1],ss[2],ndems)
dhigh = fltarr(ss[1],ss[2],ndems)

lowpass = fltarr(ss[1],ss[2],ndems)
medpass = fltarr(ss[1],ss[2],ndems)
highpass = fltarr(ss[1],ss[2],ndems)


tarray = mash['temp_array']
demstr = mash['dem']


for i=0,ndems-1 do begin

	;choose 3 temperature ranges

	dlow[*,*,i] = total(demstr[i].dem[*,*,1:3],3)		;5.6 5.8
	dmed[*,*,i] = total(demstr[i].dem[*,*,6:9],3)		;6.1 6.4
	dhigh[*,*,i] = total(demstr[i].dem[*,*,12:15],3)	;6.7 7.0

	lowpass[*,*,i] = dlow[*,*,i] ge 20. ;25
	medpass[*,*,i] = dmed[*,*,i] ge 25
	highpass[*,*,i]  = dhigh[*,*,i] ge 250

endfor


;=========================================
if keyword_set(hbinning) then begin
;bin high temperature filter by 3 time steps

for i=2,ndems-1 do begin
	highpass[*,*,i]=(total(highpass[*,*,i-2:i],3) ge 1)
endfor

endif

;========== use a filter where high and low temperatures are present together =============

fs = size(lowpass)

k = replicate(1,3,3)
ks = [[0,1,0],[1,1,1],[0,1,0]]
rfilter = fltarr(ss[1],ss[2],fs[3])

;for each exposure ===============


for t = 0,fs[3]-1 do begin
	rmask = highpass[*,*,t] or lowpass[*,*,t]

	lmask = label_region(rmask)

;	;remove some of the very small regions

	for i=0,max(lmask)-1 do begin
		ww = where(lmask eq i, sregion)
		if sregion lt 2 then rmask[ww] = 0
	endfor

	;try scricter filters and more blowing up of the regions
	rfilter[*,*,t] = dilate(dilate(dilate(dilate(rmask,ks),ks),ks),ks)
endfor

mash['dem filter'] = rfilter  ; use mash['dem mask'] in mossvar.pro
;stop
end

