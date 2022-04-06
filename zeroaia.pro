function zeroaia,dcube,errors,sig,crossreturn=crossreturn

;use zerocrossings on AIA data

;returns zcube - binary mask where crossings happen

sc = size(dcube)

crossmap = fltarr(sc[1],sc[2])
zcube = fltarr(sc[1],sc[2],sc[3])

sdtime = fltarr(sc[1],sc[2])


thrs = sig

;as in paolas code
crosscube = fltarr(sc[1],sc[2],sc[3])
crosscube[where(dcube gt thrs*errors)]=1
crosscube(where(dcube lt -thrs*errors))=-1
crosscube[where(dcube lt thrs*errors and dcube gt -thrs*errors)]=0


for i=0,sc[1]-1 do begin
	for j=0,sc[2]-1 do begin
		dint = reform(crosscube[i,j,*])

			aa = zerocrossings(dint,n_crossings=rawcross)

			ab = where(aa gt 0,n_crossings)

		crossmap[i,j] = n_crossings
		
		;array of positions where there are crossings, binary map
		dint = dint*0.
		dint[aa] = 1.
		dint[0] = 0.  ;fix for dodgy 1 IDL puts at the start (why?) and cant do zero on first data point anyway
		zcube[i,j,*] = dint
	endfor
endfor


if keyword_set(crossreturn) then str = {crossmap:crossmap, crosscube:crosscube, zcube:zcube}

if ~keyword_set(crossreturn) then str = {crossmap:crossmap, zcube:zcube}

return,str

end