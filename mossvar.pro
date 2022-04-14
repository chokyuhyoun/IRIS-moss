PRO mossvar,areg, tag, mash, dembinning, sigma, makecube=makecube, mossfilter=mossfilter, $
demfilter=demfilter, variability=variability, rundem=rundem, movie=movie, netfilter=netfilter, $
fexviii=fexviii,cleanflare=cleanflare,loopfilter=loopfilter, cnetfilter=cnetfilter, filter1700=filter1700,skip1700=skip1700, $
cube_data=cube_data,wavecorr=wavecorr, read=read, manual_data=manual_data, filebychannel=filebychannel, timecorr=timecorr,$
paper_data=paper_data,demrestore=demrestore,localcube=localcube,local_movie=local_movie, $
sdo_files=sdo_files

;tclosest : nearest time to AIA 193 A
;
;/timecorr : align AIA 193 A time series data
;/savecorr : align all AIA data based on 193A
;/netfilter : find network region (1600 A 80 DN)
;/cnetfilter : similar with /netfilter but use 1700 A
;/mossfilter : create moss area mask (193A > 1250 DN)
;/filter1700 : remove larger flare or filament brightening (1700A/1600A > 0.2)
;/fexviii : Fe XVIII image (94A - 211A/120. - 171A/450.)
;/loopfilter : find hot regions (94A > 5 DN)
;/rundem : DEM on each time step


;April 1st - updating correlation step to use tr_get_disp (should be more stable)
;March 29th - cleaned up version
;May 24th - version now with correlation on data that removes rotation, and performs correlation between 171 and 193
;May 24th - version to get data from SDO/IRIS cubes
;Aug 6th - minor cleaning up - creating new version with IRIS check

;areg - active region - folder name
;tag  - label for dataset - can be a time stamp or nothing (e.g tag='')

;channels required for proceedure
channels = ['94','131','171','193','211','335','1600','1700']
n_channels = n_elements(channels)

empty = 0.


;##############################################################
if keyword_set(read) then begin

;create ordered hash structure to put everything in - works like a python dictionary
;access with hash['key'] --- check keys with hash.Keys()

	mash = orderedhash('AR',areg)

if keyword_set(cube_data) then begin
	;restore data from premade AIA/IRIS cubes

	if keyword_set(localcube) then begin
		sdo_loc = areg+'/SDO'
		allfiles = find_files('*aia_l2*',areg+'/SDO')
	endif else begin
		allfiles = file_search(areg, '*aia_l2*')
	endelse

	;Create hash structure with data and index for each wavelength

	if n_elements(allfiles) ge 7 then begin

		for i=0,n_channels-1 do begin
;			sdofile = find_files('*aia_l2*'+channels[i]+'.fits',sdo_loc)
      sdofile = file_search(areg, '*aia_l2*'+channels[i]+'.fits')
			print,sdofile

			read_iris_l2,sdofile,index,data,/keep_null
			dstr = create_struct('index',index,'data',data)
			wavtag = 'aia_'+channels[i]
			mash[wavtag] = dstr
		endfor
		mash['status'] = 'data ok'
	endif else begin
		mash['status'] = 'no data'
	endelse
endif

;--------------------------

if keyword_set(manual_data) then begin
	;to restore data in save files created from ar_server_pro.pro

	for i=0,n_channels-1 do begin

		if ~keyword_set(filebychannel) then datafile = areg+'/cube_'+channels[i]+'_'+tag+'.sav'
		if keyword_set(filebychannel) then datafile = areg+'/'+channels[i]+'/cube_'+channels[i]+'_'+tag+'.sav'
		print,datafile
		restore,datafile
		dstr = create_struct('index',index,'data',data)
		wavtag = 'aia_'+channels[i]
		mash[wavtag] = dstr
	endfor
endif

if keyword_set(paper_data) then begin
	;to restore data prepeared by DG for 2019 paper - no one else needs to use this...

	for i=0,n_channels-1 do begin

		;a fix because 1700 is missing for these steps
		cstring = channels[i]
		if tag eq '66' or tag eq '67' or tag eq '68' and cstring eq '1700' then begin
			cstring = '1600'
			print,'skip 1700'
		endif

		datafile = areg+'/'+cstring+'/cube_'+cstring+'_'+tag+'.genx'
		print,datafile
		restgen,file=datafile,struct=gstr
		dstr = create_struct('index',gstr.index,'data',gstr.map.data)
		wavtag = 'aia_'+channels[i]
		mash[wavtag] = dstr
	endfor
endif

if mash['status'] ne 'no data' then begin

dkeys = mash.Keys()
datakeys = dkeys[1:n_channels]

print,datakeys

;do exposure time correction
for i=0,n_channels-1 do begin
	struct_working = mash[datakeys[i]]
	etime = struct_working.index.exptime
	for t=0,n_elements(etime)-1 do struct_working.data[*,*,t] = struct_working.data[*,*,t]/etime[t]
	mash[datakeys[i]] = struct_working
endfor

mash['DN/s'] = 'done'
print,'Converted to DN/s'

;a few empty variables for checks
mash['dem mask done'] = 0
mash['carbon mask done'] = 0
mash['wavecorr done'] = 0


;=======================

;useful for all

ss = size(mash['aia_193'].data)

atime94 = anytim(mash['aia_94'].index.date_obs)
atime131 = anytim(mash['aia_131'].index.date_obs)
atime171 = anytim(mash['aia_171'].index.date_obs)
atime193 = anytim(mash['aia_193'].index.date_obs)
atime211 = anytim(mash['aia_211'].index.date_obs)
atime335 = anytim(mash['aia_335'].index.date_obs)
atime1600 = anytim(mash['aia_1600'].index.date_obs)
atime1700 = anytim(mash['aia_1700'].index.date_obs)

dkeys = mash.Keys()
datakeys = dkeys[1:n_channels]

tlength = fltarr(n_channels)
for i=0,n_channels-1 do tlength[i] = n_elements(mash[datakeys[i]].index.date_obs)


tclosest = fltarr(ss[3],7)

;nearest times to 193
for t=0,ss[3]-1 do begin
	tm = min(abs(atime94-atime193[t]),t94)
		tclosest[t,0] = t94
	tm = min(abs(atime131-atime193[t]),t131)
		tclosest[t,1] = t131
	tm = min(abs(atime171-atime193[t]),t171)
		tclosest[t,2] = t171
	tm = min(abs(atime211-atime193[t]),t211)
		tclosest[t,3] = t211
	tm = min(abs(atime335-atime193[t]),t335)
		tclosest[t,4] = t335
	tm = min(abs(atime1600-atime193[t]),t1600)
		tclosest[t,5] = t1600
	tm = min(abs(atime1700-atime193[t]),t1700)
		tclosest[t,6] = t1700
endfor


mash['tclosest'] = tclosest

print,'found nearest times'
endif
;=======================

endif
;##############################################################

if mash['status'] ne 'no data' then ss = size(mash['aia_193'].data)

;##############################################################
if keyword_set(timecorr) then begin
;do a spatial correlation on the time series
;only required if not using AIA/IRIS cubes or if they have a correlation issue

a193 = mash['aia_193'].data

cc193 = tr_get_disp(a193)

for c = 0,5 do begin
	cstring = 'aia_'+channels[c]
	datastruct = mash[cstring]
	din = datastruct.data

	din = shift_img(din,cc193)

	datastruct.data = din
	mash[cstring] = datastruct
endfor

a1600 = mash['aia_1600'].data

;make a mask for the correlation
c1600 = a1600 gt 50.
cc1600 = tr_get_disp(c1600)

for c = 6,7 do begin
	cstring = 'aia_'+channels[c]
	datastruct = mash[cstring]
	din = datastruct.data

	din = shift_img(din,cc1600)

	datastruct.data = din
	mash[cstring] = datastruct
endfor



endif
;##############################################################



;correlation between EUV channels relative to 193

;use this step carefully as it seems to overshift the data occasionally when there are very bright transient features
;probably best to do the time correlation before this step

if keyword_set(wavecorr) then begin


s94 = mash['aia_94']
s131 = mash['aia_131']
s171 = mash['aia_171']
s193 = mash['aia_193']
s211 = mash['aia_211']
s335 = mash['aia_335']

a94 = s94.data
a131 = s131.data
a171= s171.data
a193 = s193.data
a211 = s211.data
a335 = s335.data


;ideally run time corr first
;aa = tr_get_disp(a193,/shift)

print,'running wavelength correlation'

for t=0,min(tlength)-1 do begin

	;use original data cube time sequence

	img = fltarr(ss[1],ss[2],6)

	img[*,*,0] = float(a94[*,*,t])
	img[*,*,1] = float(a131[*,*,t])
	img[*,*,2] = float(a171[*,*,t])
	img[*,*,3] = float(a193[*,*,t])
	img[*,*,4] = float(a211[*,*,t])
	img[*,*,5] = float(a335[*,*,t])

	pre_img = img
	aimg = fltarr(ss[1],ss[2],7)  ;leave space for alignment choice in first row
	aimg[*,*,0] = img[*,*,3]     ;using 193 as the key
	aimg[*,*,1:*] = img
	cc = tr_get_disp(aimg)

	;look for dodgy shifts
	badcorr = where(cc gt 2.,nb)
	if nb ge 1 then begin cc[badcorr] = 0.
		print,'ignoring bad shift for t = ',t
		endif
	simg = shift_img(aimg,cc)

	img = simg[*,*,1:*]
	corr_out = cc[*,1:*]

	a94[*,*,t] = img[*,*,0]
	a131[*,*,t] = img[*,*,1]
	a171[*,*,t] = img[*,*,2]
	a193[*,*,t] = img[*,*,3]
	a211[*,*,t] = img[*,*,4]
	a335[*,*,t] = img[*,*,5]


endfor

;put the correlated data back in the hash structure
s94.data = a94
s131.data = a131
s171.data = a171
s193.data = a193
s211.data = a211
s335.data = a335


;wave correlated data
mash['aia_94'] = s94
mash['aia_131'] = s131
mash['aia_171'] = s171
mash['aia_193'] = s193
mash['aia_211'] = s211
mash['aia_335'] = s335

mash['wavecorr done'] = 1

print,'finished wave correlation'


endif

;##############################################################
;use 1600 to locate network regions

if keyword_set(netfilter) then begin

ss = size(mash['aia_193'].data)

tclosest = mash['tclosest']

network_filter = intarr(ss[1],ss[2],ss[3])


for t = 0,ss[3]-1 do begin

	;closest 1600 to 193 step
	t1600 = tclosest[t,5]

	;================================
	;FIND MOSS REGIONS IN 1600

	DNchromo = 80.
	k = replicate(1,3,3)


	network_map = (mash['aia_1600'].data[*,*,t1600]) ge DNchromo
	network_map = dilate(dilate(network_map,k),k)

	network_filter[*,*,t] = network_map
endfor

mash['network mask'] = network_filter
print,'network mask created'

endif

;##############################################################
;use 1700 to locate network regions

if keyword_set(cnetfilter) then begin

;version using 1700

network_filter = intarr(ss[1],ss[2],ss[3])
tclosest = mash['tclosest']

a1700 = mash['aia_1700'].data

for t = 0,ss[3]-1 do begin

	;nearest times to 193
	t1700 = tclosest[t,6]

	;================================
	;FIND MOSS REGIONS IN 1700

	DNchromo = 1200.
	k = replicate(1,3,3)

	network_map = (a1700[*,*,t1700]) ge DNchromo
	network_map = dilate(dilate(network_map,k),k)

	network_filter[*,*,t] = network_map
endfor

mash['network mask'] = network_filter
print,'network mask 1700 created'
endif

;##############################################################
;use 193 to locate moss regions

if keyword_set(mossfilter) then begin

	DNmoss = 1250.

	moss_filter = (mash['aia_193'].data) ge DNmoss

	mash['moss mask'] = moss_filter

	print,'moss mask created'
endif

;##############################################################

if keyword_set(filter1700) then begin

;filter that uses the ratio of 1600/1700 to remove larger flare or filament brightening signatures

;upper ratio limit
rcut = 0.20

a1600 = mash['aia_1600'].data
a1700 = mash['aia_1700'].data

sso = size(mash['aia_1600'].data)

atime193 = anytim(mash['aia_193'].index.date_obs)
atime1600 = anytim(mash['aia_1600'].index.date_obs)
atime1700 = anytim(mash['aia_1700'].index.date_obs)

ratio = fltarr(sso[1],sso[2],sso[3])

tcorr = intarr(sso[3])

;nearest times to 193
tclosest = mash['tclosest']
t1600_193 = tclosest[*,5]
t1700_193 = tclosest[*,6]


for t = 0,sso[3]-1 do begin
	;1700 image closest to 1600
	tm = min(abs(atime1700-atime1600[t]),t1700)

	tcorr[t] = t1700

	;make sure that the images are within 30 seconds of each other
	if tm lt 30. then ratio[*,*,t] = float(a1600[*,*,t])/float(a1700[*,*,t1700])
endfor

precarbon = intarr(sso[1],sso[2],sso[3])
tbin = 4


if max(tcorr) ge sso[3] then begin tend = sso[3] - 1
	endif else begin
	tend = max(tcorr)
	endelse

for t = 0+tbin,tend do begin

	;binning by a few timesteps - tbin
	;mask is where the ratio is above rcut

	rmask = total((ratio[*,*,t-tbin:t] ge rcut),3) ge 1

	;dilate it a little
	k = replicate(1,3,3)
	rmask = dilate(dilate(rmask,k),k)

	;remove smallest regions
	ll = label_region(rmask)
	for i=0,max(ll) do begin
		ww = where(ll eq i,npix)
		if npix lt 150 then ll[ww]=0
	endfor

	precarbon[*,*,t] = ll ge 1
endfor

carbon = intarr(ss[1],ss[2],ss[3])

for t = 0,ss[3]-1 do begin
	;need closest 1600 time to 193
	t193 = tclosest[t,5]
	carbon[*,*,t] = precarbon[*,*,t193]
endfor

mash['carbon mask'] = carbon
mash['carbon mask done'] = 1
print,'Carbon mask created'
endif



;##############################################################
;Make Fe XVIII image

if keyword_set(fexviii) then begin

fe18 = fltarr(ss[1],ss[2],ss[3])

a94 = mash['aia_94'].data
a193 = mash['aia_193'].data
a211 = mash['aia_211'].data
a171 = mash['aia_171'].data
tclosest = mash['tclosest']


for t = 0,ss[3]-1 do begin

	;nearest times to 193 done here since it was not already in wavecorr step
	t94 = tclosest[t,0]
	t171 = tclosest[t,2]
	t211 = tclosest[t,3]

	;giulios trick
	fe18[*,*,t] = a94[*,*,t94] - (a211[*,*,t211]/120.) - (a171[*,*,t171]/450.)

endfor

mash['fe18'] = fe18
print,'Fe XVIII mask created'

endif



;##############################################################
;use 94 or Fe18 to locate hot regions

if keyword_set(loopfilter) then begin

	loop_filter = bytarr(ss[1],ss[2],ss[3])

	k = replicate(1,3,3)

	DNloop = 5.
	;DNloopwas 10 for 2018 paper
	minloopsize = 100.

	loops = mash['fe18'] ge DNloop
	loops = dilate(dilate(dilate(loops,k),k),k)  ;fill out the loops a little

	for t = 0,ss[3]-1 do begin

		looplabel = label_region(loops[*,*,t])

		;remove the smallest loops, or things that probably arent loops

		for i = 0,max(looplabel)-1 do begin
			ww = where(looplabel eq i, lsize)
			if lsize lt minloopsize then looplabel[ww] = 0
		endfor

		;COLLECT smoothed loop regions
		loop_filter[*,*,t] = looplabel ge 1
	endfor

mash['loop mask'] = loop_filter
print,'Loop mask created'

endif



;##############################################################

;run DEM on each time step


;compile:
;ar_dem_timestamp.pro
;dem_blos_imgplot.pro
;imgplot.pro

if keyword_set(rundem) then begin

temp_array = 5.5+findgen(21)*0.1

;t0 should be set to current data time
t0 = mash['aia_193'].index[0].date_obs

aia_sparse_em_init, use_lgtaxis = temp_array, timedepend=t0, bases_sigmas=[0.0,0.1,0.2,0.6]

mash['temp_array'] = temp_array


file_mkdir,areg+'/dem/'

dem4moss_cube,mash,dembinning,areg+'/dem/',/makedem,/save

endif

;##############################################################

;run DEM filters, remove transient brightenings

if keyword_set(demfilter) then begin

dfiles = find_files('emcube_'+tag+'.sav',areg+'/dem')

if keyword_set(demrestore) then begin
	restore,dfiles
	mash['dem'] = demstr
	;this should be in the array made in rundem
	mash['temp_array'] = 5.5+findgen(21)*0.1

endif

dmax = ss[3]-1

dem_filters_cube, mash, dembinning, bfilter

ds = size(bfilter)
bsize = dembinning

;Reshape the binned dem filter that the time base matches that of 193
if bsize gt 1 then begin
	dfilter = bytarr(ss[1],ss[2],ss[3])
	for i=0,ds[3]-1 do begin

		d = i*bsize
		dp = d+bsize-1

		if dp ge dmax then begin dfilter[*,*,d:-1] = rebin(bfilter[*,*,i],ss[1],ss[2],dmax-d+1)
		endif else begin
			dfilter[*,*,d:dp] = rebin(bfilter[*,*,i],ss[1],ss[2],bsize)
		endelse

	endfor
endif

mash['dem mask'] = dfilter
mash['dem mask done'] = 1
print,'DEM mask created'

endif


;##############################################################
if keyword_set(variability) then begin

a171 = mash['aia_171'].data
a193 = mash['aia_193'].data
tclosest = mash['tclosest']

ss = size(a193)

;difference image
;=== Start with 193 data

dcube = fltarr(ss[1],ss[2],ss[3])

dcube[*,*,1:*] = (a193[*,*,1:-1])-(a193[*,*,0:-2])
dcube[*,*,0] = dcube[*,*,1]

aia_error,a193,3,error

;for testing
dcube_193 = dcube

;=== perform zero crossing
ncube193 = zeroAIA(dcube,error,sigma,/crossreturn)

;=== 171 channel

se = size(mash['aia_171'].data)


;=== backup correlation for 193 and 171 incase wavecorr has issues

if mash['wavecorr done'] eq 0 then begin

	print,'running 171 to 193 wavecorr'
	shifted_171 = fltarr(se[1],se[2],se[3])

	for t=0,ss[3]-1 do begin

		t171 = tclosest[t,2]		;closest 171 to 193

		;-----test 171 to 193 correlation

		blinkarray = fltarr(ss[1],ss[2],2)
		blinkarray[*,*,0] = a193[*,*,t]
		blinkarray[*,*,1] = a171[*,*,t171]

		cc = tr_get_disp(blinkarray)

		;look for dodgy shifts
		badcorr = where(cc gt 2.,nb)
		if nb ge 1. then begin cc[badcorr] = 0.
			;print,'ignoring bad shift for t = ',t
			endif
		blinkarray = shift_img(blinkarray,cc)

		shifted_171[*,*,t171] = blinkarray[*,*,1]
		;print,cc
	endfor

	a171 = shifted_171

endif
;================



ecube = fltarr(se[1],se[2],se[3])

ecube[*,*,1:*] = (a171[*,*,1:-1])-(a171[*,*,0:-2])
ecube[*,*,0] = ecube[*,*,1]

aia_error,a171,2,error


;=== perform zero crossing
ncube171 = zeroAIA(ecube,error,sigma,/crossreturn)


;=== dilate 171 detections
k=[[0,1,0],[1,1,1],[0,1,0]]

;fill out a171
dil_map = dilate(dilate(dilate(ncube171.zcube,k),k),k)

zmatch = bytarr(ss[1],ss[2],ss[3])

;make 171 to nearest 193 check
for t=0,ss[3]-1 do begin
	t171 = tclosest[t,2]
	zmatch[*,*,t] = dil_map[*,*,t171] AND ncube193.zcube[*,*,t]
endfor


;---------------------------------------------------------------------------
;remove large areas of zero crossings - probably flare related
if keyword_set(cleanflare) then begin

sz = size(zmatch)

for t = 0,sz[3]-1 do begin

	lmask = label_region(zmatch[*,*,t])
	rmask = zmatch[*,*,t]

	;remove some of the big regions

	for i=0,max(lmask)-1 do begin
		ww = where(lmask eq i, sregion)
		if sregion gt 15 then rmask[ww] = 0
	endfor

	zmatch[*,*,t] = rmask

endfor

endif

;---------------------------------------------------------------------------


;add check to see if dem fitler has been run
;if mash['dem mask'] gt 0 then zmatch = zmatch - (zmatch AND outputs.dfilter)

;if n_dimensions(mash['carbon mask']) gt 1 then begin
;	zmatch = zmatch - (zmatch AND outputs.carbon)
;	print,'1600/1700 filter done'
;endif
plot_image,total(zmatch,3) gt 1


if mash['dem mask done'] eq 1 then begin
	dmask = mash['dem mask']
	zmatch = zmatch - (zmatch AND dmask)
	print,'used dem mask'
endif

if mash['carbon mask done'] eq 1 then begin
	cmask = mash['carbon mask']
	zmatch = zmatch - (zmatch AND cmask)
	print,'used carbon mask'
endif

netmask = mash['network mask']

zmatch = byte(zmatch AND mash['network mask'] AND mash['moss mask'])
zmatch_d = byte(zmatch AND mash['network mask'] AND mash['moss mask'] AND mash['loop mask'])

mash['zmatch no loop'] = zmatch
mash['zmatch'] = zmatch_d


plot_image,total(zmatch,3) gt 1


endif


;##############################################################

if keyword_set(movie) then begin
	file_mkdir,areg+'/movies/'
	id = mash['aia_193'].index[0].obsid
	movie_quick_cube,mash,areg+'/movies/ar'+id+'.mp4',0,0,0,0,0,0,/comp

	;movie_quick_cube,mash,ar'+tag+'.mp4',0,0,0,0,0,0,/comp
endif

if keyword_set(local_movie) then begin
	print,'Making Movie'
	file_mkdir,'movies/'
	id = mash['aia_193'].index[0].obsid
	movie_quick_cube,mash,'movies/ar'+id+'.mp4',0,0,0,0,0,0,/comp
	print,'Finished Movie'
endif


;##############################################################

;could return a structure instead of the hash by using .... someoutvariable = mash.tostruct()
end
