a1600 = mash['aia_1600'].data
a1600index = mash['aia_1600'].index

time1600 = mash['aia_1600'].index.date_obs

tiris  = 100
tm = min( abs( anytim(time1600)-anytim(index_sji[tiris].date_obs)) ,taia)


read_iris_l2, irisfiles[0], index_sji, data_sji

index2map,index_sji[100],data_sji[*,*,100],imap
index2map,a1600index[taia],a1600[*,*,taia],amap

;for finding the average iris and aia maps, we need the average coordinates
;find the first and last maps

index2map,index_sji[0],data_sji[*,*,0],imap0
index2map,index_sji[-1],data_sji[*,*,-1],imap1

index2map,a1600index[0],a1600[*,*,0],amap0
index2map,a1600index[-1],a1600[*,*,1],amap1


xfac = amap.dx/imap.dx
yfac = amap.dy/imap.dy


jsize = size(imap0[0].data)
sxl = jsize[1]
syl = jsize[2]

dsize = size(amap0.data)
xl = dsize[1]
yl = dsize[2]
tl = dsize[3]

mean_ixc = (imap0.xc + imap1.xc)/2.
mean_iyc = (imap0.yc + imap1.yc)/2.
mean_axc = (amap0.xc + amap1.xc)/2.
mean_ayc = (amap0.yc + amap1.yc)/2.

mean_idx = (imap0.dx + imap1.dx)/2.
mean_idy = (imap0.dy + imap1.dy)/2.
mean_adx = (amap0.dx + amap1.dx)/2.
mean_ady = (amap0.dy + amap1.dy)/2.


;create x and y solar position axes (IRIS)
ixleft = mean_ixc - ((sxl/2.)*mean_idx)
ixpos = ixleft + ixran
iyleft = mean_iyc - ((syl/2.)*mean_idy)
iypos = iyleft + iyran

;create x and y solar position axes (AIA)
xleft = mean_axc - ((xl/2.)*mean_adx)
xpos = xleft + xran
yleft = mean_ayc - ((yl/2.)*mean_ady)
ypos = yleft + yran

;find the edges of the SJI window on AIA
mm1 = min(abs(xpos-ixpos[0]),mleft)
mm2 = min(abs(ypos-iypos[0]),mbottom)
imleft = mleft*xfac
imbottom = mbottom*yfac

;do correlation
idata = total(data_sji,3)
adata = total(a1600,3)

;adata = amap.data
;idata = imap.data

xfac = amap.dx/imap.dx
yfac = amap.dy/imap.dy

ardata = congrid(adata,xl*xfac, yl*yfac)
arcut = ardata[imleft:(sxl-1)+imleft,imbottom:(syl-1)+imbottom]

ss = size(idata)
carray = fltarr(ss[1],ss[2],2)

carray[*,*,0] = arcut
carray[*,*,1] = idata

a2i=tr_get_disp(carray)

;shifting the map
new_imap = shift_map(imap,a2i[0,1]/xfac,a2i[1,1]/yfac)
