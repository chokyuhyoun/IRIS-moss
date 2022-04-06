PRO align_iris_aia, iris_index, iris_data, aia_index, aia_data, xshift,yshift, shift_check=shift_check

;tool to align prepared AIA CUBES with corresponding IRIS CUBES
;returns X and Y shift coeffecients using the time averaged intensity for each cube


;for finding the average iris and aia maps, we need the average coordinates
;find the first and last maps

index2map,iris_index[0],iris_data[*,*,0],imap0
index2map,iris_index[-1],iris_data[*,*,-1],imap1

index2map,aia_index[0],aia_data[*,*,0],amap0
index2map,aia_index[-1],aia_data[*,*,1],amap1

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

xfac = mean_adx/mean_idx
yfac = mean_ady/mean_idy

;x and y arrays (un-centered)
ixran = findgen(sxl)*mean_idx
iyran = findgen(syl)*mean_idy

;x and y ranges
xran = findgen(xl)*mean_adx
yran = findgen(yl)*mean_ady


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
idata = total(iris_data,3)
adata = total(aia_data,3)

ardata = congrid(adata,xl*xfac, yl*yfac)
arcut = ardata[imleft:(sxl-1)+imleft,imbottom:(syl-1)+imbottom]

ss = size(idata)
carray = fltarr(ss[1],ss[2],2)

carray[*,*,0] = arcut
carray[*,*,1] = idata

a2i = tr_get_disp(carray)

;convert shift back to AIA pixel scale
xshift = a2i[0,1]/xfac
yshift = a2i[1,1]/yfac

if keyword_set(shift_check) then begin
    ;shifting the map
    afinal = amap0
    afinal.data = adata
    afinal.dx = mean_adx
    afinal.dy = mean_ady
    afinal.xc = mean_axc
    afinal.yc = mean_ayc

    ifinal = imap0
    ifinal.data = idata
    ifinal.dx = mean_idx
    ifinal.dy = mean_idy
    ifinal.xc = mean_ixc
    ifinal.yc = mean_iyc

    plot_map,amap0
    plot_map,imap0,/cont,/over
    plot_map,shift_map(imap0,xshift,yshift),/cont,/over
endif

return

end
