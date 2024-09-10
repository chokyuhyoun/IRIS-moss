PRO map_fov, map, xfov, yfov

;INPUT
;takes a map structure
;OUTPUT
;returns the X and Y grid in arcsecs

msize = size(map.data)
xl = msize[1]
yl = msize[2]

xran = findgen(xl)*map.dx
yran = findgen(yl)*map.dy

;create x and y solar position axes
xleft = map.xc - ((xl/2.)*map.dx)
xfov = xleft + xran

yleft = map.yc - ((yl/2.)*map.dy)
yfov = yleft + yran

return

end
