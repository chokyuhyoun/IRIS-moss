dir = '/Users/khcho/Desktop/IRIS-moss-main/'
cd, dir
out_dir = dir+'paper_fig/'

file_94 = file_search(dir, '*lev1.94A*.fits')
file_193 = file_search(dir, '*lev1.193A*.fits')
file_cont = file_search(dir, '*hmi.ic*.fits')
file_mag = file_search(dir, '*hmi.m_*.fits')
files = [[file_cont], [file_mag], [file_94], [file_193]]
wave = [6173, 1, 94, 193]
wave_title = ['HMI Int. ', 'HMI B$_{los}$ ', 'AIA 94$\AA$', 'AIA 193$\AA$']
dr = [[1d4, 6d4], [-1d3, 1d3], [0, 50], [100, 3d3]]
ar_no = ['12415', '12473', '12524', '12529']

w01_sz = [8d2, 10d2]
mar = [1d2, 5d1]
ygap = 8.5d1
p01_sz = (w01_sz[0]-total(mar))/n_elements(ar_no)
hfov = 200
cen = [[427, -422], [295, -305], [-174, 360], [-40, 265]]

w01 = window(dim = w01_sz)
im01 = objarr(4, 4)
rr = 2048
for i=0, n_elements(file_cont)-1 do begin ; wavelength
  for j=0, n_elements(ar_no)-1 do begin   ; ar no
    read_sdo, files[j, i], index, data, /sil
    if i le 1 then data = rotate(temporary(data), 2)
    xp = (findgen(index.naxis1)-index.crpix1)*index.cdelt1 + index.crval1
    yp = (findgen(index.naxis2)-index.crpix2)*index.cdelt2 + index.crval2
    pos = [mar[0]+p01_sz*i,     w01_sz[1]-4.5d1-ygap*j-p01_sz*(j+1), $
           mar[0]+p01_sz*(i+1), w01_sz[1]-4.5d1-ygap*j-p01_sz*j]
    data = rebin(temporary(data), rr, rr)
    xp = rebin(temporary(xp), rr)
    yp = rebin(temporary(yp), rr)
    im01[i, j] = image_kh(data, xp, yp, /current, /dev, pos=pos, $
                    rgb_table = aia_rgb_table(wave[i]), $
                    min=dr[0, i], max=dr[1, i], $
                    xr=cen[0, j]+hfov*[-1, 1], yr=cen[1, j]+hfov*[-1, 1], $ 
                    xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)', $
                    xminor=1, yminor=1, $
                    font_size = 12)
    im01[i, j].axes[1].showtext = (i eq 0) ? 1 : 0
    t01 = text(pos[0]+10, pos[3]+5, wave_title[i]+'!c'+strmid(index.date_obs, 0, 19), $
               /dev, font_size=10, $
               align=0, vertical_align=0)
    if i eq 0 then begin
      t02 = text(pos[0]+10, pos[3]-10, 'AR '+ar_no[j], $
                 /dev, font_size=13, $
                 align=0, vertical_align=1)
    endif
  endfor
endfor
w01.save, out_dir+'ar_image.pdf', resol=200

end