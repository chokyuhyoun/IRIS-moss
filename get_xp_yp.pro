pro get_xp_yp, index, xp, yp
  xp = (findgen(index.naxis1) - 0.5*(index.naxis1-1))*index.cdelt1 + index.crval1
  yp = (findgen(index.naxis2) - 0.5*(index.naxis2-1))*index.cdelt2 + index.crval2
end

