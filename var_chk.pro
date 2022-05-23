pro var_chk, curve, sji_index, i_0, the_peak_val, the_peak_pos, fwhm
  t_lim_sji = n_elements(curve)/2
  dum = histogram(curve, loc=xbin, binsize=100)
  i_0 = xbin[(where(dum eq max(dum)))[0]]
  d_curve = curve[1:*] - curve[0:-2]
  peak_pos = where(sgn(d_curve[1:*]) - sgn(d_curve[0:-2]) eq -2) + 1
  exc_peak_val = curve[peak_pos+1]-i_0
  real_peak_pos = peak_pos[where(exc_peak_val gt 3.*sqrt(i_0) and $
                           peak_pos gt t_lim_sji-24./sji_index[0].cdelt3, count)]
  if count eq 0 then return   ; No peak 
  dist0 = min(abs(real_peak_pos - t_lim_sji), dum_ind)
  the_peak_pos = real_peak_pos[dum_ind]
  the_peak_val = curve[the_peak_pos]
  hm = 0.5*(the_peak_val - i_0)

  left_curve = reverse(curve[0:the_peak_pos])
  if total(left_curve lt hm) eq 0 then return
  inip = min(where(left_curve lt hm))
  left_pos = interpol(findgen(inip+1), left_curve[0:inip], hm)
  
  right_curve = curve[the_peak_pos:*]
  if total(right_curve lt hm) eq 0 then return
  endp = min(where(right_curve lt hm))
  right_pos = interpol(findgen(endp+1), right_curve[0:endp], hm)
  fwhm = (right_pos + left_pos)*sji_index[0].cdelt3
  if fwhm gt 60. then fwhm = !null
end