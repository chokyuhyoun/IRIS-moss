pro var_chk, curve0, dt, fwhm, i_0, the_peak_ind, init_ind_diff, end_ind_diff
  curve = curve0
  curve[where(curve lt 0., /null)] = 0.
  t_lim_sji = n_elements(curve)/2   ; t_lim_sji = center of curve = reference time for sg 
  binsize = 50
  dum = histogram(curve, loc=xbin, binsize=binsize)
  i_0 = xbin[(where(dum eq max(dum)))[0]] + 0.5*binsize
  d_curve = curve[1:*] - curve[0:-2]
  peak_ind = where(sgn(d_curve[1:*]) - sgn(d_curve[0:-2]) eq -2) + 1
  exc_peak_val = curve[peak_ind]-i_0
  the_peak_ind = peak_ind[where(exc_peak_val gt 3.*stddev(curve-i_0, /nan) and $
                                 abs(peak_ind-t_lim_sji) le 30./dt, $
                                 count)]
;  stop
  if count eq 0 then return   ; No peak
  the_peak_ind = the_peak_ind[0]
  the_peak_val = curve[the_peak_ind]
  hm = 0.5*(the_peak_val - i_0) + i_0

  left_curve = reverse(curve[0:the_peak_ind])
  if total(left_curve lt hm) eq 0 then return
  inip = min(where(left_curve lt hm))
  left_ind = interpol(findgen(inip+1), left_curve[0:inip], hm)

  right_curve = curve[the_peak_ind:*]
  if total(right_curve lt hm) eq 0 then return
  endp = min(where(right_curve lt hm))
  right_ind = interpol(findgen(endp+1), right_curve[0:endp], hm)
  if (t_lim_sji lt floor(the_peak_ind-left_ind)) or $
     (t_lim_sji gt ceil(the_peak_ind+right_ind)) then return
  init_ind_diff = floor(the_peak_ind-left_ind*2.) - t_lim_sji
  end_ind_diff = ceil(the_peak_ind+right_ind*2.) - t_lim_sji
  fwhm = (right_ind + left_ind)*dt
;  print, fwhm
;  stop
  if fwhm gt 60. then fwhm = !null
  
end