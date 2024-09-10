sav_files = file_search('/Users/khcho/Desktop/IRIS-moss-main/simulation/rh_output_kh', $
                        'sp_si4_*.save')
si_cen = 1402.77d0                        
for i=0, n_elements(sav_files)-1 do begin
;  i = 10
  restore, sav_files[i]
  sav_dir = file_dirname(sav_files[i])
  sim_si4_res = !null
  dlam = mean(wvl[1:*]-wvl[0:-2])
  nwvl = wvl[0]+findgen(n_elements(wvl))*dlam
  wavep_ind = where(nwvl gt 1402.2d0 and nwvl lt 1403.2d0)
  wv0 = nwvl[wavep_ind]
  sp_time_si4 = sp_time[wavep_ind, *]/2./4. ; spectral binning = 2 / spatial binning = 1 / exposure time = 4 s
  for j=0, 50 do begin
    res = si_iv_fit(wvl[wavep_ind], sp_time_si4[*, j], den_cal=0, spec_bin = 2)
    sim_si4_res = [sim_si4_res, res]
  endfor
  time0 = time
;  stop
  save, sp_time_si4, wv0, time0, sim_si4_res, $
        filename=sav_dir+'/si_iv_fit_res.sav'
endfor
                        
end                        