testa_2020 = ['2014-02-23T23:23', $
              '2014-04-10T02:46', $
              '2014-09-17T14:45', $
              '2014-09-17T17:20', $
              '2014-09-18T08:08', $ 
              '2015-01-29T18:29', $
              '2015-11-11T02:47', $
              '2015-11-12T01:38', $
              '2015-12-24T15:21', $
              '2016-01-29T06:23', $
              '2013-11-09T12:15']   ; for Testa et al. (2014)
testa_2020_sec = anytim(testa_2020)
testa_2020_end = anytim(testa_2020_sec+600, /ccsds)

spec_bin_arr = !null
spat_bin_arr = !null
for i=0, n_elements(testa_2020)-1 do begin
  hcr=ssw_hcr_query('https://www.lmsal.com/hek/hcr?cmd=search-events3&outputformat=json&'+$
                    'startTime='+testa_2020[i]+'&stopTime='+testa_2020_end[i]+'&hasData=true&hideMostLimbScans=true&limit=200')
  
  umodes = (strsplit(hcr[0].umodes, ',', /extract))[0]
  strput, umodes, '.', strlen(umodes)-5
  iris_dir = file_dirname(umodes)
  f = file_search(iris_dir, '*raster*.fits')
  dd = iris_obj(f[0])
  line_id = dd->getline_id()
  Si_IV_ind = (where(strmatch(line_id, '*Si IV*'), /null))[-1]
  spat_bin = dd->binning_region('FUV')
  spec_bin = dd->binning_spectral(si_iv_ind)
  spat_bin_arr = [spat_bin_arr, spat_bin]
  spec_bin_arr = [spec_bin_arr, spec_bin]
  print, spat_bin, spec_bin
  obj_destroy, dd
  hcr = 0
endfor
end