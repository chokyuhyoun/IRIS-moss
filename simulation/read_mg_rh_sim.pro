cd, '/Users/khcho/Desktop/IRIS-moss-main/'

save_files = file_search('/Users/khcho/Desktop/IRIS-moss-main/simulation/rh_output_kh/', $
  'output_ray.ncdf')
cd, current=current0


;w01 = window(dim=[8d2, 8d2])
;p01 = plot([0, 0], [-20, 20], ':', /current, xr=[-20, 20], yr=[-20, 20], $
;            xtitle='v$_{D, Mg II h3}$ (km $s^{-1}$)', ytitle='v$_{D, Mg II k3}$ (km $s^{-1}$)')
;p02 = plot(20*[-1, 1], [0, 0], ':', over=p01)
;p03 = plot(20*[-1, 1], 20*[-1, 1], ':', over=p01)
;p04 = plot(20*[1, -1], 20*[-1, 1], ':', over=p01)

for i=0, n_elements(save_files)-1 do begin       
  i = 11              
  save1=save_files[i]
  print, save1
  
  ;; constants & factors we'll need to convert from RH units to DN, change the ones in bold to synthesise Mg II k/h line
  h =6.63*1.0d-34 ;; m^2 * kg * s^-1
  cc = 3.0d8     ;;m * s^-1
  wvl_ref = 2796.3509493 ; 2803.5297192 for mgii h ; Angstrom
  wvl_m = wvl_ref*1d-10
  joules_to_ph = h*cc/(wvl_m); joules/ph
  energy = (h*cc)/wvl_m/1.d-7 ; erg/ph
  phperdn = 18. ;for NUV
  iresp = iris_get_response()
  i1 = where(iresp.lambda ge 279.63509493d0) 
  aresp_mg2 = (iresp.AREA_SG)[i1[0],1] ;cm^2 EA at SiIV
  sr_resp = aresp_mg2/(1.496d13)^2.
  apix = ((0.166*720d5)*(0.33*720d5))
  fwhm = 0.05054  & spec_iris=0.02546 & wvl_int=[2795,2805]; wvl_int=[ 2802.3, 2804.6,] for Mgii h
  
  dt = 1; 1s= time interval of nanoflare simulations
  
  inter_rh_spectra, save1,sp1, wvl, spec,time1, fwhm=fwhm,spec_iris=spec_iris, wvl_int=wvl_int, dt=dt;, /air_to_vacuum
  vel=(wvl-wvl_ref)/wvl_ref*3.e5
  
  ;; convert from RH units to DN s-1 pix-1
  sp1=sp1*1e3*cc*1e10/(wvl_ref)^2.
  ; J s-1 m-2 Hz-1 sr-1 *1.e7(joules to erg)*1.e-4(m-2 to cm-2)*cc/(lambda^2.)-> erg s-1 cm-2 A-1 sr-1 ;
  sp1=sp1*(1./energy)*sr_resp*apix*spec_iris/phperdn 
  ;  erg s-1 cm-2 A-1 sr-1* ph/erg*sr*cm2*A*DN/ph*=DN s-1 pix-1
  ;  
  sim_mg2_res = !null
  for ii=0, 40 do begin
    res0 = mg_ii_fit(wvl, sp1[*, ii])
    sim_mg2_res = [[sim_mg2_res], [res0]]
  endfor
  stop
  sp_time_mg2 = sp1[*, 40]
  cd, file_dirname(save1)
  cd, '..'
  cd, current=current
  wv = wvl
  save, sp_time_mg2, wv, sim_mg2_res, filename=current+'/mg_ii_fit_res.sav'
  
;  p05 = plot(res[1, *].v_d_3, res[0, *].v_d_3, '.r', sym_size=2, over=p01)
;  w01.save, current+'mg_ii_fit_res.png', resol=200
;  p05.delete
  cd, current0
endfor
cd, '..'

end