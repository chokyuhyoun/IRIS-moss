; fwhm=0.0286 & spec_iris=0.01298 & wvl_int=[1334., 1335.] for cii
; fwhm=0.05054  & spec_iris=0.02546 & wvl_int=[ 2802.3, 2804.6,] for
; Mgii h
; fwhm=0.05054  & spec_iris=0.02546 & wvl_int=[2797,2800 ] for Mgii
; triplet
; fwhm=0.05054  & spec_iris=0.02546 & wvl_int=[2795.5,2797 ] & dt=0.5

;
; infile is the output_ray.ncdf file inside the output folder from RH
;out_spectra is the output synthetic spectra (2D array as a function of wavelength wvl_iris and time) for Mg IIk or h.
 
pro inter_rh_spectra, infile, out_spectra, wvl_iris, spec,time, fwhm=fwhm,spec_iris=spec_iris, wvl_int=wvl_int, air_to_vacuum=air_to_vacuum, dt=dt
 
; IRIS PSF


; first psf and then interpolation
  
  wvl_int0=wvl_int[0] & wvl_int1=wvl_int[1]
  intensity=read_ncdf_var(infile, 'intensity')
  wave=read_ncdf_var(infile, 'wavelength')
  sz=size(intensity)
  time=indgen(sz[3])*dt
for j=0, sz[3]-1 do begin
 ind=where(intensity[*,0,j] ge 9.96921e+36) 
if ind[0] ne -1 then intensity[ind, 0,j]=0
endfor

 
;; stop
if keyword_set(air_to_vacuum) then wave_vacuum=airtovacuum(wave)*10 else wave_vacuum=wave*10

  ind=where(wave_vacuum lt wvl_int1  and wave_vacuum gt wvl_int0)
  wvl0=double(wave_vacuum[ind])
  dat0=double(intensity[ind, *,*])
  spec=min(wvl0[1:*]-wvl0[0:*])
  ;spec_iris=0.053; IRIS spectral lambda
  
  nwvl=round((wvl_int1-wvl_int0)/spec)
  nwvl_iris=round((wvl_int1-wvl_int0)/spec_iris)
  wvl=double(wvl_int0+indgen(nwvl+1)*spec) 
  wvl_iris=double(wvl_int0+indgen(nwvl_iris+1)*spec_iris) 
  delta_lambda=wvl0*0.
  delta_lambda[1:*]=(wvl0[1:*]-wvl0[0:*]) & delta_lambda[0]=delta_lambda[1]
  delta_nu=3d18*delta_lambda/wvl0^2.
 delta_lambda1=wvl*0.
 delta_lambda1[1:*]=(wvl[1:*]-wvl[0:*]) & delta_lambda1[0]=delta_lambda1[1]
  delta_nu1=3d18*delta_lambda1/wvl^2.
  
delta_lambda_iris=wvl_iris*0.
 delta_lambda_iris[1:*]=(wvl_iris[1:*]-wvl_iris[0:*]) & delta_lambda_iris[0]=delta_lambda_iris[1]
  delta_nu_iris=3d18*delta_lambda_iris/wvl_iris^2.

  nt=(size(dat0))[3]
  dat=fltarr(n_elements(wvl_iris), 1,nt)
  for ii=0, nt-1 do begin
;  da=interpol(dat0[*,0,ii], wvl0, wvl, /spline); interpolate on uniform grid with the minimum rh grid 
;  da_convol=gaussfold(wvl, da, fwhm) ; convolve with IRIS PSF
;    ; re-interpolate on iris_grid
; da_convol_interp=interpol(da_convol, wvl, wvl_iris, /spline);
;  dat[*,0,ii]=da_convol_interp
    
    da = interpol(dat0[*, 0, ii], wvl0, wvl_iris, /spline)
    dat[*, 0, ii] = gaussfold(wvl_iris, da, fwhm)
  endfor
 out_spectra=reform(dat)

;stop
end

