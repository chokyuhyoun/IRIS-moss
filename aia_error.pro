PRO aia_error,datacube,chanid,error

ss = size(datacube)
error = fltarr(ss[1],ss[2],ss[3])

;hacked version of aia_bp_estimate_error by DG 22 Jun 2017



	;t - time indice of the datacube

	errtable = AIA_BP_READ_ERROR_TABLE()
	dnperpht = errtable[chanid].dnperpht



for t=0,ss[3]-1 do begin

	counts = reform(datacube[*,*,t])

	n_sample = 1
	sumfactor = SQRT(n_sample) / n_sample

	nullcount = fltarr(ss[1],ss[2])


	numphot = counts / dnperpht
	stdevphot = SQRT(numphot)
	shotnoise = (stdevphot * dnperpht)
	shotnoise = nullcount + shotnoise / SQRT(n_sample)

	darknoise = nullcount + 0.18

	readnoise = nullcount + 1.15 * sumfactor	

	quantnoise = nullcount + 0.288819 * sumfactor


	compressratio = errtable[chanid].compress
	compressnoise = (shotnoise / compressratio) > 0.288819
		lowcounts = WHERE(counts lt 25, numlow)	;	Linear portion of lookup table
			if numlow gt 0 then compressnoise[lowcounts] = 0.
				compressnoise = compressnoise * sumfactor
			

	sigmadn = SQRT(shotnoise^2. + darknoise^2. + readnoise^2. + quantnoise^2. + compressnoise^2.)
	snr = counts / sigmadn

	error[*,*,t] = sigmadn

endfor


end