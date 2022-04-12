PRO find_iris_moss_events3, hcr,outdir,eout,plotevent=plotevent,savesample=savesample,movie=movie

;returns a list of frames where moss variability events are found near the IRIS slit using the AIA MOSSVAR code

;INPUTS
;HCR - takes the HCR structure from an iris data search
;can be used with multiple HCR results or a single HCR structure
;OUTDIR - where to save the output structures


;OUTPUTS
;EOUT - the output structure - can be returned after running one HCR file - useful for testing
;structure containing IRIS Spectrograph and SJI X and Y positions of detected events with the following info

;SG_TIME - spectrograph time
;CLOSEST_SJI2EVENT - index for SJI file containing closest SJI image to event
;SJI_WAVE - wavelength of this SJI file
;AIA_STEP - time step of AIA CUBE at event
;AIA_TIME - time in UT of above
;SJI_TIME_INDEX - time step of closest SJI exposure in correct SJI file
;IRIS_SG_XPOS - spectrograph x position
;IRIS_SJI_XPOS - SJI x position
;IRIS_SG_YPOS - spectrograph y position
;IRIS_SG_XPIXEL - spectrograph x pixel index
;IRIS_SJI_XPIXEL - SJI x pixel index
;IRIS_SG_SCAN - spectrograph raster file number

;HEADER - the IRIS hcr header used to process this data
;ID - raster obs ID and date code
;SJI_FILES and AIAFILES - files to locate the cubes
;SJI_MAP - SJI map for each dectection
;EVENT_MAP - event detections per frame from the AIA MOSSVAR routine
;X/YSHIFT_IRIS2AIA - x and y offsets used to correct alignment
;OBS - obstype
;SCAN_WIDTH - number of exposures per raster scan or SNS length
;SCAN_REPEATS - number of raster scans (1 for SNS)

;OPTIONS
;PLOTEVENT - creates a plot for every event detection showing the AIA events and select event position on the slit
;SAVESAMPLE - records the SJI and AIA EVENT map structures for each event
;MOVIE - saves the AIA event movie from the MOSSVAR code - see its own instructions

;EXAMPLES
;to return just the mossvar result

;UPDATES
;27/11/19 - cleaned up outputs - added check to remove events from outside SJI positive data area
;25/11/19 - moving SJI collection to before event loop
;25/11/19 - bug fixes - SJI plotter working but checking if it has correct index - also removing old variables
;21/11/19 - Adding full SJI plotting mode after event detection
;20/11/19 - Check for SJI 1400 to use in correlation - if not found use the shortest wavelength available (usually 1330)
;14/11/19 - RASTER version
;8/11/19 - lite version to prepare for RASTER mode
;7/11/19 - Added AIA to IRIS image corellation check.
;Takes the total image of each cube and finds the offet
;AIA coordinates are automatically updated for use in finding the slit position

nresults = n_elements(hcr)

slit_x = fltarr(nresults)
date = strarr(nresults)
obs_short = strarr(nresults)
udir = strarr(nresults)
number_mossframes = fltarr(nresults)

;how far to search around slit
slit_offset = 1.0

;start HCR loop
for i=0,nresults-1 do begin
    slit_x[i]=hcr[i].xcen
    date[i] = hcr[i].starttime
    obs_short[i] = hcr[i].iris_obsshort
    uu = hcr[i].umodes
    us = str_sep(uu,'/iris_')
    udir[i] = us[0]
    aiafiles = us[0]+'/aia'
    hsep = str_sep(us[0],'/')
    hname = hsep[-1]

    ;SETUP USEFUL STRUCTURE TO STORE IRIS DATA CUBES
    readme = 'IRIS SJI data and more'
    storeplot = orderedhash('stores SJI maps for plotting later',readme)

    ;GET IRIS SG FILES
    sg_files = find_files('*raster*',us[0])
    if sg_files eq '' then mash['status'] = 'no sg data'
    nraster_scans = n_elements(sg_files)


    ;RUN MOSSVAR ON AIA CUBE TO DETECT EVENTS
    mossvar,udir[i],'',mash,4,4,/read,/cube_data
    if mash['status'] eq 'data ok' then begin

        mossvar,udir[i],'',mash,4,4,/moss,/cnet,/fexviii,/loop,/filter1700
        if keyword_set(movie) then mossvar,udir[i],'',mash,4,4,/variability,/cleanflare,/local_movie
        if ~keyword_set(movie) then mossvar,udir[i],'',mash,4,4,/variability,/cleanflare

        ;cube of counted moss events
        zcube = mash['zmatch no loop']

        index193 = mash['aia_193'].index
        data193 = mash['aia_193'].data
        time193 = mash['aia_193'].index.date_obs

        index2map, index193, data193 ,map193

        ;take what's needed out of the hash structure
        aia1600index = mash['aia_1600'].index
        aia1600data =  mash['aia_1600'].data

        ;clean up memory
        delvar,mash

        ;AIA ZERO CROSSING CUBE
        zmap = map193
        zmap.data = zcube
        print,'Event detection completed'
        print,'created zcube map from MOSSVAR'

        dsize = size(zmap.data)
        xl = dsize[1]
        yl = dsize[2]
        tl = dsize[3]

        ;get SJI index
        sji_files = find_files('*SJI*',udir[i])
        storeplot['sji_files'] = sji_files
        num_sji = n_elements(sji_files)
        sji_wins = strarr(num_sji)

        ;number of sji exposures in each wavelength
        key_sji = intarr(num_sji+1)

        sji_timeline = 'empty'

        ;STORE THE SJI DATA AND OBSERVATION TIMES
        for ss = 0,num_sji-1 do begin
            read_iris_l2, sji_files[ss], sindex, sdata
            sji_wins[ss] = sindex[0].TDESC1

            dstr = create_struct('sindex',sindex,'sdata',sdata)
            storeplot['dstr_'+sji_wins[ss]] = dstr

            storeplot['time_'+sji_wins[ss]] = sindex.date_obs
                ;concantenate obstimes into one vector
            sji_timeline = [sji_timeline,sindex.date_obs]

            ;get the start index of the wavelength
            key_sji[ss+1] = key_sji[ss] + n_elements(sindex.date_obs)
        endfor

        sji_timeline = sji_timeline[1:*]
        ;change this to the END index of each wavelength
        key_sji = key_sji-1
        key_sji[0] = 0


        ;COALIGNMENT BETWEEN AIA AND IRIS, CHECK PERFORMED HERE - provides correlation coefficients

        ;USE 1400 but pick the next closest if 1400 is unavailable
        scheck = where(sji_wins eq 'SJI_1400')
        if scheck ne -1 then begin
            print,'correlation with '+sji_wins[scheck]
            align_iris_aia, storeplot['dstr_SJI_1400'].sindex, storeplot['dstr_SJI_1400'].sdata, aia1600index, aia1600data ,xcoeff,ycoeff

        endif else begin
            print,'no SJI 1400 found reading '+sji_wins[0]
            align_iris_aia, storeplot['dstr_'+sji_wins[0]].sindex, storeplot['dstr_'+sji_wins[0]].sdata, aia1600index, aia1600data ,xcoeff,ycoeff
        endelse

        delvar,aia1600index
        delvar,aia1600data

        ;modify AIA EVENT MAP coordinates
        ;caution if comparing to original AIA cube - this modifies the AIA cube
        zmap.xc = zmap.xc - xcoeff
        zmap.yc = zmap.yc - ycoeff


        ;GET IRIS SG DATA
        ;sg_files = find_files('*raster*',us[0])
        ;nraster_scans = n_elements(sg_files)

        ;READ the first spectrograph file to check if raster or sit and stare
        read_iris_l2, sg_files[0], index_sg

        if index_sg[0].cdelt3 eq 0. then obstype = 'sit and stare'
        if index_sg[0].cdelt3 gt 0. then obstype = 'raster'

        print,'observation is a ',obstype

        nraster_steps = index_sg[0].nrasterp
        rstep_size = index_sg[0].steps_av

        rlen = n_elements(index_sg.date_obs)

        ;total SG exposures in this HCR entry
        tot_exposures = rlen * nraster_scans

        ;OUTPUTS
        ;define vector to store event positions
        iris_sg_xpos = 0
        iris_sji_xpos = 0
        iris_sg_ypos = 0
        iris_sg_xpixel = 0
        iris_sji_xpixel = 0
        iris_sg_ypixel = 0
        tsji_closest = 0

        sji_time_index = 0
        sji_time = 'empty'
        iris_sji_wave = 'empty'
        iris_sg_scan = 0
        aia_step = 0
        sg_time = 'empty'
        aia_time = 'empty'

        expcount = 0

        ;read all of the slit solar x positions from the sg files
        sg_slit = fltarr(nraster_scans,rlen)

        ;create dummy maps to make a save array
        index2map,storeplot['dstr_'+sji_wins[0]].sindex[0], storeplot['dstr_'+sji_wins[0]].sdata[*,*,0], sji_save
        map_save = zmap[0]

        index_sji0 = storeplot['dstr_'+sji_wins[0]].sindex[0]

        xfac = zmap[0].dx/sji_save[0].dx
        yfac = zmap[0].dy/sji_save[0].dy


        ;===========================================================

        for rr = 0,nraster_scans-1 do begin    ;start raster loop here <<<
        print,'searching SG scan ',rr

        read_iris_l2, sg_files[rr], index_sg
        d = iris_obj(sg_files[rr])
        xx = d->getxpos()
        sg_slit[rr,*] = xx

         for tiris = 0,rlen-1 do begin
            ;go through every time step in one raster, or entire time range for SNS
            ;then check the slit position against moss event count

            ;sg index from start
            tfromstart = rlen*rr + tiris

            ;NEAREST 193 TIME TO IRIS SG EXPOSURE
            tm = min( abs( anytim(time193)-anytim(index_sg[tiris].date_obs)) ,taia)

            ;NEAREST SJI TO SG EXPOSURE
            tn = min( abs(anytim(sji_timeline)-anytim(index_sg[tiris].date_obs)) ,tline)


            ;do this horrible complicated thing to find the correct SJI cube reference index
            diffs = tline-key_sji
            tt = where(diffs gt 0)
            mint = min(diffs[tt],sjiloc)

            ;corresponding time index on selected SJI cube
            if sjiloc eq 0 then sline = tline
            if sjiloc ge 1 then sline = tline-key_sji[sjiloc]-1

            index_sji = storeplot['dstr_'+sji_wins[sjiloc]].sindex

            ;SJI X SCALE and SLIT POSITION ON SJI IMAGE
            xdim_pix = index_sji[sline].naxis1
            xcen = index_sji[sline].xcen
            xfov = index_sji[sline].fovx
            xscale = index_sji[sline].cdelt1
            xslitpix = index_sji[sline].sltpx1ix    ;proper slit pixel on SJI
            irisxrange = findgen(xdim_pix)*xscale+(xcen-(xfov/2.))

            ;store the event sji x position
            sji_slit = irisxrange[xslitpix]

            ;create x and y solar position axes (AIA)
            map_fov, zmap[taia], xpos, ypos

            ;find pixels some +- of the slit
            ;left and right limits in solar xy
            leftedge = sg_slit[rr,tiris]-slit_offset
            rightedge = sg_slit[rr,tiris]+slit_offset

            ;top and bottom edges of the slit on the AIA window (solar xy)
            topedge = index_sg[tiris].ycen+(index_sg[tiris].fovy/2.)
            bottomedge = index_sg[tiris].ycen-(index_sg[tiris].fovy/2.)

            ;slice of this on the AIA pixel scale
            ww = where(xpos ge leftedge and xpos le rightedge)
            yy = where(ypos ge bottomedge and ypos le topedge)
            yyo = where(ypos lt bottomedge or ypos gt topedge)

            zimage = reform(zcube[*,*,taia])

            ;events within slit window
            zslice = zimage[ww,*]
            zslice = zslice[*,yy]

            ;Full AIA size image of events in slit cutout
            zcut = zimage*0
            zcut[ww,*] = zimage[ww,*]
            zcut[*,yyo] = 0

            ;total event pixel count per frame in the slit window
            etotal = total(zslice)

            ;GO FIND IRIS SG yposition IF event number is 1 or more
            if etotal ge 1 then begin

                ;count exposures with any event
                expcount = expcount + 1

                print,'TIRIS ',tiris
                print,'TAIA ',taia
                print,'TSJI ',sline

                ;here - use the proper SJI cube

                ;make SJI map with correct SJI
                ;this can be slow when there are a lot of events to create the map for
                ;ideally this should be done outside of the loop but it needs to select the correct data cube...

                ;tempdata = storeplot['dstr_'+sji_wins[sjiloc]].sdata
                index2map, index_sji[sline], storeplot['dstr_'+sji_wins[sjiloc]].sdata[*,*,sline], sj_frame

                sj_data = sj_frame.data
                zdata = zmap[taia].data

                ;create x and y solar position axes (IRIS)
                map_fov, sj_frame, ixpos, iypos

                ;find the edges of the SJI window on AIA
                mm1 = min(abs(xpos-ixpos[0]),mleft)
                mm2 = min(abs(ypos-iypos[0]),mbottom)
                imleft = mleft*xfac
                imbottom = mbottom*yfac

                ;AIA event map rescaled to SJI pixel size
                jsize = size(sj_data)
                sxl = jsize[1]
                syl = jsize[2]

                zdata_resize = congrid(zcut,xl*xfac, yl*yfac)
                zdata_cut = zdata_resize[imleft:(sxl-1)+imleft,imbottom:(syl-1)+imbottom]

                ;check only where SJI has data - removes events falling outside SJI data area
                poscheck = sj_data gt -200
                zdata_cut = zdata_cut AND poscheck

                eventareas = (zdata_cut ge 1)   ;make sure this is a binary map

                ;label_region will not create an event area for detections on the edge of the window
                ;skip these as there is unlikely good IRIS data at the very edge?

                ;now split multiple detections along the slit into induvidual events
                ;make a binary mask with uniquely numbered regions
                lareas = label_region(eventareas)

                ;EVENT LOOP
                ;IF THERE ARE EVENTS ON THE SLIT GO INTO EVENT BY EVENT LOOP
                if max(lareas) ge 1 then begin

                print,'event map has ',max(lareas),' events'
                     for a = 1,max(lareas) do begin
                        print,'event ',a

                        aia_step = [aia_step, taia]

                        sg_time = [sg_time, index_sg[tiris].date_obs]
                        aia_time = [aia_time, time193[taia]]

                        sji_time = [sji_time, storeplot['dstr_'+sji_wins[sjiloc]].sindex[sline].date_obs]

                        ;make a map for each counted event

                        eventarea = lareas eq a
                        single_map = sj_frame
                        single_map.data = eventarea

                        ;compress into slit vector and get iris y slit positions
                        eventareay = total(eventarea,1)
                        wwevent = where(eventareay ge 1)

                        ;event y centre
                        eventy_index = mean(wwevent)

                        ;SAVE X AND Y EVENT POSITIONS
                        iris_sg_xpos = [iris_sg_xpos, sg_slit[rr,tiris]]
                        iris_sg_xpixel = [iris_sg_xpixel, tiris]   ;this will need to be the raster step number

                        iris_sji_xpos = [iris_sji_xpos, sji_slit]  ;check this one
                        iris_sji_xpixel = [iris_sji_xpixel, xslitpix]

                        iris_sg_ypos = [iris_sg_ypos, iypos[eventy_index]]
                        iris_sg_ypixel = [iris_sg_ypixel, eventy_index]

                        iris_sg_scan = [iris_sg_scan, rr]
                        iris_sji_wave = [iris_sji_wave, sji_wins[sjiloc]]
                        tsji_closest = [tsji_closest, sjiloc]
                        sji_time_index = [sji_time_index, sline]

                        if keyword_set(savesample) then begin
                            ;saving the maps
                            sji_save = [sji_save, sj_frame]
                            map_save = [map_save, zmap[taia]]
                        endif ;SAVESAMPLE

                        if keyword_set(plotevent) then begin
                            aia_lct,wave='193',/load
                            plot_map,sj_frame,/log

                            oplot,[leftedge,leftedge],[iypos[0],iypos[-1]]
                            oplot,[rightedge,rightedge],[iypos[0],iypos[-1]]

                            oplot,[sji_slit],[iypos[eventy_index]],psym=2,thick=2
                            oplot,[sg_slit[rr,tiris]],[iypos[eventy_index]],psym=2,color=120,thick=2

                            loadct,1
                            plot_map,zmap[taia],/cont,/over,color=130
                            loadct,0
                            pause
                        endif ;PLOTEVENT

                    endfor ;A
                endif  ;extra event check

            endif ;EVENT TOTAL GE 1
        endfor ;IRIS cube loop

    endfor ;raster file loop


        ;clean initial 0s from counting vectors
        if n_elements(sg_time) gt 1 then begin
            sg_time = sg_time[1:*]
            sji_time = sji_time[1:*]
            aia_step = aia_step[1:*]
            aia_time = aia_time[1:*]
            iris_sg_xpos = iris_sg_xpos[1:*]
            iris_sji_xpos = iris_sji_xpos[1:*]
            iris_sg_ypos = iris_sg_ypos[1:*]
            iris_sg_xpixel = iris_sg_xpixel[1:*]
            iris_sji_xpixel = iris_sji_xpixel[1:*]
            iris_sg_ypixel = iris_sg_ypixel[1:*]

            iris_sg_scan = iris_sg_scan[1:*]
            iris_sji_wave = iris_sji_wave[1:*]
            tsji_closest = tsji_closest[1:*]
            sji_time_index = sji_time_index[1:*]


            if keyword_set(savesample) then sji_save = sji_save[1:*]
            if keyword_set(savesample) then map_save = map_save[1:*]
        endif

        ;EVENT SEARCH IS NOW COMPLETE (still working inside HCR loop)


        ;summary of exposures with events in current cube
        perc = (float(expcount)/float(tot_exposures))*100.
        print,perc,'% of IRIS exposures have AIA event detections'


        if ~keyword_set(savesample) then begin
            eout = {sg_time:sg_time, closest_sji2event:tsji_closest, sji_wave:iris_sji_wave,  $
                    sji_time_index:sji_time_index, sji_time:sji_time, aia_step:aia_step, aia_time:aia_time, $
                    iris_sg_xpos:iris_sg_xpos, iris_sji_xpos:iris_sji_xpos, iris_sg_ypos:iris_sg_ypos, $
                    iris_sg_xpixel:iris_sg_xpixel, iris_sji_xpixel:iris_sji_xpixel, iris_sg_ypixel:iris_sg_ypixel, $
                    iris_sg_scan:iris_sg_scan, $
                    header:hcr[i], id:hname, sji_files:storeplot['sji_files'], aiafiles:aiafiles, sg_files:sg_files, $
                    xshift_iris2aia:xcoeff, yshift_iris2aia:ycoeff, obs:obstype, scan_width:rlen, scan_repeats:nraster_scans}
            save, eout, filename=outdir+'/events_hcr_'+hname+'.sav'
        endif


        if keyword_set(savesample) then begin
            eout = {sg_time:sg_time, closest_sji2event:tsji_closest, sji_wave:iris_sji_wave, $
                    aia_step:aia_step, aia_time:aia_time, sji_time_index:sji_time_index, $
                    iris_sg_xpos:iris_sg_xpos, iris_sji_xpos:iris_sji_xpos, iris_sg_ypos:iris_sg_ypos, $
                    iris_sg_xpixel:iris_sg_xpixel, iris_sji_xpixel:iris_sji_xpixel, iris_sg_ypixel:iris_sg_ypixel, $
                    iris_sg_scan:iris_sg_scan, $
                    header:hcr[i], id:hname, sji_files:storeplot['sji_files'], aiafiles:aiafiles, sg_files:sg_files, $
                    event_map:map_save, sji_map:sji_save, $
                    xshift_iris2aia:xcoeff, yshift_iris2aia:ycoeff, obs:obstype, scan_width:rlen, scan_repeats:nraster_scans}
            save, eout, filename=outdir+'/events_hcr_'+hname+'.sav'
        endif

        ;remove store of SJI maps before the next HCR file
        delvar, storeplot

    endif  ;data check
endfor ; MAIN HCR ENTRY LOOP
;=========================================


return
end
