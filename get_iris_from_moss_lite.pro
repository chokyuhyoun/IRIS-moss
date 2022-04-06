PRO get_iris_from_moss_lite,hcr,outdir,eout,mash_out,checkplot=checkplot,savesample=savesample,movie=movie,lightcurve=lightcurve

;returns a list of frames where moss events are found near the IRIS slit


;INPUT
;takes the HCR structure from an iris data search
;can be used with multiple HCR results or a single HCR structure

;OUTPUT
;structure containing time indices of the AIA and SJI cubes where events are found
;AIA_EVENT_INDICES
;IRIS_EVENT_INDICES
;HEADER - the IRIS hcr header used to process this data
;LEFTLIMIT, RIGHTLIMIT, XCEN_SJI - extent of the search  window per frame
;SJI_MAP - Slitjaw map for each dectection
;EVENT_MAP - event detections per frame from the AIA MOSSVAR routine

;EXAMPLES
;to return just the mossvar result

;UPDATES
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
slit_offset = 1.5

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

    ;setup output structure

    readme = 'moss events detected in IRIS'
    output = orderedhash('what is this',readme)

    ;create moss events zmatch cube

    ;>>> for later, turn the following into its own function and allow for local director or /irisa choice

    mossvar,udir[i],'',mash,4,4,/read,/cube_data
    if mash['status'] eq 'data ok' then begin

        mossvar,udir[i],'',mash,4,4,/moss,/cnet,/fexviii,/loop,/filter1700
        if keyword_set(movie) then mossvar,udir[i],'',mash,4,4,/variability,/cleanflare,/local_movie
        if ~keyword_set(movie) then mossvar,udir[i],'',mash,4,4,/variability,/cleanflare

        ;cube of counted moss events
        zcube = mash['zmatch no loop']

        index193 = mash['aia_193'].index
        data193 = mash['aia_193'].data

        index2map, index193, data193 ,map193

        ;AIA ZERO CROSSING CUBE
        zmap = map193
        zmap.data = zcube
        print,'created zcube map from MOSSVAR'

        dsize = size(zmap.data)
        xl = dsize[1]
        yl = dsize[2]
        tl = dsize[3]

        ;get SJI index
        irisfiles = find_files('*SJI*',udir[i])
        output['irisfiles'] = irisfiles

        ;READ IRIS CUBE
        read_iris_l2, irisfiles[0], index_sji, data_sji

        ;create dummy maps to make a save array
        index2map,index_sji[0],data_sji[*,*,0],sji_save
        event_map_save = zmap[0]

        xfac = zmap[0].dx/sji_save[0].dx
        yfac = zmap[0].dy/sji_save[0].dy

        jsize = size(data_sji)
        sxl = jsize[1]
        syl = jsize[2]
        stl = jsize[3]

        ;OUTPUTS
        ;define vector to store event positions
        event_iris_xpos = 0
        event_iris_ypos = 0
        event_irisx_pixel = 0
        event_irisy_pixel = 0

        event_iris_step = 0
        event_aia_step = 0
        event_iris_time = 'empty'
        event_aia_time = 'empty'
        ;===========================================================

        ;COALIGNMENT check here - provides correlation coefficients
        align_iris_aia, index_sji, data_sji, mash['aia_1600'].index, mash['aia_1600'].data,xcoeff,ycoeff

        ;modify AIA EVENT MAP coordinates
        ;caution if comparing to original AIA cube - this is a destructive change
        zmap.xc = zmap.xc - xcoeff
        zmap.yc = zmap.yc - ycoeff

        slen = n_elements(index_sji.date_obs)

        print,'read iris SJI file'

        ;empty array for IRIS slit position
        sji_slit = fltarr(n_elements(index_sji.xcen))

        ;for slit position per event
        sji_slit_xpos = 0

        time193 = mash['aia_193'].index.date_obs

        ;to collect total event number per hcr entry
        zf = fltarr(tl)

        ;count events per aia frame
        zcount = fltarr(tl)

         for tiris = 0,slen-1 do begin
            ;go through every time step in one raster, check the slit position against moss event count

            ;SJI X SCALE and IRIS SLIT POSITION
            xdim_pix = index_sji[tiris].naxis1
            xcen = index_sji[tiris].xcen
            xfov = index_sji[tiris].fovx
            xscale = index_sji[tiris].cdelt1
            xslitpix = index_sji[tiris].sltpx1ix    ;proper slit pixel
            irisxrange = findgen(xdim_pix)*xscale+(xcen-(xfov/2.))

            sji_slit[tiris] = irisxrange[xslitpix]

            ;NEAREST 193 TIME TO IRIS
            tm = min( abs( anytim(time193)-anytim(index_sji[tiris].date_obs)) ,taia)

            ;create x and y solar position axes (AIA)
            map_fov, zmap[taia], xpos, ypos

            ;find pixels some +- of the slit
            ;left and right limits in solar xy
            leftedge = sji_slit[tiris]-slit_offset
            rightedge = sji_slit[tiris]+slit_offset

            ;top and bottom edges of the slit on the AIA window (solar xy)
            topedge = index_sji[tiris].ycen+(index_sji[tiris].fovy/2.)
            bottomedge = index_sji[tiris].ycen-(index_sji[tiris].fovy/2.)

            ;slice of this on the AIA pixel scale
            ww = where(xpos ge leftedge and xpos le rightedge)
            yy = where(ypos ge bottomedge and ypos le topedge)
            yyo = where(ypos lt bottomedge or ypos gt topedge)

            zimage = reform(zcube[*,*,taia])

            ;events within slit window
            zslice = zimage[ww,*]
            zslice = zslice[*,yy]

            ;AIA image sized cutout of events in slit
            zcut = zimage*0
            zcut[ww,*] = zimage[ww,*]
            zcut[*,yyo] = 0

            ;total event count per frame in the slit window
            etotal = total(zslice)
            zcount[taia] = etotal
            zf[taia] = total(total(zcube[*,*,taia],1),1)

            ;IRIS yposition search IF event number is 1 or more
            if etotal ge 1 then begin

                index2map,index_sji[tiris],data_sji[*,*,tiris],sj_frame

                sj_data = sj_frame.data
                zdata = zmap[taia].data

                ;use if enlarging the event map to account for correllation error
                ;switch off for now... correalation is now provided earlier
                ;k=[[0,1,0],[1,1,1],[0,1,0]]
                ;zdata = dilate(zdata,k)

                ;create x and y solar position axes (IRIS)
                map_fov, sj_frame, ixpos, iypos

                ;find the edges of the SJI window on AIA
                mm1 = min(abs(xpos-ixpos[0]),mleft)
                mm2 = min(abs(ypos-iypos[0]),mbottom)
                imleft = mleft*xfac
                imbottom = mbottom*yfac

                ;AIA event map rescaled to SJI pixel size
                zdata_resize = congrid(zcut,xl*xfac, yl*yfac)
                zdata_cut = zdata_resize[imleft:(sxl-1)+imleft,imbottom:(syl-1)+imbottom]

                eventareas = (zdata_cut ge 1)   ;make sure this is a binary map
                scut = eventareas * sj_data

                ;label_region will not create an area for detections on the edge of the window
                ;we can skip these as there is unlikely IRIS data at the very edge
                ;or later extend the window?

                ;split multiple detections into induvidual events
                ;make a binary mask with uniquely numbered regions
                lareas = label_region(eventareas)

                ;EVENT LOOP
                ;IF THERE ARE EVENTS ON THE SLIT GO INTO EVENT BY EVENT LOOP
                if max(lareas) ge 1 then begin
                print,tiris
                print,'event map has ',max(lareas),' events'
                    for a = 1,max(lareas) do begin

                        event_iris_step = [event_iris_step, tiris]
                        event_aia_step = [event_aia_step, taia]

                        event_iris_time = [event_iris_time, index_sji[tiris].date_obs]
                        event_aia_time = [event_aia_time, time193[taia]]

                        ;make a map for each counted event

                        eventarea = lareas eq a
                        single_event_map = sj_frame
                        single_event_map.data = eventarea

                        ;compress into slit vector and get iris y slit positions
                        eventareay = total(eventarea,1)
                        wwevent = where(eventareay ge 1)

                        ;event y centre
                        irisy_index = mean(wwevent)

                        ;SAVE X AND Y EVENT POSITIONS
                        event_iris_xpos = [event_iris_xpos, sji_slit[tiris]]
                        event_iris_ypos = [event_iris_ypos, iypos[irisy_index]]

                        event_irisx_pixel = [event_irisx_pixel, index_sji[tiris].sltpx1ix]
                        event_irisy_pixel = [event_irisy_pixel, irisy_index]

                        if keyword_set(savesample) then begin

                            ;for saving the maps
                            sji_save = [sji_save, sj_frame]
                            event_map_save = [event_map_save, zmap[taia]]

                            if keyword_set(checkplot) then begin
                                aia_lct,wave='193',/load
                                plot_map,sj_frame,dmax=1000

                                oplot,[leftedge,leftedge],[iypos[0],iypos[-1]]
                                oplot,[rightedge,rightedge],[iypos[0],iypos[-1]]

                                oplot,[sji_slit[tiris]],[iypos[irisy_index]],psym=2
                                loadct,1
                                plot_map,zmap[taia],/cont,/over,color=130
                                loadct,0
                                pause
                            endif ;CHECKPLOT
                        endif ;SAVESAMPLE

                    endfor ;A
                endif  ;extra event check

            event_iris_step = event_iris_step[1:*]
            event_iris_time = event_iris_time[1:*]
            event_aia_step = event_aia_step[1:*]
            event_aia_time = event_aia_time[1:*]
            event_iris_xpos = event_iris_xpos[1:*]
            event_iris_ypos = event_iris_ypos[1:*]
            event_irisx_pixel = event_irisx_pixel[1:*]
            event_irisy_pixel = event_irisy_pixel[1:*]

            sji_map = sji_map[1:*]
            event_map_save = event_map_save[1,*]

            endif ;EVENT TOTAL GE 1

        endfor ;IRIS cube loop



        ;some extra stuff - check this later
        ;check for frames with anomalous event spikes and remove those from the count
        znonzero = where(zf gt 0.)
        gmean = mean(zf(znonzero))
        flag = where(zf gt 10.*gmean, nflag)
        if nflag gt 0 then zcount[flag] = 0


        tmax = max(zcount,mtime)

        ;count frames with event numbers found in the slit window above a threshold
        eventcut = 3.
        event_aia_indices = where(zcount gt eventcut,nframes)

        ;summary of event number in current cube
        perc = (float(nframes)/float(tl))*100.
        print,perc,'% of frames with detections'


        if ~keyword_set(savesample) then begin
            eout = {event_iris_step:event_iris_step, event_iris_time:event_iris_time,$
                    event_aia_step:event_aia_step, event_aia_time:event_aia_time,$
                    event_iris_xpos:event_iris_xpos, event_iris_ypos:event_iris_ypos, $
                    event_irisx_pixel:event_irisx_pixel, event_irisy_pixel:event_irisy_pixe, $
                    header:hcr[i], event_id:hname, irisfiles:output['irisfiles'], aiafiles:aiafiles}
            save, eout, filename=outdir+'/events_hcr_'+hname+'.sav'
        endif


        if keyword_set(savesample) then begin
            eout = {event_iris_step:event_iris_step, event_iris_time:event_iris_time,$
                    event_aia_step:event_aia_step, event_aia_time:event_aia_time,$
                    event_iris_xpos:event_iris_xpos, event_iris_ypos:event_iris_ypos, $
                    event_irisx_pixel:event_irisx_pixel, event_irisy_pixel:event_irisy_pixel, $
                    header:hcr[i], event_id:hname, irisfiles:output['irisfiles'], aiafiles:aiafiles, $
                    event_map:event_map_save, sji_map:sji_save}
            save, eout, filename=outdir+'/events_hcr_'+hname+'.sav'
        endif


    endif  ;data check
endfor ; MAIN HCR ENTRY LOOP
;=========================================


return
end
