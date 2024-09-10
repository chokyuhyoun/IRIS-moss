PRO get_iris_from_mossevent,hcr,outdir,eout,mash_out,checkplot=checkplot,savesample=savesample,movie=movie,lightcurve=lightcurve

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

nresults = n_elements(hcr)

slit_x = fltarr(nresults)
date = strarr(nresults)
obs_short = strarr(nresults)
udir = strarr(nresults)
number_mossframes = fltarr(nresults)

;how far to search around slit
slit_offset = 2.5

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

        zmap = map193
        zmap.data = zcube
        print,'created zcube map'

        dsize = size(zmap.data)
        xl = dsize[1]
        yl = dsize[2]
        tl = dsize[3]

        zcount = fltarr(tl)

        ;get SJI index
        irisfiles = find_files('*SJI*',udir[i])
            output['irisfiles'] = irisfiles

        if keyword_set(savesample) or keyword_set(lightcurve) then begin
            ;saving sji samples:
            read_iris_l2, irisfiles[0], index_sji, data_sji

            endif else begin

            ;only working with event counts:
            read_iris_l2, irisfiles[0], index_sji
        endelse

        slen = n_elements(index_sji.date_obs)

        print,'read iris SJI file'

        ;time dependent IRIS slit position array for cube
        sji_slit = fltarr(n_elements(index_sji.xcen))

        ;for slit position per event
        sji_slit_xpos = 0

        time193 = mash['aia_193'].index.date_obs

        ;x and y ranges
        xran = findgen(xl)*zmap[0].dx
        yran = findgen(yl)*zmap[0].dy

        ;to collect total event number per frame
        zf = fltarr(tl)

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

            ;nearest aia 193 time to iris
            tm = min( abs( anytim(time193)-anytim(index_sji[tiris].date_obs)) ,taia)

            ;create x and y solar position axes (AIA)
            xleft = zmap[taia].xc - ((xl/2.)*zmap[taia].dx)
            xpos = xleft + xran

            yleft = zmap[taia].yc - ((yl/2.)*zmap[taia].dy)
            ypos = yleft + yran


            ;find pixels some +- of the slit
            leftedge = sji_slit[tiris]-slit_offset
            rightedge = sji_slit[tiris]+slit_offset

            topedge = index_sji[tiris].ycen+(index_sji[tiris].fovy/2.)
            bottomedge = index_sji[tiris].ycen-(index_sji[tiris].fovy/2.)

            ;slice on the AIA pixel scale
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
            zcount[taia] = total(zslice)
            zf[taia] = total(total(zcube[*,*,taia],1),1)

        endfor

        ;check for frames with anomalous event spikes and remove those from the count
        znonzero = where(zf gt 0.)
        gmean = mean(zf(znonzero))
        flag = where(zf gt 10.*gmean, nflag)
        if nflag gt 0 then zcount[flag] = 0


        ;make an example plot for the most significant event

        tmax = max(zcount,mtime)

        ;count frames with event numbers found in the slit window above a threshold
        eventcut = 3.
        aia_event_indices = where(zcount gt eventcut,nframes)

        perc = (float(nframes)/float(tl))*100.
        print,perc,'% of frames with detections'

        ;=========================================

        ;find corresponding SJI frames to AIA event detections
        iris_event_indices = 0

        ;in case of no events default output is the first slit position of the data
        leftlimit = sji_slit[0]-slit_offset
        rightlimit = sji_slit[0]+slit_offset
        sji_slit_xpos = sji_slit[0]

        output['leftlimit'] = leftlimit
        output['rightlimit'] = rightlimit
        output['sji_slit_xpos'] = sji_slit_xpos

        output['aia_event_indices'] = aia_event_indices
        output['iris_event_indices'] = iris_event_indices

        if nframes gt 0 then begin
            iris_event_indices = intarr(nframes)
            leftlimit = fltarr(nframes)
            rightlimit = fltarr(nframes)
            zsave = 'empty'
            map_sji = 'empty'

            for j=0,nframes-1 do begin
                mtime = aia_event_indices[j]
                tm = min( abs( anytim(index_sji.date_obs)-anytim(time193[mtime])) ,tiris)
                iris_event_indices[j] = tiris

                leftlimit[j] = sji_slit[tiris]-slit_offset
                rightlimit[j] = sji_slit[tiris]+slit_offset

                topedge = index_sji[tiris].ycen+(index_sji[tiris].fovy/2.)
                bottomedge = index_sji[tiris].ycen-(index_sji[tiris].fovy/2.)

            endfor

            sji_slit_xpos = sji_slit[iris_event_indices]

            if keyword_set(savesample) then begin
                ;create a map structure for SJI frames where there is a detection
                index2map,index_sji[0],data_sji[*,*,0],map_sji

                zsave = zmap[0]

                for j=0,nframes-1 do begin
                    mtime = aia_event_indices[j]
                    tiris = iris_event_indices[j]

                    index2map,index_sji[tiris],data_sji[*,*,tiris],sj_frame
                    map_sji = [map_sji, sj_frame]
                    zsave = [zsave, zmap[mtime]]

                    if keyword_set(checkplot) then begin
                        aia_lct,wave='193',/load
                        plot_map,sj_frame,dmax=1000

                        oplot,[leftlimit[j],leftlimit[j]],[ypos[0],ypos[-1]]
                        oplot,[rightlimit[j],rightlimit[j]],[ypos[0],ypos[-1]]

                        loadct,1
                        plot_map,zmap[mtime],/cont,/over,color=130
                        pause
                    endif
                endfor

                map_sji = map_sji[1:*]
                zsave = zsave[1:*]

                output['sji_map'] = map_sji
                output['event_map'] = zsave
            endif

            ;======== check the lightcurve of brightest events ============

            if keyword_set(lightcurve) then begin

                ;create a map structure for SJI frames where there is a detection
                index2map,index_sji[0],data_sji[*,*,0],map_sji

                zsave = zmap[0]
                ;keep a map structure at the SJI scale for collecting induvidual events (may be more than one per time step)
                larea_map = map_sji
                larea_map.id = 'AIA EVENT MAP AT SJI PIXEL SCALE'

                ;somewhere to keep iris time indice for each event
                tiris_by_event = 0
                taia_by_event = 0

                jsize = size(map_sji[0].data)
                sxl = jsize[1]
                syl = jsize[2]

                ;x and y arrays (un-centered)
                ixran = findgen(sxl)*map_sji[0].dx
                iyran = findgen(syl)*map_sji[0].dy

                ;rescale factor for AIA to SJI
                xfac = zsave[0].dx/map_sji.dx
                yfac = zsave[0].dy/map_sji.dy


                for j=0,nframes-1 do begin

                    tiris = iris_event_indices[j]
                    taia = aia_event_indices[j]

                    index2map,index_sji[tiris],data_sji[*,*,tiris],sj_frame
                    map_sji = [map_sji, sj_frame]
                    zsave = [zsave, zmap[taia]]

                    sj_data = sj_frame.data
                    zdata = zmap[taia].data

                    ;use if enlarging the event map to account for correllation error
                    k=[[0,1,0],[1,1,1],[0,1,0]]
                    zdata = dilate(zdata,k)

                    ;create x and y solar position axes (IRIS)
                    ixleft = sj_frame.xc - ((sxl/2.)*sj_frame.dx)
                    ixpos = ixleft + ixran
                    iyleft = sj_frame.yc - ((syl/2.)*sj_frame.dy)
                    iypos = iyleft + iyran

                    ;create x and y solar position axes (AIA)
                    xleft = zmap[taia].xc - ((xl/2.)*zmap[taia].dx)
                    xpos = xleft + xran
                    yleft = zmap[taia].yc - ((yl/2.)*zmap[taia].dy)
                    ypos = yleft + yran

                    ;find pixels some +- of the slit
                    leftedge = sji_slit[tiris]-slit_offset
                    rightedge = sji_slit[tiris]+slit_offset

                    topedge = index_sji[tiris].ycen+(index_sji[tiris].fovy/2.)
                    bottomedge = index_sji[tiris].ycen-(index_sji[tiris].fovy/2.)

                    ;slice on the AIA pixel scale
                    ww = where(xpos ge leftedge and xpos le rightedge)
                    yy = where(ypos ge bottomedge and ypos le topedge)
                    yyo = where(ypos lt bottomedge or ypos gt topedge)

                    ;AIA image sized cutout of events in slit
                    zcut = zdata*0
                    zcut[ww,*] = zdata[ww,*]
                    zcut[*,yyo] = 0

                    ;find the edges of the SJI window on AIA
                    mm1 = min(abs(xpos-ixpos[0]),mleft)
                    mm2 = min(abs(ypos-iypos[0]),mbottom)
                    imleft = mleft*xfac
                    imbottom = mbottom*yfac

                    ;AIA event map rescaled to SJI pixel size
                    zdata_resize = congrid(zcut,xl*xfac, yl*yfac)
                    zdata_cut = zdata_resize[imleft:(sxl-1)+imleft,imbottom:(syl-1)+imbottom]

                    eventareas = (zdata_cut ge 1)
                    scut = eventareas * sj_data

                    ;make lightcurves for each event - sometimes multiple per image
                    ;label_region will not create an area for detections on the edge of the window
                    ;we can skip these as there is unlikely IRIS data at the very edge


                    lareas = label_region(eventareas)
                    if max(lareas) gt 0 then begin

                        tcurves = fltarr(20,max(lareas))

                        ;= ISOLATE EVENTS INTO SEPERATE MAPS
                        for a = 1,max(lareas) do begin
                            ;print,a
                            tiris_by_event = [tiris_by_event,tiris]
                            taia_by_event = [taia_by_event,taia]
                            ;make a map for each counted event
                            single_event_map = sj_frame
                            single_event_map.data = lareas eq a

                            larea_map = [larea_map, single_event_map]
                        endfor

                    endif else begin
                        tiris_by_event = [tiris_by_event,tiris]
                        taia_by_event = [taia_by_event,taia]

                        single_event_map = sj_frame
                        single_event_map.data = single_event_map.data * 0

                        larea_map = [larea_map, single_event_map]
                    endelse

                endfor  ;EVENT MAP LOOP END

                map_sji = map_sji[1:*]
                zsave = zsave[1:*]
                event_map_selected = larea_map[1:*]
                tiris_by_event = tiris_by_event[1:*]
                taia_by_event = taia_by_event[1:*]

                ncurves = n_elements(event_map_selected)
                tcurves = fltarr(20,ncurves)
                tdiff = fltarr(ncurves)

                for j = 0,ncurves-1 do begin
                    tiris = tiris_by_event[j]
                    print,tiris
                    tcurve = fltarr(20)
                    ;set so that area is centered at a=9 around the event indice

                    low = tiris-9
                    high = tiris+10

                    if slen-1-high lt 0 then high = slen-1
                    if low lt 0 then low = 0
                    print,tiris
                    print,low,'   ',high
                    count = -1
                    for tt = low,high do begin
                        count = count+1
                        tcurves[count,j] = total( (event_map_selected[j].data eq 1) * data_sji[*,*,tt])
                    endfor

                    nz = where(tcurves[*,j] gt 0)
                    tdiff[j] = max(tcurves[nz,j]) - min(tcurves[nz,j])

                endfor  ;ISOLATED EVENT LOOP

                mm = max(tdiff,maxd)
                print,mm,'   ',maxd

                wmp = where(iris_event_indices eq tiris_by_event[maxd])

                output['best_curve'] = tcurves[*,maxd]
                output['best_event_iris_index'] = tiris_by_event[maxd]
                output['tcurves'] = tcurves

                output['best_event_map'] = event_map_selected[maxd]
                output['best_sji_map'] = map_sji[wmp[0]]
                output['event_map_selected'] = event_map_selected

                output['aia_by_event'] = taia_by_event
                output['iris_by_event'] = tiris_by_event

            endif  ;END LIGHTCURVE

            ;==============================================================

            output['aia_event_indices'] = aia_event_indices
            output['iris_event_indices'] = iris_event_indices

            output['leftlimit'] = leftlimit
            output['rightlimit'] = rightlimit
            output['sji_slit_xpos'] = sji_slit_xpos


            if keyword_set(savesample) and ~keyword_set(lightcurve) then begin
                ;save structures
                eout = {aia_event_indices:output['aia_event_indices'], iris_event_indices:output['iris_event_indices'],$
                        leftlimit:output['leftlimit'], rightlimit:output['rightlimit'], xcen_sji:output['sji_slit_xpos'], header:hcr[i], event_id:hname, $
                        irisfiles:output['irisfiles'], aiafiles:aiafiles, event_map:output['event_map'], sji_map:output['sji_map']}

                save, eout, filename=outdir+'/events_hcr_'+hname+'.sav'
            endif

            if keyword_set(lightcurve) and ~keyword_set(savesample) then begin
                ;save structures
                eout = {aia_event_indices:output['aia_event_indices'], iris_event_indices:output['iris_event_indices'],$
                        leftlimit:output['leftlimit'], rightlimit:output['rightlimit'],xcen_sji:output['sji_slit_xpos'], header:hcr[i], event_id:hname, irisfiles:output['irisfiles'], $
                        aiafiles:aiafiles, best_event_map:output['best_event_map'], best_curve:output['best_curve'], best_iris_index:output['best_event_iris_index'], best_sji_map:output['best_sji_map'],$
                        aia_by_event:output['aia_by_event'], iris_by_event:output['iris_by_event']}
                save, eout, filename=outdir+'/events_hcr_'+hname+'.sav'
            endif

            if keyword_set(lightcurve) and keyword_set(savesample) then begin
                ;save structures
                eout = {aia_event_indices:output['aia_event_indices'], iris_event_indices:output['iris_event_indices'],$
                        leftlimit:output['leftlimit'], rightlimit:output['rightlimit'],xcen_sji:output['sji_slit_xpos'], header:hcr[i], event_id:hname, irisfiles:output['irisfiles'], $
                        aiafiles:aiafiles, event_map:output['event_map'], sji_map:output['sji_map'], $
                        best_event_map:output['best_event_map'], best_curve:output['best_curve'], best_iris_index:output['best_event_iris_index'], best_sji_map:output['best_sji_map'],$
                        aia_by_event:output['aia_by_event'], iris_by_event:output['iris_by_event']}
                save, eout, filename=outdir+'/events_hcr_'+hname+'.sav'
            endif

        endif else begin   ;end of nframes gt 0
        ;save structures
        eout = {aia_event_indices:0, iris_event_indices:output['iris_event_indices'],$
                leftlimit:output['leftlimit'], rightlimit:output['rightlimit'],xcen_sji:output['sji_slit_xpos'], header:hcr[i], event_id:hname, irisfiles:output['irisfiles']}
        save, eout, filename=outdir+'/events_hcr_'+hname+'.sav'
        endelse

    number_mossframes[i] = nframes


    endif ;END OF 'data ok'



;MAIN HEADER FILE LOOP
endfor

out = number_mossframes
mash_out = mash

return

end
