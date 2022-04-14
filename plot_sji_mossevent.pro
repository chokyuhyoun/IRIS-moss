PRO plot_sji_mossevent, eout

slit_offset = 1.0


for i = 0,n_elements(eout.sg_time)-1 do begin

    aia_lct,wave='193',/load
    plot_map,eout.sji_map[i],/log
    map_fov, eout.sji_map[i], ixpos, iypos

    leftedge = eout.iris_sg_xpos[i]-slit_offset
    rightedge = eout.iris_sg_xpos[i]+slit_offset

    oplot,[leftedge,leftedge],[iypos[0],iypos[-1]]
    oplot,[rightedge,rightedge],[iypos[0],iypos[-1]]

    oplot,[eout.iris_sg_xpos[i]],[eout.iris_sg_ypos[I]],psym=2,color=120,thick=2

    loadct,1, /sil
    plot_map,eout.event_map[i],/cont,/over,color=130
    loadct,0, /sil
    pause
endfor

return
end
