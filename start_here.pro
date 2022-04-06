;here is the standard method to run MOSSVAR.PRO

obsdir = '/irisa/data/level2/2014/08/06/20140806_013430_3800011454/'

mossvar,obsdir,'',mash,4,4,/read,/cube_data
print,mash['status']

;---creates some masks to narrow down the event areas
mossvar,obsdir,'',mash,4,4,/moss,/cnet,/fexviii,/loop,/filter1700

;---locates the actual brightenings and makes a movie in your local directory under $pwd/movies
mossvar,obsdir,'',mash,4,4,/variability,/cleanflare,/local_movie

;---show all of the detections
plot_image, total(mash['zmatch no loop'],3)


;some notes
;you can check the contents of output by doing IDL> mash.keys()
;access any part of this by doing IDL> help,/str,mash['variable here']

;'zmatch' contains the detection cube synced to the 193 channel
;'zmatch no loop' is the same but without the hot loop restriction


end
