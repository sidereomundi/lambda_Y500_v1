pro do_plots,dirout=dirout,wait=wait

if not is_def(wait) then wait =0.
setcolors
basedir='/big/users/saro/LRG/SVA1/GOLD/CATALOGUES/MCMC/CH'
if not is_def(dirout) then begin
    dirout = 1
    while file_exists(basedir+strtrim(string(dirout),1)) eq 1 do dirout++
    dirout--
endif
dirout = basedir+strtrim(string(dirout),1)+'/'

print,dirout

restore,dirout+'mcmc_chains.sav',/ver

plots_mcmc,outputvalue=outputvalue,outputll=outputll,param=param,input=input

stop
nparam = n_elements(param.value)
bestf =  fltarr(nparam)
sig1 =  fltarr(nparam)
for i =0, nparam -1 do begin
    histogauss68,outputvalue[i,*],kj,xtitle=param.name[i]    
    print,param.name[i] + ' marginalized constraints: ',kj[1], '+-',kj[2]
    bestf[i] = kj[1]
    sig1[i] = kj[2]
    wait,wait
endfor



stop
end



pro plots_mcmc,outputvalue=outputvalue,outputll=outputll,param=param,input=input
setcolors

nparam = n_elements(param.name)
!p.multi=[0,nparam,nparam,0,0]
;erase & multiplot,/default 
erase & multiplot,[nparam,nparam]
!P.charsize=sqrt(nparam)


for i=0,nparam -1 do begin
    for j=0,nparam-1 do begin
        print,i,j
        xrange = minmax(outputvalue(j,*))
        yrange = minmax(outputvalue(i,*))
        if j eq 0 then ytit = param.name(i) else ytit = ''
        if i eq nparam-1 then xtit = param.name(j) else xtit = ''
        if i eq j then begin
            setcolors
            autohist, outputvalue(i,*),xx,yy,/noplot
            yy = smooth(yy/max(yy),2)*(yrange[1] - yrange[0]) + yrange[0]
            plot,xx,yy,xtit=xtit,ytit=ytit,xrange=xrange,yrange=yrange,col=mycolor(0)
        endif
        if j gt i then begin
            plot,[0],xstyle=4,ystyle=4,/nodata  
            setcolors
        endif
        if j lt i then begin
            setcolors
            do_map_weight,outputvalue[j,*],outputvalue[i,*],$
              outputll,100,100,minmax(outputvalue[j,*]),minmax(outputvalue[i,*]),map,xarr,yarr,levs=levs
            map = -map
            map(where(map eq 0)) = min(map(where(map gt 0)))
            map = smooth(alog(map),2)
            map -= min(map)
;            plot,[0],xtit=xtit,ytit=ytit,xrange=xrange,yrange=yrange,/xs,/ys,col=mycolor(0),/nodata
;            do_map,map,xarr,yarr,xtit=xtit,ytit=ytit,nobar=1,/rev
            CONTOUR, map/MAX(MAP),xarr,yarr, C_COLORS = [mycolor(35),mycolor(4)],levels=[0.05,0.32],xtit=xtit,ytit=ytit,/fill
            setcolors
            if is_def(input) then oplot,input[j:j],input[i:i],psym=2,thick=4,col=mycolor(0),symsize=4
        endif
        multiplot
    endfor
endfor
stop

multiplot,/default 

end
