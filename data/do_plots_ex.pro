pro do_plots_ex
setcolors
rdfloat,'out.txt',id,m_deb,P_M_deb,y_obs,P_y_obs,y_pred,P_y_pred,p_y_interp,integrand

lg = greek('lambda',/app)
s5 = '!D500!N'
msun = '10!U14!N M!D'+sunsymbol()+'!N'

plot,m_deb,P_M_deb,xtit='M'+s5+'['+msun+']',ytit= 'P(M'+s5+'|'+lg+')'
stop

plot ,y_pred,P_y_pred,xtit='Y'+s5+'[arcmin!u2!n]',ytit='P(Y'+s5+'|'+lg+')',yrange= minmax(P_y_obs)*1.2,/ys
oplot,y_obs,P_y_obs,col=mycolor(4)

legend,['Y'+s5+'!d-predicted!N','Y'+s5+'!d-measured!N'],col=[mycolor(0),mycolor(4)],linestyle=[0,0],box=0,/right


end
