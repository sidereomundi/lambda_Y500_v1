pro read_data,data,zmax=zmax
  
data_struct={data,clus_id:0L,$    ; /** ID of the cluster */                                       
             field_id:0L,$   ; /** ID of field */                                            
             redshift:0D,$   ; /** redshift of the cluster */                                   
             lambda:0D,$     ; /** lambda observable */                                         
             dlambda:0D,$  ; /** lambda observable error */                                         
             y0: dblarr(30),$ ;                                                                    
             dy0: dblarr(30),$ ;                                                                    
             filter: dblarr(30)$ ;                                                                    
             }

datafile = '/home/moon/saro/LAMBDA_MASS_CALIBRATION/SB/DEBUG/Y500/data/data_for_likelihood_lib.sav'
restore,datafile

field_id_match = unique(fields,/sort)
Ncluster = n_elements(fields)
field_id = LONARR(Ncluster)

for i=0,Ncluster-1 do begin
    tmp = where(field_id_match eq fields[i])
    field_id[i] = tmp[0]
endfor

data = replicate(data_struct,Ncluster)

data.clus_id = cat_subsample.MEM_MATCH_ID
data.field_id = field_id
data.redshift = cat_subsample.ZRED
data.lambda = cat_subsample.LAMBDA_CHISQ
data.dlambda = cat_subsample.LAMBDA_CHISQ_E
filter = dindgen(30)*0.25+0.5

for i=0,Ncluster -1 do begin
    data[i].y0  = y0[i,*]
    data[i].dy0 = dy0[i,*]
    data[i].filter = filter
endfor

ra = cat_subsample.ra
dec = cat_subsample.dec
ir = where((ra ge 50 and ra le 100 ) or (ra ge 300 and dec le -50),nir)
data = data[ir]

i = where(finite(data.dy0[0]) eq 1,ni)
data = data[i]
if is_def(zmax) then data = data(where(data.redshift le zmax) )

end




pro calculate_p_Y500_M500,omegal=omegal,omegam=omegam,h0=h0,DY500=Dy500

if not is_def(omegam) then omegam= 0.292000D
if not is_def(omegal) then Omegal = 0.708000D
if not is_def(h0) then h0 = 68.6D/100
if not is_def(DY500) then DY500 = 0.075*alog(10)

datafile = '/home/moon/saro/LAMBDA_MASS_CALIBRATION/SB/DEBUG/Y500/data/data_for_likelihood_lib.sav'                          
restore,datafile,/ver
n = 1.e3
nyarr = 1.e3
marr = (dindgen(n)/(n-1)*(1.e16-1.e13)) +1.e13
z= cat_subsample.zred 
n_cluster= n_elements(z)
read_y0_to_ycyl,y0_to_ycyl           

Mpc_over_arcmin = zang(1.e3,z,H0 =h0*100, Omega_m =omegam, Lambda0=omegal)/60.

p_M= dblarr(n_cluster,n)

for i=0, n_cluster-1 do begin
    percentuale,i,n_cluster-1,pp
    r500_arcmin = zang(mcrit_to_rcrit(marr,z=z[i], h0=h0, omegal=omegal,omegam=omegam),$
                       z[i],H0 =h0*100, Omega_m =omegam, Lambda0=omegal,silent=1)/60.
    y500factor =  INTERPOL( y0_to_ycyl.Y0_TO_Y500CYL, y0_to_ycyl.R500, r500_arcmin) 
    y0interpol = INTERPOL(  y0[i,*],filter, r500_arcmin) 
    dy0interpol = INTERPOL(dy0[i,*],filter, r500_arcmin) 
    
    y500_measured  = y500factor*y0interpol
    dy500_measured = y500factor*dy0interpol
    
    y500_cyl_predict = m500c_to_y500cyl(marr,z[i],h0=h0,omegam=omegam,omegal=omegal)*Mpc_over_arcmin[i]^2
    

    for j=0,n-1 do begin
    max_y = max(y500_measured[j] + 7*dy500_measured[j]) > y500_cyl_predict[j]*(1.+Dy500)*7
    min_y = min(y500_measured[j] - 7*dy500_measured[j])
    y500_arr = findgen(nyarr)/(nyarr -1)*(max_y - min_y) + min_y
    p_m(i,j)=integral(y500_arr,$
                      normal(y500_arr,y500_measured[j],dy500_measured[j])*$
                      lognormal(y500_arr,y500_cyl_predict[j],Dy500))
    endfor
endfor

m500 = marr
p_y500_given_m500 = p_m
save,filename ='p_y500_given_m500.sav',m500,p_y500_given_m500

stop



end
