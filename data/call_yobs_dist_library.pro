pro get_Yobs_dist,lambda_par=lambda_par

Profiler, /SYSTEM & Profiler
read_data,data
;data[0].lambda=20.
;data= data[0:4]
;data= data[0]
Ncluster = N_elements(data.lambda)
read_params,params
I_H0 = 12
if is_def(lambda_par) then begin
    params.theta[17] = lambda_par[0]
    params.theta[18] = lambda_par[1]
    params.theta[19] = lambda_par[2]
    params.theta[20] = lambda_par[3]
endif

Ngrid=50
                                                                                                                                  
file = '~/SPT/arnaud/y0_to_Y500_arnaud.sav' 
restore,file,/ver 

y0_to_ycyl = {R500:dblarr(146),y0_to_y500cyl:dblarr(146)}
y0_to_ycyl.R500 = R500_ARNAUD 
y0_to_ycyl.y0_to_y500cyl = Y0_TO_Y500_CYL 

read_massfunction,mass_func=mass_func,h0 = params.theta[I_H0]
lk = call_yobs_dist_library(Ncluster,Ngrid,params,mass_func,data,y0_to_ycyl)
 Profiler, /REPORT

stop
end

function call_yobs_dist_library,Ncluster,Ngrid,params,mass_func,data,y0_to_Y500cyl

soname = '/home/moon/saro/LAMBDA_MASS_CALIBRATION/SB/DEBUG/Y500/trunk/Lambdacalib.so'
Yobs_dist = DBLARR(Ncluster,Ngrid)
;print, struct[0].str
;help, Ncluster,mass_func,data,params 


;help,Ncluster,params.N,params,mass_func,data,y0_to_Y500cyl,loglikelihood
tmp=call_external(soname, 'call_Obtain_Yobs_dist',Ncluster,params.N,Ngrid,params,mass_func,data,y0_to_Y500cyl,Yobs_dist)
;print,loglikelihood 
return, loglikelihood
end
