function call_likelihood_library,Ncluster,params,mass_func,data,y0_to_Y500cyl
soname = '../trunk/Lambdacalib.so'
loglikelihood  = 0D
;print, struct[0].str
;help, Ncluster,mass_func,data,params 

;help,Ncluster,params.N,params,mass_func,data,y0_to_Y500cyl,loglikelihood
tmp=call_external(soname, 'call_likelihood_library', Ncluster,params.N,params,mass_func,data,y0_to_Y500cyl,loglikelihood)
;print,loglikelihood 
return, loglikelihood
end


pro test,lambda_par=lambda_par

journal,'my_journal.pro'
read_data,data
;data= data[246]
;data= data[2]

Ncluster = N_elements(data.lambda)
read_params,params
I_H0 = 12
Profiler, /SYSTEM & Profiler

if is_def(lambda_par) then begin
    params.theta[17] = lambda_par[0]
    params.theta[18] = lambda_par[1]
    params.theta[19] = lambda_par[2]
    params.theta[20] = lambda_par[3]
endif
                                                                                                                                  
file = '~/SPT/arnaud/y0_to_Y500_arnaud.sav' 
restore,file,/ver 

; y0_to_ycyl = {R500:dblarr(146),y0_to_y500cyl:dblarr(146)}
; y0_to_ycyl.R500 = R500_ARNAUD 
; y0_to_ycyl.y0_to_y500cyl = Y0_TO_Y500_CYL 

read_massfunction,mass_func=mass_func,h0 = params.theta[I_H0]
read_y0_to_ycyl,y0_to_ycyl

y = create_struct(['MF','y0_to_ycyl'],mass_func,y0_to_ycyl)

mass_func = {MF,inherits MASS_FUNC}
mass_func = y.MF

y0_to_ycyl = {y0_to_ycyl,inherits Y0_TO_YCYL}
y0_to_ycyl = y.y0_to_ycyl


read_massfunction,mass_func=mass_func,h0 = params.theta[I_H0]
lk = call_likelihood_library(Ncluster,params,mass_func,data,y0_to_ycyl)
print, 'LK: ',lk
 Profiler, /REPORT
journal
stop
end



pro debugt,lambda_par=lambda_par

read_data,data
;data[0].lambda=20.
;data= data[0:4]
;data= data[40]
data= data[0]
Ncluster = N_elements(data.lambda)
read_params,params
I_H0 = 12
if is_def(lambda_par) then begin
    params.theta[17] = lambda_par[0]
    params.theta[18] = lambda_par[1]
    params.theta[19] = lambda_par[2]
    params.theta[20] = lambda_par[3]
endif
                                                                                                                                  
file = '~/SPT/arnaud/y0_to_Y500_arnaud.sav' 
restore,file,/ver 

y0_to_ycyl = {R500:dblarr(146),y0_to_y500cyl:dblarr(146)}
y0_to_ycyl.R500 = R500_ARNAUD 
y0_to_ycyl.y0_to_y500cyl = Y0_TO_Y500_CYL 

read_massfunction,mass_func=mass_func,h0 = params.theta[I_H0]
lk = call_likelihood_library(Ncluster,params,mass_func,data,y0_to_ycyl)

stop
end



pro get_Yobs_dist,lambda_par=lambda_par

read_data,data
;data[0].lambda=20.
;data= data[0:4]
data= data[0]
Ncluster = N_elements(data.lambda)
read_params,params
I_H0 = 12
if is_def(lambda_par) then begin
    params.theta[17] = lambda_par[0]
    params.theta[18] = lambda_par[1]
    params.theta[19] = lambda_par[2]
    params.theta[20] = lambda_par[3]
endif

Ngrid=50L
                                                                                                                                  
file = '~/SPT/arnaud/y0_to_Y500_arnaud.sav' 
restore,file,/ver 

y0_to_ycyl = {R500:dblarr(146),y0_to_y500cyl:dblarr(146)}
y0_to_ycyl.R500 = R500_ARNAUD 
y0_to_ycyl.y0_to_y500cyl = Y0_TO_Y500_CYL 

read_massfunction,mass_func=mass_func,h0 = params.theta[I_H0]
lk = call_yobs_dist_library(Ncluster,Ngrid,params,mass_func,data,y0_to_ycyl)

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
