pro read_params,params
Nparam = 23L
param = dblarr(Nparam)

I_SZ_A = 0
I_SZ_B = 1
I_SZ_C = 2
I_SZ_D = 3

I_S_A = 4
I_S_B = 5
I_S_C = 6
I_S_D = 7

I_S_D_0 = 17
I_S_D_N = 18

I_X_A = 8
I_X_B = 9
I_X_C = 10
I_X_D = 11

I_H0 = 12
I_OM = 13
I_OL = 14
I_W0 = 15
I_WA = 16

I_lambda_A = 17
I_lambda_B = 18
I_lambda_C = 19
I_lambda_D = 20

                                                                                                                                  
I_Y500cyl_D = 21                                                                                                             
I_Nfilt     = 22                                                                                                                 


param[I_SZ_A]=4.62211D
param[I_SZ_B]=1.36578D
param[I_SZ_C]=.465867D
param[I_SZ_D]=.169876D

param[I_S_A]=1
param[I_S_B]=1
param[I_S_C]=1
param[I_S_D]=1

param[I_S_D_0]=1
param[I_S_D_N]=1

param[I_X_A]=1
param[I_X_B]=1
param[I_X_C]=1
param[I_X_D]=1

param[I_H0]= 0.7;.7162D  
param[I_OM]= 0.3;.255D   
param[I_OL]= 1.D - param[I_OM]     
param[I_W0]= -1.D    
param[I_WA]= 0.D     

param[I_lambda_A]=20.562211D
param[I_lambda_B]=1.116578D  
param[I_lambda_C]=1.D  
param[I_lambda_D]=0.169876D  

param[I_Y500cyl_D]=0.172694
param[I_Nfilt]=30

params = {N:Nparam, theta:param}

end
