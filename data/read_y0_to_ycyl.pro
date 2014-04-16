pro read_y0_to_ycyl,y0_to_ycyl
file = '~/SPT/arnaud/y0_to_Y500_arnaud.sav' 
restore,file,/ver 

y0_to_ycyl = {y0_to_ycyl,R500:dblarr(146),y0_to_y500cyl:dblarr(146)}
y0_to_ycyl.R500 = R500_ARNAUD 
y0_to_ycyl.y0_to_y500cyl = Y0_TO_Y500_CYL 


end
