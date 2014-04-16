; SHOW MASS FUNCTION
PRO read_massfunction,mass_func = mass_func,H0=H0
;mass_function=mass_function,mass_array=mass_array,z_array=z_array,H0=H0

  if not is_def(h0) then h0 = 0.7162

  nx = 500L                     ;,xarr = mass
  ny = 51L                      ;yarr = redshift
  NAME = '/home/moon/saro/LAMBDA_MASS_CALIBRATION/SB/data/massfunction.dat'
  MASS_FUNCTION = DBLARR(ny,nx)
  OPENR,1,NAME
  READF,1,MASS_FUNCTION
  CLOSE,1

  MASS_FUNCTION = TRanspose(MASS_FUNCTION)

  MASSMIN = 13.5                                                                                                                        
  MASSMAX = 16.                                                                                                                         
  ZMIN = 0.025                                                                                                                           
  ZMAX = 2.525                                                                                                                          
  LN10 =  2.30258509299404568402
  dm = (MASSMAX-MASSMIN)/(NX-1.)
  dz = (ZMAX-ZMIN)/(NY-1.)
  Z_ARRAY = ZMAX-dz*DINDGEN(ny)

  mass_array = 10.D^(MASSMIN+dm*dindgen(nx))/(1.e14 * h0)

  for i=0L, nx-1 do MASS_FUNCTION[i,*]/=alog(10.*mass_array[i])
  
  mass_func= {mass_func,Nx:nx,$
              Ny:ny,$
              x:dblarr(Nx),$
              y:dblarr(Ny),$
              z:dblarr(Nx,Ny)}
  
  mass_func.x = mass_array
  mass_func.y = reverse(Z_ARRAY)
  mass_func.z = (MASS_FUNCTION)


  ;; NAME = 'mf_x.dat'
  ;; Z = FLTARR(1,nz)
  ;; OPENR,1,NAME
  ;; READF,1,Z
  ;; CLOSE,1
  
  ;; NAME = 'mf_y.dat'
  ;; M = FLTARR(1,nm)
  ;; OPENR,1,NAME
  ;; READF,1,M
  ;; CLOSE,1
  ;; M=REFORM(10^M)
  ;; Z = REFORM(z)
END

