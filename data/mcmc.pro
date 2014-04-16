;MCMC-IDL v1.0 by Ankur Desai (Dec 2010), based on ml-metro5.c by Bill Sacks and
;George Hurtt and Rob Braswell
;see mcmc_example below for a test case and mcmc for usage information
;Contact: Ankur R Desai, desai@aos.wisc.edu


FUNCTION lambda_mass_relation,x,param,_EXTRA=ex
; Returns lambda given x = [0,mass/median(mass)] x =
; [1,E(redshift)/E(median z)] and x = [2,dmass/median(mass)] and 
; param[0] = Alambda, param[1] = Blambda and param[2] = Clambda
  lambda = param[0]*((x[0,*])^(param[1]))*(x[1,*])^(param[2])
  return,lambda
END

FUNCTION mass_calib_likelihood,x,y,param,model,valid,dat=dat,_EXTRA=ex,modout=modout

data = x



mass_func = {MF,inherits MASS_FUNC}
mass_func = y.MF

y0_to_ycyl = {y0_to_ycyl,inherits Y0_TO_YCYL}
y0_to_ycyl = y.y0_to_ycyl

Ncluster = N_elements(data.lambda)
read_params,params

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

params.theta[I_lambda_A] = 10^(param[0])
params.theta[I_lambda_B] = param[1]
params.theta[I_lambda_C] = param[2]
params.theta[I_lambda_D] = param[3]

; params.theta[I_SZ_A] = value[4]
; params.theta[I_SZ_B] = value[5]
; params.theta[I_SZ_C] = value[6]
; params.theta[I_SZ_D] = value[7]

loglike = -2.D*call_likelihood_library(Ncluster,params,mass_func,data,y0_to_ycyl) 

; openu,1,'params_chains.txt',/app
; printf,1,param[0],param[1],param[2] ,param[3]
; close,1

openu,1,'ll_chains.txt',/app
printf,1,param[0],param[1],param[2] ,param[3],loglike
close,1

if finite(loglike) eq 0 then loglike=1.e31

return,loglike

END


;----MAIN PROGRAM----
PRO mcmc,x,y,param,$
         numatonce=numatonce,random_start=random_start,numchains=numchains,numspinups=numspinups,$
         iter=iter,valid_frac=valid_frac,validdata=validdata,model=model,likelihood=likelihood,$
         outputll=outputll,outputvalue=outputvalue,outputy=outputy,_EXTRA = ex,fast=fast,medium=medium,ranacc=ranacc,quiet=quiet,superfast=superfast

print,model,likelihood

;MCMC parameter estimator
;based on sipnet
;Note: this can work for multiple data types (likelihood function
;would do the math)
;However, this is a single location version
;to make work at multiple locs, would need to modify likelihood function
;also no aggregation is done here
;SEE BELOW FOR AN EXAMPLE FOR HOW TO USE IN IDL

;---REQUIRED INPUTS---
;x is the input data for the model (no requirements on format except
;for what likelihood and model expect)
;y is the output data to compare (same as above)

;param is a structure  with the following properties
;param.name is an arry of name
;param.value is parameter value (initial guess) array
;param.max max value array
;param.min min value array
;param.knob is knob array (not used in this version, can be all zero)
;param.changeable is whether it should be fixed (0) or estimated (1)

;---OPTIONAL KEYWORDS---
;valid is % of datapoints in each interval that are valid (default is 100%)
;valid_frac is min % of datapoints to accept (default is 50%)

;model is the name of the model function (string) - default is "mcmc_testmodel"
;likelihood is the name of the likelihood function - it should call
;   model and compute the likelihood - default is "mcmc_likelihood"
;_EXTRA can be used to pass extra keywords to likelihood or model

;numatonce - how often to check for chain convergence (default is 10000)
;random_start - start at init guess if 0, else randomly within prior (default is 0)
;numchains - number of chains (default is 10)
;numspinups - how much burn in on final iteration (default is 125000)
;iter - max number of iterations both for chains and final (default is 375000)
;ranacc - tuning parameter for % of "worse" likelihoods to accept (default is -1) 
;/fast, /medium, and /superfast are keywords with different default settings for the above settings - all are faster than default (see below)
;/quiet turns off all printing of messages

;---OUTPUTS---
;outputll is the likelihood of output values, best value is first
;outputvalues are accepted parameter values (same format as param.value)
;outputy is the model run with the best outputvalues parameter set (same format as y)

  A_STAR = 0.4 ;(target rate)
  THRESH = 0.02 ;(+/- target rate)
  DEC = 0.99
  INC = DEC ^ ( (A_STAR-1)/A_STAR)
;  add_fraction = 0.5
  add_fraction = 1.0

;fast mode
  IF n_elements(fast) NE 0 THEN BEGIN
    numatonce = 1000l
    random_start = 1l
    numchains = 3l
    numspinups = 5000l
    iter = 10000l
  ENDIF 

;medium mode
  if n_elements(medium) ne 0 then begin
    numatonce = 10000l
    random_start = 1l
    numchains = 6l
    numspinups = 70000l
    iter = 150000l
    thresh = 0.025
  endif  

;super fast model
  IF n_elements(superfast) NE 0 THEN BEGIN
    numatonce = 1000l
    random_start = 0l
    numchains = 1l
    numspinups = 2000l
    iter = 5000l
  ENDIF 

;set defaults
  IF n_elements(numatonce) EQ 0 THEN numatonce = 10000l
  IF n_elements(random_start) EQ 0 THEN random_start = 0l
  IF n_elements(numchains) EQ 0 THEN numchains = 10l
  IF n_elements(numspinups) EQ 0 THEN numspinups = 125000l
  IF n_elements(iter) EQ 0 THEN iter = 375000l
  IF n_elements(ranacc) EQ 0 THEN ranacc = -1.0 ;set to -5.0 for compelx cost func
  IF n_elements(validdata) EQ 0 THEN BEGIN 
    validdata = y
    validdata[*] = 0
    gy = where(finite(y),ngy)
    IF ngy GT 0 THEN validdata[gy] = 1
  ENDIF 
  IF n_elements(valid_frac) EQ 0 THEN valid_frac = 0.5
  val = validdata GE valid_frac

  IF n_elements(model) EQ 0 THEN model = 'mcmc_testmodel'
  IF n_elements(likelihood) EQ 0 THEN likelihood = 'mcmc_likelihood'

  max = double(param.max)
  min = double(param.min)
  range = max-min
  change = where(param.changeable EQ 1,nchange)
   
  verybestll = double(-1e31)

;build some chains
  FOR c = 0,numchains-1 DO BEGIN

    IF ~keyword_set(quiet) THEN print,'Starting Chain ',c+1

;reset values
    value = double(param.value)
    IF (random_start EQ 1) AND (c GT 0) THEN value[change] = min[change] + (range[change] * randomu(systime(/sec),nchange))
;    knob = double(param.knob)
    knob = replicate(add_fraction,n_elements(range))
    converged = 0
    steps = 0l
    seed = systime(/sec)
    oldvalue = value
    bestvalue = value
    ll = (-1.0) * call_FUNCTION(likelihood,x,y,value,model,val,_extra=ex)
    bestll = ll
    ll_old = bestll

;go through the chain until convergence
    WHILE (converged EQ 0) && (steps LT iter) DO BEGIN 
      ichgs = change[long(randomu(seed,numatonce,/double)*nchange)]
      tune = randomu(seed,numatonce,/double)-0.5
      ran_accept = ranacc * randomu(seed,numatonce,gamma=1,/double)
;-5.0
      yes = 0l

      FOR k = 0l,numatonce-1l DO BEGIN
;randomly pick a parameter to change
        accept = 1
        ichg = ichgs[k]
        oldval = value[ichg] 
        newval = (knob[ichg] * range[ichg] * tune[k])+oldval
        IF (newval GT max[ichg]) OR (newval LT min[ichg]) THEN accept = 0

;run the model and calculate the likelihood
        IF accept EQ 1 THEN BEGIN
          value[ichg] = newval
          ll = (-1.0) * call_FUNCTION(likelihood,x,y,value,model,val,_extra=ex)
          IF (ll LE ll_old) && (ran_accept[k] GE (ll-ll_old)) THEN accept = 0
          IF (accept EQ 1) && (ll GT bestll) THEN BEGIN 
              print, 'New best: ',reform[value,ll]
              print, 'Old best: ',reform[bestvalue,bestll]
            bestll = ll
            bestvalue = value
          ENDIF
        ENDIF

;keep track of accepted parameter sets, tune knob
        IF accept EQ 1 THEN BEGIN 
          ll_old = ll
          yes++
          knob[ichg]*=INC
        ENDIF ELSE BEGIN
          value[ichg] = oldval
          knob[ichg] = (knob[ichg]*DEC)>(1e-9)
        ENDELSE

      ENDFOR

;check for convergence of this chain
      steps+=numatonce
      IF ~keyword_set(quiet) THEN print,'  iteration ',steps,' accept ',float(yes)/numatonce,' llmax ',bestll
;      IF float(yes)/numatonce GE a_star THEN BEGIN 
      IF abs(float(yes)/numatonce - a_star) LT thresh THEN BEGIN
        converged = 1
        IF ~keyword_set(quiet) THEN print,'  Chain ',c+1,' Converged LL: ',ll
        IF ~keyword_set(quiet) THEN print,'  Values: ',value
        IF bestll GE verybestll THEN BEGIN
          IF ~keyword_set(quiet) THEN print,'    And it is the best chain so far!'
          verybestll = bestll
          verybestvalue = bestvalue
          endvalue = value
          endll = ll
          endknob = knob
        ENDIF
      ENDIF ELSE BEGIN
        IF ~keyword_set(quiet) THEN print,'  Chain ',c+1,' not yet converged'
        yes = 0l
      ENDELSE 
 
    ENDWHILE 

    IF converged EQ 0 THEN IF ~keyword_set(quiet) THEN print,'  Chain ',c+1,' did not converge'

  ENDFOR

;start at end of best chain

  IF n_elements(endvalue) EQ 0 THEN BEGIN
    IF ~keyword_set(quiet) THEN print,'No chains converged, try changing mcmc iterations'
    IF ~keyword_set(quiet) THEN print,'Starting from best value'
    endvalue = bestvalue
    endll = bestll
    endknob = knob
  ENDIF 

  value = endvalue
  ll_old = endll
  knob = endknob
  seed = systime(/sec)
  bestvalue = value
  bestll = endll

  ichgs = change[long(randomu(seed,iter,/double)*nchange)]
  tune = randomu(seed,iter,/double)-0.5
  ran_accept = ranacc * randomu(seed,iter,gamma=1,/double)
  yes = 0l
  yes2 = yes

  outputll = fltarr(iter)
  outputvalue = fltarr(n_elements(value),iter)  

  FOR k = 0l,iter-1l DO BEGIN 
    IF k MOD numatonce EQ 0 THEN IF ~keyword_set(quiet) THEN print,'Final iteration ',k,' accepted ',yes2,' saved ',yes

;randomly pick a parameter to change
    accept = 1
    ichg = ichgs[k]
    oldval = value[ichg]
    newval = (knob[ichg] * range[ichg] * tune[k])+oldval
    IF (newval GT max[ichg]) OR (newval LT min[ichg]) THEN accept = 0

;run the model and calculate the likelihood
    IF accept EQ 1 THEN BEGIN
      value[ichg] = newval
      ll = (-1.0) * call_FUNCTION(likelihood,x,y,value,model,val,_extra=ex)
      IF (ll LE ll_old) && (ran_accept[k] GE (ll-ll_old)) THEN accept = 0
      IF (accept EQ 1) && (ll GT bestll) THEN BEGIN 
        bestll = ll
        bestvalue = value
;        IF bestll GE verybestll THEN BEGIN
;          verybestll = bestll
;          verybestvalue = bestvalue
;        ENDIF 
      ENDIF
    ENDIF

;output values
    IF accept EQ 1 THEN BEGIN
      yes2++
      ll_old = ll
;if past numspinups, then start saving vals
      IF k GE numspinups THEN BEGIN 
        outputll[yes] = ll
        outputvalue[*,yes] = value
        yes++
      ENDIF 
    ENDIF ELSE BEGIN
      value[ichg] = oldval
    ENDELSE 

  ENDFOR

;create the history of accepted values
  IF yes EQ 0 THEN BEGIN 
    outputvalue = bestvalue
    outputll = bestll
  ENDIF ELSE BEGIN 
    outputvalue = [[bestvalue],[outputvalue[*,0:(yes-1)]]]
    outputll = [bestll,outputll[0:(yes-1)]]
    srt = reverse(sort(outputll))
    outputvalue = outputvalue[*,srt]
    outputll = outputll[srt]
  ENDELSE 

;output values if outputy is there
  IF arg_present(outputy) THEN BEGIN
    dummy = call_FUNCTION(likelihood,x,y,outputvalue[*,0],model,val,_extra=ex,modout=outputy)
  ENDIF 

  IF ~keyword_set(quiet) THEN print,'MCMC complete '
  IF ~keyword_set(quiet) THEN print,'Best LL: ',outputll[0]
  IF ~keyword_set(quiet) THEN print,'Values: ',outputvalue[*,0]

END

;----EXAMPLE CODE----

;----EXAMPLE CODE----
PRO mcmc_lambda,scalecut=scalecut,fast=fast,medium=medium,npoints=npoints,$
                input=input,zmax=zmax,wait=wait
if not is_def(scalecut) then scalecut =0
if not is_def(zmax) then zmax =0.7
if not is_def(lambdacut) then lambdacut =20
if not is_def(wait) then wait =0

basedir='/big/users/saro/LRG/SVA1/GOLD/CATALOGUES/MCMC/CH'
dirout = 1
while file_exists(basedir+strtrim(string(dirout),1)) eq 1 do dirout++
dirout = basedir+strtrim(string(dirout),1)+'/'

spawn,'rm -f ll_chains.txt'

openw,4,'log.txt'
printf,4,'Job started at '+systime()

;setcolors
omegam = 0.3
omegal = 0.7
M500str = 'M!D500!N [M!D'+sunsymbol()+'!N]'
lamg = greek('lambda',/append)
alg = 'A!D'+lamg+'!N'
blg = 'B!D'+lamg+'!N'
clg = 'C!D'+lamg+'!N'
dlg = 'D!D'+lamg+'!N'

; clus_cat = mrdfits('/big/users/saro/LRG/SVA1/GOLD/CATALOGUES/sva1_gold_1.0_run_redmapper_v5.9_lgt5_catalog.fit',1)
; restore,'/big/users/saro/LRG/SVA1/GOLD/CATALOGUES/SPT_values_clusters.sav'

; relyable = where(clus_cat.lambda_chisq/clus_cat.scaleval ge scalecut and clus_cat.zred le zmax $
;                  and clus_cat.lambda_chisq ge lambdacut, nii)

; clus_cat = clus_cat(relyable)
; XI       = XI(relyable)       
; Y0       = Y0(relyable)      
; SIGMA_Y0 = SIGMA_Y0(relyable)
; FIELD    = FIELD(relyable)   

; lambda = clus_cat.lambda_chisq
; dlambda = clus_cat.lambda_chisq_e
; zred = clus_cat.zred
I_H0 = 12      
read_data,data,zmax=0.7
;data = data(where(data.redshift gt 0.3))
;data = data[0]
read_params,params
read_massfunction,mass_func=mass_func,h0 = params.theta[I_H0]
read_y0_to_ycyl,y0_to_ycyl

x = data
y = create_struct(['MF','y0_to_ycyl'],mass_func,y0_to_ycyl)

;The input parameter file for mcmc
param = { name :       [alg,          blg,  clg,    dlg], $
          value :      [alog10(25) ,  1.1,    0,    0.4], $
          max :        [alog10(300),   3,    5,      1], $
          min :        [alog10(1),    0.1,   -5,   0.05], $
          knob :       [1.0,          1.0,  1.0,    1.0], $
          changeable : [1  ,            1,   1,      1] }

model = 'lambda_mass_relation'

;Call MCMC
mcmc,x,y,param,outputll=outputll,outputvalue=outputvalue,$
  model=model,likelihood='mass_calib_likelihood',validdata=1,/superfast

print,'Done..'
printf,4,'Finished at '+systime()
close,4
spawn,'mkdir '+ dirout

spawn,'mv log.txt '+dirout
spawn,'mv ll_chains.txt '+dirout
save,filename=dirout+'mcmc_chains.sav',alg,blg,clg,dlg,params,param,outputll,outputvalue,data

;plots_mcmc,outputvalue=outputvalue,outputll=outputll,param=param,input=input

stop

nparam = n_elements(param.value)
bestf =  fltarr(nparam)
sig1 =  fltarr(nparam)
for i =0, nparam -1 do begin
    histogauss68,outputvalue[i,*],kj,xtitle=param.name[i]    
    print,param.name[i] + ' marginalized constraints: ',kj[1], '+-',kj[2]
    bestf[i] = kj[1]
    sig1[i] = kj[2]
endfor
;print,'Tot clusters, Zmed and Mmed: ',n_elements(cross),zmed,mmed

;Plot outputs

best = float(outputvalue[*,0])
stop

END
