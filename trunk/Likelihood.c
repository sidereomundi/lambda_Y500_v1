#include <omp.h>
#include "utils.h"
#include "Data.h"
#include "Likelihood.h"
#include "ScalingRelation.h"
#include "Statistics.h"
#include "Function.h"
#include "linalg.h"
#include "Random.h"
#include "Cosmology.h"
#include "Integration.h"
#include "Interpolation.h"
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

// Calculate integral (P(Y500obs|mass)*P(mass|lambda)dmass)
double Get_likelihood_for_one_cluster(Parameters_struct* params,
				Function_1D* vect_p_m_lambda_debiased,
				double redshift, 
				double* y0,
				double* dy0,
				double* filter,
				y0_to_Y500cyl_structure* y0_to_Y500cyl,
				int n,double* mass_array)
{
  // Get parameter values
  double DY500 = params->theta[I_Y500cyl_D];
  double h0 = params->theta[I_H0];
  double Om = params->theta[I_OM];
  double Ol = params->theta[I_OL];
  double w0 = params->theta[I_W0];
  double wa = params->theta[I_WA];

  double Mpc_over_arcmin=0.;

  int i, j, Ngrid;
  Ngrid = vect_p_m_lambda_debiased->N;
  //  int n = mass_func->Nx;
  int Nfilt = params->theta[I_Nfilt];                                                                                                  
  /* printf("%d\n",n);  */
  /* for (i=0;i<n;i++)printf("%f\n",mass_array[i]); */
  /* return 0;  */

  double* R500c_vect = (double*) malloc(n * sizeof(double)); 
  double* Y500cyl_vect = (double*) malloc(n * sizeof(double)); 
  double* y0vect = (double*) malloc(n * sizeof(double));
  double* dy0vect = (double*) malloc(n * sizeof(double));
  double* y0_to_Y500cyl_factor = (double*) malloc(n * sizeof(double));      
  
  double mindy0 = 1.e6,miny0y500=1.e6,Y500cyl_0=0.;

  double* drdm = (double*) malloc(n * sizeof(double));
  double* dmdY500cyl = (double*) malloc(n * sizeof(double));

  double* central_Y500 = (double*) malloc(n * sizeof(double));      
  double* width_of_gaussian = (double*) malloc(n * sizeof(double));      

  double* Y500cyl_max = (double*) malloc(n * sizeof(double));      
  double* Y500cyl_min = (double*) malloc(n * sizeof(double));      
  double* Y500_arr = (double*) malloc(n * sizeof(double));      

  Function_1D* p_Y500cyl_mass = (Function_1D*) malloc(sizeof(Function_1D));
  p_Y500cyl_mass->N = n;
  p_Y500cyl_mass->x = (double*) malloc(n * sizeof (double));
  p_Y500cyl_mass->y = (double*) malloc(n * sizeof (double));

  Function_1D* p_Y500cyl_mass_interp = (Function_1D*) malloc(sizeof(Function_1D));
  p_Y500cyl_mass_interp->N = Ngrid;
  p_Y500cyl_mass_interp->x = (double*) malloc(Ngrid * sizeof (double));
  p_Y500cyl_mass_interp->y = (double*) malloc(Ngrid * sizeof (double));


  // Model prediction
  Function_1D* lognormal_Y500cyl = (Function_1D*) malloc(sizeof(Function_1D));
  lognormal_Y500cyl->N = n;
  lognormal_Y500cyl->x = (double*) malloc(n * sizeof (double));
  lognormal_Y500cyl->y = (double*) malloc(n * sizeof (double));

  // Measured value
  Function_1D* normal_Y500cyl = (Function_1D*) malloc(sizeof(Function_1D));
  normal_Y500cyl->N = n;
  normal_Y500cyl->x = (double*) malloc(n * sizeof (double));
  normal_Y500cyl->y = (double*) malloc(n * sizeof (double));

  double* integrand = (double*) malloc(n * sizeof(double));      
  double* integrand_grid = (double*) malloc(Ngrid * sizeof(double));      
  double probability,loglikelihood;

  Mpc_over_arcmin =  SebCosmolib_AngDiamDist(redshift, Om, Ol, w0)/(60.*180./PI);

  for (i=0;i<n;i++) 
    {
      R500c_vect[i] = R500_in_arcmin_given_M500(mass_array[i],redshift,h0,Om,Ol,w0,wa);
      Y500cyl_vect[i] = m500c2Y500_cyl(mass_array[i]*1.e14, redshift, h0, Om, Ol, w0, wa)/Mpc_over_arcmin/Mpc_over_arcmin;
      //      R500c_vect[i] = R500_in_arcmin_given_M500(vect_p_m_lambda_debiased->x[i],redshift,h0,Om,Ol,w0,wa);
      //      Y500cyl_vect[i] = m500c2Y500_cyl(vect_p_m_lambda_debiased->x[i]*1.e14, redshift, h0, Om, Ol, w0, wa)/Mpc_over_arcmin/Mpc_over_arcmin;
    }      
  
  Interp1D_Vector(filter,y0,Nfilt,R500c_vect,y0vect,n);   
  Interp1D_Vector(filter,dy0,Nfilt,R500c_vect,dy0vect,n);   
  Interp1D_Vector(y0_to_Y500cyl->r500,y0_to_Y500cyl->y0_to_Y500cyl,146,R500c_vect,y0_to_Y500cyl_factor,n);   
  
  
  for(i = 0 ; i < n; i++)
    {
      if(abs(dy0vect[i])<mindy0)mindy0=abs(dy0vect[i]); 
      if(abs(y0_to_Y500cyl_factor[i])<miny0y500)miny0y500=abs(y0_to_Y500cyl_factor[i]); 
    }


  for(i = 0 ; i < n; i++)
    {
      //      if (i < n-1) drdm[i] = (R500c_vect[i+1]-R500c_vect[i])/(vect_p_m_lambda_debiased->x[i+1]-vect_p_m_lambda_debiased->x[i]);
      //      if (i < n-1) dmdY500cyl[i] = ((vect_p_m_lambda_debiased->x[i+1]-vect_p_m_lambda_debiased->x[i])/(Y500cyl_vect[i+1]-Y500cyl_vect[i]));
      
      if(dy0vect[i]<mindy0)dy0vect[i]=mindy0;
      if(y0_to_Y500cyl_factor[i]<miny0y500)y0_to_Y500cyl_factor[i]=miny0y500;
      
      central_Y500[i] = y0vect[i]*y0_to_Y500cyl_factor[i];
      width_of_gaussian[i] = dy0vect[i]*y0_to_Y500cyl_factor[i];

      Y500cyl_max[i] = central_Y500[i] + 5*width_of_gaussian[i];
      Y500cyl_min[i] = central_Y500[i] - 5*width_of_gaussian[i];
      
      if(Y500cyl_min[i]>0.)Y500cyl_min[i]=1.e-10;
      if(Y500cyl_max[i]<(Y500cyl_vect[i]*(1.+DY500)*5))Y500cyl_max[i]=Y500cyl_vect[i]*(1.+DY500)*5;
      //printf("%f %f %f %f %f\n",Y500cyl_min[i],Y500cyl_max[i],central_Y500[i],width_of_gaussian[i],Y500cyl_vect[i]);
    }
  
  //  drdm[n-1] = 0.;
  //  dmdY500cyl[n-1] = 0.;
    
  for (i=0;i<n;i++)
    {
      p_Y500cyl_mass->x[i] = 0;
      p_Y500cyl_mass->y[i] = 0;
      for (j=0;j<n;j++)
	{
	  Y500_arr[j] = Y500cyl_min[i] + j*(Y500cyl_max[i] - Y500cyl_min[i])/(n - 1.);
	  Y500cyl_0 = Y500cyl_vect[i];

	  lognormal_Y500cyl->x[j] = Y500_arr[j];
	  normal_Y500cyl->x[j] = Y500_arr[j];
	  
	  lognormal_Y500cyl->y[j] = 0;
	  normal_Y500cyl->y[j] = 0;
	  if(Y500cyl_0>0)
	    if(Y500_arr[j]>0)
	      lognormal_Y500cyl->y[j] = lognormal(Y500_arr[j], Y500cyl_0, DY500);
	  normal_Y500cyl->y[j] = normal(Y500_arr[j],central_Y500[i],width_of_gaussian[i]);
	}

      //      Normalize_1D(lognormal_Y500cyl);
      //      Normalize_1D(normal_Y500cyl);
      
      for (j=0;j<n;j++)
	{
	  integrand[j] = lognormal_Y500cyl->y[j]*normal_Y500cyl->y[j];
	  /* printf("%d %f %f %f %f\n",i,Y500_arr[j],lognormal_Y500cyl->y[j],normal_Y500cyl->y[j],integrand[j]); */
	}
      //      p_Y500cyl_mass->x[i]=vect_p_m_lambda_debiased->x[i];
      p_Y500cyl_mass->x[i]=mass_array[i];
      p_Y500cyl_mass->y[i]=IntTrap_varstep(Y500_arr,integrand, n, n);
    }    

  //  Normalize_1D(p_Y500cyl_mass);      
  /* for (i=0;i<n;i++)printf("%f %f \n",p_Y500cyl_mass->x[i],p_Y500cyl_mass->y[i]); */
  /* return 0; */

//  Now interpolate on vect_p_m_lambda_debiased->x

  for (i=0;i<Ngrid;i++)
    p_Y500cyl_mass_interp->x[i]= vect_p_m_lambda_debiased->x[i];  
  //printf("%d %f\n",Ngrid,vect_p_m_lambda_debiased->x[i]);//printf("%d %f %f\n",Ngrid,p_Y500cyl_mass_interp->x[i],vect_p_m_lambda_debiased->x[i]);
  Interp1D_Vector(p_Y500cyl_mass->x, p_Y500cyl_mass->y, n, p_Y500cyl_mass_interp->x, p_Y500cyl_mass_interp->y, Ngrid);    

  //Now integrate P(Y|M) P(M|lambda_obs)
  for (i=0;i<Ngrid;i++)
    integrand_grid[i] = p_Y500cyl_mass_interp->y[i]*vect_p_m_lambda_debiased->y[i];  
  
  probability=IntTrap_varstep(vect_p_m_lambda_debiased->x,integrand_grid, Ngrid, Ngrid);

  if((probability<=0.) || (probability!=probability))
    probability = 1.e-10;

  loglikelihood = log(probability);
  
  free(R500c_vect);
  free(Y500cyl_vect);
  free(y0vect);
  free(dy0vect);
  free(y0_to_Y500cyl_factor);
  free(drdm);
  free(dmdY500cyl);
  free(central_Y500);
  free(width_of_gaussian); 
  free(Y500cyl_max);
  free(Y500cyl_min); 
  free(Y500_arr);
  Free_Function_1D(lognormal_Y500cyl);
  Free_Function_1D(normal_Y500cyl);
  Free_Function_1D(p_Y500cyl_mass);
  Free_Function_1D(p_Y500cyl_mass_interp);
  free(integrand);
  free(integrand_grid);
  return loglikelihood;
}

/* This is the main program
	In the following, the underlying value is lambda_0, the measured value is lambda_obs
	* P(lambda_obs|lambda_0)
	* P(lambda_0|lambda_obs) (not yet debiased)
	* P(M|lambda_obs)_debiased
	* P(Y500cyl|M) */
double Likelihood_Lambdacalib(Parameters_struct* params, 
			      int Ncluster, 
			      Data_structure* data,
			      Function_2D* mass_func,
			      y0_to_Y500cyl_structure* y0_to_Y500cyl){
  
  int i,j;
  FILE *fp;       
  
  double loglikelihood=0.,loglikelihood_one_cluster[Ncluster];

  // Variables for the different distributions
  Function_1D** vect_p_lambda0_lambdaobs;
  Function_1D** vect_p_m_lambda_debiased;
  vect_p_lambda0_lambdaobs = (Function_1D**) malloc(Ncluster * sizeof(Function_1D*));
  vect_p_m_lambda_debiased = (Function_1D**) malloc(Ncluster * sizeof(Function_1D*));  
	
  double** plan_lambda_m;

  int Nfilt = params->theta[I_Nfilt];
    
  // lambda scaling parameters
  double Alambda = params->theta[I_lambda_A];
  double Blambda = params->theta[I_lambda_B];
  double Clambda = params->theta[I_lambda_C];
  double Dlambda = params->theta[I_lambda_D];
  
  // Cosmology parameters
  
  double h0 = params->theta[I_H0];
  double Om = params->theta[I_OM];
  double Ol = params->theta[I_OL];
  double w0 = params->theta[I_W0];
  double wa = params->theta[I_WA];
 

  // lambda grid in lambda_vect
  double lambda_max = 300;
  double lambda_min = 1;
  int Nlambda = 50;
  double lambda0 = 1.;
  double dlambda = 10.;
  double lambda_vect[Nlambda];
  double R500c_vect[Nlambda];

#pragma omp parallel for
  for(i = 0; i < Nlambda; i++)
    {
      lambda_vect[i] = pow(10,log10(lambda_min)+i*(log10(lambda_max)-log10(lambda_min))/(Nlambda-1));
    }

  // Calculate 2D array with lambda-M relation containing Dlambda scatter and uncertainty
  plan_lambda_m = Get_array_lambda_M(lambda_vect, Nlambda, Dlambda);

  // Get P(lambda_0|lambda_obs) (not yet debiased)
#pragma omp parallel for
  for(i = 0; i < Ncluster; i++)
    {
      vect_p_lambda0_lambdaobs[i] = Get_P_lambda0_lambdaobs(plan_lambda_m, lambda_vect, Nlambda, data[i].lambda,data[i].dlambda);
    }
  
#pragma omp parallel for
  for(i = 0; i < Nlambda; i++)
    free(plan_lambda_m[i]);
  free(plan_lambda_m);
  
  /* Get P(M|lambda)_debiased: P(lambda_0|lambda_obs) * mass function
     finally normalize it */
#pragma omp parallel for
    for(i = 0; i < Ncluster; i++)
    {	      	  
      vect_p_m_lambda_debiased[i] = Get_P_M_lambda(params, vect_p_lambda0_lambdaobs[i], mass_func, data[i].redshift);
      Normalize_1D(vect_p_m_lambda_debiased[i]);      
/*       for (j=0;j<Nlambda;j++){printf(" %f %f %f %f \n",vect_p_lambda0_lambdaobs[i]->x[j],vect_p_lambda0_lambdaobs[i]->y[j], */
/* 				     vect_p_m_lambda_debiased[i]->x[j],vect_p_m_lambda_debiased[i]->y[j]);} */
    }
  
    int N_mass_array = 30;
    double max_mass = log10(3.e15);
    double min_mass = log10(5.e13);
    double mass_array[N_mass_array];
    for(i=0;i<N_mass_array;i++)
      { 
	mass_array[i] = pow(10,min_mass + i*(max_mass - min_mass)/(N_mass_array -1.))/1.e14;
	//	printf("%f\n", mass_array[i]);
      }

#pragma omp parallel for
    for(i = 0; i < Ncluster; i++)
      loglikelihood_one_cluster[i]=  Get_likelihood_for_one_cluster(params,vect_p_m_lambda_debiased[i],data[i].redshift,data[i].y0,data[i].dy0,data[i].filter,y0_to_Y500cyl,N_mass_array,mass_array);

/*     for(i = 0; i < Ncluster; i++) */
/*       printf("%d %f\n",i,loglikelihood_one_cluster[i]); */
    
    for(i = 0; i < Ncluster; i++)
      loglikelihood+=loglikelihood_one_cluster[i];

    // Free everything
#pragma omp parallel for
    for(i = 0; i < Ncluster; i++)
      {
      Free_Function_1D(vect_p_lambda0_lambdaobs[i]);
      Free_Function_1D(vect_p_m_lambda_debiased[i]);
    }
  free(vect_p_lambda0_lambdaobs);
  free(vect_p_m_lambda_debiased);
  return loglikelihood;
}


int call_likelihood_library (int      argc, void *   argv[])
{
  Parameters_struct* params; 
  int Ncluster,Nparams,Nfilt;
  Data_structure* data;
  Function_2D* mass_func;
  int i,j;
  double * loglikelihood;
  y0_to_Y500cyl_structure* y0_to_Y500cyl; 

  Ncluster=*(int *) argv[0];
  Nparams=*(int *) argv[1];
  
  struct  Function_Param{
    int N;
    double THETA[Nparams];
  } *PF_IDL;
  
  PF_IDL=((struct Function_Param *)argv[2]);
  params = malloc(sizeof(Parameters_struct));
  /* number of parameters in the likelihood */
  params->N = PF_IDL[0].N;
  params->theta = malloc(params->N*sizeof(double));
  //#pragma omp parallel for
  for (i=0;i<PF_IDL[0].N;i++)
    params->theta[i] = PF_IDL[0].THETA[i];

  Nfilt = params->theta[I_Nfilt];

  struct  Function_MF{
    int Nx;
    int Ny;
    double X[500];
    double Y[51];
    double Z[51][500];
  } *MF_IDL;

  MF_IDL=((struct Function_MF *)argv[3]);
  mass_func = malloc(sizeof(Function_2D));
  mass_func->Nx = MF_IDL[0].Nx; /* 1st dimension is mass */
  mass_func->Ny = MF_IDL[0].Ny; /* 2nd dimension is redshift */
  //  printf("%d %d %d %d\n ",MF_IDL[0].Nx,mass_func->Nx,MF_IDL[0].Ny,mass_func->Ny);
  mass_func->x = malloc(mass_func->Nx*sizeof(double));
  mass_func->y = malloc(mass_func->Ny*sizeof(double));
  mass_func->z = malloc(mass_func->Nx*sizeof(double*));
  for(i = 0; i < mass_func->Nx; ++i) 
    mass_func->z[i] = malloc(mass_func->Ny*sizeof(double));

  for (i=0;i<MF_IDL[0].Nx;i++)
    mass_func->x[i] = MF_IDL[0].X[i];

  for (i=0;i<MF_IDL[0].Ny;i++)
    mass_func->y[i] = MF_IDL[0].Y[i];
  
  for(i = 0; i < mass_func->Nx; i++)
    for(j = 0; j < mass_func->Ny; j++)
      { mass_func->z[i][j] = MF_IDL[0].Z[j][i];
	//	printf("%d %d %lf \n",i,j,MF_IDL[0].Z[j][i]*1.e10);
      }

  struct  Function_Data{
    int     CLUS_ID;    
    int     FIELD_ID;       
    double   REDSHIFT;       
    double   LAMBDA;         
    double   dLAMBDA;      
    double   y0[Nfilt];             
    double   dy0[Nfilt];          
    double   filter[Nfilt];          
  } *DF_IDL;
  

  DF_IDL=((struct Function_Data *)argv[4]);
  //  for (i=0;i<Nfilt;i++){printf("y0: %e dy0: %e filter %f \n",DF_IDL[0].y0[i],DF_IDL[0].dy0[i],DF_IDL[0].filter[i]);};    
  //  for (id=0;id<Ncluster;id++){printf("%d: lambda %f\n",id,DF_IDL[id].LAMBDA);};
  data = malloc(Ncluster * sizeof (Data_structure));
  
  //#pragma omp parallel for
  for (i=0;i<Ncluster;i++)
    {
      data[i].clus_id    =    DF_IDL[i].CLUS_ID;
      data[i].field_id   =    DF_IDL[i].FIELD_ID;
      data[i].redshift   =    DF_IDL[i].REDSHIFT;
      data[i].lambda	 =    DF_IDL[i].LAMBDA;	 
      data[i].dlambda    =    DF_IDL[i].dLAMBDA;
      for(j=0;j<Nfilt;j++)
	{
	  data[i].y0[j]     =    DF_IDL[i].y0[j];
	  data[i].dy0[j]    =    DF_IDL[i].dy0[j];
	  data[i].filter[j] =    DF_IDL[i].filter[j];
	}
    }
  
  y0_to_Y500cyl = ((struct y0_to_Y500cyl_structure *)argv[5]);

  loglikelihood = (double *) argv[6];
  
  /*   for (i=0;i<10;i++){printf("lambda: %f\n",data[i].lambda);};     */
  //  printf("clus_id %d \n",data[0].clus_id);
  //  for (i=0;i<Nfilt;i++){printf("y0: %e dy0: %e filter %f \n",data[0].y0[i],data[0].dy0[i],data[0].filter[i]);};    
  //  for (i=0;i<Ny500interp;i++){printf("r500: %f y0y500: %f\n",y0_to_Y500cyl->r500[i],y0_to_Y500cyl->y0_to_Y500cyl[i]);};    

  *loglikelihood = Likelihood_Lambdacalib(params,Ncluster,data,mass_func,y0_to_Y500cyl);
  
  free(params->theta);  
  free(params);        
  free(mass_func->x);                                                                                                        
  free(mass_func->y);                                                                                                        
  for(i = 0; i < mass_func->Nx; ++i)                                                                                         
    free(mass_func->z[i]);                                                                                             
  free(mass_func->z);                                                                                                        
  free(mass_func);  
  free(data);

  return 0;
}

