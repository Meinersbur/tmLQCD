/***********************************************************************
 *
 * Copyright (C) 2006,2007,2008 Thomas Chiarappa, Carsten Urbach
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "linsolve.h"
#include "linalg_eo.h"
#include "start.h"
#include "tm_operators.h"
#include "Nondegenerate_Matrix.h"
#include "phmc.h"
#include "chebyshev_polynomial_nd.h"

#include "bgq/bgq_spinorfield.h"



#define PI 3.141592653589793

double func(double u, double exponent){
  return pow(u,exponent);
}


void chebyshev_coefs(double aa, double bb, double c[], int n, double exponent){
  int k,j;
  double fac,bpa,bma,*f;
  double inv_n;


  inv_n=1./(double)n;
  f=calloc(n,sizeof(double));/*vector(0,n-1);*/
  if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
    printf("PHMC: chebyshev_polynomial\n");
    printf("PHMC: n= %d inv_n=%e \n",n,inv_n);
    printf("PHMC: allocation !!!\n");
  }
  fflush(stdout);
  bma=0.5*(bb-aa);
  bpa=0.5*(bb+aa);
  for (k=0;k<n;k++) {
    double y=cos(PI*(k+0.5)*inv_n);
    f[k]=func(y*bma+bpa,exponent);
  }
  fac=2.0*inv_n;
  for (j=0;j<n;j++) {
    double sum=0.0;
    for (k=0;k<n;k++)
      sum += f[k]*cos(PI*j*(k+0.5)*inv_n);
    c[j]=fac*sum;
  }
  free(f);


}
#undef PI


/****************************************************************************  
 *
 * computation of, despite of the name, (Q Q^dagger) on a vector
 *   by using the chebyshev approximation for the function ()^1/4
 * subtraction of low-lying eigenvalues is not yet implemented for this
 *
 **************************************************************************/


void QdaggerQ_poly(spinor *R_s, spinor *R_c, double *c, int n, 
                   spinor *S_s, spinor *S_c){

  int j;
  double fact1, fact2, temp1, temp2, temp3, temp4;

  spinor *svs_=NULL, *svs=NULL, *ds_=NULL, *ds=NULL, *dds_=NULL, *dds=NULL, 
         *auxs_=NULL, *auxs=NULL, *aux2s_=NULL, *aux2s=NULL, *aux3s_=NULL, 
         *aux3s=NULL;
  spinor *svc_=NULL, *svc=NULL, *dc_=NULL, *dc=NULL, *ddc_=NULL, 
         *ddc=NULL, *auxc_=NULL, *auxc=NULL, *aux2c_=NULL, *aux2c=NULL, 
         *aux3c_=NULL, *aux3c=NULL;


#if 0
   svs_  = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
   svs   = (spinor *)(((unsigned long int)(svs_)+ALIGN_BASE)&~ALIGN_BASE);
   ds_   = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
   ds    = (spinor *)(((unsigned long int)(ds_)+ALIGN_BASE)&~ALIGN_BASE);
   dds_  = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
   dds   = (spinor *)(((unsigned long int)(dds_)+ALIGN_BASE)&~ALIGN_BASE);
   auxs_ = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
   auxs  = (spinor *)(((unsigned long int)(auxs_)+ALIGN_BASE)&~ALIGN_BASE);
   aux2s_= calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
   aux2s = (spinor *)(((unsigned long int)(aux2s_)+ALIGN_BASE)&~ALIGN_BASE);
   aux3s_= calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
   aux3s = (spinor *)(((unsigned long int)(aux3s_)+ALIGN_BASE)&~ALIGN_BASE);
   svc_  = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
   svc   = (spinor *)(((unsigned long int)(svc_)+ALIGN_BASE)&~ALIGN_BASE);
   dc_   = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
   dc    = (spinor *)(((unsigned long int)(dc_)+ALIGN_BASE)&~ALIGN_BASE);
   ddc_  = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
   ddc   = (spinor *)(((unsigned long int)(ddc_)+ALIGN_BASE)&~ALIGN_BASE);
   auxc_ = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
   auxc  = (spinor *)(((unsigned long int)(auxc_)+ALIGN_BASE)&~ALIGN_BASE);
   aux2c_= calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
   aux2c = (spinor *)(((unsigned long int)(aux2c_)+ALIGN_BASE)&~ALIGN_BASE);
   aux3c_= calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
   aux3c = (spinor *)(((unsigned long int)(aux3c_)+ALIGN_BASE)&~ALIGN_BASE);
#else
   spinor *mem = malloc_aligned(12 * VOLUMEPLUSRAND/2 * sizeof(spinor), BGQ_ALIGNMENT_L2);
   svs = &mem[0 * VOLUMEPLUSRAND/2];
   ds = &mem[1 * VOLUMEPLUSRAND/2];
   dds = &mem[2 * VOLUMEPLUSRAND/2];
   auxs = &mem[3 * VOLUMEPLUSRAND/2];
   aux2s = &mem[4 * VOLUMEPLUSRAND/2];
   aux3s = &mem[5 * VOLUMEPLUSRAND/2];
   svc = &mem[6 * VOLUMEPLUSRAND/2];
   dc = &mem[7 * VOLUMEPLUSRAND/2];
   ddc = &mem[8 * VOLUMEPLUSRAND/2];
   auxc = &mem[9 * VOLUMEPLUSRAND/2];
   aux2c = &mem[10 * VOLUMEPLUSRAND/2];
   aux3c = &mem[11 * VOLUMEPLUSRAND/2];
   bgq_weylfield_collection *col = bgq_spinorfields_allocate(12, mem, VOLUMEPLUSRAND/2);
#endif


   fact1=4/(phmc_cheb_evmax-phmc_cheb_evmin);
   fact2=-2*(phmc_cheb_evmax+phmc_cheb_evmin)/(phmc_cheb_evmax-phmc_cheb_evmin);

   zero_spinor_field(&ds[0],VOLUME/2);
   spinorfield_setOddness(&ds[0], 0);
   zero_spinor_field(&dds[0],VOLUME/2); 
   spinorfield_setOddness(&dds[0], 0);
   zero_spinor_field(&dc[0],VOLUME/2);
   spinorfield_setOddness(&dc[0], 1);
   zero_spinor_field(&ddc[0],VOLUME/2); 
   spinorfield_setOddness(&ddc[0], 1);

   /*   sub_low_ev(&aux3[0], &S[0]);  */
   assign(&aux3s[0], &S_s[0],VOLUME/2);  
   assign(&aux3c[0], &S_c[0],VOLUME/2);  
   
   /*  Use the Clenshaw's recursion for the Chebysheff polynomial */
   for (j=n-1; j>=1; j--) {
     assign(&svs[0],&ds[0],VOLUME/2);
     assign(&svc[0],&dc[0],VOLUME/2); 
       
     /*     
     if ( (j%10) == 0 ) {
  	 sub_low_ev(&aux[0], &d[0]);
     }
     else { */
     assign(&auxs[0], &ds[0], VOLUME/2);
     assign(&auxc[0], &dc[0], VOLUME/2);
     /*   } */  


     Q_Qdagger_ND(&R_s[0], &R_c[0], &auxs[0], &auxc[0]);

     temp1=-1.0;
     temp2=c[j];
     assign_mul_add_mul_add_mul_add_mul_r(&ds[0] , &R_s[0], &dds[0], &aux3s[0], fact2, fact1, temp1, temp2,VOLUME/2);
     assign_mul_add_mul_add_mul_add_mul_r(&dc[0] , &R_c[0], &ddc[0], &aux3c[0], fact2, fact1, temp1, temp2,VOLUME/2);
     assign(&dds[0], &svs[0],VOLUME/2);
     assign(&ddc[0], &svc[0],VOLUME/2);

   }
     
   /*     sub_low_ev(&R[0],&d[0]);  */ 
   assign(&R_s[0], &ds[0],VOLUME/2);  
   assign(&R_c[0], &dc[0],VOLUME/2);  


   Q_Qdagger_ND(&auxs[0], &auxc[0], &R_s[0], &R_c[0]);

   temp1=-1.0;
   temp2=c[0]/2;
   temp3=fact1/2;
   temp4=fact2/2;
   assign_mul_add_mul_add_mul_add_mul_r(&auxs[0], &ds[0], &dds[0], &aux3s[0], temp3, temp4, temp1, temp2,VOLUME/2);
   assign_mul_add_mul_add_mul_add_mul_r(&auxc[0], &dc[0], &ddc[0], &aux3c[0], temp3, temp4, temp1, temp2,VOLUME/2);
   assign(&R_s[0], &auxs[0],VOLUME/2);
   assign(&R_c[0], &auxc[0],VOLUME/2);
     
   /*     addproj_q_invsqrt(&R[0], &S[0]); */
    
   /*
#ifndef _SOLVER_OUTPUT
     if(g_proc_id == g_stdio_proc){
       printf("Order of Chebysheff approximation = %d\n",j); 
       fflush( stdout);};
#endif
   */

    
   bgq_spinorfields_free(col);
   free(mem);
}
  


double cheb_eval(int M, double *c, double s){

  double d=0,dd=0, sv, z, z2, res;
  int j;

  z = (2.0*s - phmc_cheb_evmin - phmc_cheb_evmax)/(double)(phmc_cheb_evmax - phmc_cheb_evmin);
  z2 = 2.0*z;

  for(j=M-1; j>=1; j--){
    sv = d;
    d = z2*d - dd + c[j];
    dd = sv;
    }

  res = z*d - dd + 0.5*c[0];

  return(res);  
}

/**************************************************************************
 *
 * The externally accessible function is
 *
 *   void degree_of_polynomial_nd(void)
 *     Computation of (QdaggerQ)^1/4
 *     by using the chebyshev approximation for the function ()^1/4  
 *
 *
 *****************************************************************************/


void degree_of_polynomial_nd(const int degree_of_p){
  int j;
  double temp, temp2;
  static int ini=0;

  double sum=0.0;

  spinor *ss=NULL, *ss_=NULL, *sc=NULL, *sc_=NULL;
  spinor *auxs=NULL, *auxs_=NULL, *auxc=NULL, *auxc_=NULL;
  spinor *aux2s=NULL, *aux2s_=NULL, *aux2c=NULL, *aux2c_=NULL;

  phmc_dop_n_cheby=degree_of_p+1;
  if(ini==0){
    phmc_dop_cheby_coef = calloc(phmc_dop_n_cheby,sizeof(double));
    ini=1;
  }


#if 0
  ss_   = calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
  auxs_ = calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
  aux2s_= calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
  sc_   = calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
  auxc_ = calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
  aux2c_= calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
  
  ss    = (spinor *)(((unsigned long int)(ss_)+ALIGN_BASE)&~ALIGN_BASE);
  auxs  = (spinor *)(((unsigned long int)(auxs_)+ALIGN_BASE)&~ALIGN_BASE);
  aux2s = (spinor *)(((unsigned long int)(aux2s_)+ALIGN_BASE)&~ALIGN_BASE);
  sc    = (spinor *)(((unsigned long int)(sc_)+ALIGN_BASE)&~ALIGN_BASE);
  auxc  = (spinor *)(((unsigned long int)(auxc_)+ALIGN_BASE)&~ALIGN_BASE);
  aux2c = (spinor *)(((unsigned long int)(aux2c_)+ALIGN_BASE)&~ALIGN_BASE);
  
#else
  spinor *mem = malloc_aligned(6 * VOLUMEPLUSRAND/2 * sizeof(spinor), BGQ_ALIGNMENT_L2);
  ss   = &mem[0 * VOLUMEPLUSRAND/2];
  auxs = &mem[1 * VOLUMEPLUSRAND/2];
  aux2s= &mem[2 * VOLUMEPLUSRAND/2];
  sc   = &mem[3 * VOLUMEPLUSRAND/2];
  auxc = &mem[4 * VOLUMEPLUSRAND/2];
  aux2c= &mem[5 * VOLUMEPLUSRAND/2];
  bgq_weylfield_collection *col = bgq_spinorfields_allocate(6, mem, VOLUMEPLUSRAND/2);
#endif
  
  
  chebyshev_coefs(phmc_cheb_evmin, phmc_cheb_evmax, phmc_dop_cheby_coef, phmc_dop_n_cheby, -0.5);

  random_spinor_field(ss,VOLUME/2, 1);
  spinorfield_setOddness(ss, 0);
  random_spinor_field(sc,VOLUME/2, 1);
  spinorfield_setOddness(sc, 1);

  if((g_proc_id == g_stdio_proc) && (g_debug_level > 0)){
    printf("NDPOLY MD Polynomial: EVmin = %e  EVmax = %e  \n", phmc_cheb_evmin, phmc_cheb_evmax);
    printf("NDPOLY MD Polynomial: the degree was set to: %d\n", phmc_dop_n_cheby);
    fflush(stdout);
  }

  /* Here we check the accuracy */
  QdaggerQ_poly(&auxs[0], &auxc[0], phmc_dop_cheby_coef, phmc_dop_n_cheby, &ss[0], &sc[0]);
  Q_Qdagger_ND(&aux2s[0], &aux2c[0], &auxs[0], &auxc[0]);
  QdaggerQ_poly(&auxs[0], &auxc[0], phmc_dop_cheby_coef, phmc_dop_n_cheby, &aux2s[0], &aux2c[0]);

  diff(&aux2s[0],&auxs[0],&ss[0],VOLUME/2);
  temp=square_norm(&aux2s[0],VOLUME/2, 1)/square_norm(&ss[0],VOLUME/2, 1)/4.0;

  diff(&aux2c[0],&auxc[0],&sc[0],VOLUME/2);
  temp2 = square_norm(&aux2c[0],VOLUME/2, 1)/square_norm(&sc[0],VOLUME/2, 1)/4.0;

  if(g_epsbar == 0.){ 
    temp2 = 0.0;
  }

  if(g_proc_id == g_stdio_proc && g_debug_level > 0){
    /* this is || (P S P - 1)X ||^2 /|| 2X ||^2 */
    /* where X is a random spinor field         */
    printf("NDPOLY MD Polynomial: relative squared accuracy in components:\n UP=%e  DN=%e \n", temp, temp2);
    /*     printf("NDPOLY: Sum remaining | c_n | = %e \n", sum); */
    fflush(stdout);
  }

  if(g_debug_level > 1) {
    temp = cheb_eval(phmc_dop_n_cheby, phmc_dop_cheby_coef, phmc_cheb_evmin);
    temp *= phmc_cheb_evmin;
    temp *= cheb_eval(phmc_dop_n_cheby, phmc_dop_cheby_coef, phmc_cheb_evmin);
    temp = 0.5*fabs(temp - 1);
    if(g_proc_id == g_stdio_proc) {
      printf("PHMC: Delta_IR at s=%f:    | P s_low P - 1 |/2 = %e \n", phmc_cheb_evmin, temp);
    }
  }
  /* RECALL THAT WE NEED AN EVEN DEGREE !!!! */

#if 0
   free(ss_);   
   free(auxs_); 
   free(aux2s_);
   free(sc_);   
   free(auxc_); 
   free(aux2c_);
#else
   bgq_spinorfields_free(col);
   free(mem);
#endif

}
