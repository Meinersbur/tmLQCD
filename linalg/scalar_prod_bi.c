/* $Id$ */

/****************************************************************************
 *
 *     Scalar product routine adapted for the bispinor case
 *
 * Author: Thomas Chiarappa
 *         Thomas.Chiarappa@mib.infn.it
 *
 ***************************************************************************/

#include <stdlib.h>
#ifdef MPI
#include <mpi.h>
#endif
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include "su3.h"
#include "scalar_prod_bi.h"


/*  <S,R>=S^* times R */
complex scalar_prod_bi(bispinor * const S, bispinor * const R, const int N){

  int ix;
  static double ks,kc,ds,tr,ts,tt;
  spinor *s,*r, *t, *u;
  complex c;
  
  /* Real Part */

  ks=0.0;
  kc=0.0;
  
  for (ix = 0; ix < N; ix++){

    s=(spinor *) &S[ix].sp_up;
    r=(spinor *) &R[ix].sp_up;
    t=(spinor *) &S[ix].sp_dn;
    u=(spinor *) &R[ix].sp_dn;

    
    ds=(*r).s0.c0.re*(*s).s0.c0.re+(*r).s0.c0.im*(*s).s0.c0.im+
       (*r).s0.c1.re*(*s).s0.c1.re+(*r).s0.c1.im*(*s).s0.c1.im+
       (*r).s0.c2.re*(*s).s0.c2.re+(*r).s0.c2.im*(*s).s0.c2.im+
       (*r).s1.c0.re*(*s).s1.c0.re+(*r).s1.c0.im*(*s).s1.c0.im+
       (*r).s1.c1.re*(*s).s1.c1.re+(*r).s1.c1.im*(*s).s1.c1.im+
       (*r).s1.c2.re*(*s).s1.c2.re+(*r).s1.c2.im*(*s).s1.c2.im+
       (*r).s2.c0.re*(*s).s2.c0.re+(*r).s2.c0.im*(*s).s2.c0.im+
       (*r).s2.c1.re*(*s).s2.c1.re+(*r).s2.c1.im*(*s).s2.c1.im+
       (*r).s2.c2.re*(*s).s2.c2.re+(*r).s2.c2.im*(*s).s2.c2.im+
       (*r).s3.c0.re*(*s).s3.c0.re+(*r).s3.c0.im*(*s).s3.c0.im+
       (*r).s3.c1.re*(*s).s3.c1.re+(*r).s3.c1.im*(*s).s3.c1.im+
       (*r).s3.c2.re*(*s).s3.c2.re+(*r).s3.c2.im*(*s).s3.c2.im+
       (*u).s0.c0.re*(*t).s0.c0.re+(*u).s0.c0.im*(*t).s0.c0.im+
       (*u).s0.c1.re*(*t).s0.c1.re+(*u).s0.c1.im*(*t).s0.c1.im+
       (*u).s0.c2.re*(*t).s0.c2.re+(*u).s0.c2.im*(*t).s0.c2.im+
       (*u).s1.c0.re*(*t).s1.c0.re+(*u).s1.c0.im*(*t).s1.c0.im+
       (*u).s1.c1.re*(*t).s1.c1.re+(*u).s1.c1.im*(*t).s1.c1.im+
       (*u).s1.c2.re*(*t).s1.c2.re+(*u).s1.c2.im*(*t).s1.c2.im+
       (*u).s2.c0.re*(*t).s2.c0.re+(*u).s2.c0.im*(*t).s2.c0.im+
       (*u).s2.c1.re*(*t).s2.c1.re+(*u).s2.c1.im*(*t).s2.c1.im+
       (*u).s2.c2.re*(*t).s2.c2.re+(*u).s2.c2.im*(*t).s2.c2.im+
       (*u).s3.c0.re*(*t).s3.c0.re+(*u).s3.c0.im*(*t).s3.c0.im+
       (*u).s3.c1.re*(*t).s3.c1.re+(*u).s3.c1.im*(*t).s3.c1.im+
       (*u).s3.c2.re*(*t).s3.c2.re+(*u).s3.c2.im*(*t).s3.c2.im;

    /* Kahan Summation */    
    tr=ds+kc;
    ts=tr+ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;
  }
  kc=ks+kc;

#if defined MPI
  MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  kc = ks;
#endif

  c.re = kc;

  /* Imaginary Part */

  ks=0.0;
  kc=0.0;
  
  for (ix=0;ix<N;ix++){

    s=(spinor *) &S[ix].sp_up;
    r=(spinor *) &R[ix].sp_up;
    t=(spinor *) &S[ix].sp_dn;
    u=(spinor *) &R[ix].sp_dn;
    
    ds=-(*r).s0.c0.re*(*s).s0.c0.im+(*r).s0.c0.im*(*s).s0.c0.re-
       (*r).s0.c1.re*(*s).s0.c1.im+(*r).s0.c1.im*(*s).s0.c1.re-
       (*r).s0.c2.re*(*s).s0.c2.im+(*r).s0.c2.im*(*s).s0.c2.re-
       (*r).s1.c0.re*(*s).s1.c0.im+(*r).s1.c0.im*(*s).s1.c0.re-
       (*r).s1.c1.re*(*s).s1.c1.im+(*r).s1.c1.im*(*s).s1.c1.re-
       (*r).s1.c2.re*(*s).s1.c2.im+(*r).s1.c2.im*(*s).s1.c2.re-
       (*r).s2.c0.re*(*s).s2.c0.im+(*r).s2.c0.im*(*s).s2.c0.re-
       (*r).s2.c1.re*(*s).s2.c1.im+(*r).s2.c1.im*(*s).s2.c1.re-
       (*r).s2.c2.re*(*s).s2.c2.im+(*r).s2.c2.im*(*s).s2.c2.re-
       (*r).s3.c0.re*(*s).s3.c0.im+(*r).s3.c0.im*(*s).s3.c0.re-
       (*r).s3.c1.re*(*s).s3.c1.im+(*r).s3.c1.im*(*s).s3.c1.re-
       (*r).s3.c2.re*(*s).s3.c2.im+(*r).s3.c2.im*(*s).s3.c2.re-
       (*u).s0.c0.re*(*t).s0.c0.im+(*u).s0.c0.im*(*t).s0.c0.re-
       (*u).s0.c1.re*(*t).s0.c1.im+(*u).s0.c1.im*(*t).s0.c1.re-
       (*u).s0.c2.re*(*t).s0.c2.im+(*u).s0.c2.im*(*t).s0.c2.re-
       (*u).s1.c0.re*(*t).s1.c0.im+(*u).s1.c0.im*(*t).s1.c0.re-
       (*u).s1.c1.re*(*t).s1.c1.im+(*u).s1.c1.im*(*t).s1.c1.re-
       (*u).s1.c2.re*(*t).s1.c2.im+(*u).s1.c2.im*(*t).s1.c2.re-
       (*u).s2.c0.re*(*t).s2.c0.im+(*u).s2.c0.im*(*t).s2.c0.re-
       (*u).s2.c1.re*(*t).s2.c1.im+(*u).s2.c1.im*(*t).s2.c1.re-
       (*u).s2.c2.re*(*t).s2.c2.im+(*u).s2.c2.im*(*t).s2.c2.re-
       (*u).s3.c0.re*(*t).s3.c0.im+(*u).s3.c0.im*(*t).s3.c0.re-
       (*u).s3.c1.re*(*t).s3.c1.im+(*u).s3.c1.im*(*t).s3.c1.re-
       (*u).s3.c2.re*(*t).s3.c2.im+(*u).s3.c2.im*(*t).s3.c2.re;
    
    tr=ds+kc;
    ts=tr+ks;
    tt=ts-ks;
    ks=ts;
    kc=tr-tt;
  }
  kc=ks+kc;


#if defined MPI
  MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  kc = ks;
#endif

  c.im = kc;
  return(c);

}
