#include "Hopping_unopt.h"

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "su3.h"
#include "sse.h"
#ifdef MPI
#  include "xchange_field.h"
#  if defined _USE_HALFSPINOR
#    include "xchange_halffield.h"
#  endif
#endif
#include "boundary.h"
#include "init_dirac_halfspinor.h"
#include "update_backward_gauge.h"

#include "geometry_eo.c"


/* 3. */
/* l output , k input*/
/* for ieo=0, k resides on  odd sites and l on even sites */
void Hopping_unopt_halfspinor(const int ieo, spinor * const l, spinor * const k){
/* total 1632 flops */
  int i,ix;
  su3 * restrict U ALIGN;
  spinor * restrict s ALIGN;
  spinor rs;
  static su3_vector psi, chi, psi2, chi2;
  halfspinor * restrict * phi ALIGN;
  halfspinor32 * restrict * phi32 ALIGN;

#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge();
  }
#endif

  if(k == l){
    printf("Error in H_psi (simple.c):\n");
    printf("Arguments k and l must be different\n");
    printf("Program aborted\n");
    exit(1);
  }
  s = k;

  if(ieo == 0) {
    U = g_gauge_field_copy[0][0];
  }
  else {
    U = g_gauge_field_copy[1][0];
  }
    phi = NBPointer[ieo];

    /**************** loop over all lattice sites ****************/
    ix=0;
    /* #pragma ivdep*/
    for(i = 0; i < (VOLUME)/2; i++){
      _vector_assign(rs.s0, (*s).s0);
      _vector_assign(rs.s1, (*s).s1);
      _vector_assign(rs.s2, (*s).s2);
      _vector_assign(rs.s3, (*s).s3);
      s++;
      /*********************** direction +0 ************************/

      _vector_add(psi, rs.s0, rs.s2);     /* 6 flops */
      _vector_add(psi2, rs.s1, rs.s3); /* 6 flops */
      _su3_multiply(chi,(*U),psi);     /* 66 flops (36 muls, 30 adds) */
      _su3_multiply(chi2,(*U),psi2);     /* 66 flops (36 muls, 30 adds) */
      _complex_times_vector((*phi[ix]).s0, ka0, chi); /* 18 flops (12 muls, 6 adds)    */
      _complex_times_vector((*phi[ix]).s1, ka0, chi2);/* 18 flops (12 muls, 6 adds) */

      U++;
      ix++;

      /*********************** direction -0 ************************/

      _vector_sub((*phi[ix]).s0, rs.s0, rs.s2);/* 6 flops */
      _vector_sub((*phi[ix]).s1, rs.s1, rs.s3);/* 6 flops */

      ix++;

      /*********************** direction +1 ************************/

      _vector_i_add(psi, rs.s0, rs.s3);/* 6 flops */
      _vector_i_add(psi2, rs.s1, rs.s2);/* 6 flops */
      _su3_multiply(chi, (*U), psi);/* 66 flops (36 muls, 30 adds) */
      _su3_multiply(chi2, (*U), psi2);/* 66 flops (36 muls, 30 adds) */
      _complex_times_vector((*phi[ix]).s0, ka1, chi);/* 18 flops (12 muls, 6 adds) */
      _complex_times_vector((*phi[ix]).s1, ka1, chi2);/* 18 flops (12 muls, 6 adds) */

      U++;
      ix++;

      /*********************** direction -1 ************************/

      _vector_i_sub((*phi[ix]).s0, rs.s0, rs.s3);/* 6 flops */
      _vector_i_sub((*phi[ix]).s1, rs.s1, rs.s2);/* 6 flops */

      ix++;
      /*********************** direction +2 ************************/

      _vector_add(psi, rs.s0, rs.s3);/* 6 flops */
      _vector_sub(psi2, rs.s1, rs.s2);/* 6 flops */
      _su3_multiply(chi,(*U),psi);/* 66 flops (36 muls, 30 adds) */
      _su3_multiply(chi2,(*U),psi2);/* 66 flops (36 muls, 30 adds) */
      _complex_times_vector((*phi[ix]).s0, ka2, chi);/* 18 flops (12 muls, 6 adds) */
      _complex_times_vector((*phi[ix]).s1, ka2, chi2);/* 18 flops (12 muls, 6 adds) */

      U++;
      ix++;

      /*********************** direction -2 ************************/

      _vector_sub((*phi[ix]).s0, rs.s0, rs.s3);/* 6 flops */
      _vector_add((*phi[ix]).s1, rs.s1, rs.s2);/* 6 flops */
      ix++;

      /*********************** direction +3 ************************/

      _vector_i_add(psi, rs.s0, rs.s2);/* 6 flops */
      _vector_i_sub(psi2, rs.s1, rs.s3);/* 6 flops */
      _su3_multiply(chi, (*U), psi);/* 66 flops (36 muls, 30 adds) */
      _su3_multiply(chi2,(*U),psi2);/* 66 flops (36 muls, 30 adds) */
      _complex_times_vector((*phi[ix]).s0, ka3, chi);/* 18 flops (12 muls, 6 adds) */
      _complex_times_vector((*phi[ix]).s1, ka3, chi2);/* 18 flops (12 muls, 6 adds) */

      U++;
      ix++;
      /*********************** direction -3 ************************/

      _vector_i_sub((*phi[ix]).s0, rs.s0, rs.s2);/* 6 flops */
      _vector_i_add((*phi[ix]).s1, rs.s1, rs.s3);/* 6 flops */

      ix++;
      /************************ end of loop ************************/
    }
#    if (defined MPI && !defined _NO_COMM)
    xchange_halffield();
#    endif
    s = l;
    phi = NBPointer[2 + ieo];
    if(ieo == 0) {
      U = g_gauge_field_copy[1][0];
    }
    else {
      U = g_gauge_field_copy[0][0];
    }

    ix = 0;
    /* #pragma ivdep */
    for(i = 0; i < (VOLUME)/2; i++){
      /*********************** direction +0 ************************/
      _vector_assign(rs.s0, (*phi[ix]).s0);
      _vector_assign(rs.s2, (*phi[ix]).s0);
      _vector_assign(rs.s1, (*phi[ix]).s1);
      _vector_assign(rs.s3, (*phi[ix]).s1);
      ix++;
      /*********************** direction -0 ************************/
      _su3_inverse_multiply(chi,(*U),(*phi[ix]).s0);/* 66 flops (36 muls, 30 adds) */
      _su3_inverse_multiply(chi2,(*U),(*phi[ix]).s1);/* 66 flops (36 muls, 30 adds) */
      _complexcjg_times_vector(psi,ka0,chi);/* 18 flops (12 muls, 6 adds) */
      _complexcjg_times_vector(psi2,ka0,chi2);/* 18 flops (12 muls, 6 adds) */
      _vector_add_assign(rs.s0, psi);/* 6 flops */
      _vector_sub_assign(rs.s2, psi);/* 6 flops */
      _vector_add_assign(rs.s1, psi2);/* 6 flops */
      _vector_sub_assign(rs.s3, psi2);/* 6 flops */
      ix++;
      U++;
      /*********************** direction +1 ************************/

      _vector_add_assign(rs.s0, (*phi[ix]).s0);/* 6 flops */
      _vector_i_sub_assign(rs.s3, (*phi[ix]).s0);/* 6 flops */

      _vector_add_assign(rs.s1, (*phi[ix]).s1);/* 6 flops */
      _vector_i_sub_assign(rs.s2, (*phi[ix]).s1);/* 6 flops */

      ix++;
      /*********************** direction -1 ************************/

      _su3_inverse_multiply(chi,(*U), (*phi[ix]).s0);/* 66 flops (36 muls, 30 adds) */
      _su3_inverse_multiply(chi2, (*U), (*phi[ix]).s1);/* 66 flops (36 muls, 30 adds) */
      _complexcjg_times_vector(psi,ka1,chi);/* 18 flops (12 muls, 6 adds) */
      _complexcjg_times_vector(psi2,ka1,chi2);/* 18 flops (12 muls, 6 adds) */
      _vector_add_assign(rs.s0, psi);/* 6 flops */
      _vector_i_add_assign(rs.s3, psi);/* 6 flops */
      _vector_add_assign(rs.s1, psi2);/* 6 flops */
      _vector_i_add_assign(rs.s2, psi2);/* 6 flops */

      U++;
      ix++;

      /*********************** direction +2 ************************/

      _vector_add_assign(rs.s0, (*phi[ix]).s0);/* 6 flops */
      _vector_add_assign(rs.s3, (*phi[ix]).s0);/* 6 flops */

      _vector_add_assign(rs.s1, (*phi[ix]).s1);/* 6 flops */
      _vector_sub_assign(rs.s2, (*phi[ix]).s1);/* 6 flops */

      ix++;
      /*********************** direction -2 ************************/

      _su3_inverse_multiply(chi,(*U), (*phi[ix]).s0);/* 66 flops (36 muls, 30 adds) */
      _su3_inverse_multiply(chi2, (*U), (*phi[ix]).s1);/* 66 flops (36 muls, 30 adds) */
      _complexcjg_times_vector(psi,ka2,chi);/* 18 flops (12 muls, 6 adds) */
      _complexcjg_times_vector(psi2,ka2,chi2);/* 18 flops (12 muls, 6 adds) */
      _vector_add_assign(rs.s0, psi);/* 6 flops */
      _vector_sub_assign(rs.s3, psi);/* 6 flops */
      _vector_add_assign(rs.s1, psi2);/* 6 flops */
      _vector_add_assign(rs.s2, psi2);/* 6 flops */

      U++;
      ix++;
      /*********************** direction +3 ************************/

      _vector_add_assign(rs.s0, (*phi[ix]).s0);/* 6 flops */
      _vector_i_sub_assign(rs.s2, (*phi[ix]).s0);/* 6 flops */

      _vector_add_assign(rs.s1, (*phi[ix]).s1);/* 6 flops */
      _vector_i_add_assign(rs.s3, (*phi[ix]).s1);/* 6 flops */

      ix++;

      /*********************** direction -3 ************************/

      _su3_inverse_multiply(chi,(*U), (*phi[ix]).s0);/* 66 flops (36 muls, 30 adds) */
      _su3_inverse_multiply(chi2, (*U), (*phi[ix]).s1);/* 66 flops (36 muls, 30 adds) */
      _complexcjg_times_vector(psi,ka3,chi);/* 18 flops (12 muls, 6 adds) */
      _complexcjg_times_vector(psi2,ka3,chi2);/* 18 flops (12 muls, 6 adds) */
      _vector_add((*s).s0, rs.s0, psi);/* 6 flops */
      _vector_i_add((*s).s2, rs.s2, psi);/* 6 flops */
      _vector_add((*s).s1, rs.s1, psi2);/* 6 flops */
      _vector_i_sub((*s).s3, rs.s3, psi2);/* 6 flops */

      U++;
      ix++;
      s++;
    }
  
}

/*###########################################################################*/


#define CYC_T(it) (it==-1 ? T+1 : it)
#define CYC_X(ix) (ix==-1 ? LX+1 : ix)
#define CYC_Y(iy) (iy==-1 ? LY+1 : iy)
#define CYC_Z(iz) (iz==-1 ? LZ+1 : iz)

spinor *ref_spinorfield(spinor * const arr,bool isEven,int it,int ix,int iy,int iz) {
    assert(isEven == (it+ix+iy+iz)%2);
    return &(arr)[g_lexic2eosub[Index(it,ix,iy,iz)]]; 
    // alterntive: arr[g_lexic2eo[g_ipt[CYC_T(it)][CYC_X(ix)][CYC_Y(iy)][CYC_Z(iz)]]
}

/* get an element of a spinorfield */
#define REF_SPINORFIELD(arr,isEven,it,ix,iy,iz) *ref_spinorfield(arr,isEven,it,ix,iy,iz)


// arr contains only even elements
#define REF_SPINORFIELD_EVEN(arr,it,ix,iy,iz) REF_SPINORFIELD(arr,true,it,ix,iy,iz)
// arr contains only odd elements
#define REF_SPINORFIELD_ODD(arr,it,ix,iy,iz) REF_SPINORFIELD(arr,false,it,ix,iy,iz)
    

/* get element of spinorfield, assume that arr points to correct half of spinorfield (either even half or odd half)*/
#define REF_SPINORFIELD_EOSUB(arr,it,ix,iy,iz) arr[g_lexic2eosub[Index(it,ix,iy,iz)]]

#define REF_GAUGEFIELD(arr,it,ix,iy,iz) arr[g_lexic2eo[Index(it,ix,iy,iz)]]


/*
idx_noneo = 0...VOLUME-1

int* g_t;       g_t: idx_noneo -> t of idx=Index(t,x,y,z); t={0..LT-1}
int* g_x;
int* g_y;
int* g_z;

int**** g_ipt;  g_ipt: t,x,y,z -> Index(t,x,y,z)
t = 0..LT+1 (-1 mapped to LT+1) or 0..LT-1 depending on PARALLELT

g_iup[ix][0] = Index(x0+1, x1, x2, x3);
g_idn[ix][0] = Index(x0-1, x1, x2, x3);
g_iup[ix][1] = Index(x0, x1+1, x2, x3);
g_idn[ix][1] = Index(x0, x1-1, x2, x3);
g_iup[ix][2] = Index(x0, x1, x2+1, x3);
g_idn[ix][2] = Index(x0, x1, x2-1, x3);
g_iup[ix][3] = Index(x0, x1, x2, x3+1);
g_idn[ix][3] = Index(x0, x1, x2, x3-1);

g_coord[ix][t=0/x=1/y=2/z=3] -> global coordinate
*/



/* 8. */
/* l output , k input*/
/* for ieo=0, k resides on odd sites and l on even sites */
/* 1608 flops */
void Hopping_unopt_fullspinor(int ieo, spinor * const l, spinor * const k){
  int ix,iy;
  int ioff,ioff2,icx,icy;
  su3 * restrict up, * restrict um;
  spinor * restrict r, * restrict sp, * restrict sm;
  spinor temp;
  su3_vector psi1, psi2, psi, chi, phi1, phi3;
  bool isEven = (ieo == 0);

  /* for parallelization */
#    if (defined MPI && !(defined _NO_COMM))
  xchange_field(k, ieo);
#    endif

  if(k == l){
    printf("Error in H_psi (simple.c):\n");
    printf("Arguments k and l must be different\n");
    printf("Program aborted\n");
    exit(1);
  }
  if(ieo == 0)
    ioff = 0;
  else
    ioff = (VOLUME+RAND)/2;
  
  ioff2 = (VOLUME+RAND)/2-ioff;
  /**************** loop over all lattice sites ****************/

  //for (icx = ioff; icx < (VOLUME/2 + ioff); icx++){
  for (int jt = 0; jt < T; jt++)
    for (int jx = 0; jx < LX; jx++)
        for (int jy = 0; jy < LY; jy++)
            for (int jz = 0; jz < LZ; jz++) {
                if ((jt+jx+jy+jz) % 2 == ieo)   // ieo == 0 = even; WARNING: even/odd based on global (not local) coordinates
                    continue;
                icx = g_lexic2eo[Index(jt,jx,jy,jz)]; /* icx = ioff..ioff+VOLUME/2-1 */

    //ix=g_eo2lexic[icx];
                ix = Index(jt,jx,jy,jz); /* ix = 0..VOLUME-1 (mod 2) */

    // Spinor to save the result
    //r=l+(icx-ioff);
                //int r_ix = g_lexic2eosub[ix]; /* r_ix = 0..VOLUME/2-1 */
                //r = &l[r_ix];
                r = &REF_SPINORFIELD(l, isEven, jt,jx,jy,jz) ;

    /*********************** direction +0 ************************/
    //iy=g_iup[ix][0]; 
                iy = Index(jt+1,jx,jy,jz);
    icy=g_lexic2eosub[iy];
                
    //sp=k+icy;
        sp = &REF_SPINORFIELD(k,isEven,jt+1,jx,jy,jz);
    //up=&g_gauge_field[ix][0];          
        up = &REF_GAUGEFIELD(g_gauge_field,jt,jx,jy,jz)[0]; // gauge from (jt,jx,jy,jz) to (jt+1,jx,jy,jz)

    _vector_add(psi,(*sp).s0,(*sp).s2); /* 6 flops */ // (su3_vector) psi = sp.s0 + sp.s2

    _su3_multiply(chi,(*up),psi); /* 66 flops (36 muls, 30 adds) */ // (su3_vector) chi = up * psi (Matrix-Vector-Multiplication)
    _complex_times_vector(psi,ka0,chi); /* 18 flops (12 muls, 6 adds) */ // (su3_vector) psi = ka_0 * chi (Scalar-Vector-Multiplcation)

    _vector_assign(temp.s0,psi); // temp.s0 = psi
    _vector_assign(temp.s2,psi); // temp.s2 = psi

    _vector_add(psi,(*sp).s1,(*sp).s3); /* 6 flops */ // (su3_vector) psi = sp.s1 + sp.s3

    _su3_multiply(chi,(*up),psi); /* 66 flops (36 muls, 30 adds) */ // (su3_vector) chi = up * psi (Matrix-Vector-Multiplication)
    _complex_times_vector(psi,ka0,chi); /* 18 flops (12 muls, 6 adds) */ // (su3_vector) psi = ka_0 * chi (Scalar-Vector-Multiplcation)

    _vector_assign(temp.s1,psi); // temp.s1 = psi
    _vector_assign(temp.s3,psi); // temp.s3 = psi

    /*********************** direction -0 ************************/

    //iy=g_idn[ix][0]; 
        iy = Index(jt-1,jx,jy,jz);    
    icy=g_lexic2eosub[iy];

    //sm=k+icy;
        sm = &REF_SPINORFIELD(k,IsEven,jt-1,jx,jy,jz);
    //um=&g_gauge_field[iy][0];
        um = &REF_GAUGEFIELD(g_gauge_field,jt-1,jx,jy,jz)[0]; // gauge from (jt-1,jx,jy,jz) to (jt,jx,jy,jz)
        // um^-1 = um^\dagger = um* (conjugate transpose; deutsch: adjungierte Matrix) = gauge from (jt,jx,jy,jz) to (jt-1,jx,jy,jz)

    _vector_sub(psi,(*sm).s0,(*sm).s2); /* 6 flops */ // (su3_vector) psi = sm.s0 + sm.s2

    _su3_inverse_multiply(chi,(*um),psi); /* 66 flops (36 muls, 30 adds) */ // (su3_vector) chi = um^-1 * psi (Matrix-Vector-Multiplication)
    _complexcjg_times_vector(psi,ka0,chi);/* 18 flops (12 muls, 6 adds) */ 

    _vector_add_assign(temp.s0,psi);/* 6 flops */
    _vector_sub_assign(temp.s2,psi);/* 6 flops */

    _vector_sub(psi,(*sm).s1,(*sm).s3);/* 6 flops */

    _su3_inverse_multiply(chi,(*um),psi);  /* 66 flops (36 muls, 30 adds) */
    _complexcjg_times_vector(psi,ka0,chi);/* 18 flops (12 muls, 6 adds) */ 

    _vector_add_assign(temp.s1,psi); /* 6 flops */
    _vector_sub_assign(temp.s3,psi);/* 6 flops */

    /*********************** direction +1 ************************/

    iy=g_iup[ix][1]; icy=g_lexic2eosub[iy];

    sp=k+icy;

#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif

    _vector_i_add(psi,(*sp).s0,(*sp).s3);/* 6 flops */

    _su3_multiply(chi,(*up),psi);/* 66 flops (36 muls, 30 adds) */ 
    _complex_times_vector(psi,ka1,chi);/* 18 flops (12 muls, 6 adds) */ 

    _vector_add_assign(temp.s0,psi);/* 6 flops */
    _vector_i_sub_assign(temp.s3,psi);/* 6 flops */

    _vector_i_add(psi,(*sp).s1,(*sp).s2);/* 6 flops */

    _su3_multiply(chi,(*up),psi);/* 66 flops (36 muls, 30 adds) */ 
    _complex_times_vector(psi,ka1,chi);/* 18 flops (12 muls, 6 adds) */ 

    _vector_add_assign(temp.s1,psi);/* 6 flops */
    _vector_i_sub_assign(temp.s2,psi);/* 6 flops */

    /*********************** direction -1 ************************/

    iy=g_idn[ix][1]; icy=g_lexic2eosub[iy];

    sm=k+icy;
#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][1];
#    else
    um=up+1;
#    endif

    _vector_i_sub(psi,(*sm).s0,(*sm).s3);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka1,chi);

    _vector_add_assign(temp.s0,psi);
    _vector_i_add_assign(temp.s3,psi);

    _vector_i_sub(psi,(*sm).s1,(*sm).s2);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka1,chi);

    _vector_add_assign(temp.s1,psi);
    _vector_i_add_assign(temp.s2,psi);

    /*********************** direction +2 ************************/

    iy=g_iup[ix][2]; icy=g_lexic2eosub[iy];

    sp=k+icy;
#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    _vector_add(psi,(*sp).s0,(*sp).s3);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka2,chi);

    _vector_add_assign(temp.s0,psi);
    _vector_add_assign(temp.s3,psi);

    _vector_sub(psi,(*sp).s1,(*sp).s2);

    _su3_multiply(chi,(*up),psi);
    _complex_times_vector(psi,ka2,chi);

    _vector_add_assign(temp.s1,psi);
    _vector_sub_assign(temp.s2,psi);


    /*********************** direction -2 ************************/

    iy=g_idn[ix][2]; icy=g_lexic2eosub[iy];

    sm=k+icy;
#    ifndef _GAUGE_COPY
    um = &g_gauge_field[iy][2];
#    else
    um = up +1;
#    endif

    _vector_sub(psi,(*sm).s0,(*sm).s3);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka2,chi);

    _vector_add_assign(temp.s0,psi);
    _vector_sub_assign(temp.s3,psi);

    _vector_add(psi,(*sm).s1,(*sm).s2);

    _su3_inverse_multiply(chi,(*um),psi);
    _complexcjg_times_vector(psi,ka2,chi);

    _vector_add_assign(temp.s1,psi);
    _vector_add_assign(temp.s2,psi);

    /*********************** direction +3 ************************/

    iy=g_iup[ix][3]; icy=g_lexic2eosub[iy];

    sp=k+icy;
#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    _vector_i_add(psi,(*sp).s0,(*sp).s2); /* 6 flops */

    _su3_multiply(chi,(*up),psi);/* 66 flops (36 muls, 30 adds) */ 
    _complex_times_vector(psi,ka3,chi);/* 18 flops (12 muls, 6 adds) */ 

    _vector_add_assign(temp.s0,psi);/* 6 flops */
    _vector_i_sub_assign(temp.s2,psi);/* 6 flops */

    _vector_i_sub(psi,(*sp).s1,(*sp).s3);/* 6 flops */

    _su3_multiply(chi,(*up),psi);/* 66 flops (36 muls, 30 adds) */ 
    _complex_times_vector(psi,ka3,chi);/* 18 flops (12 muls, 6 adds) */ 

    _vector_add_assign(temp.s1,psi);/* 6 flops */
    _vector_i_add_assign(temp.s3,psi);/* 6 flops */

    /*********************** direction -3 ************************/

    iy=g_idn[ix][3]; icy=g_lexic2eosub[iy];

    sm=k+icy;
#    ifndef _GAUGE_COPY
    um = &g_gauge_field[iy][3];
#    else
    um = up+1;
#    endif

    _vector_i_sub(psi,(*sm).s0,(*sm).s2);/* 6 flops */

    _su3_inverse_multiply(chi,(*um),psi);/* 66 flops (36 muls, 30 adds) */ 
    _complexcjg_times_vector(psi,ka3,chi);/* 18 flops (12 muls, 6 adds) */ 

    _vector_add((*r).s0, temp.s0, psi);/* 6 flops */
    _vector_i_add((*r).s2, temp.s2, psi);/* 6 flops */

    _vector_i_add(psi,(*sm).s1,(*sm).s3);/* 6 flops */

    _su3_inverse_multiply(chi,(*um),psi);/* 66 flops (36 muls, 30 adds) */ 
    _complexcjg_times_vector(psi,ka3,chi);/* 18 flops (12 muls, 6 adds) */ 

    _vector_add((*r).s1, temp.s1, psi);/* 6 flops */
    _vector_i_sub((*r).s3, temp.s3, psi);/* 6 flops */
    /************************ end of loop ************************/
  }
}
