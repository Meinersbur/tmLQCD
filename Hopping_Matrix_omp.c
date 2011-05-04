#include <omp.h>

#if (defined _USE_HALFSPINOR && defined BGL && defined XLC && !defined NOBGP && defined OMP)
/**********************************
 *
 * Blue Gene/L Version
 *
 * Author: Carsten Urbach
 *
 **********************************/
#undef HopVerMsg
#define HopVerMsg printf("Hopping_Matrix half_spinor BGL OpenMP edition");

/* 2. */
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k){
  su3 * restrict Ustart ALIGN;
  halfspinor * restrict * phi ALIGN;
#pragma disjoint(*Ustart, *phi, *l, *l)

#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge();
  }
#endif

  __alignx(16, l);
  __alignx(16, k);

    __alignx(16, HalfSpinor);
    /* We will run through the source vector now */
    /* instead of the solution vector            */

    _prefetch_spinor(k);

    /* s contains the source vector */

    if(ieo == 0) {
      Ustart = g_gauge_field_copy[0][0];
    }
    else {
      Ustart = g_gauge_field_copy[1][0];
    }
    phi = NBPointer[ieo];

    _prefetch_su3(Ustart);
    /**************** loop over all lattice sites ******************/
    #pragma omp for
    for(int i = 0; i < (VOLUME)/2; i++){
          int ix;
  su3 * restrict U ALIGN;
  spinor * restrict s ALIGN;
  halfspinor32 * restrict * phi32 ALIGN;
  /* We have 32 registers available */
  double _Complex reg00, reg01, reg02, reg03, reg04, reg05;
  double _Complex reg10, reg11, reg12, reg13, reg14, reg15;
  /* For the gauge field, reuse the first three!*/
  double _Complex u00, u01, u02, u10, u11, u12;
  double _Complex reg20, reg21;
  /* The following contains the result spinor (12 regs) */
  double _Complex rs00, rs01, rs02, rs10, rs11, rs12, rs20, rs21, rs22, rs30, rs31, rs32;
  #pragma disjoint(*Ustart, *phi, *l, *l, *U, *s)

        ix = 0 + i*8;
        s = k + i;
        U = Ustart + i*4;

      _prefetch_halfspinor(phi[ix+4]);
      _bgl_load_rs0((*s).s0);
      _bgl_load_rs1((*s).s1);
      _bgl_load_rs2((*s).s2);
      _bgl_load_rs3((*s).s3);
      s++;
      _prefetch_spinor(s);
      /*********************** direction +0 ************************/
      _prefetch_su3(U+1);
      _bgl_vector_add_rs2_to_rs0_reg0();
      _bgl_vector_add_rs3_to_rs1_reg1();

      _bgl_su3_multiply_double((*U));
      _bgl_vector_cmplx_mul_double(ka0);
      /* result is now in regx0, regx1, regx2 , x=0,1 */

      _bgl_store_reg0_up((*phi[ix]).s0);
      _bgl_store_reg1_up((*phi[ix]).s1);
      U++;
      ix++;


      /*********************** direction -0 ************************/
      _prefetch_halfspinor(phi[ix+4]);
      _bgl_vector_sub_rs2_from_rs0_reg0();
      _bgl_vector_sub_rs3_from_rs1_reg1();

      _bgl_store_reg0((*phi[ix]).s0);
      _bgl_store_reg1((*phi[ix]).s1);
      ix++;


      /*********************** direction +1 ************************/
      _prefetch_halfspinor(phi[ix+4]);
      _prefetch_su3(U+1);
      _bgl_vector_i_mul_add_rs3_to_rs0_reg0();
      _bgl_vector_i_mul_add_rs2_to_rs1_reg1();

      _bgl_su3_multiply_double((*U));
      _bgl_vector_cmplx_mul_double(ka1);

      _bgl_store_reg0_up((*phi[ix]).s0);
      _bgl_store_reg1_up((*phi[ix]).s1);
      ix++;
      U++;

      /*********************** direction -1 ************************/
      _prefetch_halfspinor(phi[ix+4]);
      _bgl_vector_i_mul_sub_rs3_from_rs0_reg0();
      _bgl_vector_i_mul_sub_rs2_from_rs1_reg1();

      _bgl_store_reg0((*phi[ix]).s0);
      _bgl_store_reg1((*phi[ix]).s1);
      ix++;


      /*********************** direction +2 ************************/
      _prefetch_halfspinor(phi[ix+4]);
      _prefetch_su3(U+1);

      _bgl_vector_add_rs3_to_rs0_reg0();
      _bgl_vector_sub_rs2_from_rs1_reg1();

      _bgl_su3_multiply_double((*U));
      _bgl_vector_cmplx_mul_double(ka2);

      _bgl_store_reg0_up((*phi[ix]).s0);
      _bgl_store_reg1_up((*phi[ix]).s1);
      ix++;
      U++;

      /*********************** direction -2 ************************/
      _prefetch_halfspinor(phi[ix+4]);
      _bgl_vector_sub_rs3_from_rs0_reg0();
      _bgl_vector_add_rs2_to_rs1_reg1();

      _bgl_store_reg0((*phi[ix]).s0);
      _bgl_store_reg1((*phi[ix]).s1);
      ix++;

      /*********************** direction +3 ************************/
      _prefetch_halfspinor(phi[ix+4]);
      _prefetch_su3(U+1);

      _bgl_vector_i_mul_add_rs2_to_rs0_reg0();
      _bgl_vector_i_mul_sub_rs3_from_rs1_reg1();

      _bgl_su3_multiply_double((*U));
      _bgl_vector_cmplx_mul_double(ka3);

      _bgl_store_reg0_up((*phi[ix]).s0);
      _bgl_store_reg1_up((*phi[ix]).s1);
      ix++;
      U++;

      /*********************** direction -3 ************************/
      _prefetch_halfspinor(phi[ix+4]);
      _bgl_vector_i_mul_sub_rs2_from_rs0_reg0();
      _bgl_vector_i_mul_add_rs3_to_rs1_reg1();

      _bgl_store_reg0((*phi[ix]).s0);
      _bgl_store_reg1((*phi[ix]).s1);
      ix++;

      /************************ end of loop ************************/
    }

#    if (defined MPI && !defined _NO_COMM)
    xchange_halffield();
#    endif

    //s = l;
    phi = NBPointer[2 + ieo];
    _prefetch_halfspinor(phi[0]);
    if(ieo == 0) {
      Ustart = g_gauge_field_copy[1][0];
    }
    else {
      Ustart = g_gauge_field_copy[0][0];
    }
    _prefetch_su3(Ustart);

    /* Now we sum up and expand to a full spinor */
    //ix = 0;
    /*   _prefetch_spinor_for_store(s); */
    #pragma omp for
    for(int i = 0; i < (VOLUME)/2; i++){
                  int ix;
  su3 * restrict U ALIGN;
  spinor * restrict s ALIGN;
  halfspinor32 * restrict * phi32 ALIGN;
  /* We have 32 registers available */
  double _Complex reg00, reg01, reg02, reg03, reg04, reg05;
  double _Complex reg10, reg11, reg12, reg13, reg14, reg15;
  /* For the gauge field, reuse the first three!*/
  double _Complex u00, u01, u02, u10, u11, u12;
  double _Complex reg20, reg21;
  /* The following contains the result spinor (12 regs) */
  double _Complex rs00, rs01, rs02, rs10, rs11, rs12, rs20, rs21, rs22, rs30, rs31, rs32;
  #pragma disjoint(*Ustart, *phi, *l, *l, *U, *s)

     ix = 0 + i*8;
     s = l;
     U = Ustart + i*4;

      /* This causes a lot of trouble, do we understand this? */
      _prefetch_halfspinor(phi[ix+3]);
      /*********************** direction +0 ************************/
      _bgl_load_rs0((*phi[ix]).s0);
      rs20 = rs00;
      rs21 = rs01;
      rs22 = rs02;
      _bgl_load_rs1((*phi[ix]).s1);
      rs30 = rs10;
      rs31 = rs11;
      rs32 = rs12;
      ix++;
      /*********************** direction -0 ************************/
      _prefetch_halfspinor(phi[ix+3]);
      _prefetch_su3(U+1);

      _bgl_load_reg0((*phi[ix]).s0);
      _bgl_load_reg1((*phi[ix]).s1);

      _bgl_su3_inverse_multiply_double((*U));
      _bgl_vector_cmplxcg_mul_double(ka0);
      /* result is in the upper parts of reg0? and reg1? */

      _bgl_add_to_rs0_reg0();
      _bgl_sub_from_rs2_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_sub_from_rs3_reg1();
      U++;
      ix++;

      /*********************** direction +1 ************************/
      _prefetch_halfspinor(phi[ix+3]);
      _bgl_load_reg0_up((*phi[ix]).s0);
      _bgl_load_reg1_up((*phi[ix]).s1);

      _bgl_add_to_rs0_reg0();
      _bgl_i_mul_sub_from_rs3_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_i_mul_sub_from_rs2_reg1();
      ix++;

      /*********************** direction -1 ************************/
      _prefetch_halfspinor(phi[ix+3]);
      _prefetch_su3(U+1);

      _bgl_load_reg0((*phi[ix]).s0);
      _bgl_load_reg1((*phi[ix]).s1);

      _bgl_su3_inverse_multiply_double((*U));
      _bgl_vector_cmplxcg_mul_double(ka1);

      _bgl_add_to_rs0_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_i_mul_add_to_rs3_reg0();
      _bgl_i_mul_add_to_rs2_reg1();
      U++;
      ix++;

      /*********************** direction +2 ************************/
      _prefetch_halfspinor(phi[ix+3]);
      _bgl_load_reg0_up((*phi[ix]).s0);
      _bgl_load_reg1_up((*phi[ix]).s1);

      _bgl_add_to_rs0_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_sub_from_rs2_reg1();
      _bgl_add_to_rs3_reg0();
      ix++;

      /*********************** direction -2 ************************/

      _prefetch_halfspinor(phi[ix+3]);
      _prefetch_su3(U+1);

      _bgl_load_reg0((*phi[ix]).s0);
      _bgl_load_reg1((*phi[ix]).s1);

      _bgl_su3_inverse_multiply_double((*U));
      _bgl_vector_cmplxcg_mul_double(ka2);

      _bgl_add_to_rs0_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_add_to_rs2_reg1();
      _bgl_sub_from_rs3_reg0();
      U++;
      ix++;

      /*********************** direction +3 ************************/

      _prefetch_halfspinor(phi[ix+3]);
      _bgl_load_reg0_up((*phi[ix]).s0);
      _bgl_load_reg1_up((*phi[ix]).s1);

      _bgl_add_to_rs0_reg0();
      _bgl_add_to_rs1_reg1();
      _bgl_i_mul_sub_from_rs2_reg0();
      _bgl_i_mul_add_to_rs3_reg1();
      ix++;

      /*********************** direction -3 ************************/
      _prefetch_spinor(s);
      _prefetch_halfspinor(phi[ix+3]);
      _prefetch_su3(U+1);

      _bgl_load_reg0((*phi[ix]).s0);
      _bgl_load_reg1((*phi[ix]).s1);

      _bgl_su3_inverse_multiply_double((*U));
      _bgl_vector_cmplxcg_mul_double(ka3);

      _bgl_add_to_rs0_reg0();
      _bgl_store_rs0((*s).s0);
      _bgl_i_mul_add_to_rs2_reg0();
      _bgl_store_rs2((*s).s2);

      _bgl_add_to_rs1_reg1();
      _bgl_store_rs1((*s).s1);
      _bgl_i_mul_sub_from_rs3_reg1();
      _bgl_store_rs3((*s).s3);

      U++;
      ix++;
      s++;
    }
}

#elif (!defined _USE_HALFSPINOR && defined BGL && defined XLC && !defined NOBGP && defined OMP)
#undef HopVerMsg
#define HopVerMsg printf("Hopping_Matrix full_spinor bgl OpenMP edition");
/**********************************
 *
 * Blue Gene/L Version
 *
 * Author: Carsten Urbach
 *
 **********************************/

/* 6. */
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k){
 int ioff,ioff2;

#pragma disjoint(*l, *k)

#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge();
  }
#endif

#    if (defined MPI && !(defined _NO_COMM))
  xchange_field(k, ieo);
#    endif

  if(ieo == 0){
    ioff = 0;
  }
  else{
    ioff = (VOLUME+RAND)/2;
  }
  ioff2 = (VOLUME+RAND)/2-ioff;

  /**************** loop over all lattice sites ******************/
#pragma omp for firstprivate(ioff, ioff2)
  for(int icx = ioff; icx < (VOLUME/2+ioff); icx++){
  int icy,icz;
  int ix,iy,iz;
  su3 * restrict up ALIGN;
  su3 * restrict um ALIGN;
  spinor * restrict sp ALIGN;
  spinor * restrict sm ALIGN;
  spinor * restrict rn ALIGN;
  /* We have 32 registers available */
  double _Complex reg00, reg01, reg02, reg03, reg04, reg05;
  double _Complex reg10, reg11, reg12, reg13, reg14, reg15;
  /* For the gauge field, reuse the first three!*/
  double _Complex u00, u01, u02, u10, u11, u12;
  double _Complex reg20, reg21;
  /* The following contains the result spinor (12 regs) */
  double _Complex rs00, rs01, rs02, rs10, rs11, rs12, rs20, rs21, rs22, rs30, rs31, rs32;

#pragma disjoint(*sp, *sm, *rn, *up, *um, *l)
  __alignx(16,l);
  __alignx(16,k);

  ix=g_eo2lexic[icx];
  iy=g_iup[ix][0];
  icy=g_lexic2eosub[iy];

  sp=k+icy;

#    if ((defined _GAUGE_COPY))
  up=&g_gauge_field_copy[icx][0];
#    else
  up=&g_gauge_field[ix][0];
#    endif

    rn=l+(icx-ioff);
    ix=g_eo2lexic[icx];
    /*********************** direction +0 ************************/
    iy=g_idn[ix][0];
    icy=g_lexic2eosub[iy];
#    if (!defined _GAUGE_COPY)
    um=&g_gauge_field[iy][0];
#    else
    um=up+1;
#    endif
    _prefetch_su3(um);
    sm=k+icy;
    _prefetch_spinor(sm);

    _bgl_load_reg0((*sp).s0);
    _bgl_load_reg1((*sp).s1);
    _bgl_load_reg0_up((*sp).s2);
    _bgl_load_reg1_up((*sp).s3);
    _bgl_vector_add_reg0();
    _bgl_vector_add_reg1();
    /* result is now in regx0, regx1, regx2 x = 0,1 */

    _bgl_su3_multiply_double((*up));
    _bgl_vector_cmplx_mul_double(ka0);
    _bgl_store_reg0_up_rs0();
    _bgl_store_reg0_up_rs2();
    _bgl_store_reg1_up_rs1();
    _bgl_store_reg1_up_rs3();

    /*********************** direction -0 ************************/

    iy=g_iup[ix][1];
    icy=g_lexic2eosub[iy];

#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    _prefetch_su3(up);
    sp=k+icy;
    _prefetch_spinor(sp);

    _bgl_load_reg0((*sm).s0);
    _bgl_load_reg1((*sm).s1);
    _bgl_load_reg0_up((*sm).s2);
    _bgl_load_reg1_up((*sm).s3);
    _bgl_vector_sub_reg0();
    _bgl_vector_sub_reg1();

    _bgl_su3_inverse_multiply_double((*um));
    _bgl_vector_cmplxcg_mul_double(ka0);

    _bgl_add_to_rs0_reg0();
    _bgl_sub_from_rs2_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_sub_from_rs3_reg1();

    /*********************** direction +1 ************************/

    iy=g_idn[ix][1];
    icy=g_lexic2eosub[iy];

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][1];
#    else
    um = up+1;
#    endif
    _prefetch_su3(um);
    sm=k+icy;
    _prefetch_spinor(sm);

    _bgl_load_reg0((*sp).s0);
    _bgl_load_reg1((*sp).s1);
    _bgl_load_reg0_up((*sp).s3);
    _bgl_load_reg1_up((*sp).s2);
    _bgl_vector_i_mul_add_reg0();
    _bgl_vector_i_mul_add_reg1();

    _bgl_su3_multiply_double((*up));
    _bgl_vector_cmplx_mul_double(ka1);

    _bgl_add_to_rs0_reg0();
    _bgl_i_mul_sub_from_rs3_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_i_mul_sub_from_rs2_reg1();

    /*********************** direction -1 ************************/

    iy=g_iup[ix][2];
    icy=g_lexic2eosub[iy];

#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    _prefetch_su3(up);
    sp=k+icy;
    _prefetch_spinor(sp);

    _bgl_load_reg0((*sm).s0);
    _bgl_load_reg1((*sm).s1);
    _bgl_load_reg0_up((*sm).s3);
    _bgl_load_reg1_up((*sm).s2);
    _bgl_vector_i_mul_sub_reg0();
    _bgl_vector_i_mul_sub_reg1();

    _bgl_su3_inverse_multiply_double((*um));
    _bgl_vector_cmplxcg_mul_double(ka1);

    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_i_mul_add_to_rs3_reg0();
    _bgl_i_mul_add_to_rs2_reg1();

    /*********************** direction +2 ************************/

    iy=g_idn[ix][2];
    icy=g_lexic2eosub[iy];

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][2];
#    else
    um= up+1;
#    endif
    _prefetch_su3(um);
    sm=k+icy;
    _prefetch_spinor(sm);

    _bgl_load_reg0((*sp).s0);
    _bgl_load_reg1((*sp).s1);
    _bgl_load_reg1_up((*sp).s2);
    _bgl_load_reg0_up((*sp).s3);
    _bgl_vector_add_reg0();
    _bgl_vector_sub_reg1();

    _bgl_su3_multiply_double((*up));
    _bgl_vector_cmplx_mul_double(ka2);

    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_sub_from_rs2_reg1();
    _bgl_add_to_rs3_reg0();


    /*********************** direction -2 ************************/

    iy=g_iup[ix][3];
    icy=g_lexic2eosub[iy];

#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    _prefetch_su3(up);
    sp=k+icy;
    _prefetch_spinor(sp);

    _bgl_load_reg0((*sm).s0);
    _bgl_load_reg1((*sm).s1);
    _bgl_load_reg1_up((*sm).s2);
    _bgl_load_reg0_up((*sm).s3);
    _bgl_vector_sub_reg0();
    _bgl_vector_add_reg1();

    _bgl_su3_inverse_multiply_double((*um));
    _bgl_vector_cmplxcg_mul_double(ka2);

    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_add_to_rs2_reg1();
    _bgl_sub_from_rs3_reg0();

    /*********************** direction +3 ************************/

    iy=g_idn[ix][3];
    icy=g_lexic2eosub[iy];

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][3];
#    else
    um=up+1;
#    endif
    _prefetch_su3(um);
    sm=k+icy;
    _prefetch_spinor(sm);

    _bgl_load_reg0((*sp).s0);
    _bgl_load_reg1((*sp).s1);
    _bgl_load_reg0_up((*sp).s2);
    _bgl_load_reg1_up((*sp).s3);
    _bgl_vector_i_mul_add_reg0();
    _bgl_vector_i_mul_sub_reg1();

    _bgl_su3_multiply_double((*up));
    _bgl_vector_cmplx_mul_double(ka3);

    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_i_mul_sub_from_rs2_reg0();
    _bgl_i_mul_add_to_rs3_reg1();

    /*********************** direction -3 ************************/

    icz=icx+1;
    if(icz==((VOLUME+RAND)/2+ioff)) icz=ioff;
    iz=g_eo2lexic[icz];
    iy=g_iup[iz][0]; icy=g_lexic2eosub[iy];



#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up=&g_gauge_field[iz][0];
#    endif
    _prefetch_su3(up);
    sp=k+icy;
    _prefetch_spinor(sp);

    _bgl_load_reg0((*sm).s0);
    _bgl_load_reg1((*sm).s1);
    _bgl_load_reg0_up((*sm).s2);
    _bgl_load_reg1_up((*sm).s3);
    _bgl_vector_i_mul_sub_reg0();
    _bgl_vector_i_mul_add_reg1();

    _bgl_su3_inverse_multiply_double((*um));
    _bgl_vector_cmplxcg_mul_double(ka3);



    _bgl_add_to_rs0_reg0();
    _bgl_store_rs0((*rn).s0);
    _bgl_i_mul_add_to_rs2_reg0();
    _bgl_store_rs2((*rn).s2);

    _bgl_add_to_rs1_reg1();
    _bgl_store_rs1((*rn).s1);
    _bgl_i_mul_sub_from_rs3_reg1();
    _bgl_store_rs3((*rn).s3);

    /************************ end of loop ************************/
  }
}


#elif (defined _USE_HALFSPINOR && defined OMP)

#undef HopVerMsg
#define HopVerMsg printf("Hopping_Matrix half_spinor generic OpenMP edition\n");
/* 3. */
/* l output , k input*/
/* for ieo=0, k resides on  odd sites and l on even sites */
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k){
    su3 * restrict Ustart ALIGN;
    halfspinor * restrict * phi ALIGN;

#ifdef XLC
#pragma disjoint(*l, *k, *Ustart, *phi)
#endif

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
  //s = k;

  if(ieo == 0) {
    Ustart = g_gauge_field_copy[0][0];
  }
  else {
    Ustart = g_gauge_field_copy[1][0];
  }

    phi = NBPointer[ieo];

    /**************** loop over all lattice sites ****************/
    /* #pragma ivdep*/
    #pragma omp for
    for(int i = 0; i < (VOLUME)/2; i++){
          int ix;
  su3 * restrict U ALIGN;
  spinor * restrict s ALIGN;
  spinor rs;
  su3_vector psi, chi, psi2, chi2;
  halfspinor32 * restrict * phi32 ALIGN;
  #ifdef XLC
#pragma disjoint(*l, *k, *U, *s, *Ustart)
#endif

        ix = 0 + i*8;
        s = k + i;
        U = Ustart + i*4;

      _vector_assign(rs.s0, (*s).s0);
      _vector_assign(rs.s1, (*s).s1);
      _vector_assign(rs.s2, (*s).s2);
      _vector_assign(rs.s3, (*s).s3);
      s++;
      /*********************** direction +0 ************************/

      _vector_add(psi, rs.s0, rs.s2);
      _vector_add(psi2, rs.s1, rs.s3);
      _su3_multiply(chi,(*U),psi);
      _su3_multiply(chi2,(*U),psi2);
      _complex_times_vector((*phi[ix]).s0, ka0, chi);
      _complex_times_vector((*phi[ix]).s1, ka0, chi2);

      U++;
      ix++;

      /*********************** direction -0 ************************/

      _vector_sub((*phi[ix]).s0, rs.s0, rs.s2);
      _vector_sub((*phi[ix]).s1, rs.s1, rs.s3);

      ix++;

      /*********************** direction +1 ************************/

      _vector_i_add(psi, rs.s0, rs.s3);
      _vector_i_add(psi2, rs.s1, rs.s2);
      _su3_multiply(chi, (*U), psi);
      _su3_multiply(chi2, (*U), psi2);
      _complex_times_vector((*phi[ix]).s0, ka1, chi);
      _complex_times_vector((*phi[ix]).s1, ka1, chi2);

      U++;
      ix++;

      /*********************** direction -1 ************************/

      _vector_i_sub((*phi[ix]).s0, rs.s0, rs.s3);
      _vector_i_sub((*phi[ix]).s1, rs.s1, rs.s2);

      ix++;
      /*********************** direction +2 ************************/

      _vector_add(psi, rs.s0, rs.s3);
      _vector_sub(psi2, rs.s1, rs.s2);
      _su3_multiply(chi,(*U),psi);
      _su3_multiply(chi2,(*U),psi2);
      _complex_times_vector((*phi[ix]).s0, ka2, chi);
      _complex_times_vector((*phi[ix]).s1, ka2, chi2);

      U++;
      ix++;

      /*********************** direction -2 ************************/

      _vector_sub((*phi[ix]).s0, rs.s0, rs.s3);
      _vector_add((*phi[ix]).s1, rs.s1, rs.s2);
      ix++;

      /*********************** direction +3 ************************/

      _vector_i_add(psi, rs.s0, rs.s2);
      _vector_i_sub(psi2, rs.s1, rs.s3);
      _su3_multiply(chi, (*U), psi);
      _su3_multiply(chi2,(*U),psi2);
      _complex_times_vector((*phi[ix]).s0, ka3, chi);
      _complex_times_vector((*phi[ix]).s1, ka3, chi2);

      U++;
      ix++;
      /*********************** direction -3 ************************/

      _vector_i_sub((*phi[ix]).s0, rs.s0, rs.s2);
      _vector_i_add((*phi[ix]).s1, rs.s1, rs.s3);

      ix++;
      /************************ end of loop ************************/
    }

#    if (defined MPI && !defined _NO_COMM)
    xchange_halffield();
#    endif
    //s = l;
    phi = NBPointer[2 + ieo];
    if(ieo == 0) {
      Ustart = g_gauge_field_copy[1][0];
    }
    else {
      Ustart = g_gauge_field_copy[0][0];
    }

    //ix = 0;
    /* #pragma ivdep */
    #pragma omp for
    for(int i = 0; i < (VOLUME)/2; i++){
                  int ix;
  su3 * restrict U ALIGN;
  spinor * restrict s ALIGN;
  spinor rs;
  su3_vector psi, chi, psi2, chi2;
  halfspinor32 * restrict * phi32 ALIGN;
  #ifdef XLC
#pragma disjoint(*l, *k, *U, *s, *Ustart)
#endif

        ix = 0 + i*8;
        s = l;
        U = Ustart + i*4;

      /*********************** direction +0 ************************/
      _vector_assign(rs.s0, (*phi[ix]).s0);
      _vector_assign(rs.s2, (*phi[ix]).s0);
      _vector_assign(rs.s1, (*phi[ix]).s1);
      _vector_assign(rs.s3, (*phi[ix]).s1);
      ix++;
      /*********************** direction -0 ************************/
      _su3_inverse_multiply(chi,(*U),(*phi[ix]).s0);
      _su3_inverse_multiply(chi2,(*U),(*phi[ix]).s1);
      _complexcjg_times_vector(psi,ka0,chi);
      _complexcjg_times_vector(psi2,ka0,chi2);
      _vector_add_assign(rs.s0, psi);
      _vector_sub_assign(rs.s2, psi);
      _vector_add_assign(rs.s1, psi2);
      _vector_sub_assign(rs.s3, psi2);
      ix++;
      U++;
      /*********************** direction +1 ************************/

      _vector_add_assign(rs.s0, (*phi[ix]).s0);
      _vector_i_sub_assign(rs.s3, (*phi[ix]).s0);

      _vector_add_assign(rs.s1, (*phi[ix]).s1);
      _vector_i_sub_assign(rs.s2, (*phi[ix]).s1);

      ix++;
      /*********************** direction -1 ************************/

      _su3_inverse_multiply(chi,(*U), (*phi[ix]).s0);
      _su3_inverse_multiply(chi2, (*U), (*phi[ix]).s1);
      _complexcjg_times_vector(psi,ka1,chi);
      _complexcjg_times_vector(psi2,ka1,chi2);
      _vector_add_assign(rs.s0, psi);
      _vector_i_add_assign(rs.s3, psi);
      _vector_add_assign(rs.s1, psi2);
      _vector_i_add_assign(rs.s2, psi2);

      U++;
      ix++;

      /*********************** direction +2 ************************/

      _vector_add_assign(rs.s0, (*phi[ix]).s0);
      _vector_add_assign(rs.s3, (*phi[ix]).s0);

      _vector_add_assign(rs.s1, (*phi[ix]).s1);
      _vector_sub_assign(rs.s2, (*phi[ix]).s1);

      ix++;
      /*********************** direction -2 ************************/

      _su3_inverse_multiply(chi,(*U), (*phi[ix]).s0);
      _su3_inverse_multiply(chi2, (*U), (*phi[ix]).s1);
      _complexcjg_times_vector(psi,ka2,chi);
      _complexcjg_times_vector(psi2,ka2,chi2);
      _vector_add_assign(rs.s0, psi);
      _vector_sub_assign(rs.s3, psi);
      _vector_add_assign(rs.s1, psi2);
      _vector_add_assign(rs.s2, psi2);

      U++;
      ix++;
      /*********************** direction +3 ************************/

      _vector_add_assign(rs.s0, (*phi[ix]).s0);
      _vector_i_sub_assign(rs.s2, (*phi[ix]).s0);

      _vector_add_assign(rs.s1, (*phi[ix]).s1);
      _vector_i_add_assign(rs.s3, (*phi[ix]).s1);

      ix++;

      /*********************** direction -3 ************************/

      _su3_inverse_multiply(chi,(*U), (*phi[ix]).s0);
      _su3_inverse_multiply(chi2, (*U), (*phi[ix]).s1);
      _complexcjg_times_vector(psi,ka3,chi);
      _complexcjg_times_vector(psi2,ka3,chi2);
      _vector_add((*s).s0, rs.s0, psi);
      _vector_i_add((*s).s2, rs.s2, psi);
      _vector_add((*s).s1, rs.s1, psi2);
      _vector_i_sub((*s).s3, rs.s3, psi2);

      U++;
      ix++;
      s++;
    }
}

#endif
