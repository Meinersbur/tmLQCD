/**********************************************************************
 *
 *
 * Copyright (C) 2012 Carsten Urbach, Bartosz Kostrzewa
 *
 * This file is based on an implementation of the Dirac operator 
 * written by Martin Luescher, modified by Martin Hasenbusch in 2002 
 * and modified and extended by Carsten Urbach from 2003-2008
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
 *
 **********************************************************************/

  int ioff;
  int * hi;
  su3 * restrict ALIGN up;
  su3 * restrict ALIGN um;
  spinor * restrict ALIGN sp;
  spinor * restrict ALIGN sm;
  spinor * restrict ALIGN rn;
  su3_vector psi_v0, chi_v0, rho_v0;
  su3_vector psi_v1, chi_v1, rho_v1;
  
#ifdef XLC
#  pragma disjoint(*sp, *sm, *rn, *up, *um, *l)
#endif
  _declare_regs();

  if(ieo == 0){
    ioff = 0;
  } 
  else{
    ioff = (VOLUME+RAND)/2;
  }

#ifndef OMP
  hi = &g_hi[16*ioff];

#  if ((defined _GAUGE_COPY))
  up=&g_gauge_field_copy[ioff][0];
#  else
  up=&g_gauge_field[(*hi)][0];
#  endif
  hi++;
  sp=k+(*hi);
  hi++;
#endif

  /**************** loop over all lattice sites ******************/
#ifdef OMP
#  pragma omp for
#endif
  for(int icx = ioff; icx < (VOLUME/2+ioff); icx++){
#ifdef OMP
    hi = &g_hi[16*icx];
#  if ((defined _GAUGE_COPY))
    up=&g_gauge_field_copy[icx][0];
#  else
    up=&g_gauge_field[(*hi)][0];
#  endif
    hi++;
    sp=k+(*hi);
    hi++;
#endif
    rn=l+(icx-ioff);
#ifdef _TM_SUB_HOP
    pn=p+(icx-ioff);
#endif

    const int ix/*lexic*/=g_eo2lexic[icx];
	const int t = g_coord[ix][0];
	const int x = g_coord[ix][1];
	const int y = g_coord[ix][2];
	const int z = g_coord[ix][3];

    /*********************** direction +t ************************/
#    if (!defined _GAUGE_COPY)
    um=&g_gauge_field[(*hi)][0]; 
#    else
    um=up+1;
#    endif
    hi++;
    sm=k+(*hi);
    hi+=2;

    _vector_add(psi_v0,sp->s0,sp->s2); // 6
    _su3_multiply(chi_v0,(*up),psi_v0); // 3*3*6 + 2*3*2
    _complex_times_vector(rho_v0,ka0,chi_v0);
    _vector_assign(temp.s0,rho_v0);
    _vector_assign(temp.s2,rho_v0);

    _vector_add(psi_v1,sp->s1,sp->s3); // 6
    _su3_multiply(chi_v1,(*up),psi_v1); // 3*3*6 + 2*3*2
    _complex_times_vector(rho_v1,ka0,chi_v1);
    _vector_assign(temp.s1,rho_v1);
    _vector_assign(temp.s3,rho_v1);

    bgq_setdesc(BGQREF_TUP_SOURCE, "BGQREF_TUP_SOURCE");
    bgq_setdesc(BGQREF_TUP, "BGQREF_TUP");
    bgq_setdesc(BGQREF_TUP_GAUGE, "BGQREF_TUP_GAUGE");
    bgq_setdesc(BGQREF_TUP_WEYL, "BGQREF_TUP_WEYL");
    bgq_setdesc(BGQREF_TUP_KAMUL, "BGQREF_TUP_KAMUL");
    bgq_setdesc(BGQREF_TUP_RECV, "BGQREF_TUP_RECV");
    bgq_setdesc(BGQREF_TUP_ACCUM, "BGQREF_TUP_ACCUM");
    bgq_setrefvalue(t,x,y,z, BGQREF_TUP_SOURCE, sp->s0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_TUP, psi_v0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_TUP_GAUGE, up->c00);
    bgq_setrefvalue(t,x,y,z, BGQREF_TUP_WEYL, chi_v0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_TUP_KAMUL, rho_v0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_TUP_RECV, rho_v0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_TUP_ACCUM, temp.s0.c0);

    /*********************** direction -t ************************/
#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    sp=k+(*hi);
    hi++;

    _vector_sub(psi_v0,sm->s0,sm->s2);
    _su3_inverse_multiply(chi_v0,(*um),psi_v0);
    _complexcjg_times_vector(rho_v0,ka0,chi_v0);
    _vector_add_assign(temp.s0,rho_v0); // 6
    _vector_sub_assign(temp.s2,rho_v0); // 6

    _vector_sub(psi_v1,sm->s1,sm->s3);
    _su3_inverse_multiply(chi_v1,(*um),psi_v1);
    _complexcjg_times_vector(rho_v1,ka0,chi_v1);
    _vector_add_assign(temp.s1,rho_v1); // 6
    _vector_sub_assign(temp.s3,rho_v1); // 6


    bgq_setdesc(BGQREF_TDOWN_SOURCE, "BGQREF_TDOWN_SOURCE");
    bgq_setdesc(BGQREF_TDOWN, "BGQREF_TDOWN");
    bgq_setdesc(BGQREF_TDOWN_GAUGE, "BGQREF_TDOWN_GAUGE");
    bgq_setdesc(BGQREF_TDOWN_WEYL, "BGQREF_TDOWN_WEYL");
    bgq_setdesc(BGQREF_TDOWN_KAMUL, "BGQREF_TDOWN_KAMUL");
    bgq_setdesc(BGQREF_TDOWN_RECV, "BGQREF_TDOWN_RECV");
    bgq_setdesc(BGQREF_TDOWN_ACCUM, "BGQREF_TDOWN_ACCUM");
    bgq_setrefvalue(t,x,y,z, BGQREF_TDOWN_SOURCE, sm->s0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_TDOWN, psi_v0.c2);
    bgq_setrefvalue(t,x,y,z, BGQREF_TDOWN_GAUGE, um->c00);
    bgq_setrefvalue(t,x,y,z, BGQREF_TDOWN_WEYL, chi_v0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_TDOWN_KAMUL, rho_v0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_TDOWN_RECV, rho_v0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_TDOWN_ACCUM, temp.s0.c0);


    /*********************** direction +x ************************/
#    ifndef _GAUGE_COPY
    um=&g_gauge_field[(*hi)][1]; 
#    else
    um = up+1;
#    endif
    hi++;
    sm=k+(*hi);
    hi+=2;

    _vector_i_add(psi_v0,sp->s0,sp->s3);
    _su3_multiply(chi_v0,(*up),psi_v0);
    _complex_times_vector(rho_v0,ka1,chi_v0);
    _vector_add_assign(temp.s0,rho_v0);
    _vector_i_sub_assign(temp.s3,rho_v0);

    _vector_i_add(psi_v1,sp->s1,sp->s2);
    _su3_multiply(chi_v1,(*up),psi_v1);
    _complex_times_vector(rho_v1,ka1,chi_v1);
    _vector_add_assign(temp.s1,rho_v1);
    _vector_i_sub_assign(temp.s2,rho_v1);

    bgq_setdesc(BGQREF_XUP_SOURCE, "BGQREF_XUP_SOURCE");
    bgq_setdesc(BGQREF_XUP, "BGQREF_XUP");
    bgq_setdesc(BGQREF_XUP_GAUGE, "BGQREF_XUP_GAUGE");
    bgq_setdesc(BGQREF_XUP_WEYL, "BGQREF_XUP_WEYL");
    bgq_setdesc(BGQREF_XUP_KAMUL, "BGQREF_XUP_KAMUL");
    bgq_setdesc(BGQREF_XUP_RECV, "BGQREF_XUP_RECV");
    bgq_setdesc(BGQREF_XUP_ACCUM, "BGQREF_XUP_ACCUM");
    bgq_setrefvalue(t,x,y,z, BGQREF_XUP_SOURCE, sp->s0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_XUP, psi_v0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_XUP_GAUGE, up->c00);
    bgq_setrefvalue(t,x,y,z, BGQREF_XUP_WEYL, chi_v0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_XUP_KAMUL, rho_v0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_XUP_RECV, rho_v0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_XUP_ACCUM, temp.s0.c0);

    /*********************** direction -x ************************/
#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    sp=k+(*hi);
    hi++;

    _vector_i_sub(psi_v0,sm->s0,sm->s3);
    _su3_inverse_multiply(chi_v0,(*um),psi_v0);
    _complexcjg_times_vector(rho_v0,ka1,chi_v0);
    _vector_add_assign(temp.s0,rho_v0);
    _vector_i_add_assign(temp.s3,rho_v0);

    _vector_i_sub(psi_v1,sm->s1,sm->s2);
    _su3_inverse_multiply(chi_v1,(*um),psi_v1);
    _complexcjg_times_vector(rho_v1,ka1,chi_v1);
    _vector_add_assign(temp.s1,rho_v1);
    _vector_i_add_assign(temp.s2,rho_v1);

    bgq_setdesc(BGQREF_XDOWN, "BGQREF_XDOWN");
    bgq_setdesc(BGQREF_XDOWN_GAUGE, "BGQREF_XDOWN_GAUGE");
    bgq_setdesc(BGQREF_XDOWN_WEYL, "BGQREF_XDOWN_WEYL");
    bgq_setdesc(BGQREF_XDOWN_KAMUL, "BGQREF_XDOWN_KAMUL");
    bgq_setdesc(BGQREF_XDOWN_RECV, "BGQREF_XDOWN_RECV");
    bgq_setdesc(BGQREF_XDOWN_ACCUM, "BGQREF_XDOWN_ACCUM");
    bgq_setrefvalue(t,x,y,z, BGQREF_XDOWN, psi_v0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_XDOWN_GAUGE, um->c00);
    bgq_setrefvalue(t,x,y,z, BGQREF_XDOWN_WEYL, chi_v0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_XDOWN_KAMUL, rho_v0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_XDOWN_RECV, rho_v0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_XDOWN_ACCUM, temp.s0.c0);

    /*********************** direction +y ************************/
#    ifndef _GAUGE_COPY
    um=&g_gauge_field[(*hi)][2]; 
#    else
    um= up+1;
#    endif
    hi++;
    sm=k+(*hi);
    hi+=2;

    _vector_add(psi_v0,sp->s0,sp->s3);
    _su3_multiply(chi_v0,(*up),psi_v0);
    _complex_times_vector(rho_v0,ka2,chi_v0);
    _vector_add_assign(temp.s0,rho_v0);
    _vector_add_assign(temp.s3,rho_v0);

    _vector_sub(psi_v1,sp->s1,sp->s2);
    _su3_multiply(chi_v1,(*up),psi_v1);
    _complex_times_vector(rho_v1,ka2,chi_v1);
    _vector_add_assign(temp.s1,rho_v1);
    _vector_sub_assign(temp.s2,rho_v1);

    bgq_setdesc(BGQREF_YUP_RECV, "BGQREF_YUP_RECV");
    bgq_setdesc(BGQREF_YUP_ACCUM, "BGQREF_YUP_ACCUM");
    bgq_setrefvalue(t,x,y,z, BGQREF_YUP_RECV, rho_v0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_YUP_ACCUM, temp.s0.c0);

    /*********************** direction -y ************************/
#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    sp=k+(*hi);
    hi++;

    _hop_y_m();

    bgq_setdesc(BGQREF_YDOWN_ACCUM, "BGQREF_YDOWN_ACCUM");
    bgq_setrefvalue(t,x,y,z, BGQREF_YDOWN_ACCUM, temp.s0.c0);

    /*********************** direction +z ************************/
#    ifndef _GAUGE_COPY
    um=&g_gauge_field[(*hi)][3]; 
#    else
    um=up+1;
#    endif
    hi++;
    sm=k+(*hi);
    hi++;

    _vector_i_add(psi_v0,sp->s0,sp->s2);
    _su3_multiply(chi_v0,(*up),psi_v0);
    _complex_times_vector(rho_v0,ka3,chi_v0);
    _vector_add_assign(temp.s0,rho_v0);
    _vector_i_sub_assign(temp.s2,rho_v0);

    _vector_i_sub(psi_v1,sp->s1,sp->s3);
    _su3_multiply(chi_v0,(*up),psi_v1);
    _complex_times_vector(rho_v1,ka3,chi_v0);
    _vector_add_assign(temp.s1,rho_v1);
    _vector_i_add_assign(temp.s3,rho_v1);

    bgq_setdesc(BGQREF_ZUP, "BGQREF_ZUP");
    bgq_setdesc(BGQREF_ZUP_GAUGE, "BGQREF_ZUP_GAUGE");
    bgq_setdesc(BGQREF_ZUP_KAMUL, "BGQREF_ZUP_KAMUL");
    bgq_setdesc(BGQREF_ZUP_RECV, "BGQREF_ZUP_RECV");
    bgq_setdesc(BGQREF_ZUP_ACCUM, "BGQREF_ZUP_ACCUM");
    bgq_setrefvalue(t,x,y,z, BGQREF_ZUP, psi_v0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_ZUP_GAUGE, up->c00);
    bgq_setrefvalue(t,x,y,z, BGQREF_ZUP_KAMUL, rho_v0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_ZUP_RECV, rho_v0.c0);
    bgq_setrefvalue(t,x,y,z, BGQREF_ZUP_ACCUM, temp.s0.c0);

    /*********************** direction -z ************************/
#ifndef OMP
#  if ((defined _GAUGE_COPY))
    up=um+1;
#  else
    up=&g_gauge_field[(*hi)][0];
#  endif
    hi++;
    sp=k+(*hi);
    hi++;
#endif
    _hop_z_m();

#ifdef _MUL_G5_CMPLX
    _hop_mul_g5_cmplx_and_store();
#elif defined _TM_SUB_HOP
    _g5_cmplx_sub_hop_and_g5store();
#else
    _store_res();
#endif

    bgq_setdesc(BGQREF_ACCUM, "BGQREF_ACCUM");
    bgq_setrefvalue(t,x,y,z, BGQREF_ACCUM, temp.s0.c0);
  }
