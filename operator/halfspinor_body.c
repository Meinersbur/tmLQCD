/**********************************************************************
 *
 * Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2012 Carsten Urbach
 *
 * This file is based on an implementation of the Dirac operator 
 * written by Martin Luescher, modified by Martin Hasenbusch in 2002 
 * this is a new version based on the aforementioned implementations
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


int ix;
su3 * restrict U ALIGN;
spinor * restrict s ALIGN;
halfspinor * restrict * phi ALIGN;
halfspinor32 * restrict * phi32 ALIGN;
_declare_hregs();

#ifdef XLC
# pragma disjoint(*l, *k)
# pragma disjoint(*k, *U)
# pragma disjoint(*l, *U)
# pragma disjoint(*U, *s)
# pragma disjoint(*k, *s)
# pragma disjoint(*l, *s)
__alignx(32, l);
__alignx(32, k);
__alignx(32, U);
__alignx(32, s);
#endif


#ifdef _KOJAK_INST
#pragma pomp inst begin(hoppingmatrix)
#endif

#ifndef OMP  
s = k;
_prefetch_spinor(s);
if(ieo == 0) {
  U = g_gauge_field_copy[0][0];
 }
 else {
   U = g_gauge_field_copy[1][0];
 }
_prefetch_su3(U);
#else
if(ieo == 0) {
  u0 = g_gauge_field_copy[0][0];
 }
 else {
   u0 = g_gauge_field_copy[1][0];
 }
#endif

#if (defined SSE2 || defined SSE3)
g_sloppy_precision = 0;
#endif
if(g_sloppy_precision == 1 && g_sloppy_precision_flag == 1) {
  phi32 = NBPointer32[ieo];

#if 0//def SPI

/* compute the surface terms to be sent first */

  surface = g_surface + ieo*SURFACE/2;
  body = g_body + ieo*BODY/2;

#ifdef OMP
#pragma omp for
#endif
  for(unsigned int i = 0; i < SURFACE/2; ++i){
    U=u0+surface[i]*4;
    s=k+surface[i];
    ix=surface[i]*8;
    
    _hop_t_p_pre32();
    U++;
    ix++;
     
    _hop_t_m_pre32();
    ix++;
     
    _hop_x_p_pre32();
    U++;
    ix++;
     
    _hop_x_m_pre32();
    ix++;
     
    _hop_y_p_pre32();
    U++;
    ix++;
     
    _hop_y_m_pre32();
    ix++;
     
    _hop_z_p_pre32();
    U++;
    ix++;
     
    _hop_z_m_pre32();
    }

#else // SPI

/* without interleaving we compute surface and body together */

#ifdef OMP
#pragma omp for
#else
  ix=0;
#endif
  for(unsigned int i = 0; i < (VOLUME)/2; i++){
#ifdef OMP
    U=u0+i*4;
    s=k+i;
    ix=i*8;
#endif
    _hop_t_p_pre32();
    U++;
    ix++;
    
    _hop_t_m_pre32();
    ix++;
    
    _hop_x_p_pre32();
    U++;
    ix++;
    
    _hop_x_m_pre32();
    ix++;
    
    _hop_y_p_pre32();
    U++;
    ix++;
    
    _hop_y_m_pre32();
    ix++;
    
    _hop_z_p_pre32();
    U++;
    ix++;
    
    _hop_z_m_pre32();
    
#ifndef OMP
    s++;
    ix++;
#endif
  }

#endif // SPI
  
#ifdef OMP
#  if 0//def SPI
#pragma omp single nowait
#  else
#pragma omp single
#  endif
  {
#endif
    
#    if (defined MPI && !defined _NO_COMM)
#      if 0//def SPI

     // Initialize the barrier, resetting the hardware.
     int rc = MUSPI_GIBarrierInit ( &GIBarrier, 0 /*comm world class route */);
     if(rc) {
       printf("MUSPI_GIBarrierInit returned rc = %d\n", rc);
       exit(__LINE__);
     }
     // reset the recv counter 
     recvCounter = totalMessageSize/2;
     global_barrier(); // make sure everybody is set recv counter

     //#pragma omp for nowait
     for (unsigned int j = 0; j < spi_num_dirs; j++) {
       descCount[ j ] =
	 msg_InjFifoInject ( injFifoHandle,
			     j,
			     &SPIDescriptors32[j]);
     }
#      else
    xchange_halffield32(); 
#      endif
#    endif
    
#ifdef OMP
  }
#endif

#if 0//def SPI

/* compute the body while the communication is in progress 
   the barrier at the end of this for loop will also catch the
   thread carrying out the "single nowait" above */
   
#ifdef OMP
#pragma omp for
#endif
  for(int i = 0; i < BODY/2; ++i)
  {
    U=u0+body[i]*4;
    s=k+body[i];
    ix=body[i]*8;
    
    _hop_t_p_pre32();
    U++;
    ix++;
     
    _hop_t_m_pre32();
    ix++;
     
    _hop_x_p_pre32();
    U++;
    ix++;
     
    _hop_x_m_pre32();
    ix++;
     
    _hop_y_p_pre32();
    U++;
    ix++;
     
    _hop_y_m_pre32();
    ix++;
     
    _hop_z_p_pre32();
    U++;
    ix++;
     
    _hop_z_m_pre32(); 
  }

  /* wait for the communication to finish */
#ifdef OMP
#pragma omp single
{
#endif
  while ( recvCounter > 0 );
  _bgq_msync();
#ifdef OMP
}
#endif

#endif // SPI

/* move on to compute the result as usual */
  
#ifndef OMP
  s = l;
  if(ieo == 0) {
    U = g_gauge_field_copy[1][0];
  }
  else {
    U = g_gauge_field_copy[0][0];
  }
#else
  if(ieo == 0) {
    u0 = g_gauge_field_copy[1][0];
  }
  else {
    u0 = g_gauge_field_copy[0][0];
  }
#endif
  
  phi32 = NBPointer32[2 + ieo];
  
#ifdef OMP
#pragma omp for
#else
  ix = 0;
#endif
  for(unsigned int i = 0; i < (VOLUME)/2; i++){
#ifdef OMP
    ix=i*8;
    s=l+i;
    U=u0+i*4;
#endif
#ifdef _TM_SUB_HOP
     pn=p+i;
#endif
    _hop_t_p_post32();
    ix++;
    
    _hop_t_m_post32();
    ix++;
    U++;
    
    _hop_x_p_post32();
    ix++;
    
    _hop_x_m_post32();
    U++;
    ix++;
    
    _hop_y_p_post32();
    ix++;
    
    _hop_y_m_post32();
    U++;
    ix++;
    
    _hop_z_p_post32();
    ix++;
    
    _hop_z_m_post32();
    
#ifdef _MUL_G5_CMPLX
    _hop_mul_g5_cmplx_and_store(s);
#elif defined _TM_SUB_HOP
     _g5_cmplx_sub_hop_and_g5store(s);
#else
    _hop_store_post(s);
#endif
    
#ifndef OMP
    U++;
    ix++;
    s++;
#endif
  }
 }
 else {

  phi = NBPointer[ieo];

#if 0//def SPI

/* compute the surface terms to be sent first */

  surface = g_surface + ieo*SURFACE/2;
  body = g_body + ieo*BODY/2;

#ifdef OMP
#pragma omp for
#endif
  for(unsigned int i = 0; i < SURFACE/2; ++i){
    U=u0+surface[i]*4;
    _prefetch_su3(U);
    s=k+surface[i];
    _prefetch_spinor(s);
    ix=surface[i]*8;
    
    _hop_t_p_pre();
    U++;
    ix++;
     
    _hop_t_m_pre();
    ix++;
     
    _hop_x_p_pre();
    U++;
    ix++;
     
    _hop_x_m_pre();
    ix++;
     
    _hop_y_p_pre();
    U++;
    ix++;
     
    _hop_y_m_pre();
    ix++;
     
    _hop_z_p_pre();
    U++;
    ix++;
     
    _hop_z_m_pre();
    }

#else // SPI
  
  /* when not interleaving we compute body and surface together */
   
#ifdef OMP
#pragma omp for
#else
   ix=0;
#endif
   for(unsigned int i = 0; i < (VOLUME)/2; i++){
#ifdef OMP
     s=k+i;
     _prefetch_spinor(s);
     ix=i*8;
     U=u0+i*4;
     _prefetch_su3(U);
#endif
     
     _hop_t_p_pre();
     U++;
     ix++;
     
     _hop_t_m_pre();
     ix++;
     
     _hop_x_p_pre();
     U++;
     ix++;
     
     _hop_x_m_pre();
     ix++;
     
     _hop_y_p_pre();
     U++;
     ix++;
     
     _hop_y_m_pre();
     ix++;
     
     _hop_z_p_pre();
     U++;
     ix++;
     
     _hop_z_m_pre();
     
#ifndef OMP
     s++;            
     ix++;
#endif
   }

#endif // SPI

/* with SPI we can easily interleave communication and computation
   so we don't wait for the communication to complete */

#ifdef OMP
# if 0//def SPI
#pragma omp single nowait
# else
#pragma omp single
# endif // SPI
   {
#endif // OMP
     
#    if (defined MPI && !defined _NO_COMM)
#      if 0//def SPI

     // Initialize the barrier, resetting the hardware.
     int rc = MUSPI_GIBarrierInit ( &GIBarrier, 0 /*comm world class route */);
     if(rc) {
       printf("MUSPI_GIBarrierInit returned rc = %d\n", rc);
       exit(__LINE__);
     }
     // reset the recv counter 
     recvCounter = totalMessageSize;
     global_barrier(); // make sure everybody is set recv counter

     //#pragma omp for nowait
     for (unsigned int j = 0; j < spi_num_dirs; j++) {
       descCount[ j ] =
	 msg_InjFifoInject ( injFifoHandle,
			     j,
			     &SPIDescriptors[j]);
     }
#      else // SPI
     xchange_halffield(); 
#      endif // SPI
#    endif
     
#ifdef OMP
   }
#endif


#if 0//def SPI

/* compute the body while the communication is in progress 
   the barrier at the end of this for loop will also catch the
   thread carrying out the "single nowait" above */
   
#ifdef OMP
#pragma omp for
#endif
  for(int i = 0; i < BODY/2; ++i)
  {
    U=u0+body[i]*4;
    _prefetch_su3(U);
    s=k+body[i];
    _prefetch_spinor(s);
    ix=body[i]*8;
    
    _hop_t_p_pre();
    U++;
    ix++;
     
    _hop_t_m_pre();
    ix++;
     
    _hop_x_p_pre();
    U++;
    ix++;
     
    _hop_x_m_pre();
    ix++;
     
    _hop_y_p_pre();
    U++;
    ix++;
     
    _hop_y_m_pre();
    ix++;
     
    _hop_z_p_pre();
    U++;
    ix++;
     
    _hop_z_m_pre(); 
  }

  /* wait for the communication to finish */
#ifdef OMP
#pragma omp single
{
#endif
  while ( recvCounter > 0 );
  _bgq_msync();
#ifdef OMP
}
#endif

#endif // SPI

/* now we can move on to compute the solution as usual */

#ifndef OMP
   s = l;
   if(ieo == 0) {
     U = g_gauge_field_copy[1][0];
   }
   else {
     U = g_gauge_field_copy[0][0];
   }
   _prefetch_su3(U);
#else
   if(ieo == 0) {
     u0 = g_gauge_field_copy[1][0];
   }
   else {
     u0 = g_gauge_field_copy[0][0];
   }
#endif
   
   phi = NBPointer[2 + ieo];
   
#ifdef OMP
#pragma omp for
#else
   ix = 0;
#endif
   /* #pragma ivdep */
   for(unsigned int i = 0; i < (VOLUME)/2; i++){
#ifdef OMP
     ix=i*8;
     U=u0+i*4;
     _prefetch_su3(U);
     s=l+i;
     _prefetch_spinor(s);
#endif
#ifdef _TM_SUB_HOP
     pn=p+i;
#endif
     _hop_t_p_post();
     ix++;
     
     _hop_t_m_post();
     ix++;
     U++;
     
     _hop_x_p_post();
     ix++;
     
     _hop_x_m_post();
     U++;
     ix++;
     
     _hop_y_p_post();
     ix++;
     
     _hop_y_m_post();
     U++;
     ix++;
     
     _hop_z_p_post();
     ix++;
     
     _hop_z_m_post();
     
#ifdef _MUL_G5_CMPLX
     _hop_mul_g5_cmplx_and_store(s);
#elif defined _TM_SUB_HOP
     _g5_cmplx_sub_hop_and_g5store(s);
#else
     _hop_store_post(s);
#endif
     
#ifndef OMP
     U++;
     ix++;
     s++;
#endif
   }
 }
#ifdef _KOJAK_INST
#pragma pomp inst end(hoppingmatrix)
#endif



