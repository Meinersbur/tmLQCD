/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
/* $Id$ */

/**********************************************************
 * 
 * exchange routines for half spinor fields
 *
 * Author: Carsten Urbach 
 *
 **********************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef MPI
# include <mpi.h>
#endif
#ifdef _USE_SHMEM
# include <mpp/shmem.h>
#endif
#include "global.h"
#if (defined XLC && defined BGL)
#  include "bgl.h"
#endif
#include "mpi_init.h"
#include "su3.h"
#include "init_dirac_halfspinor.h"
#include "xchange_halffield.h"

/* GG */
static int ivisit = 0;

/* GG */
//extern halfspinor * GalfSpinor ALIGN;

#if (!defined _INDEX_INDEP_GEOM)
#if (defined _USE_SHMEM && defined _USE_HALFSPINOR)
# include <mpp/shmem.h>
/* 1. */
void xchange_halffield() {

#ifdef _KOJAK_INST
#pragma pomp inst begin(xchangehalf)
#endif
#  ifdef MPI

  shmem_barrier_all();

  shmem_double_put((double*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ/2), 
		   (double*)(HalfSpinor + 4*VOLUME),
                   (LX*LY*LZ*6), g_nb_t_up);
  shmem_double_put((double*)(HalfSpinor + 4*VOLUME + RAND/2), 
		   (double*)(HalfSpinor + 4*VOLUME + LX*LY*LZ/2),
                   (LX*LY*LZ*6), g_nb_t_dn);

#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  shmem_double_put((double*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ/2), 
		   (double*)(HalfSpinor + 4*VOLUME + LX*LY*LZ),
                   (T*LY*LZ*6), g_nb_x_up);
  shmem_double_put((double*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ),
		   (double*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ/2),
                   (T*LY*LZ*6), g_nb_x_dn);

#    endif

#    if (defined PARALLELXYT || defined PARALLELXYZT)
  shmem_double_put((double*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2),
		   (double*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ),
                   (T*LX*LZ*6), g_nb_y_up);
  shmem_double_put((double*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ),
		   (double*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2),
                   (T*LX*LZ*6), g_nb_y_dn);

#    endif

#    if (defined PARALLELXYZT)
  shmem_double_put((double*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2),
		   (double*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ),
                   (T*LX*LY*6), g_nb_z_up);
  shmem_double_put((double*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ),
		   (double*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2),
                   (T*LX*LY*6), g_nb_z_dn);

#    endif

  shmem_barrier_all();
#  endif /* MPI */
  return;
#ifdef _KOJAK_INST
#pragma pomp inst end(xchangehalf)
#endif
}

#elif (defined _USE_HALFSPINOR && defined _PERSISTENT)

MPI_Request prequests[16];

/* 2. */
void init_xchange_halffield() {

#  ifdef MPI

#  ifdef PARALLELT
  int reqcount = 4;
#  elif defined PARALLELXT
  int reqcount = 8;
#  elif defined PARALLELXYT
  int reqcount = 12;
#  elif defined PARALLELXYZT
  int x0=0, x1=0, x2=0, ix=0;
  int reqcount = 16;
#  endif
#  if (defined XLC && defined BGL)
  __alignx(16, HalfSpinor);
#  endif

  /* send the data to the neighbour on the right in t direction */
  /* recieve the data from the neighbour on the left in t direction */
  MPI_Send_init((void*)(HalfSpinor + 4*VOLUME), LX*LY*LZ*12/2, MPI_DOUBLE, 
		g_nb_t_up, 81, g_cart_grid, &prequests[0]);
  
  MPI_Recv_init((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ/2), LX*LY*LZ*12/2, MPI_DOUBLE, 
		g_nb_t_dn, 81, g_cart_grid, &prequests[1]);

  /* send the data to the neighbour on the left in t direction */
  /* recieve the data from the neighbour on the right in t direction */
  MPI_Send_init((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ/2), LX*LY*LZ*12/2, MPI_DOUBLE, 
		g_nb_t_dn, 82, g_cart_grid, &prequests[2]);

  MPI_Recv_init((void*)(HalfSpinor + 4*VOLUME + RAND/2), LX*LY*LZ*12/2, MPI_DOUBLE, 
		g_nb_t_up, 82, g_cart_grid, &prequests[3]);

#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  MPI_Send_init((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ), T*LY*LZ*12/2, MPI_DOUBLE, 
		g_nb_x_up, 91, g_cart_grid, &prequests[4]);

  MPI_Recv_init((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ/2), T*LY*LZ*12/2, MPI_DOUBLE,
		g_nb_x_dn, 91, g_cart_grid, &prequests[5]);

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */  
  MPI_Send_init((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ/2), T*LY*LZ*12/2, MPI_DOUBLE,
		g_nb_x_dn, 92, g_cart_grid, &prequests[6]);

  MPI_Recv_init((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ), T*LY*LZ*12/2, MPI_DOUBLE,
		g_nb_x_up, 92, g_cart_grid, &prequests[7]);
#    endif
    
#    if (defined PARALLELXYT || defined PARALLELXYZT)
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
  MPI_Send_init((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ), T*LX*LZ*12/2, MPI_DOUBLE, 
		g_nb_y_up, 101, g_cart_grid, &prequests[8]);

  MPI_Recv_init((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2), T*LX*LZ*12/2, MPI_DOUBLE, 
		g_nb_y_dn, 101, g_cart_grid, &prequests[9]);

    /* send the data to the neighbour on the leftt in y direction */
    /* recieve the data from the neighbour on the right in y direction */
  MPI_Send_init((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2), T*LX*LZ*12/2, MPI_DOUBLE, 
		g_nb_y_dn, 102, g_cart_grid, &prequests[10]);

  MPI_Recv_init((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ), T*LX*LZ*12/2, MPI_DOUBLE, 
		g_nb_y_up, 102, g_cart_grid, &prequests[11]);
#    endif
    
#    if (defined PARALLELXYZT)
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  MPI_Send_init((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ), 
		T*LX*LY*12/2, MPI_DOUBLE, g_nb_z_up, 503, g_cart_grid, &prequests[12]); 

  MPI_Recv_init((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2), 
		T*LX*LY*12/2, MPI_DOUBLE, g_nb_z_dn, 503, g_cart_grid, &prequests[13]); 

  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Send_init((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2), 
		12*T*LX*LY/2, MPI_DOUBLE, g_nb_z_dn, 504, g_cart_grid, &prequests[14]); 

  MPI_Recv_init((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ), 
		T*LX*LY*12/2, MPI_DOUBLE, g_nb_z_up, 504, g_cart_grid, &prequests[15]); 
#  endif
#  endif /* MPI */
  return;
}

/* 3. */
void xchange_halffield() {
#  ifdef MPI

  MPI_Status status[16];
#    ifdef PARALLELT
  int reqcount = 4;
#    elif defined PARALLELXT
  int reqcount = 8;
#    elif defined PARALLELXYT
  int reqcount = 12;
#    elif defined PARALLELXYZT
  int x0=0, x1=0, x2=0, ix=0;
  int reqcount = 16;
#    endif
#    if (defined XLC && defined BGL)
  __alignx(16, HalfSpinor);
#    endif
  MPI_Startall(reqcount, prequests);

  MPI_Waitall(reqcount, prequests, status); 
#  endif /* MPI */
  return;
}

#else /* def _USE_HALFSPINOR && (_USE_SHMEM || _PERSISTENT) */ 

/* 4. */
void xchange_halffield() {

#  ifdef MPI

  /* GG */
  //int buf_size, ilg_max = 1000000, ilg;
  int ilg = -1, g_mon = -1, bufsH = buf_size;
#define TRANSTST 0

  MPI_Request requests[16];
  MPI_Status status[16];
#  ifdef PARALLELT
  int reqcount = 4;
#  elif defined PARALLELXT
  int reqcount = 8;
#  elif defined PARALLELXYT || defined PARALLELXYZ
  int reqcount = 12;
#  elif defined PARALLELXYZT
  int reqcount = 16;
#  endif
#  if (defined XLC && defined BGL)
  __alignx(16, HalfSpinor);
#  endif

#ifdef _KOJAK_INST
#pragma pomp inst begin(xchangehalf)
#endif

#if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  /* send the data to the neighbour on the right in t direction */
  /* recieve the data from the neighbour on the left in t direction */
  MPI_Isend((void*)(HalfSpinor + 4*VOLUME), LX*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_t_up, 81, g_cart_grid, &requests[0]);

  MPI_Irecv((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ/2), LX*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_t_dn, 81, g_cart_grid, &requests[1]);

  /* send the data to the neighbour on the left in t direction */
  /* recieve the data from the neighbour on the right in t direction */
  MPI_Isend((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ/2), LX*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_t_dn, 82, g_cart_grid, &requests[2]);

  MPI_Irecv((void*)(HalfSpinor + 4*VOLUME + RAND/2), LX*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_t_up, 82, g_cart_grid, &requests[3]);
#endif

#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXYZ)

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  MPI_Isend((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ), T*LY*LZ*12/2, MPI_DOUBLE, 
	    g_nb_x_up, 91, g_cart_grid, &requests[4]);

  MPI_Irecv((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ/2), T*LY*LZ*12/2, MPI_DOUBLE,
	    g_nb_x_dn, 91, g_cart_grid, &requests[5]);

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */  
  MPI_Isend((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ/2), T*LY*LZ*12/2, MPI_DOUBLE,
 	    g_nb_x_dn, 92, g_cart_grid, &requests[6]);

  MPI_Irecv((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ), T*LY*LZ*12/2, MPI_DOUBLE,
 	    g_nb_x_up, 92, g_cart_grid, &requests[7]);
#    endif
    
#    if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXYZ)
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
  MPI_Isend((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ), T*LX*LZ*12/2, MPI_DOUBLE, 
	    g_nb_y_up, 101, g_cart_grid, &requests[8]);

  MPI_Irecv((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2), T*LX*LZ*12/2, MPI_DOUBLE, 
	    g_nb_y_dn, 101, g_cart_grid, &requests[9]);

    /* send the data to the neighbour on the leftt in y direction */
    /* recieve the data from the neighbour on the right in y direction */
  MPI_Isend((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2), T*LX*LZ*12/2, MPI_DOUBLE, 
	    g_nb_y_dn, 102, g_cart_grid, &requests[10]);

  MPI_Irecv((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ), T*LX*LZ*12/2, MPI_DOUBLE, 
	    g_nb_y_up, 102, g_cart_grid, &requests[11]);
#    endif
    
#    if (defined PARALLELXYZT || defined PARALLELXYZ)
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */

  /* GG */
#if TRANSTST
  if ( ivisit == 7 ) {
#ifdef ORIG
    g_mon = 0;
#else
    g_mon = source_node;
#endif
    if( g_debug_level > 0 && (g_proc_id == g_mon) ) {
      printf("xchange_halffield ggmonS isend %d g_nb_z_up %d <= g_proc_id %d Word % 23.16e \n", 503, g_nb_z_up, g_proc_id,
	     (*(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ)).s0.c0.re );
    }
    if( g_debug_level > 0 && (g_proc_id == g_mon) ) {
      printf("xchange_halffield ggmonS isend %d g_nb_z_dn %d <= g_proc_id %d Word % 23.16e \n", 504, g_nb_z_dn, g_proc_id,
	     (*(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2)).s0.c0.re );
      GalfSpinorSrc[0] = HalfSpinor[4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2];
    }
  }
#endif

  MPI_Isend((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ), 
	    T*LX*LY*12/2, MPI_DOUBLE, g_nb_z_up, 503, g_cart_grid, &requests[12]); 

  MPI_Irecv((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2), 
	    T*LX*LY*12/2, MPI_DOUBLE, g_nb_z_dn, 503, g_cart_grid, &requests[13]); 

  /* GG */
#if TRANSTST
  if ( ivisit == 7 ) {
    if(g_debug_level > 0 && g_proc_id == g_stdio_proc) {
      printf("xchange_halffield ggmonA isend %d g_nb_z_up %d <= g_proc_id %d <= g_nb_z_dn %d ivisit %d size %d Y g_nb_y_up %d g_nb_y_dn %d X g_nb_x_up %d g_nb_x_dn %d T g_nb_t_up %d g_nb_t_dn %d\n", 
	     503, g_nb_z_up, g_proc_id, g_nb_z_dn, ivisit, T*LX*LY*12/2,
	     g_nb_y_up, g_nb_y_dn, g_nb_x_up, g_nb_x_dn, g_nb_t_up, g_nb_t_dn
	     );
    }
  }
#endif

  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Isend((void*)(HalfSpinor + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2), 
	    12*T*LX*LY/2, MPI_DOUBLE, g_nb_z_dn, 504, g_cart_grid, &requests[14]); 

  MPI_Irecv((void*)(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ), 
	    T*LX*LY*12/2, MPI_DOUBLE, g_nb_z_up, 504, g_cart_grid, &requests[15]); 
#    endif

  MPI_Waitall(reqcount, requests, status); 

  /* GG */
#if TRANSTST
  if ( ivisit == 7 ) {
#ifdef ORIG
    g_mon = 0;
#else
    g_mon = target_node;
#endif
    if( g_debug_level > 0 && (g_proc_id == g_mon) ) {
      printf("xchange_halffield ggmonR isend %d g_nb_z_dn %d => g_proc_id %d Word % 23.16e \n", 503, g_nb_z_dn, g_proc_id,
	     (*(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2)).s0.c0.re );
    }
    if( g_debug_level > 0 && (g_proc_id == g_mon) ) {
      printf("xchange_halffield ggmonR isend %d g_nb_z_up %d => g_proc_id %d Word % 23.16e \n", 504, g_nb_z_up, g_proc_id,
	     (*(HalfSpinor + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ)).s0.c0.re );
      GalfSpinorSrc[0] = HalfSpinor[4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ];
    }
  }
  /* Beginning of IB transfer testing loop */
  /* buf_size, ilg_max, source_node, target_node are
     input options in the main program (curr.: invert.c)
     buf_size is a number of halfspinor in input,
     but the bufsW printed value is in number of MPI_DOUBLE words 
     while bufsB is in bytes */
  if ( ivisit == 7 ) {
    // Filling up the test buffer ...
    for (ilg=1; ilg<buf_size; ilg++) {
      //GalfSpinorSrc[ilg] = HalfSpinor[4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ];
      GalfSpinorSrc[ilg] = GalfSpinorSrc[0];
    }
    // Computing the true size in number of MPI_DOUBLE words ...
    //buf_size = 1024*sizeof(halfspinor)/8;
    buf_size *= sizeof(halfspinor)/8;
    // Looping on the source_node side ...
    if( g_debug_level > 0 && (g_proc_id == source_node ) ) {
      for (ilg=0; ilg<ilg_max; ilg++) {
	MPI_Isend((void*)(GalfSpinorSrc), buf_size, MPI_DOUBLE, target_node, 503, g_cart_grid, &requests[12]);
	MPI_Irecv((void*)(GalfSpinorDst), buf_size, MPI_DOUBLE, target_node, 504, g_cart_grid, &requests[13]);
	MPI_Waitall(2, &requests[12], status); 
      }
      printf("xchange_halffield ggmonL bufsW %d bufsB %d ilgm %d srcn %d tgtn %d isend %d \n", 
	     buf_size, buf_size*sizeof(MPI_DOUBLE), ilg_max, source_node, target_node, 503); fflush(stdout);
      //printf("xchange_halffield ggmonSrc Gsrc % 23.16e Gdst % 23.16e isend %d \n", GalfSpinorSrc[0].s0.c0.re, GalfSpinorDst[0].s0.c0.re, 503); fflush(stdout);
      printf("xchange_halffield ggmonSrc Gsrc % 23.16e Gdst % 23.16e isend %d \n", GalfSpinorSrc[bufsH-1].s1.c2.im, GalfSpinorDst[bufsH-1].s1.c2.im, 503); fflush(stdout);
    }
    // Looping on the target_node side ...
    if( g_debug_level > 0 && (g_proc_id == target_node) ) {
      for (ilg=0; ilg<ilg_max; ilg++) {
	MPI_Isend((void*)(GalfSpinorSrc), buf_size, MPI_DOUBLE, source_node, 504, g_cart_grid, &requests[14]);
	MPI_Irecv((void*)(GalfSpinorDst), buf_size, MPI_DOUBLE, source_node, 503, g_cart_grid, &requests[15]);
	//printf("xchange_halffield ggmonL ilg %d isend %d \n", ilg, 504); fflush(stdout);
	MPI_Waitall(2, &requests[14], status); 
      }
      //printf("xchange_halffield ggmonDst Gsrc % 23.16e Gdst % 23.16e isend %d \n", GalfSpinorSrc[0].s0.c0.re, GalfSpinorDst[0].s0.c0.re, 504); fflush(stdout);
      printf("xchange_halffield ggmonDst Gsrc % 23.16e Gdst % 23.16e isend %d \n", GalfSpinorSrc[bufsH-1].s1.c2.im, GalfSpinorDst[bufsH-1].s1.c2.im, 504); fflush(stdout);
    }
  }
    ivisit++;
#endif

#  endif
  return;

#ifdef _KOJAK_INST
#pragma pomp inst end(xchangehalf)
#endif
}

#endif /* def _USE_HALFSPINOR && (_USE_SHMEM || _PERSISTENT) */ 

#if (defined _USE_SHMEM && defined _USE_HALFSPINOR)
# include <mpp/shmem.h>

/* 32-1. */
void xchange_halffield32() {

#ifdef _KOJAK_INST
#pragma pomp inst begin(xchangehalf32)
#endif
#  ifdef MPI

  shmem_barrier_all();

  shmem_float_put((float*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ/2), 
		   (float*)(HalfSpinor32 + 4*VOLUME),
                   (LX*LY*LZ*6), g_nb_t_up);
  shmem_float_put((float*)(HalfSpinor32 + 4*VOLUME + RAND/2), 
		   (float*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ/2),
                   (LX*LY*LZ*6), g_nb_t_dn);

#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  shmem_float_put((float*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ/2), 
		   (float*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ),
                   (T*LY*LZ*6), g_nb_x_up);
  shmem_float_put((float*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ),
		   (float*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ/2),
                   (T*LY*LZ*6), g_nb_x_dn);

#    endif

#    if (defined PARALLELXYT || defined PARALLELXYZT)
  shmem_float_put((float*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2),
		   (float*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ),
                   (T*LX*LZ*6), g_nb_y_up);
  shmem_float_put((float*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ),
		   (float*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2),
                   (T*LX*LZ*6), g_nb_y_dn);

#    endif

#    if (defined PARALLELXYZT)
  shmem_float_put((float*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2),
		   (float*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ),
                   (T*LX*LY*6), g_nb_z_up);
  shmem_float_put((float*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ),
		   (float*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2),
                   (T*LX*LY*6), g_nb_z_dn);

#    endif

  shmem_barrier_all();
#  endif /* MPI */
  return;
#ifdef _KOJAK_INST
#pragma pomp inst end(xchangehalf32)
#endif
}

#else /* (defined _USE_SHMEM && defined _USE_HALFSPINOR) */

/* 32-2. */
void xchange_halffield32() {

#  ifdef MPI

  MPI_Request requests[16];
  MPI_Status status[16];
#  ifdef PARALLELT
  int reqcount = 4;
#  elif defined PARALLELXT
  int reqcount = 8;
#  elif defined PARALLELXYT || defined PARALLELXYZ
  int reqcount = 12;
#  elif defined PARALLELXYZT
  int reqcount = 16;
#  endif
#ifdef _KOJAK_INST
#pragma pomp inst begin(xchangehalf32)
#endif
#  if (defined XLC && defined BGL)
  __alignx(16, HalfSpinor32);
#  endif

#if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)

  /* send the data to the neighbour on the right in t direction */
  /* recieve the data from the neighbour on the left in t direction */
  MPI_Isend((void*)(HalfSpinor32 + 4*VOLUME), LX*LY*LZ*12/2, MPI_FLOAT, 
	    g_nb_t_up, 81, g_cart_grid, &requests[0]);

  MPI_Irecv((void*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ/2), LX*LY*LZ*12/2, MPI_FLOAT, 
	    g_nb_t_dn, 81, g_cart_grid, &requests[1]);

  /* send the data to the neighbour on the left in t direction */
  /* recieve the data from the neighbour on the right in t direction */
  MPI_Isend((void*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ/2), LX*LY*LZ*12/2, MPI_FLOAT, 
	    g_nb_t_dn, 82, g_cart_grid, &requests[2]);

  MPI_Irecv((void*)(HalfSpinor32 + 4*VOLUME + RAND/2), LX*LY*LZ*12/2, MPI_FLOAT, 
	    g_nb_t_up, 82, g_cart_grid, &requests[3]);
#endif


#    if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXYZ)

  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */
  MPI_Isend((void*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ), T*LY*LZ*12/2, MPI_FLOAT, 
	    g_nb_x_up, 91, g_cart_grid, &requests[4]);

  MPI_Irecv((void*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ/2), T*LY*LZ*12/2, MPI_FLOAT,
	    g_nb_x_dn, 91, g_cart_grid, &requests[5]);

  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */  
  MPI_Isend((void*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ/2), T*LY*LZ*12/2, MPI_FLOAT,
 	    g_nb_x_dn, 92, g_cart_grid, &requests[6]);

  MPI_Irecv((void*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ), T*LY*LZ*12/2, MPI_FLOAT,
 	    g_nb_x_up, 92, g_cart_grid, &requests[7]);
#    endif
    
#    if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXYZ)
    /* send the data to the neighbour on the right in y direction */
    /* recieve the data from the neighbour on the left in y direction */
  MPI_Isend((void*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ), T*LX*LZ*12/2, MPI_FLOAT, 
	    g_nb_y_up, 101, g_cart_grid, &requests[8]);

  MPI_Irecv((void*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2), T*LX*LZ*12/2, MPI_FLOAT, 
	    g_nb_y_dn, 101, g_cart_grid, &requests[9]);

    /* send the data to the neighbour on the leftt in y direction */
    /* recieve the data from the neighbour on the right in y direction */
  MPI_Isend((void*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ/2), T*LX*LZ*12/2, MPI_FLOAT, 
	    g_nb_y_dn, 102, g_cart_grid, &requests[10]);

  MPI_Irecv((void*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ), T*LX*LZ*12/2, MPI_FLOAT, 
	    g_nb_y_up, 102, g_cart_grid, &requests[11]);
#    endif
    
#    if (defined PARALLELXYZT || defined PARALLELXYZ)
  /* send the data to the neighbour on the right in z direction */
  /* recieve the data from the neighbour on the left in z direction */
  MPI_Isend((void*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ), 
	    T*LX*LY*12/2, MPI_FLOAT, g_nb_z_up, 503, g_cart_grid, &requests[12]); 

  MPI_Irecv((void*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2), 
	    T*LX*LY*12/2, MPI_FLOAT, g_nb_z_dn, 503, g_cart_grid, &requests[13]); 

  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Isend((void*)(HalfSpinor32 + 4*VOLUME + LX*LY*LZ + T*LY*LZ + T*LX*LZ + T*LX*LY/2), 
	    12*T*LX*LY/2, MPI_FLOAT, g_nb_z_dn, 504, g_cart_grid, &requests[14]); 

  MPI_Irecv((void*)(HalfSpinor32 + 4*VOLUME + RAND/2 + LX*LY*LZ + T*LY*LZ + T*LX*LZ), 
	    T*LX*LY*12/2, MPI_FLOAT, g_nb_z_up, 504, g_cart_grid, &requests[15]); 
#    endif

  MPI_Waitall(reqcount, requests, status); 
#  endif /* MPI */
  return;
#ifdef _KOJAK_INST
#pragma pomp inst end(xchangehalf32)
#endif
}
#endif /* (defined _USE_SHMEM && defined _USE_HALFSPINOR) */
#endif /* (!defined _INDEX_INDEP_GEOM) */

static char const rcsid[] = "$Id$";
