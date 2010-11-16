/***********************************************************************
 * $Id$ 
 *
 * Copyright (C) 2008 Albert Deuzeman, Siebren Reker, Carsten Urbach
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
#include <config.h>
#endif
#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "global.h"
#include "su3.h"
#include "complex.h"
#include "start.h"
#include "ranlxs.h"
#include "D_psi.h"
/* GG */
#include "read_input.h"
#include <mcheck.h>

#include "poly_precon.h"
#include "Msap.h"
#include "gmres_precon.h"
#include "linalg_eo.h"
#include "gram-schmidt.h"
#include "lu_solve.h"
#include "block.h"
#include "little_D.h"
#include "gcr4complex.h"
#include "boundary.h"
#include "generate_dfl_subspace.h"

int init_little_dfl_subspace(const int N_s);

spinor ** dfl_fields = NULL;
static spinor * _dfl_fields = NULL;
complex ** little_dfl_fields = NULL;
static complex *_little_dfl_fields = NULL;
static int init_subspace = 0;
static int init_little_subspace = 0;

static void random_fields(const int Ns) {
  
  int i, j, ix;
  float r,s[24];
  double *t;
  
  r=(float)(1.0/sqrt(24.0*(double)(VOLUME)));
  
  for (i=0;i<Ns;i++) {
    t=(double*)(dfl_fields[i]);
    for (ix = 0; ix < VOLUME; ix++){
      ranlxs(s,24);
      for (j = 0; j < 24; j++) {
	(*t)=(double)(r*(s[j]-0.5f));
	t+=1;
      }
    }
  }  
  return;
}

int generate_dfl_subspace(const int Ns, const int N) {
  int i, j, blk, vpr = VOLUMEPLUSRAND*sizeof(spinor)/sizeof(complex), 
    vol = VOLUME*sizeof(spinor)/sizeof(complex);
  double nrm, e = 0.3, d = 1.1, atime, etime;
  complex s;
  complex * work = NULL;
  spinor * r, * q, * p;

#ifdef MPI
  atime = MPI_Wtime();
#else
  atime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif
  work = (complex*)malloc(2*9*Ns*sizeof(complex));
  /* GG */
/*   if ( errno ) { */
  if ( work == NULL ) {
    printf(" Not enough memory in generate_dfl_subspace ! \n"); fflush(stdout);
#ifdef MPI
    MPI_Finalize();
#endif
    exit(69);
  }

  if (init_subspace == 0) {
    i = 0;
    i = init_dfl_subspace(Ns);
    /* GG */
    if ( i ) {
      printf(" Not enough memory in init_dfl_subspace called by generate_dfl_subspace ! \n"); fflush(stdout);
#ifdef MPI
      MPI_Finalize();
#endif
      exit(69);
    }
  }

  if (init_little_subspace == 0) {
    i = 0;
    i = init_little_dfl_subspace(Ns);
    /* GG */
    if ( i ) {
      printf(" Not enough memory in init_little_dfl_subspace called by generate_dfl_subspace ! \n"); fflush(stdout);
#ifdef MPI
      MPI_Finalize();
#endif
      exit(69);
    }
  }

  random_fields(Ns);

  if(g_debug_level > 4) {
    for(e = 0.; e < 1.; e=e+0.05) {
      random_spinor_field(dfl_fields[0], N, 0);
      nrm = sqrt(square_norm(dfl_fields[0], N, 1));
      mul_r(dfl_fields[0], 1./nrm, dfl_fields[0], N);
      d = 1.1;
/*       gmres_precon(g_spinor_field[DUM_SOLVER], dfl_fields[0], 20, 1, 1.e-20, 0, N, &D_psi); */
      poly_nonherm_precon(g_spinor_field[DUM_SOLVER], dfl_fields[0], e, d, 30, N);
      D_psi(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
      diff(g_spinor_field[DUM_SOLVER], g_spinor_field[DUM_SOLVER+1], dfl_fields[0], N);
      nrm = square_norm(g_spinor_field[DUM_SOLVER], N, 1);
      if(g_proc_id == 0) {
        printf(" e= %f d= %f nrm = %1.5e\n", e, d, nrm);
      }
    }
    d = 1.1;
    e=0.3;
  }

  /* GG */
  //mtrace();

  boundary(0.1586);
  for(i = 0; i < Ns; i++) {
    ModifiedGS((complex*)dfl_fields[i], vol, i, (complex*)dfl_fields[0], vpr);
    nrm = sqrt(square_norm(dfl_fields[i], N, 1));
    mul_r(dfl_fields[i], 1./nrm, dfl_fields[i], N);

    /* GG */
/*      for(j = 0; j < 2; j++) {  */
    for(j = 0; j < gilbert_loop_parameter; j++) { 

      g_sloppy_precision = 1;
/*        Msap(g_spinor_field[0], dfl_fields[i], 4); */
    /* GG */
      
      poly_nonherm_precon(g_spinor_field[0], dfl_fields[i], e, d, gilbert_poly_parameter, N);
      
      /*
      gmres_precon(g_spinor_field[DUM_SOLVER], dfl_fields[i], gilbert_poly_parameter, 1, 1.e-14, 0, N, &D_psi);
      */
/*       gmres_precon(g_spinor_field[DUM_SOLVER], dfl_fields[i], 20, 1, 1.e-20, 0, N, &D_psi); */
      g_sloppy_precision = 0;
      ModifiedGS((complex*)g_spinor_field[0], vol, i, (complex*)dfl_fields[0], vpr);
      nrm = sqrt(square_norm(g_spinor_field[0], N, 1));
      mul_r(dfl_fields[i], 1./nrm, g_spinor_field[0], N);
/* GG */
/*
      if ( ((j+1)%8) == 0 ) {
	if(g_debug_level > -1) {
	  D_psi(g_spinor_field[DUM_SOLVER], dfl_fields[i]);
	  nrm = sqrt(square_norm(g_spinor_field[DUM_SOLVER], N, 1));
	  if(g_proc_id == 0) {
	    printf(" proc %d loopind %d ||D psi_%d||/||psi_%d|| = %1.5e\n", g_proc_id, j+1, i, i, nrm); fflush(stdout); 
	  }
	}
      }
*/
    }
    /* test quality */
    if(g_debug_level > -1) {
      D_psi(g_spinor_field[DUM_SOLVER], dfl_fields[i]);
      nrm = sqrt(square_norm(g_spinor_field[DUM_SOLVER], N, 1));
      if(g_proc_id == 0) {
	printf(" ||D psi_%d||/||psi_%d|| = %1.5e\n", i, i, nrm); 
      }
    }
  }
  g_sloppy_precision = 0;
  boundary(g_kappa);
  if(g_debug_level > 4) {
    for(i = 0; i < Ns; i++) {
      for(j = 0; j < Ns; j++) {
	s = scalar_prod(dfl_fields[i], dfl_fields[j], N, 1);
	if(g_proc_id == 0) {
	  printf("<%d, %d> = %1.3e +i %1.3e\n", i, j, s.re, s.im);
	}
      }
    }
  }
  for (i = 0; i < Ns; i++) { 
    /* add it to the basis */
    split_global_field(block_list[0].basis[i], block_list[1].basis[i], dfl_fields[i]);
  }

  /* perform local orthonormalization */
  block_orthonormalize(block_list);
  block_orthonormalize(block_list+1);

  dfl_subspace_updated = 1;

  for(j = 0; j < Ns; j++) {
    for(i = 0; i < 2*9*Ns; i++) {
      _complex_zero(little_dfl_fields[j][i]);
      _complex_zero(work[i]);
    }
  }

  /* compute the little little basis */
  r = g_spinor_field[DUM_SOLVER];
  q = g_spinor_field[DUM_SOLVER+1];
  
  for(i = 0; i < Ns; i++) {
    split_global_field(r, q,  dfl_fields[i]);
    /* now take the local scalar products */
    for(j = 0; j < Ns; j++) {
      p = r;
      for(blk = 0; blk < 2; blk++) {
	if(blk == 0) p = r;
	else p = q;
	little_dfl_fields[i][j + blk*Ns] = scalar_prod(block_list[blk].basis[j], p, block_list[0].volume, 0);
      }
    }
  }
  
  /* orthonormalise */
  for(i = 0; i < Ns; i++) {
    for (j = 0; j < i; j++) {
      s = lscalar_prod(little_dfl_fields[j], little_dfl_fields[i], 2*Ns, 1);
      lassign_diff_mul(little_dfl_fields[i], little_dfl_fields[j], s, 2*Ns);
    }
    s.re = lsquare_norm(little_dfl_fields[i], 2*Ns, 1);
    lmul_r(little_dfl_fields[i], 1./s.re, little_dfl_fields[i], 2*Ns);
  }
  if(g_debug_level > 4) {
    for(i = 0; i < Ns; i++) {
      for(j = 0; j < Ns; j++) {
	s = lscalar_prod(little_dfl_fields[i], little_dfl_fields[j], 2*Ns, 1);
	if(g_proc_id == 0) {
	  printf("<%d, %d> = %1.3e +i %1.3e\n", i, j, s.re, s.im);
	}
      }
    }
  }
  
  for(i = 0; i < Ns; i++) {
    little_D(work, little_dfl_fields[i]);
    for(j = 0; j < Ns; j++) {
      little_A[i * Ns + j]  = lscalar_prod(little_dfl_fields[j], work, 2*Ns, 1);
      if(g_proc_id == 0 && g_debug_level > 4) {
	printf("%1.3e %1.3ei, ", little_A[i * Ns + j].re, little_A[i * Ns + j].im);
      }
    }
    if(g_proc_id == 0 && g_debug_level > 4) printf("\n");
  }
  if(g_proc_id == 0 && g_debug_level > 4) printf("\n");
  /* the precision in the inversion is not yet satisfactory! */
  LUInvert(Ns, little_A, Ns);
  /* inverse of little little D now in little_A */

#ifdef MPI
  etime = MPI_Wtime();
#else
  etime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif
  if(g_proc_id == 0) {
    printf("time for subspace generation %1.3e s\n", etime-atime);
    fflush(stdout);
  }

  free_dfl_subspace();
  free(work);
  return(0);
}

int generate_dfl_subspace_free(const int Ns, const int N) {
  int i,j, vpr = VOLUMEPLUSRAND*sizeof(spinor)/sizeof(complex), 
    vol = VOLUME*sizeof(spinor)/sizeof(complex);
  double nrm;
  complex s;

  if (init_subspace == 0) {
    i = 0;
    i = init_dfl_subspace(Ns);

  /* GG */
  if ( i ) {
    printf(" Not enough memory in init_dfl_subspace called by generate_dfl_subspace_free ! \n"); fflush(stdout);
#ifdef MPI
    MPI_Finalize();
#endif
    exit(69);
  }
  }

  for(i = 0; i < 12; i++) {
    constant_spinor_field(dfl_fields[i], i, N);
    ModifiedGS((complex*)dfl_fields[i], vol, i, (complex*)dfl_fields[0], vpr);
    nrm = sqrt(square_norm(dfl_fields[i], N, 1));
    mul_r(dfl_fields[i], 1./nrm, dfl_fields[i], N);

    /* test quality */
    if(g_debug_level > -1) {
      D_psi(g_spinor_field[DUM_SOLVER], dfl_fields[i]);
      nrm = sqrt(square_norm(g_spinor_field[DUM_SOLVER], N, 1));
      if(g_proc_id == 0) {
	printf(" ||D psi_%d||/||psi_%d|| = %1.5e\n", i, i, nrm); 
      }
    }
  }

  if(g_debug_level > 4) {
    for(i = 0; i < 12; i++) {
      for(j = 0; j < 12; j++) {
	s = scalar_prod(dfl_fields[i], dfl_fields[j], N, 1);
	if(g_proc_id == 0) {
	  printf("<%d, %d> = %1.3e +i %1.3e\n", i, j, s.re, s.im);
	}
      }
    }
  }

  return(0);
}

int init_little_dfl_subspace(const int N_s) {
  int i;
  if(init_little_subspace == 0) {
    if((void*)(_little_dfl_fields = (complex*)calloc((N_s)*2*9*N_s+4, sizeof(complex))) == NULL) {
      return(1);
    }
    if((void*)(little_dfl_fields = (complex**)calloc(N_s, sizeof(complex*))) == NULL) {
      return(1);
    }
#if ( defined SSE || defined SSE2 || defined SSE3)
    little_dfl_fields[0] = (complex*)(((unsigned long int)(_little_dfl_fields)+ALIGN_BASE)&~ALIGN_BASE);
#else
    little_dfl_fields[0] = _little_dfl_fields;
#endif
    for (i = 1; i < N_s; i++) {
      little_dfl_fields[i] = little_dfl_fields[i-1] + 2*9*N_s;
    }
    if((void*)(little_A = (complex*)calloc(N_s*N_s, sizeof(complex))) == NULL) {
      return(1);
    }
    init_little_subspace = 1;
  }
  return(0);
}

int init_dfl_subspace(const int N_s) {
  int i;
  init_subspace = 1;
  if((void*)(_dfl_fields = calloc((N_s)*VOLUMEPLUSRAND+1, sizeof(spinor))) == NULL) {
    return(1);
  }
  if ((void*)(dfl_fields = calloc((N_s), sizeof(spinor *))) == NULL) {
    return(1);
  }
#if ( defined SSE || defined SSE2 || defined SSE3)
  dfl_fields[0] = (spinor*)(((unsigned long int)(_dfl_fields)+ALIGN_BASE)&~ALIGN_BASE);
#else
  dfl_fields[0] = _dfl_fields;
#endif
  for (i = 1; i < N_s; ++i) {
    dfl_fields[i] = dfl_fields[i-1] + VOLUMEPLUSRAND;
  }
  return(0);
}

int free_dfl_subspace() {
  if(init_subspace == 1) {
    free(dfl_fields);
    free(_dfl_fields);
    init_subspace = 0;
  }
  return (0);
}
