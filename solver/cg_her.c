/***********************************************************************
 * $Id$
 *
 * Copyright (C) 2001 Martin Hasenbusch
 *               2003 Thomas Chiarappa
 *               2002,2003,2004,2005,2010 Carsten Urbach
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
 * $Id$
 *
 * File: cg_her.c
 *
 * CG solver for hermitian f only!
 *
 * The externally accessible functions are
 *
 *
 *   int cg(spinor * const P, spinor * const Q, double m, const int subtract_ev)
 *     CG solver
 *
 * input:
 *   Q: source
 * inout:
 *   P: initial guess and result
 *
 *
 **************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "su3.h"
#include "linalg_eo.h"
#include "start.h"
#include "solver/matrix_mult_typedef.h"
#include "sub_low_ev.h"
#include "poly_precon.h"
#include "cg_her.h"

int cg_her(spinor * const P, spinor * const Q, const int max_iter,
	   double eps_sq, const int rel_prec, const int N, matrix_mult f) {

  static double normsq,pro,err,alpha_cg,beta_cg,squarenorm;
  int iteration;
  int save_sloppy = g_sloppy_precision;
  double atime, etime, flops;

  /* GG */
  double mflops_mpi, mflops_local;

  /* initialize residue r and search vector p */
#ifdef MPI
  atime = MPI_Wtime();
#else
  atime = ((double)clock())/((double)(CLOCKS_PER_SEC));
#endif
  squarenorm = square_norm(Q, N, 1);

  f(g_spinor_field[DUM_SOLVER], P);

  diff(g_spinor_field[DUM_SOLVER+1], Q, g_spinor_field[DUM_SOLVER], N);
  assign(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1], N);
  normsq=square_norm(g_spinor_field[DUM_SOLVER+1], N, 1);

  /* main loop */
  for(iteration = 1; iteration <= max_iter; iteration++) {
    f(g_spinor_field[DUM_SOLVER], g_spinor_field[DUM_SOLVER+2]);
    pro = scalar_prod_r(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER], N, 1);

    /*  Compute alpha_cg(i+1)   */
    alpha_cg = normsq / pro;

    /*  Compute x_(i+1) = x_i + alpha_cg(i+1) p_i    */
    assign_add_mul_r(P, g_spinor_field[DUM_SOLVER+2], alpha_cg, N);

    /*  Compute r_(i+1) = r_i - alpha_cg(i+1) Qp_i   */
    assign_mul_add_r(g_spinor_field[DUM_SOLVER], -alpha_cg, g_spinor_field[DUM_SOLVER+1], N);

    /* Check whether the precision is reached ... */
    err=square_norm(g_spinor_field[DUM_SOLVER], N, 1);

    if(g_proc_id == g_stdio_proc && g_debug_level > 1) {

            /* GG */
      etime = MPI_Wtime();
      //MK printf("cg_her: %d\t%g\t%e\n",iteration,err, etime-atime); fflush( stdout);
      //printf("CG: iterations: %d res^2 %e\n", iteration, err); fflush(stdout);
    }

    if (((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*squarenorm) && (rel_prec == 1))) {
      break;
    }
#ifdef _USE_HALFSPINOR
    if(((err*err <= eps_sq) && (rel_prec == 0)) || ((err*err <= eps_sq*squarenorm) && (rel_prec == 1))) {
      g_sloppy_precision = 1;
      if(g_debug_level > 2 && g_proc_id == g_stdio_proc) {
	//MK printf("sloppy precision on\n"); fflush( stdout);
      }
    }
#endif

    /* Compute beta_cg(i+1)
       Compute p_(i+1) = r_i+1 + beta_(i+1) p_i     */
    beta_cg = err / normsq;
    assign_mul_add_r(g_spinor_field[DUM_SOLVER+2], beta_cg, g_spinor_field[DUM_SOLVER], N);
    assign(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER], N);
    normsq = err;
  }
#ifdef MPI
  etime = MPI_Wtime();
#else
  etime = ((double)clock())/((double)(CLOCKS_PER_SEC));
#endif
  g_sloppy_precision = save_sloppy;
  /* 2 A + 2 Nc Ns + N_Count ( 2 A + 10 Nc Ns ) */
  /* 2*1320.0 because the linalg is over VOLUME/2 */
  flops = (2*(2*1320.0+2*3*4) + 2*3*4 + iteration*(2.*(2*1320.0+2*3*4) + 10*3*4))*N/1.0e6f;

  /* GG */
  mflops_local = flops/(etime-atime);
  mflops_mpi = mflops_local;
#ifdef MPI
  MPI_Reduce(&mflops_local, &mflops_mpi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(g_proc_id==0 && g_debug_level > 0) {
    printf("CGgg: flopcount: t/s: %1.4e mflops_local: %.1f mflops_global: %.1f\n",
	   etime-atime, mflops_local, mflops_mpi); fflush(stdout);
  }
#endif

  if(g_debug_level > 0 && g_proc_id == 0) {
    printf("cg_her: iter: %d eps_sq: %1.4e t/s: %1.4e\n", iteration, eps_sq, etime-atime);
    printf("cg_her: flopcount (for tmWilson with even/odd only): t/s: %1.4e mflops_local: %.1f mflops: %.1f\n",
	   etime-atime, flops/(etime-atime), g_nproc*flops/(etime-atime));
  }
  if(iteration > max_iter) return(-1);
  return(iteration);
}

static char const rcsid[] = "$Id$";








