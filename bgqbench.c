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
/*******************************************************************************
*
* Benchmark program for the even-odd preconditioned Wilson-Dirac operator
*
*
*******************************************************************************/

#define MAIN_PROGRAM
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#if (defined BGL && !defined BGP)
#  include <rts.h>
#endif
#ifdef MPI
# include <mpi.h>
# ifdef HAVE_LIBLEMON
#  include <io/params.h>
#  include <io/gauge.h>
# endif
#endif
#ifdef OMP
# include <omp.h>
# include "init_omp_accumulators.h"
#endif
#include "gettime.h"
#include "su3.h"
#include "su3adj.h"
#include "ranlxd.h"
#include "geometry_eo.h"
#include "read_input.h"
#include "start.h"
#include "boundary.h"
#include "Hopping_Matrix.h"
#include "Hopping_Matrix_nocom.h"
#include "tm_operators.h"
#include "global.h"
#include "xchange.h"
#include "init_gauge_field.h"
#include "init_geometry_indices.h"
#include "init_spinor_field.h"
#include "init_moment_field.h"
#include "init_dirac_halfspinor.h"
#include "test/check_geometry.h"
#include "xchange_halffield.h"
#include "D_psi.h"
#include "phmc.h"
#include "mpi_init.h"


/* BEGIN MK */
#include "bgq/bgq_field.h"
#include "bgq/bgq_HoppingMatrix.h"
#include "bgq/bgq_qpx.h"
#include "bgq/bgq_utils.h"
#include <omp.h>
#include "bgq/mypapi.h"

#ifdef XLC
#include <l1p/pprefetch.h>
#include <l1p/sprefetch.h>
#endif

#define COUNTOF(arr) (sizeof(arr) / sizeof((arr)[0]))

static void exec_bench();
static void exec_fakebench();
static void check_correctness_double(bool nocom);
static void check_correctness_float();
static void exec_table_asm(bool sloppiness, bgq_hmflags additional_opts, int select);

typedef bgq_spinorfield bgq_spinorfield_double;
typedef bgq_spinorfield bgq_spinorfield_float;

typedef bgq_gaugefield bgq_gaugefield_double;
typedef bgq_gaugefield bgq_gaugefield_float;

typedef bgq_gaugesite bgq_gaugesite_double;
typedef bgq_gaugesite bgq_gaugesite_float;

typedef void (*hm_func_double)(bool isOdd, bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bgq_hmflags opts);
typedef void (*hm_func_float)(bool isOdd, bgq_spinorfield_float targetfield, bgq_spinorfield_float spinorfield, bgq_gaugefield_float gaugefield, bgq_hmflags opts);
/* END MK */



#ifdef PARALLELT
#  define SLICE (LX*LY*LZ/2)
#elif defined PARALLELXT
#  define SLICE ((LX*LY*LZ/2)+(T*LY*LZ/2))
#elif defined PARALLELXYT
#  define SLICE ((LX*LY*LZ/2)+(T*LY*LZ/2) + (T*LX*LZ/2))
#elif defined PARALLELXYZT
#  define SLICE ((LX*LY*LZ/2)+(T*LY*LZ/2) + (T*LX*LZ/2) + (T*LX*LY/2))
#elif defined PARALLELX
#  define SLICE ((LY*LZ*T/2))
#elif defined PARALLELXY
#  define SLICE ((LY*LZ*T/2) + (LX*LZ*T/2))
#elif defined PARALLELXYZ
#  define SLICE ((LY*LZ*T/2) + (LX*LZ*T/2) + (LX*LY*T/2))
#endif


int check_xchange();

int main(int argc,char *argv[])
{
  int j,j_max,k,k_max = 1;
#ifdef HAVE_LIBLEMON
  paramsXlfInfo *xlfInfo;
#endif
  int status = 0;

  static double t1,t2,dt,sdt,dts,qdt,sqdt;
  double antioptaway=0.0;

#ifdef MPI
  static double dt2;

  DUM_DERI = 6;
  DUM_SOLVER = DUM_DERI+2;
  DUM_MATRIX = DUM_SOLVER+6;
  NO_OF_SPINORFIELDS = DUM_MATRIX+2;

#  ifdef OMP
  int mpi_thread_provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &mpi_thread_provided);
#  else
  MPI_Init(&argc, &argv);
#  endif
  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);

#else
  g_proc_id = 0;
#endif

  g_rgi_C1 = 1.;

    /* Read the input file */
  if((status = read_input("benchmark.input")) != 0) {
    fprintf(stderr, "Could not find input file: benchmark.input\nAborting...\n");
    exit(-1);
  }

#ifdef OMP
  if(omp_num_threads > 0)
  {
     omp_set_num_threads(omp_num_threads);
  }
  else {
    if( g_proc_id == 0 )
      printf("# No value provided for OmpNumThreads, running in single-threaded mode!\n");

    omp_num_threads = 1;
    omp_set_num_threads(omp_num_threads);
  }

  init_omp_accumulators(omp_num_threads);
#endif

  tmlqcd_mpi_init(argc, argv);

  if(g_proc_id==0) {
#ifdef SSE
    printf("# The code was compiled with SSE instructions\n");
#endif
#ifdef SSE2
    printf("# The code was compiled with SSE2 instructions\n");
#endif
#ifdef SSE3
    printf("# The code was compiled with SSE3 instructions\n");
#endif
#ifdef P4
    printf("# The code was compiled for Pentium4\n");
#endif
#ifdef OPTERON
    printf("# The code was compiled for AMD Opteron\n");
#endif
#ifdef _GAUGE_COPY
    printf("# The code was compiled with -D_GAUGE_COPY\n");
#endif
#ifdef BGL
    printf("# The code was compiled for Blue Gene/L\n");
#endif
#ifdef BGP
    printf("# The code was compiled for Blue Gene/P\n");
#endif
#ifdef _USE_HALFSPINOR
    printf("# The code was compiled with -D_USE_HALFSPINOR\n");
#endif
#ifdef _USE_SHMEM
    printf("# The code was compiled with -D_USE_SHMEM\n");
#  ifdef _PERSISTENT
    printf("# The code was compiled for persistent MPI calls (halfspinor only)\n");
#  endif
#endif
#ifdef MPI
#  ifdef _NON_BLOCKING
    printf("# The code was compiled for non-blocking MPI calls (spinor and gauge)\n");
#  endif
#endif
    printf("\n");
    fflush(stdout);
  }


#ifdef _GAUGE_COPY
  init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 1);
#else
  init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 0);
#endif
  init_geometry_indices(VOLUMEPLUSRAND + g_dbw2rand);

  if(even_odd_flag) {
    j = init_spinor_field(VOLUMEPLUSRAND/2, 2*k_max+1);
  }
  else {
    j = init_spinor_field(VOLUMEPLUSRAND, 2*k_max);
  }

  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for spinor fields! Aborting...\n");
    exit(0);
  }
  j = init_moment_field(VOLUME, VOLUMEPLUSRAND + g_dbw2rand);
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for moment fields! Aborting...\n");
    exit(0);
  }

  if(g_proc_id == 0) {
    fprintf(stdout,"# The number of processes is %d \n",g_nproc);
    printf("# The lattice size is %d x %d x %d x %d\n",
	   (int)(T*g_nproc_t), (int)(LX*g_nproc_x), (int)(LY*g_nproc_y), (int)(g_nproc_z*LZ));
    printf("# The local lattice size is %d x %d x %d x %d\n",
	   (int)(T), (int)(LX), (int)(LY),(int) LZ);
    if(even_odd_flag) {
      printf("# benchmarking the even/odd preconditioned Dirac operator\n");
    }
    else {
      printf("# benchmarking the standard Dirac operator\n");
    }
    fflush(stdout);
  }

  /* define the geometry */
  geometry();
  /* define the boundary conditions for the fermion fields */
  boundary(g_kappa);

#ifdef _USE_HALFSPINOR
  j = init_dirac_halfspinor();
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for halfspinor fields! Aborting...\n");
    exit(0);
  }
  if(g_sloppy_precision_flag == 1) {
    g_sloppy_precision = 1;
    j = init_dirac_halfspinor32();
    if ( j!= 0) {
      fprintf(stderr, "Not enough memory for 32-Bit halfspinor fields! Aborting...\n");
      exit(0);
    }
  }
#  if (defined _PERSISTENT)
  init_xchange_halffield();
#  endif
#endif

  status = check_geometry();
  if (status != 0) {
    fprintf(stderr, "Checking of geometry failed. Unable to proceed.\nAborting....\n");
    exit(1);
  }
#if (defined MPI && !(defined _USE_SHMEM))
  check_xchange();
#endif

  start_ranlux(1, 123456);
  random_gauge_field(reproduce_randomnumber_flag);

#ifdef MPI
  /*For parallelization: exchange the gaugefield */
  xchange_gauge(g_gauge_field);
#endif

  /* BEGIN MK */
  assert(even_odd_flag);


	//bgq_init_gaugefield_allprec();
	//bgq_init_spinorfields_allprec(2 * k_max + 1, 0);
	//bgq_hm_init_allprec();

	//bgq_update_backward_gauge();


#ifdef _GAUGE_COPY
	//update_backward_gauge();
#endif

	bool done = false;
#pragma omp parallel for schedule(static) ordered firstprivate(done)
	for (int j = 0; j < omp_get_max_threads(); j+=1)
	{
		if (g_proc_id == 0 && !done) {
#pragma omp ordered
			{
				printf("MK Here is omp_id %2d, cpu_id %2d running on SMT-Thread %d of Core %2d\n", omp_get_thread_num(), Kernel_ProcessorID(), Kernel_ProcessorThreadID(), Kernel_ProcessorCoreID());
				done = true;
			}
		}
	}


	for (int k = 0; k < k_max; k+=1) {
		/*initialize the pseudo-fermion fields*/
		random_spinor_field(g_spinor_field[k], VOLUME / 2, 0);
		//bgq_transfer_spinorfield_double(true, g_spinorfields_double[k], g_spinor_field[k]);

		//double compare_transfer = bgq_spinorfield_compare_double(true, g_spinorfields_double[k], g_spinor_field[k], false);
		//assert(compare_transfer == 0 /* fields should be bit-identical*/);
	}


#if !BGQ_REPLACE
	//check_correctness_double(true);
	//check_correctness_double(false);
	//check_correctness_float();
#endif

	exec_table_asm(false, 0, 0);
	exec_table_asm(false, 0, 1);
	exec_table_asm(false, 0, 2);

	//exec_fakebench();
	exec_bench();
  /* END MK */



    while(sdt < 30.) {
#ifdef MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      t1 = gettime();
      antioptaway=0.0;
      for (j=0;j<j_max;j++) {
        for (k=0;k<k_max;k++) {
          Hopping_Matrix(0, g_spinor_field[k+k_max], g_spinor_field[k]);
          Hopping_Matrix(1, g_spinor_field[2*k_max], g_spinor_field[k+k_max]);
          antioptaway+=creal(g_spinor_field[2*k_max][0].s0.c0);
        }
      }
      t2 = gettime();
      dt = t2-t1;
#ifdef MPI
      MPI_Allreduce (&dt, &sdt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      sdt = dt;
#endif
      qdt=dt*dt;
#ifdef MPI
      MPI_Allreduce (&qdt, &sqdt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      sqdt = qdt;
#endif
      sdt=sdt/((double)g_nproc);
      sqdt=sqrt(sqdt/g_nproc-sdt*sdt);
      j_max*=2;
    }
    j_max=j_max/2;
    dts=dt;
    sdt=1.0e6f*sdt/((double)(k_max*j_max*(VOLUME)));
    sqdt=1.0e6f*sqdt/((double)(k_max*j_max*(VOLUME)));

    if(g_proc_id==0) {
      printf("# The following result is just to make sure that the calculation is not optimized away: %e\n", antioptaway);
      printf("# Total compute time %e sec, variance of the time %e sec. (%d iterations).\n", sdt, sqdt, j_max);
      printf("# Communication switched on:\n# (%d Mflops [%d bit arithmetic])\n", (int)(1608.0f/sdt),(int)sizeof(spinor)/3);
#ifdef OMP
      printf("# Mflops per OpenMP thread ~ %d\n",(int)(1608.0f/(omp_num_threads*sdt)));
#endif
      printf("\n");
      fflush(stdout);
    }

#ifdef MPI
    /* isolated computation */
    t1 = gettime();
    antioptaway=0.0;
    for (j=0;j<j_max;j++) {
      for (k=0;k<k_max;k++) {
        Hopping_Matrix_nocom(0, g_spinor_field[k+k_max], g_spinor_field[k]);
        Hopping_Matrix_nocom(1, g_spinor_field[2*k_max], g_spinor_field[k+k_max]);
        antioptaway += creal(g_spinor_field[2*k_max][0].s0.c0);
      }
    }
    t2 = gettime();
    dt2 = t2-t1;
    /* compute the bandwidth */
    dt=dts-dt2;
    MPI_Allreduce (&dt, &sdt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sdt=sdt/((double)g_nproc);
    MPI_Allreduce (&dt2, &dt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dt=dt/((double)g_nproc);
    dt=1.0e6f*dt/((double)(k_max*j_max*(VOLUME)));
    if(g_proc_id==0) {
      printf("# The following result is printed just to make sure that the calculation is not optimized away: %e\n",antioptaway);
      printf("# Communication switched off: \n# (%d Mflops [%d bit arithmetic])\n", (int)(1608.0f/dt),(int)sizeof(spinor)/3);
#ifdef OMP
      printf("# Mflops per OpenMP thread ~ %d\n",(int)(1608.0f/(omp_num_threads*dt)));
#endif
      printf("\n");
      fflush(stdout);
    }
    sdt=sdt/((double)k_max);
    sdt=sdt/((double)j_max);
    sdt=sdt/((double)(2*SLICE));
    if(g_proc_id==0) {
      printf("# The size of the package is %d bytes.\n",(SLICE)*192);
#ifdef _USE_HALFSPINOR
      printf("# The bandwidth is %5.2f + %5.2f MB/sec\n", 192./sdt/1024/1024, 192./sdt/1024./1024);
#else
      printf("# The bandwidth is %5.2f + %5.2f MB/sec\n", 2.*192./sdt/1024/1024, 2.*192./sdt/1024./1024);
#endif
    }
#endif
    fflush(stdout);


#ifdef HAVE_LIBLEMON
  if(g_proc_id==0) {
    printf("# Performing parallel IO test ...\n");
  }
  xlfInfo = construct_paramsXlfInfo(0.5, 0);
  write_gauge_field( "conf.test", 64, xlfInfo);
  free(xlfInfo);
  if(g_proc_id==0) {
    printf("# done ...\n");
  }
#endif


#ifdef MPI
  MPI_Finalize();
#endif
#ifdef OMP
  free_omp_accumulators();
#endif
  free_gauge_field();
  free_geometry_indices();
  free_spinor_field();
  free_moment_field();
  return(0);
}














complexdouble nooptaway;
#if BGQ_PREFETCH_LIST
L1P_Pattern_t *(g_patternhandle[64]);
#endif






typedef struct {
	double avgtime;
	double localrmstime;
	double globalrmstime;

	double lups;
	double flops;

	double error;

	mypapi_counters counters;
	bgq_hmflags opts;
} benchstat;




static void Hopping_Matrix_switch(const int ieo, spinor * const l, spinor * const k, bool nocom) {
	if (nocom) {
		Hopping_Matrix_nocom(ieo, l, k);
	} else {
		Hopping_Matrix(ieo, l, k);
	}
}


inline void check_correctness_double(bool nocom) {
	master_print("MK Checking double precision%s correctness...\n", nocom ? " (no communication)" : "");

//#pragma omp parallel
	{
		int k = 0;
		int k_max = 1;
		//bgq_hmflags hmflags = hm_nocom * nocom | hm_nooverlap;
		double compare_even_before=0;
		//double compare_even=0;

//#pragma omp master
		{
			//bgq_transfer_spinorfield_double(true, g_spinorfields_double[k], g_spinor_field[k]);
			//compare_even_before = bgq_spinorfield_compare_double(true, g_spinorfields_double[k], g_spinor_field[k], false);
			assert(compare_even_before == 0 /* should be bitwise identical */);

			//bgq_initbgqref();
		}
//#pragma omp barrier
		//bgq_HoppingMatrix_double(false, g_spinorfields_double[k + k_max], g_spinorfields_double[k], g_gaugefield_double, hmflags);
//#pragma omp master
		{
			//master_print("MK HM_orig start\n");
			Hopping_Matrix_switch(0, g_spinor_field[k + k_max], g_spinor_field[k], nocom);
			//master_print("MK HM_orig end\n");
			//bgq_savebgqref();
			//__asm__("int3");
			//compare_even = bgq_spinorfield_compare_double(false, g_spinorfields_double[k + k_max], g_spinor_field[k + k_max], false);
			//assert(compare_even < 0.001);
		}
//#pragma omp barrier
		//bgq_HoppingMatrix_double(true, g_spinorfields_double[2 * k_max], g_spinorfields_double[k + k_max], g_gaugefield_double, hmflags);
//#pragma omp master
		{
			Hopping_Matrix_switch(1, g_spinor_field[k + 2 * k_max], g_spinor_field[k + k_max], nocom);
			//double compare_odd = bgq_spinorfield_compare_double(true, g_spinorfields_double[k + 2 * k_max], g_spinor_field[k + 2 * k_max], false);
			//assert(compare_odd < 0.001);

			//master_print("Numerical instability between double precision implementations: even %e, odd %e\n", compare_even, compare_odd);
		}
	}
}

inline void check_correctness_float() {
	master_print("MK Checking float precision correctness...\n");

//#pragma omp parallel
	{
		int k = 0;
		int k_max = 1;
		//bgq_hmflags hmflags = 0;
		double compare_even_before=0;
		//double compare_even;

//#pragma omp master
		{
			//bgq_transfer_spinorfield_float(true, g_spinorfields_float[k], g_spinor_field[k]);
			//compare_even_before = bgq_spinorfield_compare_float(true, g_spinorfields_float[k], g_spinor_field[k], false);
			assert(compare_even_before == 0 /* should be bitwise identical */);

			//bgq_initbgqref();
		}
//#pragma omp barrier
		//bgq_HoppingMatrix_float(false, g_spinorfields_float[k + k_max], g_spinorfields_float[k], g_gaugefield_float, hmflags);
//#pragma omp master
		{
			Hopping_Matrix(0, g_spinor_field[k + k_max], g_spinor_field[k]);
			//bgq_savebgqref();
			//compare_even = bgq_spinorfield_compare_float(false, g_spinorfields_float[k + k_max], g_spinor_field[k + k_max], false);
			//assert(compare_even < 0.001);
		}
//#pragma omp barrier
		//bgq_HoppingMatrix_float(true, g_spinorfields_float[2 * k_max], g_spinorfields_float[k + k_max], g_gaugefield_float, hmflags);
//#pragma omp master
		{
			Hopping_Matrix(1, g_spinor_field[k + 2 * k_max], g_spinor_field[k + k_max]);
			//double compare_odd = bgq_spinorfield_compare_float(true, g_spinorfields_float[k + 2 * k_max], g_spinor_field[k + 2 * k_max], false);
			//assert(compare_odd < 0.001);

			//master_print("Numerical instability between float precision implementations: even %g, odd %g\n", compare_even, compare_odd);
		}
	}
}


static double runcheck(bool sloppyprec, bgq_hmflags hmflags, int k_max) {
	const int k = 0;
	hmflags = hmflags & ~hm_nokamul;
	double result;
	double compare_even_before=0;
	double compare_even=0;
	double compare_odd=0;

	if (!sloppyprec) {
//#pragma omp master
		{
		//memset(g_spinorfields_double[k + k_max], 0, sizeof(bgq_spinorsite_double) * VOLUME_SITES);
		//memset(g_spinorfields_double[k + 2*k_max], 0, sizeof(bgq_spinorsite_double) * VOLUME_SITES);

		//bgq_transfer_spinorfield_double(true, g_spinorfields_double[k], g_spinor_field[k]);
		//compare_even_before = bgq_spinorfield_compare_double(true, g_spinorfields_double[k], g_spinor_field[k], false);
		}
		//bgq_HoppingMatrix_double(false, g_spinorfields_double[k + k_max], g_spinorfields_double[k], g_gaugefield_double, hmflags);
//#pragma omp master
		{
		Hopping_Matrix_switch(0, g_spinor_field[k + k_max], g_spinor_field[k], hmflags & hm_nocom);
		//compare_even = bgq_spinorfield_compare_double(false, g_spinorfields_double[k + k_max], g_spinor_field[k + k_max], true);
		}
		//bgq_HoppingMatrix_double(true, g_spinorfields_double[2 * k_max], g_spinorfields_double[k + k_max], g_gaugefield_double, hmflags);
//#pragma omp master
		{
		Hopping_Matrix_switch(1, g_spinor_field[k + 2*k_max], g_spinor_field[k + k_max], hmflags & hm_nocom);
		//compare_odd = bgq_spinorfield_compare_double(true, g_spinorfields_double[k + 2*k_max], g_spinor_field[k + 2*k_max], true);

		result = max(max(compare_even_before, compare_even), compare_odd);
		}
	} else {
//#pragma omp master
		{
		//memset(g_spinorfields_float[k + k_max], 0, sizeof(bgq_spinorsite_float) * VOLUME_SITES);
		//memset(g_spinorfields_float[k + 2*k_max], 0, sizeof(bgq_spinorsite_float) * VOLUME_SITES);

		//bgq_transfer_spinorfield_float(true, g_spinorfields_float[k], g_spinor_field[k]);
		//compare_even_before = bgq_spinorfield_compare_float(true, g_spinorfields_float[k], g_spinor_field[k], true);
		}
		//bgq_HoppingMatrix_float(false, g_spinorfields_float[k + k_max], g_spinorfields_float[k], g_gaugefield_float, hmflags);
//#pragma omp master
		{
		Hopping_Matrix_switch(0, g_spinor_field[k + k_max], g_spinor_field[k], hmflags & hm_nocom);
		//compare_even = bgq_spinorfield_compare_float(false, g_spinorfields_float[k + k_max], g_spinor_field[k + k_max], true);
		}
		//bgq_HoppingMatrix_float(true, g_spinorfields_float[2 * k_max], g_spinorfields_float[k + k_max], g_gaugefield_float, hmflags);
//#pragma omp master
		{
		Hopping_Matrix_switch(1, g_spinor_field[k + 2*k_max], g_spinor_field[k + k_max], hmflags & hm_nocom);
		//compare_odd = bgq_spinorfield_compare_float(true, g_spinorfields_float[k + 2*k_max], g_spinor_field[k + 2*k_max], true);

		result = max(max(compare_even_before, compare_even), compare_odd);
		}
	}

	return result;
}


static benchstat runbench(int k_max, int j_max, bool sloppyprec, int ompthreads, bgq_hmflags opts, hm_func_double hm_double, hm_func_float hm_float) {
	//const bool nocom = opts & hm_nocom;
	//const bool nooverlap = opts & hm_nooverlap;
	const bool nokamul = opts & hm_nokamul;
	//const bool noprefetchlist = opts & hm_noprefetchlist;
	//const bool noprefetchstream = opts & hm_noprefetchstream;
	//const bool noprefetchexplicit = opts & hm_noprefetchexplicit;
	const bool noweylsend = opts & hm_noweylsend;
	const bool nobody = opts & hm_nobody;
	const bool nosurface = opts & hm_nosurface;
	//const bool experimental = opts & hm_experimental;
	//bgq_hmflags implicitprefetch = opts & (hm_prefetchimplicitdisable | hm_prefetchimplicitoptimistic | hm_prefetchimplicitconfirmed);

	// Setup options
#ifdef XLC
	L1P_StreamPolicy_t pol;
	switch (implicitprefetch) {
	case hm_prefetchimplicitdisable:
		pol = noprefetchstream ? L1P_stream_disable : L1P_confirmed_or_dcbt/*No option to selectively disable implicit stream*/;
		break;
	case hm_prefetchimplicitoptimistic:
		pol = L1P_stream_optimistic;
		break;
	case hm_prefetchimplicitconfirmed:
		pol = noprefetchstream ? L1P_stream_confirmed : L1P_confirmed_or_dcbt;
		break;
	default:
		// Default setting
		pol = L1P_confirmed_or_dcbt;
		break;
	}


	omp_set_num_threads(ompthreads);
#pragma omp parallel
	{
		L1P_CHECK(L1P_SetStreamPolicy(pol));
	}

#pragma omp parallel
	{
		// Test whether it persists between parallel sections
		//L1P_StreamPolicy_t getpol = 0;
		//L1P_CHECK(L1P_GetStreamPolicy(&getpol));
		//if (getpol != pol)
		//	fprintf(stderr, "MK StreamPolicy not accepted\n");

		//L1P_CHECK(L1P_SetStreamDepth());
		//L1P_CHECK(L1P_SetStreamTotalDepth());

		//L1P_GetStreamDepth
		//L1P_GetStreamTotalDepth


		// Peter Boyle's setting
		// Note: L1P_CFG_PF_USR_pf_stream_establish_enable is never set in L1P_SetStreamPolicy
		//uint64_t *addr = ((uint64_t*)(Kernel_L1pBaseAddress() + L1P_CFG_PF_USR_ADJUST));
		//*addr |=  L1P_CFG_PF_USR_pf_stream_est_on_dcbt | L1P_CFG_PF_USR_pf_stream_optimistic | L1P_CFG_PF_USR_pf_stream_prefetch_enable | L1P_CFG_PF_USR_pf_stream_establish_enable; // Enable everything???
	}
#endif

#if BGQ_REPLACE
	// Cannot perform test, reference implementation has no been linked
	double err = 0;
#else
	// Error checking
	double err = runcheck(sloppyprec, opts, k_max);
#endif


	double localsumtime = 0;
	double localsumsqtime = 0;
	mypapi_counters counters;
	counters.init = false;
	int iterations = 1 + j_max;
	if (iterations < 1 + MYPAPI_SETS)
		iterations = 1 + MYPAPI_SETS;

	for (int j = 0; j < iterations; j += 1) {
		//master_print("Starting iteration %d of %d\n", j+1, iterations);
		bool isWarmup = (j == 0);
		//bool isLast = (j == iterations-1);
		bool isPapi = !isWarmup && (j >= iterations - MYPAPI_SETS);
		bool isJMax = !isWarmup && (j >= iterations - j_max);

		if (isPapi) {
			mypapi_start(j - (iterations - MYPAPI_SETS));
		}
		double start_time = MPI_Wtime();
		if (sloppyprec) {
			for (int k = 0; k < 1; k += 1) {
				//hm_float(false, g_spinorfields_float[k + k_max], g_spinorfields_float[k], g_gaugefield_float, opts);
				//hm_float(true, g_spinorfields_float[2 * k_max], g_spinorfields_float[k + k_max], g_gaugefield_float, opts);
			}
		} else {
			for (int k = 0; k < 1; k += 1) {
				//hm_double(false, g_spinorfields_double[k + k_max], g_spinorfields_double[k], g_gaugefield_double, opts);
				//hm_double(true, g_spinorfields_double[2 * k_max], g_spinorfields_double[k + k_max], g_gaugefield_double, opts);
			}
		}
		double end_time = MPI_Wtime();
		mypapi_counters curcounters;
		if (isPapi) {
			curcounters = mypapi_stop();
			counters = mypapi_merge_counters(&counters, &curcounters);
		}

		double runtime = end_time - start_time;
		if (isJMax) {
			localsumtime += runtime;
			localsumsqtime += sqr(runtime);
		}
	}


	double localavgtime = localsumtime / j_max;
	double localavgsqtime = sqr(localavgtime);
	double localrmstime = sqrt((localsumsqtime / j_max) - localavgsqtime);

	double localtime[] = { localavgtime, localavgsqtime, localrmstime };
	double sumreduce[3] = { -1, -1, -1 };
	MPI_Allreduce(&localtime, &sumreduce, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	double sumtime = sumreduce[0];
	double sumsqtime = sumreduce[1];
	double sumrmstime = sumreduce[2];

	double avgtime = sumtime / g_nproc;
	double avglocalrms = sumrmstime / g_nproc;
	double rmstime = sqrt((sumsqtime / g_nproc) - sqr(avgtime));

	int lups = LOCAL_LT * LOCAL_LX * LOCAL_LY * LOCAL_LZ;
	int lups_body=0;// = BODY_ZLINES * PHYSICAL_LP*PHYSICAL_LK * LOCAL_LZ;
	int lups_surface=0;// = SURFACE_ZLINES * PHYSICAL_LP*PHYSICAL_LK * LOCAL_LZ;
	assert(lups == lups_body + lups_surface);
	assert(lups == VOLUME);

	int flops_per_lup_body = 0;
	if (!nobody) {
		flops_per_lup_body += 1320;
		if (!nokamul)
			flops_per_lup_body += 8 * 2 * 3 * 6;
	}

	int flops_per_lup_surface = 0;
	if (!noweylsend)
		flops_per_lup_surface += 6 * 2 * 2;
	if (!nosurface) {
		flops_per_lup_surface += 1320 - (6*2*2);
		if (!nokamul)
			flops_per_lup_surface += 8 * 2 * 3 * 6;
	}
	double flops = ((double)flops_per_lup_body * (double)lups_body) + ((double)flops_per_lup_surface * (double)lups_surface);


	benchstat result;
	result.avgtime = avgtime;
	result.localrmstime = avglocalrms;
	result.globalrmstime = rmstime;
	result.lups = lups;
	result.flops = flops / avgtime;
	result.error = err;
	result.counters = counters;
	result.opts = opts;
	return result;
}



static bool kamuls[] = { false, true };
static char *kamuls_desc[] = { "dslash", "kamul" };

static bool sloppinessess[] = { false, true };
static char *sloppinessess_desc[] = { "double", "float" };

static int omp_threads[] = { 1, 2, 4, 8, 16, 32, 33, 48, 56, 64 };
static char *omp_threads_desc[] = { "1", "2", "4", "8", "16", "32", "33", "48", "56", "64" };


static bgq_hmflags flags[] = {
//		               hm_noprefetchlist | hm_noprefetchstream | hm_noprefetchexplicit,
//		hm_nooverlap | hm_noprefetchlist | hm_noprefetchstream | hm_noprefetchexplicit,
		hm_nocom | hm_noprefetchlist | hm_noprefetchstream | hm_noprefetchexplicit,
//		hm_fixedoddness | hm_nocom | hm_noprefetchlist | hm_noprefetchstream | hm_noprefetchexplicit,
		hm_nocom | hm_nosurface | hm_noweylsend | hm_noprefetchlist | hm_noprefetchstream | hm_noprefetchexplicit,
		hm_nocom | hm_nobody | hm_noprefetchlist | hm_noprefetchstream | hm_noprefetchexplicit,
		hm_nocom | hm_noweylsend | hm_noprefetchlist | hm_noprefetchstream | hm_noprefetchexplicit,
		hm_nocom | hm_noweylsend | hm_noprefetchlist | hm_noprefetchstream | hm_noprefetchexplicit | hm_prefetchimplicitdisable,
		hm_nocom | hm_noweylsend | hm_noprefetchlist | hm_noprefetchstream | hm_noprefetchexplicit | hm_prefetchimplicitconfirmed,
		hm_nocom | hm_noweylsend | hm_noprefetchlist | hm_noprefetchstream | hm_noprefetchexplicit | hm_prefetchimplicitoptimistic,
//		hm_nocom | hm_noweylsend | hm_noprefetchlist                       | hm_noprefetchexplicit | hm_prefetchimplicitdisable,
		hm_nocom | hm_noweylsend | hm_noprefetchlist | hm_noprefetchstream                         | hm_prefetchimplicitdisable,
//		hm_nocom | hm_noweylsend                     | hm_noprefetchstream | hm_noprefetchexplicit | hm_prefetchimplicitdisable
		hm_nocom | hm_noweylsend                     | hm_noprefetchstream | hm_noprefetchexplicit | hm_prefetchimplicitdisable | hm_l1pnonstoprecord,
		hm_nocom | hm_noweylsend                     | hm_noprefetchstream | hm_noprefetchexplicit | hm_prefetchimplicitdisable | hm_experimental
		};
static char* flags_desc[] = {
//		"Com async",
//		"Com sync",
		"Com off",
//		"^fixedodd",
		"bodyonly",
		"surfaceonly",
		"pf def",
		"pf disable",
		"pf confirmed",
		"pf optimistic",
//		"pf stream",
		"pf explicit",
//		"pf list",
		"pf list nonstop",
		"pf list exp"
		};


static void print_repeat(const char * const str, const int count) {
	if (g_proc_id == 0) {
		for (int i = 0; i < count; i += 1) {
			printf("%s", str);
		}
	}
}


#define CELLWIDTH 15
#define SCELLWIDTH TOSTRING(CELLWIDTH)


static void print_stats(benchstat stats[COUNTOF(flags)], int bodySites, int surfaceSites, int haloSites) {
#if PAPI
	int threads = omp_get_num_threads();

	for (mypapi_interpretations j = 0; j < __pi_COUNT; j+=1) {
		printf("%10s|", "");
		char *desc = NULL;

		for (int i3 = 0; i3 < COUNTOF(flags); i3 += 1) {
			char str[80];
			str[0] = '\0';
			benchstat *stat = &stats[i3];
			bgq_hmflags opts = stat->opts;

			//uint64_t bodySites = BODY_SITES*4;
			//uint64_t surfaceSites = SURFACE_SITES*4;
			uint64_t sites = bodySites + surfaceSites;

			double nCycles = stats[i3].counters.native[PEVT_CYCLES];
			double nCoreCycles = stats[i3].counters.corecycles;
			double nNodeCycles = stats[i3].counters.nodecycles;
			double nInstructions = stats[i3].counters.native[PEVT_INST_ALL];
			double nStores = stats[i3].counters.native[PEVT_LSU_COMMIT_STS];
			double nL1IStalls = stats[i3].counters.native[PEVT_LSU_COMMIT_STS];
			double nL1IBuffEmpty = stats[i3].counters.native[PEVT_IU_IBUFF_EMPTY_CYC];

			double nIS1Stalls = stats[i3].counters.native[PEVT_IU_IS1_STALL_CYC];
			double nIS2Stalls = stats[i3].counters.native[PEVT_IU_IS2_STALL_CYC];

			double nCachableLoads = stats[i3].counters.native[PEVT_LSU_COMMIT_CACHEABLE_LDS];
			double nL1Misses = stats[i3].counters.native[PEVT_LSU_COMMIT_LD_MISSES];
			double nL1Hits = nCachableLoads - nL1Misses;

			double nL1PMisses = stats[i3].counters.native[PEVT_L1P_BAS_MISS];
			double nL1PHits = stats[i3].counters.native[PEVT_L1P_BAS_HIT];
			double nL1PAccesses = nL1PHits + nL1PMisses;

			double nL2Misses = stats[i3].counters.native[PEVT_L2_MISSES];
			double nL2Hits = stats[i3].counters.native[PEVT_L2_HITS];
			double nL2Accesses = nL2Misses + nL2Hits;

			double nDcbtHits = stats[i3].counters.native[PEVT_LSU_COMMIT_DCBT_HITS];
			double nDcbtMisses = stats[i3].counters.native[PEVT_LSU_COMMIT_DCBT_MISSES];
			double nDcbtAccesses = nDcbtHits + nDcbtMisses;

			double nXUInstr = stats[i3].counters.native[PEVT_INST_XU_ALL];
			double nAXUInstr = stats[i3].counters.native[PEVT_INST_QFPU_ALL];
			double nXUAXUInstr = nXUInstr + nAXUInstr;

			double nNecessaryInstr = 0;
			if (!(opts & hm_nobody))
				nNecessaryInstr += bodySites * (240/*QFMA*/ + 180/*QMUL+QADD*/ + 180/*LD+ST*/)/2;
			if (!(opts & hm_noweylsend))
				nNecessaryInstr += haloSites * (2*3*2/*QMUL+QADD*/ + 4*3/*LD*/ + 2*3/*ST*/)/2;
			if (!(opts & hm_nosurface))
				nNecessaryInstr += surfaceSites * (240/*QFMA*/ + 180/*QMUL+QADD*/ + 180/*LD+ST*/ - 2*3*2/*QMUL+QADD*/ - 4*3/*LD*/ + 2*3/*LD*/)/2;
			if (!(opts & hm_nokamul))
				nNecessaryInstr += sites * (8*2*3*1/*QFMA*/ + 8*2*3*1/*QMUL*/)/2;

			uint64_t nL1PListStarted = stats[i3].counters.native[PEVT_L1P_LIST_STARTED];
			uint64_t nL1PListAbandoned= stats[i3].counters.native[PEVT_L1P_LIST_ABANDON];
			uint64_t nL1PListMismatch= stats[i3].counters.native[PEVT_L1P_LIST_MISMATCH];
			uint64_t nL1PListSkips = stats[i3].counters.native[PEVT_L1P_LIST_SKIP];
			uint64_t nL1PListOverruns = stats[i3].counters.native[PEVT_L1P_LIST_CMP_OVRUN_PREFCH];

			double nL1PLatePrefetchStalls = stats[i3].counters.native[PEVT_L1P_BAS_LU_STALL_LIST_RD_CYC];

			uint64_t nStreamDetectedStreams = stats[i3].counters.native[PEVT_L1P_STRM_STRM_ESTB];
			double nL1PSteamUnusedLines = stats[i3].counters.native[PEVT_L1P_STRM_EVICT_UNUSED];
			double nL1PStreamPartiallyUsedLines = stats[i3].counters.native[PEVT_L1P_STRM_EVICT_PART_USED];
			double nL1PStreamLines = stats[i3].counters.native[PEVT_L1P_STRM_LINE_ESTB];
			double nL1PStreamHits = stats[i3].counters.native[PEVT_L1P_STRM_HIT_LIST];

			switch (j) {
			case pi_detstreams:
				desc = "Detected streams";
				snprintf(str, sizeof(str), "%llu", nStreamDetectedStreams);
				break;
			case pi_l1pstreamunusedlines:
				desc = "Unused (partially) lines";
				snprintf(str, sizeof(str), "%.2f%% (%.2f%%)", 100.0 * nL1PSteamUnusedLines / nL1PStreamLines, 100.0 * nL1PStreamPartiallyUsedLines / nL1PStreamLines);
				break;
			case pi_l1pstreamhitinl1p:
				desc = "Loads that hit in L1P stream";
				snprintf(str, sizeof(str), "%.2f %%", 100.0 * nL1PStreamHits / nCachableLoads);
				break;
			case pi_cpi:
				desc = "Cycles per instruction (Thread)";
				snprintf(str, sizeof(str), "%.3f", nCycles / nInstructions);
				break;
			case pi_corecpi:
				desc = "Cycles per instruction (Core)";
				snprintf(str, sizeof(str), "%.3f", nCoreCycles / nInstructions);
				break;
			case pi_l1istalls:
				desc = "Empty instr buffer";
				snprintf(str, sizeof(str), "%.2f %%", nL1IBuffEmpty / nCycles);
				break;
			case pi_is1stalls:
				desc = "IS1 Stalls (dependency)";
				snprintf(str, sizeof(str), "%.2f %%", 100 * nIS1Stalls / nCycles);
				break;
			case pi_is2stalls:
				desc = "IS2 Stalls (func unit)";
				snprintf(str, sizeof(str), "%.2f %%", 100 * nIS2Stalls / nCycles);
				break;
			case pi_hitinl1:
				desc = "Loads that hit in L1";
				snprintf(str, sizeof(str), "%.2f %%",  100 * nL1Hits / nCachableLoads);
				break;
			case pi_l1phitrate:
				desc = "L1P hit rate";
				snprintf(str, sizeof(str), "%.2f %%",  100 * nL1PHits / nL1PAccesses);
				break;
			case pi_overhead:
				desc = "Instr overhead";
				snprintf(str, sizeof(str), "%.2f %%", 100 * (nInstructions - nNecessaryInstr) / nInstructions);
				break;
			//case pi_hitinl1p:
			//	desc = "Loads that hit in L1P";
			//	snprintf(str, sizeof(str), "%f %%" ,  100 * nL1PHits / nCachableLoads);
			//	break;
			case pi_l2hitrate:
				desc = "L2 hit rate";
				snprintf(str, sizeof(str), "%.2f %%", 100 * nL2Hits / nL2Accesses);
				break;
			case pi_dcbthitrate:
				desc = "dcbt hit rate";
				snprintf(str, sizeof(str), "%.2f %%", 100 * nDcbtHits / nDcbtAccesses);
				break;
			case pi_axufraction:
				desc = "FXU instrs";
				snprintf(str, sizeof(str), "%.2f %%", 100 * nAXUInstr / nXUAXUInstr);
				break;
			case pi_l1pliststarted:
				desc = "List prefetch started";
				snprintf(str, sizeof(str), "%llu", nL1PListStarted);
				break;
			case pi_l1plistabandoned:
				desc = "List prefetch abandoned";
				snprintf(str, sizeof(str), "%llu", nL1PListAbandoned);
				break;
			case pi_l1plistmismatch:
				desc = "List prefetch mismatch";
				snprintf(str, sizeof(str), "%llu", nL1PListMismatch);
				break;
			case pi_l1plistskips:
				desc = "List prefetch skip";
				snprintf(str, sizeof(str), "%llu", nL1PListSkips);
				break;
			case pi_l1plistoverruns:
				desc = "List prefetch overrun";
				snprintf(str, sizeof(str), "%llu", nL1PListOverruns);
				break;
			case pi_l1plistlatestalls:
				desc = "Stalls list prefetch behind";
				snprintf(str, sizeof(str), "%.2f", nL1PLatePrefetchStalls / nCoreCycles);
				break;
			default:
				continue;
			}

			printf("%"SCELLWIDTH"s|", str);
		}

		printf(" %s\n", desc);
	}
#endif
}


static void exec_table(bool sloppiness, hm_func_double hm_double, hm_func_float hm_float, bgq_hmflags additional_opts) {
	benchstat excerpt;

	if (g_proc_id == 0)
		printf("%10s|", "");
	for (int i3 = 0; i3 < COUNTOF(flags); i3 += 1) {
		if (g_proc_id == 0)
			printf("%-"SCELLWIDTH"s|", flags_desc[i3]);
	}
	if (g_proc_id == 0)
		printf("\n");
	print_repeat("-", 10 + 1 + (CELLWIDTH + 1) * COUNTOF(flags));
	if (g_proc_id == 0)
		printf("\n");
	for (int i2 = 0; i2 < COUNTOF(omp_threads); i2 += 1) {
		int threads = omp_threads[i2];

		if (g_proc_id == 0)
			printf("%-10s|", omp_threads_desc[i2]);

		benchstat stats[COUNTOF(flags)];
		for (int i3 = 0; i3 < COUNTOF(flags); i3 += 1) {
			bgq_hmflags hmflags = flags[i3];
			hmflags = hmflags | additional_opts;

			benchstat result = runbench(1, 3, sloppiness, threads, hmflags, hm_double, hm_float);
			stats[i3] = result;

			if (threads == 64 && i3 == 2) {
				excerpt = result;
			}

			char str[80] = { 0 };
			snprintf(str, sizeof(str), "%.0f mflop/s%s", result.flops / MEGA, (result.error > 0.001) ? "X" : "");
			if (g_proc_id == 0)
				printf("%"SCELLWIDTH"s|", str);
			if (g_proc_id == 0)
				fflush(stdout);
		}
		if (g_proc_id == 0)
			printf("\n");

		if (g_proc_id == 0) {
			//print_stats(stats, BODY_SITES*4, SURFACE_SITES*4, HALO_SITES*4);
			print_stats(stats, 0,0,0);
		}

		print_repeat("-", 10 + 1 + (CELLWIDTH + 1) * COUNTOF(flags));
		if (g_proc_id == 0)
			printf("\n");
	}

	if (g_proc_id == 0) {
		printf("Hardware counter excerpt (64 threads, nocom):\n");
		mypapi_print_counters(&excerpt.counters);
	}
	if (g_proc_id == 0)
		printf("\n");
}


static void exec_bench() {
	print_repeat("\n", 2);

	for (int i0 = 0; i0 < COUNTOF(kamuls); i0 += 1) {
		//bool kamul = kamuls[i0];
		if (g_proc_id == 0) printf("Benchmark: %s\n", kamuls_desc[i0]);

		for (int i1 = 0; i1 < COUNTOF(sloppinessess); i1 += 1) {
			//bool sloppiness = sloppinessess[i1];

			if (g_proc_id == 0) printf("Benchmarking precision: %s\n", sloppinessess_desc[i1]);
			//exec_table(sloppiness, bgq_HoppingMatrix_double, bgq_HoppingMatrix_float, !kamul * hm_nokamul);
		}

		if (g_proc_id == 0) printf("\n");
	}

	if (g_proc_id == 0) printf("Benchmark done\n");
}


#define DCBT_DIFF 1

#if BGQ_PREFETCH_LIST
L1P_Pattern_t *(*fake_sequential_patterns)[2/*even/odd*/][64/*For each thread*/][6/*total threads*/] = NULL;
#endif



static void hm_fake_sequential(bool isOdd, bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bgq_hmflags opts) {
#if 0
	const bool noprefetchlist = opts & hm_noprefetchlist;
	const bool noprefetchstream = opts & hm_noprefetchstream;
	const bool noprefetchexplicit = opts & hm_noprefetchexplicit;

	bgq_vector4double_decl(q);
	bgq_cconst(q, 5.1, 0);

#pragma omp parallel firstprivate(bgq_vars(q))
	{
#if BGQ_PREFETCH_LIST
		if (!noprefetchlist) {
			if (!fake_sequential_patterns) {
				fake_sequential_patterns = malloc(sqr(g_num_total_spinorfields_double) * sizeof(*fake_sequential_patterns));
				memset(fake_sequential_patterns, 0, sqr(g_num_total_spinorfields_double) * sizeof(*fake_sequential_patterns));
			}

			int targetindex = bgq_spinorfield_find_index_double(targetfield);
			int spinorindex = bgq_spinorfield_find_index_double(spinorfield);
			int totThreads = omp_get_num_threads();
			L1P_Pattern_t **patternhandle = &fake_sequential_patterns[targetindex + spinorindex][isOdd][Kernel_ProcessorID()][ilog(totThreads)];
			bool l1p_first = false;
			if (!*patternhandle) {
				//const uint64_t prefsize = 10000/*how to choose this?*/;
				const uint64_t prefsize = 2/*safety*/ * 2*(2 + 2) * VOLUME_SITES / totThreads;
				int retval = L1P_GetPattern(patternhandle);
				if (retval == L1P_NOTCONFIGURED) {
					L1P_CHECK(L1P_PatternConfigure(prefsize));
					L1P_CHECK(L1P_GetPattern(patternhandle));
				} else {
					L1P_CHECK(retval);
					L1P_CHECK(L1P_AllocatePattern(prefsize, patternhandle));
				}
				l1p_first = true;
			}
			assert(*patternhandle);

			// Don't have acces to L1P_CFG_PF_USR (segfault)
			//L1P_CHECK(L1P_PatternSetEnable(1));

			L1P_CHECK(L1P_SetPattern(*patternhandle));


			L1P_CHECK(L1P_PatternStart(true));

			// Don't have acces to L1P_CFG_PF_USR (segfault)
			//int pattern_enable = 0;
			//L1P_CHECK(L1P_PatternGetEnable(&pattern_enable));
			//if (!pattern_enable)
			//	printf("error: List prefetcher not enabled\n");
		}
#endif

		if (!noprefetchstream) {
			bgq_prefetch_forward(&spinorfield[0]);
		}

#pragma omp for schedule(static)
		for (int i = 0; i < VOLUME_SITES; i += 1) {
			if (!noprefetchexplicit) {
				bgq_su3_spinor_prefetch_double(&spinorfield[i+DCBT_DIFF]);
			}
			if (!noprefetchstream) {
				bgq_prefetch_forward(&spinorfield[i]);
			}

			bgq_su3_spinor_decl(spinor);
			bgq_su3_spinor_load_double(spinor, &spinorfield[i]);

			// Do some transformation
			bgq_su3_spinor_decl(result);
			bgq_su3_cvmul(result_v0, q, spinor_v0);
			bgq_su3_cvmul(result_v1, q, spinor_v1);
			bgq_su3_cvmul(result_v2, q, spinor_v2);
			bgq_su3_cvmul(result_v3, q, spinor_v3);

			if (!noprefetchexplicit) {
				bgq_su3_spinor_zeroload_double(&targetfield[i]);
			}
			bgq_su3_spinor_store_double(&targetfield[i], result);
		}

#if BGQ_PREFETCH_LIST
		if (!noprefetchlist) {
			L1P_CHECK(L1P_PatternStop());
		}
#endif
	}
#endif
}

#if BGQ_PREFETCH_LIST
L1P_Pattern_t *(*fake_random_patterns)[2/*even/odd*/][64/*For each thread*/][6/*total threads*/] = NULL;
#endif

static void hm_fake_random(bool isOdd, bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bgq_hmflags opts) {
#if 0
	const bool noprefetchlist = opts & hm_noprefetchlist;
	const bool noprefetchstream = opts & hm_noprefetchstream;
	const bool noprefetchexplicit = opts & hm_noprefetchexplicit;

	bgq_vector4double_decl(q);
	bgq_cconst(q, 5.1, 0);

#pragma omp parallel firstprivate(bgq_vars(q))
	{
#if BGQ_PREFETCH_LIST
		if (!noprefetchlist) {
			if (!fake_random_patterns) {
				fake_random_patterns = malloc(sqr(g_num_total_spinorfields_double) * sizeof(*fake_random_patterns));
				memset(fake_random_patterns, 0, sqr(g_num_total_spinorfields_double) * sizeof(*fake_random_patterns));
			}

			int targetindex = bgq_spinorfield_find_index_double(targetfield);
			int spinorindex = bgq_spinorfield_find_index_double(spinorfield);
			int totThreads = omp_get_num_threads();
			L1P_Pattern_t **patternhandle = &fake_random_patterns[targetindex + spinorindex][isOdd][Kernel_ProcessorID()][ilog(totThreads)];
			bool l1p_first = false;
			if (!*patternhandle) {
				//const uint64_t prefsize = 10000/*how to choose this?*/;
				const uint64_t prefsize = 2 * /*safety*/2 * (2 + 2) * VOLUME_SITES / totThreads;
				int retval = L1P_GetPattern(patternhandle);
				if (retval == L1P_NOTCONFIGURED) {
					L1P_CHECK(L1P_PatternConfigure(prefsize));
					L1P_CHECK(L1P_GetPattern(patternhandle));
				} else {
					L1P_CHECK(retval);
					L1P_CHECK(L1P_AllocatePattern(prefsize, patternhandle));
				}
				l1p_first = true;
			}
			assert(*patternhandle);

			// Don't have acces to L1P_CFG_PF_USR (segfault)
			//L1P_CHECK(L1P_PatternSetEnable(1));

			L1P_CHECK(L1P_SetPattern(*patternhandle));

			// From Peter Boye's bfm
			 //uint64_t *pf_usr = (uint64_t *)L1P_CFG_PF_USR;
			 //*pf_usr |=  L1P_CFG_PF_USR_pf_stream_est_on_dcbt | L1P_CFG_PF_USR_pf_stream_optimistic | L1P_CFG_PF_USR_pf_stream_prefetch_enable | L1P_CFG_PF_USR_pf_stream_establish_enable;


			L1P_CHECK(L1P_PatternStart(l1p_first));

			// Don't have acces to L1P_CFG_PF_USR (segfault)
			//int pattern_enable = 0;
			//L1P_CHECK(L1P_PatternGetEnable(&pattern_enable));
			//if (!pattern_enable)
			//	printf("error: List prefetcher not enabled\n");
		}
#endif


#pragma omp for schedule(static)
		for (int j = 0; j < VOLUME_SITES; j += 1) {
			if (!noprefetchexplicit) {
				const long long lj = j + DCBT_DIFF;
				const int k = (7 * lj + lj * lj) % VOLUME_SITES;

				bgq_su3_spinor_prefetch_double(&spinorfield[k]);
			}

			// Some pseudo-random access
			const long long lj = j;
			const int i = (7 * lj + lj * lj) % VOLUME_SITES;
			if (!noprefetchstream) {
				bgq_prefetch_forward(&spinorfield[i]);
			}

			bgq_su3_spinor_decl(spinor);
			bgq_su3_spinor_load_double(spinor, &spinorfield[i]);

			// Do some transformation
			bgq_su3_spinor_decl(result);
			bgq_su3_cvmul(result_v0, q, spinor_v0);
			bgq_su3_cvmul(result_v1, q, spinor_v1);
			bgq_su3_cvmul(result_v2, q, spinor_v2);
			bgq_su3_cvmul(result_v3, q, spinor_v3);

			if (!noprefetchexplicit) {
				bgq_su3_spinor_zeroload_double(&targetfield[i]);
			}
			bgq_su3_spinor_store_double(&targetfield[i], result);
		}

#if BGQ_PREFETCH_LIST
		if (!noprefetchlist) {
			L1P_CHECK(L1P_PatternStop());
		}
#endif
	}
#endif
}


inline void exec_fakebench() {
	if (g_proc_id == 0) printf("Sequential access:\n");
	exec_table(false, hm_fake_sequential, NULL, 0);
	if (g_proc_id == 0) printf("\n");

	if (g_proc_id == 0) printf("Random access:\n");
	exec_table(false, hm_fake_random, NULL, 0);
	if (g_proc_id == 0) printf("\n");
}





void bgq_HoppingMatrix_asm(bool isOdd, bgq_spinorfield_double restrict targetfield, bgq_spinorfield_double restrict spinorfield, bgq_gaugefield_double restrict gaugefield, bgq_hmflags opts);
void bgq_HoppingMatrix_asm_parallel(bool isOdd, bgq_spinorfield_double restrict targetfield, bgq_spinorfield_double restrict spinorfield, bgq_gaugefield_double restrict gaugefield, bgq_hmflags opts, int tid, int tot_threads);


static void runbench_exec_sequential(int tid, int tot_threads) {
	double *baseaddr_read = 0;// = (double*)g_spinorfields_double[0];
	double *baseaddr_write = 0;// = (double*)g_spinorfields_double[1];

	size_t count = PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LY*PHYSICAL_LZ;
	size_t threadcount = (count+tot_threads-1)/tot_threads;
	assert(count==threadcount*tot_threads);
	for (size_t counter = 0; counter < threadcount; counter += 1) {
		const size_t txyz = tid*threadcount + counter;
		if (txyz >= count)
			break;

		const size_t i = txyz;
		//const size_t i_next = txyz+1;

		double *addr_read = baseaddr_read + i*2*4*3*2;
		double *addr_write = baseaddr_write + i*2*4*3*2;

		bgq_su3_spinor_decl(spinor);
		bgq_su3_spinor_load_double(spinor, addr_read);

		bgq_su3_spinor_store_double(addr_write, spinor);
	}
}

static void runbench_exec_random(int tid, int tot_threads) {
	double *baseaddr_read=0;// = (double*)g_spinorfields_double[0];
	double *baseaddr_write=0;// = (double*)g_spinorfields_double[1];

	size_t count = PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LY*PHYSICAL_LZ;
	size_t threadcount = (count+tot_threads-1)/tot_threads;
	assert(count==threadcount*tot_threads);
	for (size_t counter = 0; counter < threadcount; counter += 1) {
		const size_t txyz = tid*threadcount + counter;
		if (txyz >= count)
			break;

		const size_t i = 7*txyz + txyz*txyz % count; // something seemingly random
		//const size_t i_next = 7*(txyz+1) + (txyz+1)*(txyz+1) % count;

		double *addr_read = baseaddr_read + i*2*4*3*2;
		double *addr_write = baseaddr_write + i*2*4*3*2;

		bgq_su3_spinor_decl(spinor);
		bgq_su3_spinor_load_double(spinor, addr_read);

		bgq_su3_spinor_store_double(addr_write, spinor);
	}
}

static void runbench_exec_asm_parallel(int k_max, int j_max, bool sloppyprec, bgq_hmflags opts, bool noprefetchlist, bool l1pnonstoprecord, bool *l1plist_first, int tid, int tot_threads, bool isWarmup, int select) {
	static bgq_gaugefield_double mem;

	if (!mem) {
		mem = malloc_aligned(8*4*VOLUME/2 * sizeof(bgq_gaugesite_double), 128 /*L2 cache line size*/);
	}

	if (sloppyprec) {
		// No impl	yet
	} else {
		for (int j = 0; j < j_max; j += 1) {
#ifdef XLC
			if (!noprefetchlist) {
				L1P_PatternStart(l1pnonstoprecord || *l1plist_first);
			}
#endif
			for (int k = 0; k < k_max; k += 1) {
				switch (select) {
				case 0:
					//bgq_HoppingMatrix_asm_parallel(false, g_spinorfields_double[k + k_max], g_spinorfields_double[k], mem, opts, tid, tot_threads);
					break;
				case 1:
					runbench_exec_sequential(tid, tot_threads);
					break;
				case 2:
					runbench_exec_random(tid, tot_threads);
					break;
				}
			}
#ifdef XLC
			if (!noprefetchlist) {
				L1P_PatternStop();
				*l1plist_first = false;
			}
#endif
		}
	}
}


static benchstat runbench_asm(int k_max, int j_max, bool sloppyprec, int ompthreads, bgq_hmflags opts, int select) {
	//const bool nocom = opts & hm_nocom;
	//const bool nooverlap = opts & hm_nooverlap;
	const bool nokamul = opts & hm_nokamul;
	bool noprefetchlist = opts & hm_noprefetchlist;
	//const bool noprefetchstream = opts & hm_noprefetchstream;
	//const bool noprefetchexplicit = opts & hm_noprefetchexplicit;
	//const bool noweylsend = opts & hm_noweylsend;
	const bool nobody = opts & hm_nobody;
	//const bool nosurface = opts & hm_nosurface;
	//const bgq_hmflags implicitprefetch = opts & (hm_prefetchimplicitdisable | hm_prefetchimplicitoptimistic | hm_prefetchimplicitconfirmed);
	const bool l1pnonstoprecord = opts & hm_l1pnonstoprecord;
	//const bool experimental = opts & hm_experimental;

	omp_set_num_threads(ompthreads);

	// Setup options
#ifdef XLC



#pragma omp parallel
	{

	}

#pragma omp parallel
	{
		// Test whether it persists between parallel sections
		//L1P_StreamPolicy_t getpol = 0;
		//L1P_CHECK(L1P_GetStreamPolicy(&getpol));
		//if (getpol != pol)
		//	fprintf(stderr, "MK StreamPolicy not accepted\n");

		//L1P_CHECK(L1P_SetStreamDepth());
		//L1P_CHECK(L1P_SetStreamTotalDepth());

		//L1P_GetStreamDepth
		//L1P_GetStreamTotalDepth


		// Peter Boyle's setting
		// Note: L1P_CFG_PF_USR_pf_stream_establish_enable is never set in L1P_SetStreamPolicy
		//uint64_t *addr = ((uint64_t*)(Kernel_L1pBaseAddress() + L1P_CFG_PF_USR_ADJUST));
		//*addr |=  L1P_CFG_PF_USR_pf_stream_est_on_dcbt | L1P_CFG_PF_USR_pf_stream_optimistic | L1P_CFG_PF_USR_pf_stream_prefetch_enable | L1P_CFG_PF_USR_pf_stream_establish_enable; // Enable everything???
	}
#endif



	int iterations = 1/*warmup*/ + MYPAPI_SETS;

	mypapi_counters counters = {0};
	counters.init = false;

	double localsumtime = 0;
	double localsumsqtime = 0;

#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		int tot_threads = omp_get_num_threads();

#ifdef XLC
		if (experimental) {
			noprefetchlist = false;
			uint64_t *addr = ((uint64_t*)(Kernel_L1pBaseAddress() + L1P_CFG_PF_USR_ADJUST));
			*addr |=  L1P_CFG_PF_USR_pf_stream_est_on_dcbt | L1P_CFG_PF_USR_pf_stream_optimistic | L1P_CFG_PF_USR_pf_stream_prefetch_enable | L1P_CFG_PF_USR_pf_stream_establish_enable; // Enable everything???
		} else {
		L1P_StreamPolicy_t pol;
		switch (implicitprefetch) {
		case hm_prefetchimplicitdisable:
			pol = noprefetchstream ? L1P_stream_disable : L1P_confirmed_or_dcbt/*No option to selectively disable implicit stream*/;
			break;
		case hm_prefetchimplicitoptimistic:
			pol = L1P_stream_optimistic;
			break;
		case hm_prefetchimplicitconfirmed:
			pol = noprefetchstream ? L1P_stream_confirmed : L1P_confirmed_or_dcbt;
			break;
		default:
			// Default setting
			pol = L1P_confirmed_or_dcbt;
			break;
		}
		L1P_CHECK(L1P_SetStreamPolicy(pol));
		}

		if (!noprefetchlist) {
			L1P_CHECK(L1P_PatternConfigure((sizeof(bgq_spinorsite_double)*9*VOLUME/4 + sizeof(bgq_gaugesite_double)*8*VOLUME/4)/64));
		}
#endif

		double threadsumtime = 0;
		double threadsumsqtime = 0;
		bool l1p_first = true;

		for (int j = 0; j < iterations; j += 1) {
			bool isWarmup = (j == 0);
			//bool isLast = (j == iterations-1);
			bool isPapi = !isWarmup && (j >= iterations - MYPAPI_SETS);

			// Barrier outside timed loop
			#pragma omp barrier

			//if (!noprefetchlist) {
			//	L1P_PatternStart(isWarmup);
			//}

			if (isPapi) {
				mypapi_start(j - (iterations - MYPAPI_SETS));
			}
			double start_time = MPI_Wtime();

			runbench_exec_asm_parallel(k_max, j_max, sloppyprec, opts, noprefetchlist, l1pnonstoprecord, &l1p_first, tid, tot_threads, isWarmup, select);

			double end_time = MPI_Wtime();
			if (isPapi) {
				mypapi_counters curcounters = mypapi_stop();
				if (tid == 0) {
					counters = mypapi_merge_counters(&counters, &curcounters);
				}
			}

			//if (!noprefetchlist) {
			//	L1P_PatternStop();
			//}


			double runtime = end_time - start_time;
			if (!isWarmup) {
				threadsumtime += runtime;
				threadsumsqtime += sqr(runtime);
			}
		}

		if (tid == 0) {
			localsumtime = threadsumtime / (iterations-1);
			localsumsqtime = threadsumsqtime / (iterations-1);
		}

#ifdef XLC
		if (!noprefetchlist) {
			L1P_CHECK(L1P_PatternUnconfigure());
		}
#endif
	}

	double localavgtime = localsumtime / j_max;
	double localavgsqtime = sqr(localavgtime);
	double localrmstime = sqrt((localsumsqtime / j_max) - localavgsqtime);

	double localtime[] = { localavgtime, localavgsqtime, localrmstime };
	double sumreduce[3] = { -1, -1, -1 };
	MPI_Allreduce(&localtime, &sumreduce, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	double sumtime = sumreduce[0];
	double sumsqtime = sumreduce[1];
	double sumrmstime = sumreduce[2];

	double avgtime = sumtime / g_nproc;
	double avglocalrms = sumrmstime / g_nproc;
	double rmstime = sqrt((sumsqtime / g_nproc) - sqr(avgtime));

	int lups = VOLUME / 2;
	int flops_per_lup_body = 0;
	if (!nobody) {
		flops_per_lup_body += 1320;
		if (!nokamul)
			flops_per_lup_body += 8 * 2 * 3 * 6;
	}
	double flops = ((double)flops_per_lup_body * (double)lups);

	benchstat result;
	result.avgtime = avgtime;
	result.localrmstime = avglocalrms;
	result.globalrmstime = rmstime;
	result.lups = lups;
	result.flops = flops / avgtime;
	result.error = 0;
	result.counters = counters;
	result.opts = opts;
	return result;
}


static void exec_table_asm(bool sloppiness, bgq_hmflags additional_opts, int select) {
	switch (select) {
	case 0:
	master_print("Benchmarking plain:\n");
	break;
	case 1:
		master_print("Benchmarking sequential:\n");
		break;
	case 2:
		master_print("Benchmarking ransdom:\n");
		break;
	}
	benchstat excerpt;

	if (g_proc_id == 0)
		printf("%10s|", "");
	for (int i3 = 0; i3 < COUNTOF(flags); i3 += 1) {
		if (g_proc_id == 0)
			printf("%-"SCELLWIDTH"s|", flags_desc[i3]);
	}
	if (g_proc_id == 0)
		printf("\n");
	print_repeat("-", 10 + 1 + (CELLWIDTH + 1) * COUNTOF(flags));
	if (g_proc_id == 0)
		printf("\n");
	for (int i2 = 0; i2 < COUNTOF(omp_threads); i2 += 1) {
		int threads = omp_threads[i2];

		if (g_proc_id == 0)
			printf("%-10s|", omp_threads_desc[i2]);

		benchstat stats[COUNTOF(flags)];
		for (int i3 = 0; i3 < COUNTOF(flags); i3 += 1) {
			bgq_hmflags hmflags = flags[i3];
			hmflags = hmflags | additional_opts;

			benchstat result = runbench_asm(1, 3, sloppiness, threads, hmflags, select);
			stats[i3] = result;

			if (threads == 64 && i3 == 2) {
				excerpt = result;
			}

			char str[80] = { 0 };
			snprintf(str, sizeof(str), "%.0f mflop/s%s", result.flops / MEGA, (result.error > 0.001) ? "X" : "");
			if (g_proc_id == 0)
				printf("%"SCELLWIDTH"s|", str);
			if (g_proc_id == 0)
				fflush(stdout);
		}
		if (g_proc_id == 0)
			printf("\n");

		if (g_proc_id == 0) {
			print_stats(stats, VOLUME/2, 0, 0);
		}

		print_repeat("-", 10 + 1 + (CELLWIDTH + 1) * COUNTOF(flags));
		if (g_proc_id == 0)
			printf("\n");
	}

	if (g_proc_id == 0) {
		printf("Hardware counter excerpt (64 threads, nocom):\n");
		mypapi_print_counters(&excerpt.counters);
	}
	if (g_proc_id == 0)
		printf("\n");
}
