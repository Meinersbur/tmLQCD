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
/*******************************************************************************
 *
 * Benchmark program for the even-odd preconditioned Wilson-Dirac operator
 *
 *
 *******************************************************************************/

#define MAIN_PROGRAM
//#define DD1_L1P_Workaround 1
#if HAVE_CONFIG_H
#include "config.h"
#endif
//#undef BGL

#include <stdlib.h>
#include <stdio.h>

#include <unistd.h>   /* GG */
#include <fcntl.h> /* GG */

/* GG */
extern FILE *fmemopen(void *__s, size_t __len, __const char *__modes) __THROW;

// MK
#include "mypapi.h"

#include <math.h>
#include <time.h>
#include <string.h>
#include "su3.h"
#include "su3adj.h"

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
#include "update_backward_gauge.h"
#include "test/check_geometry.h"
#include "xchange_halffield.h"
#include "D_psi.h"
#include "phmc.h"
#include "mpi_init.h"

#ifdef BGQ
#include "bgq/bgq_field.h"
#include "bgq/bgq_HoppingMatrix.h"
#include "bgq/bgq.h"
#endif
#include <omp.h>


//#if BGQ_PREFETCH_LIST
#include <l1p/pprefetch.h>
//#endif
#ifdef XLC
#include <l1p/sprefetch.h>
#endif

#define COUNTOF(arr) (sizeof(arr) / sizeof((arr)[0]))

/* GG */

void usage()
{
	fprintf(stdout, "Benchmark for Wilson twisted mass QCD\n");
	/*   fprintf(stdout, "Version %s \n\n", PACKAGE_VERSION); */
	/*   fprintf(stdout, "Please send bug reports to %s\n", PACKAGE_BUGREPORT); */
	fprintf(stdout, "Usage:   benchmark [options]\n");
	fprintf(stdout, "Options: [-f input-filename]\n");
	fprintf(stdout, "         [-v] verbose output (default no)\n");
	fprintf(stdout, "         [-h|-? this help]\n");
	exit(0);
}

/* #define BACKTRACE yes */
#if defined BACKTRACE
void
print_trace (void)
{
	/* GG */
	if ( g_proc_id )
	return;

	void *array[20];
	size_t size;
	char **strings;
	size_t i;

	size = backtrace (array, 20);
	strings = backtrace_symbols (array, size);

	//printf ("Obtained %zd stack frames.\n", size);

	for (i = 0; i < size; i++)
	if ( ! strstr(strings[i], "../invert [") && ! strstr(strings[i], "../invert(print_t") && ! strstr(strings[i], "/lib64/libc") )
	printf ("%s==", strings[i]);
	/*           printf ("%s\n", strings[i]); */

	free (strings);
}
#else
void print_trace(void) {
}
#endif

/* GG */
double gatime, getime;

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



static void exec_bench();
static void exec_fakebench();
static void check_correctness_double(bool nocom);
static void check_correctness_float();


complexdouble nooptaway;
L1P_Pattern_t *(g_patternhandle[64]);

static void listprefetch_test(bool l1plist) {
	bgq_spinorfield_double spinorfield = g_spinorfields_double[0];

	//mypapi_start();

//#pragma omp parallel
	{
#if BGQ_PREFETCH_LIST
		bool l1p_first = false;
		if (l1plist) {
			int totThreads = omp_get_num_threads();
			const int bitsp[] = {};
			L1P_Pattern_t **patternhandle = &g_patternhandle[Kernel_ProcessorID()];
			if (!*patternhandle) {
				const uint64_t prefsize = 1000000/*how to choose this?*/;
				int retval = L1P_GetPattern(patternhandle);
				if (retval == L1P_NOTCONFIGURED) {
					L1P_CHECK(L1P_PatternConfigure(prefsize));
					L1P_CHECK(L1P_GetPattern(patternhandle));
				} else {
					L1P_CHECK(retval);
					L1P_CHECK(L1P_AllocatePattern(prefsize,patternhandle));
				}
				l1p_first = true;
			}
			assert(*patternhandle);

			L1P_CHECK(L1P_SetPattern(*patternhandle));
			L1P_CHECK(L1P_PatternStart(!l1p_first));
	}
#endif

	complexdouble sum = 0;
	for (int i = 0; i < VOLUME; i+=1) {
		// Some pseudo-random access
		long long li = i;
		int j = (7*li + li*li) % VOLUME;

		bgq_spinorsite_double *site = &spinorfield[j];
		complexdouble val = site->s[i % 4][i % 3][i % 2];
		sum += val;
	}
	nooptaway = sum;

#if BGQ_PREFETCH_LIST
	if (l1plist) {
	uint64_t fetch_depth;
	uint64_t generate_depth;
	L1P_PatternGetCurrentDepth(&fetch_depth, &generate_depth);
	L1P_PatternStop();

	L1P_Status_t st;
	L1P_PatternStatus(&st);
	master_print("L1P_LIST: l1pfirst=%d maxed=%d abandoned=%d finished=%d endoflist=%d fetch_depth=%d generate_depth=%d\n", l1p_first, st.s.maximum, st.s.abandoned, st.s.finished, st.s.endoflist, (int)fetch_depth, (int)generate_depth);
	}
#endif

	//mypapi_stop();
	}
}


int main(int argc, char *argv[])
{
	int j;
	int k_max = 1;
#ifdef HAVE_LIBLEMON
	paramsXlfInfo *xlfInfo;
#endif

	/* GG */
	FILE* yyingg = NULL;
	void* yybufgg = NULL;
	int yyfd;
	char *input_filename = NULL;
	int c;
#ifdef MPI


	DUM_DERI = 6;
	DUM_SOLVER = DUM_DERI + 2;
	DUM_MATRIX = DUM_SOLVER + 6;
	NO_OF_SPINORFIELDS = DUM_MATRIX + 2;

	verbose = 0;

	//MPI_Init(&argc, &argv);
	int provided_threadlevel;
	MPI_CHECK(MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided_threadlevel));

	/* GG */
	MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id));
	master_print("provided threadlevel=%d\n", provided_threadlevel);
#endif
	g_rgi_C1 = 1.;

	/* GG modified ... */
	while ((c = getopt(argc, argv, "vh?f:")) != -1) {
		switch (c) {
		case 'f':
			input_filename = calloc(200, sizeof(char));
			strcpy(input_filename, optarg);
			break;
		case 'v':
			verbose = 1;
			break;

		case 'h':
			case '?':
			default:
			usage();
			break;
		}
	}

	/* Read the input file */
	/*   read_input("benchmark.input"); */
	if (input_filename == NULL ) {
		fprintf(stderr, "MK You forgot to specify a benchamrk.input\n");
		strcpy(input_filename, "benchmark.input");
	}

	/* GG */
	/*   strcpy(input_filename, "benchmark.input"); */
#define MPIO yes
#ifndef MPIO
	read_input(input_filename);
#else
	/* New read style with MPI_Bcast */
	yybufgg = (void *) malloc(8192*sizeof(char));
	yyingg = (FILE*) malloc(sizeof(FILE*));
	if (g_proc_id == 0) {
		yyfd = open(input_filename, O_RDONLY);
		if (!yyfd) {
			fprintf(stderr, "MK Cannot open file %s\n", input_filename);
			exit(2);
		}
		read(yyfd, yybufgg, 8192);
		close(yyfd);
	}
	MPI_Bcast(yybufgg, 8192, MPI_CHAR, 0, MPI_COMM_WORLD );
	yyingg = fmemopen(yybufgg, strlen(yybufgg), "r");
	read_input_fh(yyingg);
#endif

	/* GG */
	if (g_proc_id)
		verbose = 0;

	tmlqcd_mpi_init(argc, argv);
	mypapi_init();

	if (g_proc_id == 0) {
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
#ifdef PREFETCH
		printf("# the code was compiled with -DPREFETCH\n");
#else
		printf("# the code was compiled without -DPREFETCH\n");
#endif
		//print_prefetch;
#ifdef OMP
		printf("# the code was compiled with -DOMP\n");
#else
		printf("# the code was compiled without -DOMP\n");
#endif
#ifdef BGP
		printf("# The code was compiled for Blue Gene/P\n");
#endif
#ifdef _USE_HALFSPINOR
		printf("# The code was compiled with -D_USE_HALFSPINOR\n");
#endif
#ifdef _USE_SHMEM
		printf("# the code was compiled with -D_USE_SHMEM\n");
#  ifdef _PERSISTENT
		printf("# the code was compiled for persistent MPI calls (halfspinor only)\n");
#  endif
#endif
#ifdef MPI
#  ifdef _NON_BLOCKING
		printf("# the code was compiled for non-blocking MPI calls (spinor and gauge)\n");
#  endif
#endif
		//printHopVerMsg();
		printf("\n");
		fflush(stdout);
	}

#ifdef _GAUGE_COPY
	init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 1);
#else
	init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 0);
#endif
	init_geometry_indices(VOLUMEPLUSRAND + g_dbw2rand);

	if (even_odd_flag) {
		j = init_spinor_field(VOLUMEPLUSRAND / 2, 2 * k_max + 1);
	}
	else {
		j = init_spinor_field(VOLUMEPLUSRAND, 2 * k_max);
	}

	if (j != 0) {
		fprintf(stderr, "Not enough memory for spinor fields! Aborting...\n");
		exit(0);
	}
	j = init_moment_field(VOLUME, VOLUMEPLUSRAND);
	if (j != 0) {
		fprintf(stderr, "Not enough memory for moment fields! Aborting...\n");
		exit(0);
	}

	if (g_proc_id == 0) {
		fprintf(stdout, "The number of processes is %d \n", g_nproc);
		printf("# The lattice size is %d x %d x %d x %d\n",
		        (int) (T * g_nproc_t), (int) (LX * g_nproc_x), (int) (LY * g_nproc_y), (int) (g_nproc_z * LZ));
		printf("# The local lattice size is %d x %d x %d x %d\n",
		        (int) (T), (int) (LX), (int) (LY), (int) LZ);
		if (even_odd_flag) {
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
	if (j != 0) {
		fprintf(stderr, "Not enough memory for halfspinor fields! Aborting...\n");
		exit(0);
	}
	if (g_sloppy_precision_flag == 1) {
		g_sloppy_precision = 1;
		j = init_dirac_halfspinor32();
		if (j != 0) {
			fprintf(stderr, "Not enough memory for 32-Bit halfspinor fields! Aborting...\n");
			exit(0);
		}
	}
#  if (defined _PERSISTENT)
	init_xchange_halffield();
#  endif
#endif

	//check_geometry();
#if (defined MPI && !(defined _USE_SHMEM))
	//check_xchange();
#endif

	start_ranlux(1, 123456);
	random_gauge_field(reproduce_randomnumber_flag);

#ifdef MPI
	/*For parallelization: exchange the gaugefield */
	xchange_gauge();
#endif

#ifdef BGQ
	bgq_init_gaugefield_allprec();
	bgq_init_spinorfields_allprec(2 * k_max + 1, 0);
	bgq_hm_init_allprec();

	bgq_update_backward_gauge();
#endif

#ifdef _GAUGE_COPY
	update_backward_gauge();
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
		random_spinor_field(g_spinor_field[k], VOLUME / 2, 0);
		bgq_transfer_spinorfield_double(true, g_spinorfields_double[k], g_spinor_field[k]);

		double compare_transfer = bgq_spinorfield_compare_double(true, g_spinorfields_double[k], g_spinor_field[k], false);
		assert(compare_transfer == 0 /* fields should be bit-identical*/);
	}

//#pragma omp parallel 
{
	for (int i = 0; i < 2; i+=1) {
		//listprefetch_test(false);
	}
	for (int i = 0; i < 3; i+=1) {
		//listprefetch_test(true);
	}
}
	exec_fakebench();

	//check_correctness_double(true);
	//check_correctness_double(false);
	//check_correctness_float();
	assert(even_odd_flag);
	exec_bench();


	/* GG */
#if 0
#ifdef HAVE_LIBLEMON
	if(g_proc_id==0) {
		printf("doing parallel IO test ...\n");
	}
	xlfInfo = construct_paramsXlfInfo(0.5, 0);
	write_lemon_gauge_field_parallel( "conf.test", 64, xlfInfo);
	free(xlfInfo);
	if(g_proc_id==0) {
		printf("done ...\n");
	}
#endif
#endif

#ifdef MPI
	MPI_Finalize();
#endif
	free_gauge_field();
	free_geometry_indices();
	free_spinor_field();
	free_moment_field();
	return (0);
}


static inline double sqr(double val) {
	return val*val;
}


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
		Hopping_Matrix_nocom(ieo,l, k);
	} else {
		Hopping_Matrix(ieo, l, k);
	}
}


static void check_correctness_double(bool nocom) {
	master_print("MK Checking double precision%s correctness...\n", nocom ? " (no communication)" : "");

//#pragma omp parallel
	{
		int k = 0;
		int k_max = 1;
		bgq_hmflags hmflags = hm_nocom * nocom | hm_nooverlap;
		double compare_even_before;
		double compare_even;

//#pragma omp master
		{
			bgq_transfer_spinorfield_double(true, g_spinorfields_double[k], g_spinor_field[k]);
			compare_even_before = bgq_spinorfield_compare_double(true, g_spinorfields_double[k], g_spinor_field[k], false);
			assert(compare_even_before == 0 /* should be bitwise identical */);

			//bgq_initbgqref();
		}
//#pragma omp barrier
		bgq_HoppingMatrix_double(false, g_spinorfields_double[k + k_max], g_spinorfields_double[k], g_gaugefield_double, hmflags);
//#pragma omp master
		{
			//master_print("MK HM_orig start\n");
			Hopping_Matrix_switch(0, g_spinor_field[k + k_max], g_spinor_field[k], nocom);
			//master_print("MK HM_orig end\n");
			//bgq_savebgqref();
			//__asm__("int3");
			compare_even = bgq_spinorfield_compare_double(false, g_spinorfields_double[k + k_max], g_spinor_field[k + k_max], false);
			assert(compare_even < 0.001);
		}
//#pragma omp barrier
		bgq_HoppingMatrix_double(true, g_spinorfields_double[2 * k_max], g_spinorfields_double[k + k_max], g_gaugefield_double, hmflags);
//#pragma omp master
		{
			Hopping_Matrix_switch(1, g_spinor_field[k + 2 * k_max], g_spinor_field[k + k_max], nocom);
			double compare_odd = bgq_spinorfield_compare_double(true, g_spinorfields_double[k + 2 * k_max], g_spinor_field[k + 2 * k_max], false);
			assert(compare_odd < 0.001);

			master_print("Numerical instability between double precision implementations: even %e, odd %e\n", compare_even, compare_odd);
		}
	}
}

static void check_correctness_float() {
	master_print("MK Checking float precision correctness...\n");

//#pragma omp parallel
	{
		int k = 0;
		int k_max = 1;
		bgq_hmflags hmflags = 0;
		double compare_even_before;
		double compare_even;

//#pragma omp master
		{
			bgq_transfer_spinorfield_float(true, g_spinorfields_float[k], g_spinor_field[k]);
			compare_even_before = bgq_spinorfield_compare_float(true, g_spinorfields_float[k], g_spinor_field[k], false);
			assert(compare_even_before == 0 /* should be bitwise identical */);

			//bgq_initbgqref();
		}
//#pragma omp barrier
		bgq_HoppingMatrix_float(false, g_spinorfields_float[k + k_max], g_spinorfields_float[k], g_gaugefield_float, hmflags);
//#pragma omp master
		{
			Hopping_Matrix(0, g_spinor_field[k + k_max], g_spinor_field[k]);
			//bgq_savebgqref();
			compare_even = bgq_spinorfield_compare_float(false, g_spinorfields_float[k + k_max], g_spinor_field[k + k_max], false);
			assert(compare_even < 0.001);
		}
//#pragma omp barrier
		bgq_HoppingMatrix_float(true, g_spinorfields_float[2 * k_max], g_spinorfields_float[k + k_max], g_gaugefield_float, hmflags);
//#pragma omp master
		{
			Hopping_Matrix(1, g_spinor_field[k + 2 * k_max], g_spinor_field[k + k_max]);
			double compare_odd = bgq_spinorfield_compare_float(true, g_spinorfields_float[k + 2 * k_max], g_spinor_field[k + 2 * k_max], false);
			assert(compare_odd < 0.001);

			master_print("Numerical instability between float precision implementations: even %g, odd %g\n", compare_even, compare_odd);
		}
	}
}


static double runcheck(bool sloppyprec, bgq_hmflags hmflags, int k_max) {
	const int k = 0;
	hmflags = hmflags & ~hm_nokamul;
	double result;
	double compare_even_before;
	double compare_even;
	double compare_odd;

	if (!sloppyprec) {
//#pragma omp master
		{
		memset(g_spinorfields_double[k + k_max], 0, sizeof(bgq_spinorsite_double) * VOLUME_SITES);
		memset(g_spinorfields_double[k + 2*k_max], 0, sizeof(bgq_spinorsite_double) * VOLUME_SITES);

		bgq_transfer_spinorfield_double(true, g_spinorfields_double[k], g_spinor_field[k]);
		compare_even_before = bgq_spinorfield_compare_double(true, g_spinorfields_double[k], g_spinor_field[k], false);
		}
		bgq_HoppingMatrix_double(false, g_spinorfields_double[k + k_max], g_spinorfields_double[k], g_gaugefield_double, hmflags);
//#pragma omp master
		{
		Hopping_Matrix_switch(0, g_spinor_field[k + k_max], g_spinor_field[k], hmflags & hm_nocom);
		compare_even = bgq_spinorfield_compare_double(false, g_spinorfields_double[k + k_max], g_spinor_field[k + k_max], true);
		}
		bgq_HoppingMatrix_double(true, g_spinorfields_double[2 * k_max], g_spinorfields_double[k + k_max], g_gaugefield_double, hmflags);
//#pragma omp master
		{
		Hopping_Matrix_switch(1, g_spinor_field[k + 2*k_max], g_spinor_field[k + k_max], hmflags & hm_nocom);
		compare_odd = bgq_spinorfield_compare_double(true, g_spinorfields_double[k + 2*k_max], g_spinor_field[k + 2*k_max], true);

		result = max(max(compare_even_before, compare_even), compare_odd);
		}
	} else {
//#pragma omp master
		{
		memset(g_spinorfields_float[k + k_max], 0, sizeof(bgq_spinorsite_float) * VOLUME_SITES);
		memset(g_spinorfields_float[k + 2*k_max], 0, sizeof(bgq_spinorsite_float) * VOLUME_SITES);

		bgq_transfer_spinorfield_float(true, g_spinorfields_float[k], g_spinor_field[k]);
		compare_even_before = bgq_spinorfield_compare_float(true, g_spinorfields_float[k], g_spinor_field[k], true);
		}
		bgq_HoppingMatrix_float(false, g_spinorfields_float[k + k_max], g_spinorfields_float[k], g_gaugefield_float, hmflags);
//#pragma omp master
		{
		Hopping_Matrix_switch(0, g_spinor_field[k + k_max], g_spinor_field[k], hmflags & hm_nocom);
		compare_even = bgq_spinorfield_compare_float(false, g_spinorfields_float[k + k_max], g_spinor_field[k + k_max], true);
		}
		bgq_HoppingMatrix_float(true, g_spinorfields_float[2 * k_max], g_spinorfields_float[k + k_max], g_gaugefield_float, hmflags);
//#pragma omp master
		{
		Hopping_Matrix_switch(1, g_spinor_field[k + 2*k_max], g_spinor_field[k + k_max], hmflags & hm_nocom);
		compare_odd = bgq_spinorfield_compare_float(true, g_spinorfields_float[k + 2*k_max], g_spinor_field[k + 2*k_max], true);

		result = max(max(compare_even_before, compare_even), compare_odd);
		}
	}

	return result;
}


static benchstat runbench(int k_max, int j_max, bool sloppyprec, int ompthreads, bgq_hmflags opts, hm_func_double hm_double, hm_func_float hm_float) {
	const bool nocom = opts & hm_nocom;
	const bool nooverlap = opts & hm_nooverlap;
	const bool nokamul = opts & hm_nokamul;
	const bool noprefetchlist = opts & hm_noprefetchlist;
	const bool noprefetchstream = opts & hm_noprefetchstream;
	const bool noprefetchexplicit = opts & hm_noprefetchexplicit;
	const bool noweylsend = opts & hm_noweylsend;
	const bool nobody = opts & hm_nobody;
	const bool nosurface = opts & hm_nosurface;
	bgq_hmflags implicitprefetch = opts & (hm_prefetchimplicitdisable | hm_prefetchimplicitoptimistic | hm_prefetchimplicitconfirmed);

	// Setup options
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

	// Error checking
	double errtmp = runcheck(sloppyprec, opts, k_max);


	double localsumtime = 0;
	double localsumsqtime = 0;
	double err = 0;
	mypapi_counters counters;
	counters.init = false;
	int iterations = 1 + j_max;
	if (iterations < 1 + MYPAPI_SETS)
		iterations = 1 + MYPAPI_SETS;

	for (int j = 0; j < iterations; j += 1) {
		//master_print("Starting iteration %d of %d\n", j+1, iterations);
		bool isWarmup = (j == 0);
		bool isLast = (j == iterations-1);
		bool isPapi = !isWarmup && (j >= iterations - MYPAPI_SETS);
		bool isJMax = !isWarmup && (j >= iterations - j_max);

		if (isPapi) {
			mypapi_start(j - (iterations - MYPAPI_SETS));
		}
		double start_time = MPI_Wtime();
		if (sloppyprec) {
			for (int k = 0; k < 1; k += 1) {
				hm_float(false, g_spinorfields_float[k + k_max], g_spinorfields_float[k], g_gaugefield_float, opts);
				hm_float(true, g_spinorfields_float[2 * k_max], g_spinorfields_float[k + k_max], g_gaugefield_float, opts);
			}
		} else {
			for (int k = 0; k < 1; k += 1) {
				hm_double(false, g_spinorfields_double[k + k_max], g_spinorfields_double[k], g_gaugefield_double, opts);
				hm_double(true, g_spinorfields_double[2 * k_max], g_spinorfields_double[k + k_max], g_gaugefield_double, opts);
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
	int lups_body = BODY_ZLINES * PHYSICAL_LP*PHYSICAL_LK * LOCAL_LZ;
	int lups_surface = SURFACE_ZLINES * PHYSICAL_LP*PHYSICAL_LK * LOCAL_LZ;
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

static int omp_threads[] = { 1, 2, 4, 8, 16, 32, 64 };
static char *omp_threads_desc[] = { "1", "2", "4", "8", "16", "32", "64" };

#if 0
static struct {
	bool nocom;
	bool comnooverlap;
	bool noweylsend;
	bool nobody;
	bool nosurface;
//} coms[] = { { false, false, false, false, false }, { false, true, false, false, false }, { true, -1, false, false, false }, { true, -1, true, false, true }, { true, -1, false, true, false } };
} coms[] = { { true, -1, false, false, false } };
//static char* com_desc[] = { "Com async", "Com sync", "Com off", "bodyonly", "surfaceonly" };
static char* com_desc[] = { "Com off" };
#endif

static bgq_hmflags flags[] = {
		0,
		hm_nooverlap,
		hm_nocom,
		hm_nocom | hm_nosurface | hm_noweylsend,
		hm_nocom | hm_nobody,
		hm_nocom | hm_noweylsend | hm_noprefetchlist | hm_noprefetchstream | hm_noprefetchexplicit | hm_prefetchimplicitdisable,
		hm_nocom | hm_noweylsend | hm_noprefetchlist | hm_noprefetchstream | hm_noprefetchexplicit | hm_prefetchimplicitconfirmed,
		hm_nocom | hm_noweylsend | hm_noprefetchlist | hm_noprefetchstream | hm_noprefetchexplicit | hm_prefetchimplicitoptimistic,
		hm_nocom | hm_noweylsend | hm_noprefetchlist                       | hm_noprefetchexplicit | hm_prefetchimplicitdisable,
		hm_nocom | hm_noweylsend | hm_noprefetchlist | hm_noprefetchstream                         | hm_prefetchimplicitdisable,
		hm_nocom | hm_noweylsend                     | hm_noprefetchstream | hm_noprefetchexplicit | hm_prefetchimplicitdisable };
static char* flags_desc[] = {
		"Com async",
		"Com sync",
		"Com off",
		"bodyonly",
		"surfaceonly",
		"pf disable",
		"pf confirmed",
		"pf optimistic",
		"pf stream",
		"pf explicit",
		"pf list"
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


static void print_stats(benchstat stats[COUNTOF(flags)]) {
	int threads = omp_get_num_threads();

	for (mypapi_interpretations j = 0; j < __pi_COUNT; j+=1) {
		printf("%10s|", "");
		char *desc = NULL;

		for (int i3 = 0; i3 < COUNTOF(flags); i3 += 1) {
			char str[80];
			str[0] = '\0';
			benchstat *stat = &stats[i3];
			bgq_hmflags opts = stat->opts;

			uint64_t bodySites = BODY_SITES*4;
			uint64_t surfaceSites = SURFACE_SITES*4;
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
				nNecessaryInstr += HALO_SITES * (2*3*2/*QMUL+QADD*/ + 4*3/*LD*/ + 2*3/*ST*/)/2;
			if (!(opts & hm_nosurface))
				nNecessaryInstr += surfaceSites * (240/*QFMA*/ + 180/*QMUL+QADD*/ + 180/*LD+ST*/ - 2*3*2/*QMUL+QADD*/ - 4*3/*LD*/ + 2*3/*LD*/)/2;
			if (!(opts & hm_nokamul))
				nNecessaryInstr += sites * 8*2*3*6;

			uint64_t nL1PListStarted = stats[i3].counters.native[PEVT_L1P_LIST_STARTED];
			uint64_t nL1PListAbandoned= stats[i3].counters.native[PEVT_L1P_LIST_ABANDON];
			uint64_t nL1PListMismatch= stats[i3].counters.native[PEVT_L1P_LIST_MISMATCH];
			uint64_t nL1PListSkips = stats[i3].counters.native[PEVT_L1P_LIST_SKIP];
			uint64_t nL1PListOverruns = stats[i3].counters.native[PEVT_L1P_LIST_CMP_OVRUN_PREFCH];

			double nL1PLatePrefetchStalls = stats[i3].counters.native[PEVT_L1P_BAS_LU_STALL_LIST_RD_CYC];


			switch (j) {
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
			print_stats(stats);
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
		bool kamul = kamuls[i0];
		if (g_proc_id == 0) printf("Benchmark: %s\n", kamuls_desc[i0]);

		for (int i1 = 0; i1 < COUNTOF(sloppinessess); i1 += 1) {
			bool sloppiness = sloppinessess[i1];

			if (g_proc_id == 0) printf("Benchmarking precision: %s\n", sloppinessess_desc[i1]);
			exec_table(sloppiness, bgq_HoppingMatrix_double, bgq_HoppingMatrix_float, !kamul * hm_nokamul);
		}

		if (g_proc_id == 0) printf("\n");
	}

	if (g_proc_id == 0) printf("Benchmark done\n");
}


#define DCBT_DIFF 1

L1P_Pattern_t *(*fake_sequential_patterns)[2/*even/odd*/][64/*For each thread*/][6/*total threads*/] = NULL;

static void hm_fake_sequential(bool isOdd, bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bgq_hmflags opts) {
	const bool noprefetchlist = opts & hm_noprefetchlist;
	const bool noprefetchstream = opts & hm_noprefetchstream;
	const bool noprefetchexplicit = opts & hm_noprefetchexplicit;

	bgq_vector4double_decl(q);
	bgq_cconst(q, 5.1, 0);

#pragma omp parallel firstprivate(bgq_vars(q))
	{
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

		if (!noprefetchlist) {
			L1P_CHECK(L1P_PatternStop());
		}
	}
}

L1P_Pattern_t *(*fake_random_patterns)[2/*even/odd*/][64/*For each thread*/][6/*total threads*/] = NULL;

static void hm_fake_random(bool isOdd, bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bgq_hmflags opts) {
	const bool noprefetchlist = opts & hm_noprefetchlist;
	const bool noprefetchstream = opts & hm_noprefetchstream;
	const bool noprefetchexplicit = opts & hm_noprefetchexplicit;

	bgq_vector4double_decl(q);
	bgq_cconst(q, 5.1, 0);

#pragma omp parallel firstprivate(bgq_vars(q))
	{

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


#pragma omp for schedule(static)
		for (int j = 0; j < VOLUME_SITES; j += 1) {
			if (!noprefetchexplicit) {
				const long long lj = j + DCBT_DIFF;
				const int i = (7 * lj + lj * lj) % VOLUME_SITES;

				bgq_su3_spinor_prefetch_double(&spinorfield[i]);
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

		if (!noprefetchlist) {
			L1P_CHECK(L1P_PatternStop());
		}
	}
}


static void exec_fakebench() {
	if (g_proc_id == 0) printf("Sequential access:\n");
	exec_table(false, hm_fake_sequential, NULL, 0);
	if (g_proc_id == 0) printf("\n");

	if (g_proc_id == 0) printf("Random access:\n");
	exec_table(false, hm_fake_random, NULL, 0);
	if (g_proc_id == 0) printf("\n");
}

