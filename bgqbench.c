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
#endif
#include <omp.h>

//#if BGQ_PREFETCH_LIST
#include <l1p/pprefetch.h>
//#endif

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
static void check_correctness_double(bool nocom);
static void check_correctness_float();


complexdouble nooptaway;
L1P_Pattern_t *(g_patternhandle[64]);

static void listprefetch_test(bool l1plist) {
	bgq_spinorfield_double spinorfield = g_spinorfields_double[0];

	mypapi_start();

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

	mypapi_stop();
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
		printHopVerMsg();
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
	bgq_init_spinorfields_allprec(2 * k_max + 1);
	bgq_hm_init_allprec();

	bgq_update_backward_gauge();
#endif

#ifdef _GAUGE_COPY
	update_backward_gauge();
#endif


	for (int k = 0; k < k_max; k+=1) {
		random_spinor_field(g_spinor_field[k], VOLUME / 2, 0);
		bgq_transfer_spinorfield_double(true, g_spinorfields_double[k], g_spinor_field[k]);

		double compare_transfer = bgq_spinorfield_compare_double(true, g_spinorfields_double[k], g_spinor_field[k], false);
		assert(compare_transfer == 0 /* fields should be bit-identical*/);
	}

//#pragma omp parallel 
{
	for (int i = 0; i < 2; i+=1) {
		listprefetch_test(false);
		}
	for (int i = 0; i < 3; i+=1) {
		listprefetch_test(true);
	}
}
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

#pragma omp parallel
	{
		int k = 0;
		int k_max = 1;
		bgq_hmflags hmflags = hm_nocom * nocom | hm_nooverlap;
		double compare_even_before;
		double compare_even;

#pragma omp master
		{
			bgq_transfer_spinorfield_double(true, g_spinorfields_double[k], g_spinor_field[k]);
			compare_even_before = bgq_spinorfield_compare_double(true, g_spinorfields_double[k], g_spinor_field[k], false);
			assert(compare_even_before == 0 /* should be bitwise identical */);

			//bgq_initbgqref();
		}
#pragma omp barrier
		bgq_HoppingMatrix_double(false, g_spinorfields_double[k + k_max], g_spinorfields_double[k], g_gaugefield_double, hmflags);
#pragma omp master
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
#pragma omp master
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

#pragma omp parallel
	{
		int k = 0;
		int k_max = 1;
		bgq_hmflags hmflags = 0;
		double compare_even_before;
		double compare_even;

#pragma omp master
		{
			bgq_transfer_spinorfield_float(true, g_spinorfields_float[k], g_spinor_field[k]);
			compare_even_before = bgq_spinorfield_compare_float(true, g_spinorfields_float[k], g_spinor_field[k], false);
			assert(compare_even_before == 0 /* should be bitwise identical */);

			//bgq_initbgqref();
		}
#pragma omp barrier
		bgq_HoppingMatrix_float(false, g_spinorfields_float[k + k_max], g_spinorfields_float[k], g_gaugefield_float, hmflags);
#pragma omp master
		{
			Hopping_Matrix(0, g_spinor_field[k + k_max], g_spinor_field[k]);
			//bgq_savebgqref();
			compare_even = bgq_spinorfield_compare_float(false, g_spinorfields_float[k + k_max], g_spinor_field[k + k_max], false);
			assert(compare_even < 0.001);
		}
//#pragma omp barrier
		bgq_HoppingMatrix_float(true, g_spinorfields_float[2 * k_max], g_spinorfields_float[k + k_max], g_gaugefield_float, hmflags);
#pragma omp master
		{
			Hopping_Matrix(1, g_spinor_field[k + 2 * k_max], g_spinor_field[k + k_max]);
			double compare_odd = bgq_spinorfield_compare_float(true, g_spinorfields_float[k + 2 * k_max], g_spinor_field[k + 2 * k_max], false);
			assert(compare_odd < 0.001);

			master_print("Numerical instability between float precision implementations: even %g, odd %g\n", compare_even, compare_odd);
		}
	}
}



static double runcheck(bool sloppyprec, bgq_hmflags hmflags) {
	int k = 0;
	int k_max = 1;
	hmflags = hmflags & ~hm_nokamul;
	double result;
	double compare_even_before;
	double compare_even;

	if (!sloppyprec) {
#pragma omp master
		{
		memset(g_spinorfields_double[k + k_max], 0, sizeof(bgq_spinorsite_double) * VOLUME_SITES);
		memset(g_spinorfields_double[k + 2*k_max], 0, sizeof(bgq_spinorsite_double) * VOLUME_SITES);

		bgq_transfer_spinorfield_double(true, g_spinorfields_double[k], g_spinor_field[k]);
		double compare_even_before = bgq_spinorfield_compare_double(true, g_spinorfields_double[k], g_spinor_field[k], false);
		}
		bgq_HoppingMatrix_double(false, g_spinorfields_double[k + k_max], g_spinorfields_double[k], g_gaugefield_double, hmflags);
#pragma omp master
		{
		Hopping_Matrix_switch(0, g_spinor_field[k + k_max], g_spinor_field[k], hmflags & hm_nocom);
		double compare_even = bgq_spinorfield_compare_double(false, g_spinorfields_double[k + k_max], g_spinor_field[k + k_max], true);
		}
		bgq_HoppingMatrix_double(true, g_spinorfields_double[2 * k_max], g_spinorfields_double[k + k_max], g_gaugefield_double, hmflags);
#pragma omp master
		{
		Hopping_Matrix_switch(1, g_spinor_field[k + 2*k_max], g_spinor_field[k + k_max], hmflags & hm_nocom);
		double compare_odd = bgq_spinorfield_compare_double(true, g_spinorfields_double[k + 2*k_max], g_spinor_field[k + 2*k_max], true);

		result = max(max(compare_even_before, compare_even), compare_odd);
		}
	} else {
#pragma omp master
		{
		memset(g_spinorfields_float[k + k_max], 0, sizeof(bgq_spinorsite_float) * VOLUME_SITES);
		memset(g_spinorfields_float[k + 2*k_max], 0, sizeof(bgq_spinorsite_float) * VOLUME_SITES);

		bgq_transfer_spinorfield_float(true, g_spinorfields_float[k], g_spinor_field[k]);
		double compare_even_before = bgq_spinorfield_compare_float(true, g_spinorfields_float[k], g_spinor_field[k], true);
		}
		bgq_HoppingMatrix_float(false, g_spinorfields_float[k + k_max], g_spinorfields_float[k], g_gaugefield_float, hmflags);
#pragma omp master
		{
		Hopping_Matrix_switch(0, g_spinor_field[k + k_max], g_spinor_field[k], hmflags & hm_nocom);
		double compare_even = bgq_spinorfield_compare_float(false, g_spinorfields_float[k + k_max], g_spinor_field[k + k_max], true);
		}
		bgq_HoppingMatrix_float(true, g_spinorfields_float[2 * k_max], g_spinorfields_float[k + k_max], g_gaugefield_float, hmflags);
#pragma omp master
		{
		Hopping_Matrix_switch(1, g_spinor_field[k + 2*k_max], g_spinor_field[k + k_max], hmflags & hm_nocom);
		double compare_odd = bgq_spinorfield_compare_float(true, g_spinorfields_float[k + 2*k_max], g_spinor_field[k + 2*k_max], true);

		result = max(max(compare_even_before, compare_even), compare_odd);
		}
	}

	if (result > 0.001) {
		int a = 0;
	}
	return result;
}

#ifdef XLC
#include <l1p/sprefetch.h>
#endif

static benchstat runbench(int k_max, int j_max, bool sloppyprec, int ompthreads, bgq_hmflags hmflags) {
	int iterations = 0;

	double localsumtime = 0;
	double localsumsqtime = 0;
	double err = 0;

	const bool nocom = hmflags & hm_nocom;
	const bool nooverlap = hmflags & hm_nooverlap;
	const bool nokamul = hmflags & hm_nokamul;
	const bool kamul = !nokamul;
	const bool prefetchlist = hmflags & hm_prefetchlist;
	const bool prefetchstream = hmflags & hm_prefetchstream;
	const bool prefetchexplicit = hmflags & hm_prefetchexplicit;
	const bool noweylsend = hmflags & hm_noweylsend;
	const bool nobody = hmflags & hm_nobody;
	const bool nosurface = hmflags & hm_nosurface;
	const bool nol1plist = hmflags & hm_nol1plist;

	omp_set_num_threads(ompthreads);

#pragma omp parallel
	{
		//L1P_CHECK(L1P_SetStreamPolicy(L1P_stream_disable));

		//L1P_StreamPolicy_t pol;
		//L1P_CHECK(L1P_GetStreamPolicy(&pol));
		//if (pol != L1P_stream_disable)
		//	master_print("MK StreamPolicy not accepted\n");

	double errtmp = runcheck(sloppyprec, hmflags);

	for (int j = 0; j < j_max; j += 1) {
		////////////////////////////////////////////////////////////////////////////////
#pragma omp master
		{
		if (j == 1)
			mypapi_start();
		}
		double start_time = MPI_Wtime();
		if (sloppyprec) {
			for (int k = 0; k < 1; k += 1) {
				bgq_HoppingMatrix_float(false, g_spinorfields_float[k + k_max], g_spinorfields_float[k], g_gaugefield_float, hmflags);
				//bgq_HoppingMatrix_float(true, g_spinorfields_float[2 * k_max], g_spinorfields_float[k + k_max], g_gaugefield_float, hmflags);
#pragma omp master
		{
				iterations += 1;
		}
			}
		} else {
			for (int k = 0; k < 1; k += 1) {
				bgq_HoppingMatrix_double(false, g_spinorfields_double[k + k_max], g_spinorfields_double[k], g_gaugefield_double, hmflags);
				//bgq_HoppingMatrix_double(true, g_spinorfields_double[2 * k_max], g_spinorfields_double[k + k_max], g_gaugefield_double, hmflags);
#pragma omp master
		{
				iterations += 1;
		}
			}
		}
		////////////////////////////////////////////////////////////////////////////////
#pragma omp master
		{
		double end_time = MPI_Wtime();
		double runtime = end_time - start_time;

		if (j == 0) {
			// Just a warmup
			iterations = 0;
		} else {
			localsumtime += runtime;
			localsumsqtime += sqr(runtime);
		}
		err = errtmp;
	}
	}

	mypapi_stop();
	}


	assert(iterations == j_max-1);
	double localavgtime = localsumtime / iterations;
	double localavgsqtime = sqr(localavgtime);
	double localrmstime = sqrt((localsumsqtime / iterations) - localavgsqtime);

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

	int flops_per_lup = 1320;
	int flops_per_lup_body = 0;
	if (!nobody) {
		flops_per_lup_body += 1320;
		if (kamul)
			flops_per_lup_body += 8 * 2 * 3 * 6;
	}

	int flops_per_lup_surface = 0;
	if (!noweylsend)
		flops_per_lup_surface += 6 * 2 * 2;
	if (!nosurface) {
		flops_per_lup_surface += 1320 - (6*2*2);
		if (kamul)
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
	return result;
}

#define COUNTOF(arr) (sizeof(arr) / sizeof((arr)[0]))

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

static bgq_hmflags flags[] = { hm_nol1plist, 0 };
static char* flags_desc[] = { "noprefetch" , "list prefetch" };


static void print_repeat(const char * const str, const int count) {
	if (g_proc_id == 0) {
		for (int i = 0; i < count; i += 1) {
			printf("%s", str);
		}
	}
}


#define CELLWIDTH 15
#define SCELLWIDTH TOSTRING(CELLWIDTH)


static void exec_bench() {
	print_repeat("\n", 2);

	for (int i0 = 0; i0 < COUNTOF(kamuls); i0 += 1) {
		bool kamul = kamuls[i0];
		if (g_proc_id == 0) printf("Benchmark: %s\n", kamuls_desc[i0]);

		for (int i1 = 0; i1 < COUNTOF(sloppinessess); i1 += 1) {
			bool sloppiness = sloppinessess[i1];

			if (g_proc_id == 0) printf("Benchmarking precision: %s\n", sloppinessess_desc[i1]);
			if (g_proc_id == 0) printf("%10s|", "");
			for (int i3 = 0; i3 < COUNTOF(flags); i3 += 1) {
				if (g_proc_id == 0) printf("%"SCELLWIDTH"s|", flags_desc[i3]);
			}
			if (g_proc_id == 0) printf("\n");
			print_repeat("-", 10 + 1 + (CELLWIDTH + 1)*COUNTOF(flags));
			if (g_proc_id == 0) printf("\n");
			for (int i2 = 0; i2 < COUNTOF(omp_threads); i2 += 1) {
				int threads = omp_threads[i2];

				if (g_proc_id == 0) printf("%10s|", omp_threads_desc[i2]);

				for (int i3 = 0; i3 < COUNTOF(flags); i3 += 1) {
					bgq_hmflags hmflags = flags[i3];
					hmflags = hmflags | (!kamul * hm_nokamul);

					benchstat result = runbench(1, 3, sloppiness, threads, hmflags);

					char str[80] = {0};
					snprintf(str, sizeof(str), "%.0f mflops/s%s" , result.flops / MEGA, (result.error > 0.001) ? "X" : "");
					if (g_proc_id == 0) printf("%"SCELLWIDTH"s|", str);
					if (g_proc_id == 0) fflush(stdout);
				}
				if (g_proc_id == 0) printf("\n");
				print_repeat("-", 10 + 1 + (CELLWIDTH + 1)*COUNTOF(flags));
				if (g_proc_id == 0) printf("\n");
			}
			if (g_proc_id == 0) printf("\n");
		}

		if (g_proc_id == 0) printf("\n");
	}

	if (g_proc_id == 0) printf("Benchmark done\n");
}
