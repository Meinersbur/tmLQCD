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

	MPI_Init(&argc, &argv);
	/* GG */
	MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);
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

	update_backward_gauge();
#endif

#ifdef _GAUGE_COPY
	update_backward_gauge();
#endif



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
} benchstat;

benchstat runbench(int k_max, int j_max, bool sloppyprec, int ompthreads, bool nocom, bool nooverlap) {
	int iterations = 0;

	double localsumtime = 0;
	double localsumsqtime = 0;

	omp_set_num_threads(ompthreads);

	bgq_hmflags hmflags = 0;

	hmflags |= nocom*hm_nocom;
	hmflags |= nooverlap*hm_nooverlap;

	for (int j = 0; j < j_max; j += 1) {
		////////////////////////////////////////////////////////////////////////////////
		double start_time = MPI_Wtime();
		if (sloppyprec) {
			for (int k = 0; k < 1; k += 1) {
				bgq_HoppingMatrix_float(false, g_spinorfields_float[k + k_max], g_spinorfields_float[k], g_gaugefield_float, hmflags);
				bgq_HoppingMatrix_float(true, g_spinorfields_float[2 * k_max], g_spinorfields_float[k + k_max], g_gaugefield_float,hmflags);
				iterations += 1;
			}
		} else {
			for (int k = 0; k < 1; k += 1) {
				bgq_HoppingMatrix_double(false, g_spinorfields_double[k + k_max], g_spinorfields_double[k], g_gaugefield_double, hmflags);
				bgq_HoppingMatrix_double(true, g_spinorfields_double[2 * k_max], g_spinorfields_double[k + k_max], g_gaugefield_double, hmflags);
				iterations += 1;
			}
		}
		////////////////////////////////////////////////////////////////////////////////
		double end_time = MPI_Wtime();
		double runtime = end_time - start_time;
		localsumtime += runtime;
		localsumsqtime += sqr(runtime);
	}

	double localavgtime = localsumtime / j_max;
	double localavgsqtime = sqr(localavgtime);
	double localrmstime = sqrt((localsumsqtime / j_max) - localavgsqtime);

	double sumtime = -1;
	double sumsqtime = -1;
	MPI_Allreduce(&sumtime, &localavgtime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce(&sumsqtime, &localavgsqtime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

	double avgtime = sumtime / g_nproc;
	double rmstime = sqrt((sumsqtime / j_max) - sqr(sumtime));

	double lups = iterations * LOCAL_LT * LOCAL_LX * LOCAL_LY * LOCAL_LZ;
	double flops_per_lup = 1320;
	double flops = lups * flops_per_lup;

	benchstat result;

	result.avgtime = avgtime;
	result.localrmstime = localrmstime;
	result.globalrmstime = rmstime;

	result.lups = lups;
	result.flops = flops / avgtime;

	return result;
}

#define COUNTOF(arr) (sizeof(arr) / sizeof((arr)[0]))

static bool sloppinessess[] = { false, true };
static char* sloppinessess_desc[] = { "double", "float" };

static int omp_threads[] = { 1, 2, 4, 8, 16, 32, 64 };
static char* omp_threads_desc[] = { "1", "2", "4", "8", "16", "32", "64" };

static struct {
	bool nocom;
	bool comnooverlap;
} coms[] = { { false, false }, { false, true }, { true, -1 } };
static char* com_desc[] = { "Com async", "Com sync", "Com off" };


static void print_repeat(const char * const str, const int count) {
	if (g_proc_id == 0) {
		for (int i = 0; i < count; i += 1) {
			printf("%s", str);
		}
	}
}

#define KILO 1000.0
#define MEGA (1000.0 * 1000.0)
#define GIGA (1000.0 * 1000.0 * 1000.0)
#define CELLWIDTH 15
#define SCELLWIDTH TOSTRING(CELLWIDTH)

static void exec_bench() { print_repeat("\n", 2);
	for (int i1 = 0; i1 < COUNTOF(sloppinessess); i1 += 1) {
		bool sloppiness = sloppinessess[i1];

		if (g_proc_id == 0) printf("Benchmarking precision: %s\n", sloppinessess_desc[i1]);
		if (g_proc_id == 0) printf("%10s|", "");
		for (int i3 = 0; i3 < COUNTOF(coms); i3 += 1) {
			if (g_proc_id == 0) printf("%"SCELLWIDTH"s|", com_desc[i3]);
		}
		if (g_proc_id == 0) printf("\n");
		print_repeat("-", CELLWIDTH + 1 + (CELLWIDTH + 1)*COUNTOF(coms));
		if (g_proc_id == 0) printf("\n");
		for (int i2 = 0; i2 < COUNTOF(omp_threads); i2 += 1) {
			int threads = omp_threads[i2];

			if (g_proc_id == 0) printf("%10s|\n", omp_threads_desc[i2]);

			for (int i3 = 0; i3 < COUNTOF(coms); i3 += 1) {
				bool nocom = coms[i3].nocom;
				bool nooverlap = coms[i3].comnooverlap;

				benchstat result = runbench(1, 2, sloppiness, threads, nocom, nooverlap);

				char str[80] = {0};
				snprintf( str, sizeof(str), "%f.0 mflops/s" , result.flops / MEGA);
				if (g_proc_id == 0) printf("%"SCELLWIDTH"s|", str);
				if (g_proc_id == 0) fflush(stdout);
			}
			if (g_proc_id == 0) printf("\n");
			print_repeat("-", CELLWIDTH + 1 + (CELLWIDTH + 1)*COUNTOF(coms));
			if (g_proc_id == 0) printf("\n");
		}
		if (g_proc_id == 0) printf("\n");
	}

	if (g_proc_id == 0) printf("Benchmark done\n");
}

