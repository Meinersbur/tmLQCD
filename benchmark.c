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
#include "bgq/bgq_field_double.h"
#include "bgq/bgq_HoppingMatrix.h"
#endif

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


double bgl_wtime() {
	return (MPI_Wtime());
}

int check_xchange();

int main(int argc, char *argv[])
{
	int j, j_max, k, k_max = 1;
#ifdef HAVE_LIBLEMON
	paramsXlfInfo *xlfInfo;
#endif

	/* GG */
	int intrig;
	FILE* yyingg = NULL;
	MPI_Status yyStatus;
	void* yybufgg = NULL;
	int yyCount;
	int yyfd;
	char *input_filename = NULL;
	int c;

	static double t1, t2, dt, sdt, dts, qdt, sqdt;
	double antioptaway = 0.0;
#ifdef MPI
	static double dt2;

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
		yyCount = read(yyfd, yybufgg, 8192);
		intrig = close(yyfd);
	}
	intrig = MPI_Bcast(yybufgg, 8192, MPI_CHAR, 0, MPI_COMM_WORLD );
	yyingg = fmemopen(yybufgg, strlen(yybufgg), "r");
	intrig = read_input_fh(yyingg);
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

	check_geometry();
#if (defined MPI && !(defined _USE_SHMEM))
	check_xchange();
#endif

	start_ranlux(1, 123456);
	random_gauge_field(reproduce_randomnumber_flag);

#ifdef MPI
	/*For parallelization: exchange the gaugefield */
	xchange_gauge();
#endif

#ifdef BGQ
  bgq_init_gaugefield();
  bgq_init_spinorfields(2 * k_max + 1);
  bgq_hm_init();

  update_backward_gauge();
#endif

#ifdef _GAUGE_COPY
	update_backward_gauge();
#endif

	if (even_odd_flag) {
		/*initialize the pseudo-fermion fields*/
		j_max = 1;
		sdt = 0.0;
		for (k = 0; k < k_max; k++) {
			random_spinor_field(g_spinor_field[k], VOLUME / 2, 0);
#if BGQ
			bgq_transfer_spinorfield_double(true, g_spinorfields_double[k], g_spinor_field[k]);
#endif
		}

		int again = 0;
		int iterations = 0;
		while (sdt < 30.) {
			if (again && !g_proc_id)
				fprintf(stderr, "MK_Running again %d-th time, sdt=%f\n", again, sdt);
			again++;
#ifdef MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif
			if (sdt > 17)
				mypapi_start();
#if 1
			t1 = bgl_wtime();
#else
			t1 = (double) clock();
#endif
			antioptaway = 0.0;
			for (j = 0; j < j_max; j++) {
				for (k = 0; k < 1; k++) {
					Hopping_Matrix(0, g_spinor_field[k + k_max], g_spinor_field[k]);
					Hopping_Matrix(0, g_spinor_field[k + k_max], g_spinor_field[k]);
					//Hopping_Matrix(1, g_spinor_field[2 * k_max], g_spinor_field[k + k_max]);
					antioptaway += g_spinor_field[2 * k_max][0].s0.c0.re;
					iterations += 1;
				}
			}
#if 1
			t2 = bgl_wtime();
			dt = t2 - t1;
#else
			t2 = (double) clock();
			dt = (t2 - t1) / ((double) (CLOCKS_PER_SEC));
#endif
			if (sdt > 17) {
				mypapi_stop();
			}
#ifdef MPI
			MPI_Allreduce(&dt, &sdt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
			sdt = dt;
#endif
			qdt = dt * dt;
#ifdef MPI
			MPI_Allreduce(&qdt, &sqdt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
			sqdt = qdt;
#endif
			sdt = sdt / ((double) g_nproc);
			sqdt = sqrt(sqdt / g_nproc - sdt * sdt);
			j_max *= 2;
		}
		j_max = j_max / 2;
		double sdtsave = sdt;
		dts = dt;
		sdt = 1.0e6 * sdt / ((double) (k_max * j_max * (VOLUME)));
		sqdt = 1.0e6 * sqdt / ((double) (k_max * j_max * (VOLUME)));

		if (g_proc_id == 0) {
			printf("MK VOLUME=%d\n", VOLUME);
			printf("MK %d full iterations in %f secs avg\n", iterations, sdtsave/(double)iterations);
			printf("Print 1 result, to make sure that the calculation is not optimized away: %e  \n", antioptaway);
			printf("total time %e secs, Variance of the time %e sec. (iterations=%d). \n", sdt, sqdt, j_max);
			printf("\n");
			printf(" (%d Mflops [%d bit arithmetic])\n",
			        (int) (1320.0 / sdt), (int) sizeof(spinor) / 3);
			printf("\n");
			fflush(stdout);
		}

		mypapi_start();
#ifdef MPI
		/* isolated computation */
#if 1
		t1 = bgl_wtime();
#else
		t1 = (double) clock();
#endif
		iterations = 0;
		antioptaway = 0.0;
		for (j = 0; j < j_max; j++) {
			for (k = 0; k < k_max; k++) {
				Hopping_Matrix_nocom(0, g_spinor_field[k + k_max], g_spinor_field[k]);
				Hopping_Matrix_nocom(1, g_spinor_field[2 * k_max], g_spinor_field[k + k_max]);
				antioptaway += g_spinor_field[2 * k_max][0].s0.c0.re;
				iterations += 1;
			}
		}
#if 1
		t2 = bgl_wtime();
		dt2 = t2 - t1;
#else
		t2 = (double) clock();
		dt2 = (t2 - t1) / ((double) (CLOCKS_PER_SEC));
#endif
		mypapi_stop();

		/* compute the bandwidth */
		double dt2_single = dt2;
		dt = dts - dt2;
		MPI_Allreduce(&dt, &sdt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
		sdt = sdt / ((double) g_nproc);
		MPI_Allreduce(&dt2, &dt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
		dt = dt / ((double) g_nproc);
		dt = 1.0e6f * dt / ((double) (k_max * j_max * (VOLUME)));
		if (g_proc_id == 0) {
			printf("MKnocom %d iterations in %f secs svg\n", iterations, sdt/iterations);
			printf("Print 1 result, to make sure that the calculation is not optimized away: %e  \n", antioptaway);
			printf("communication switched off \n");
			printf("total time %e sec, Variance of the time %e sec. (iterations=%d). \n", sdt, sqdt, j_max);
			printf(" (%d Mflops [%d bit arithmetic])\n",
			        (int) (1320.0f / dt), (int) sizeof(spinor) / 3);
			printf("MK_ (Control: %.2f secs, %d Mflops)\n", dt2_single, (int) (1320.0 / (1.0e6 * dt2_single / (double) (k_max * j_max * VOLUME)))),
			        printf("\n");
			fflush(stdout);
		}
		sdt = sdt / ((double) k_max);
		sdt = sdt / ((double) j_max);
		sdt = sdt / ((double) (2 * SLICE));
		if (g_proc_id == 0) {
			printf("The size of the package is %d Byte \n", (SLICE) * 192);
#ifdef _USE_HALFSPINOR
			printf("The bandwidth is %5.2f + %5.2f   MB/sec\n",
			        192. / sdt / 1024 / 1024, 192. / sdt / 1024. / 1024);
#else
			printf("The bandwidth is %5.2f + %5.2f   MB/sec\n",
					2.*192./sdt/1024/1024, 2.*192./sdt/1024./1024);
#endif
			fflush(stdout);
		}
#endif
		fflush(stdout);
	}
	else {
		fprintf(stderr, "MK !even_odd_flag\n");

		/* the non even/odd case now */
		/*initialize the pseudo-fermion fields*/
		j_max = 1;
		sdt = 0.;
		for (k = 0; k < k_max; k++) {
			random_spinor_field(g_spinor_field[k], VOLUME, 0);
		}

		while (sdt < 3.) {
#ifdef MPI
			MPI_Barrier(MPI_COMM_WORLD );
#endif
#if 1
			t1 = bgl_wtime();
#else
			t1 = (double) clock();
#endif
			for (j = 0; j < j_max; j++) {
				for (k = 0; k < k_max; k++) {
					D_psi(g_spinor_field[k + k_max], g_spinor_field[k]);
					antioptaway += g_spinor_field[k + k_max][0].s0.c0.re;
				}
			}
#if 1
			t2 = bgl_wtime();
			dt = t2 - t1;
#else
			t2 = (double) clock();
			dt = (t2 - t1) / ((double) (CLOCKS_PER_SEC));
#endif
#ifdef MPI
			MPI_Allreduce(&dt, &sdt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#else
			sdt = dt;
#endif
			qdt = dt * dt;
#ifdef MPI
			MPI_Allreduce(&qdt, &sqdt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#else
			sqdt = qdt;
#endif
			sdt = sdt / ((double) g_nproc);
			sqdt = sqrt(sqdt / g_nproc - sdt * sdt);
			j_max *= 2;
		}
		j_max = j_max / 2;
		dts = dt;
		sdt = 1.0e6f * sdt / ((double) (k_max * j_max * (VOLUME)));
		sqdt = 1.0e6f * sqdt / ((double) (k_max * j_max * (VOLUME)));

		if (g_proc_id == 0) {
			printf("Print 1 result, to make sure that the calculation is not optimized away: %e  \n", antioptaway);
			printf("total time %e sec, Variance of the time %e sec. (iterations=%d). \n", sdt, sqdt, j_max);
			printf("\n");
			printf(" (%d Mflops [%d bit arithmetic])\n",
			        (int) (1392.0f / sdt), (int) sizeof(spinor) / 3);
			printf("\n");
			fflush(stdout);
		}
	}

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
