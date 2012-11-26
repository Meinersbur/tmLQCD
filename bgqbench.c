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
#include "config.h"
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
#include "bgq/bgq_gaugefield.h"
#include "bgq/bgq_spinorfield.h"
#include "bgq/bgq_HoppingMatrix.h"
#include "bgq/bgq_qpx.h"
#include "bgq/bgq_utils.h"
#include "bgq/bgq_dispatch.h"
#include "bgq/bgq_comm.h"
#include <omp.h>
#include "bgq/mypapi.h"
#include <getopt.h>
#define __USE_GNU
#include <fenv.h>

#ifdef XLC
#include <l1p/pprefetch.h>
#include <l1p/sprefetch.h>
#endif

//typedef bgq_spinorfield bgq_spinorfield_double;
//typedef bgq_spinorfield bgq_spinorfield_float;

//typedef bgq_gaugefield bgq_gaugefield_double;
//typedef bgq_gaugefield bgq_gaugefield_float;

//typedef bgq_gaugesite bgq_gaugesite_double;
//typedef bgq_gaugesite bgq_gaugesite_float;

typedef void (*benchfunc_t)(bgq_hmflags flags, int k, int k_max);
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

typedef struct {
	double avgtime;
	double localrmstime;
	double globalrmstime;
	double totcycles;
	double localavgflop;

	ucoord sites_body;
	ucoord sites_surface;
	ucoord sites;

	ucoord lup_body;
	ucoord lup_surface;
	ucoord lup;
	//double flops;

	double error;

	mypapi_counters counters;
	bgq_hmflags opts;
	double avgovhtime;
} benchstat;

typedef struct {
	int j_max;
	int k_max;
	benchfunc_t benchfunc;
	bgq_hmflags opts;
	benchstat result;
} master_args;

static ucoord flop_per_bodysite(bgq_hmflags opts) {
	ucoord result = 0;

	if (!(opts & hm_nobody)) {
		// Compute weyl
		result += 8/*dirs*/ * (2 * 3)/*cmplx per weyl*/ * 2/*flops*/;

		// Su3 M*V
		result += 8/*dirs*/ * 2/*su3vec per weyl*/ * (9*6 + 6*2)/*flop per su3 mv-mul*/;

		// Accummulate spinor
		// Assuming readWeyllayout:
		result += 7/*dirs*/ * (4 * 3)/*cmplx per spinor*/ * 2/*flops accumm*/;
		assert(result == 1320);

		if (!(opts & hm_nokamul)) {
			result += 6/*flops cmplx mul*/* (2 * 3)/*cmplx per weyl*/* 8/*dirs*/;
		}
	}

	return result;
}


static ucoord flop_per_surfacesite(bgq_hmflags opts) {
	ucoord result = 0;

	if (!(opts & hm_nodistribute)) {
		// Compute weyl
		result += 8/*dirs*/ * (2 * 3)/*cmplx per weyl*/ * 2/*flops*/;

		// Su3 M*V
		result += 8/*dirs*/ * 2/*su3vec per weyl*/ * (9*6 + 6*2)/*flop per su3 mv-mul*/;

		// Accummulate spinor
		// Assuming readWeyllayout:
		result += 7/*dirs*/ * (4 * 3)/*cmplx per spinor*/ * 2/*flops accumm*/;
		assert(result == 1320);

		if (!(opts & hm_nokamul)) {
			result += 6/*flops cmpl mul*/* (2 * 3)/*cmpl per weyl*/* 8/*dirs*/;
		}
	}

	return result;
}


static uint64_t compute_flop(bgq_hmflags opts, uint64_t lup_body, uint64_t lup_surface) {
	uint64_t flop_body = flop_per_bodysite(opts) * lup_body;
	uint64_t flop_surface = flop_per_surfacesite(opts) * lup_surface;
	return flop_body+flop_surface;
}


static void benchmark_setup_worker(void *argptr, size_t tid, size_t threads) {
	mypapi_init();

#ifndef NDEBUG
	//feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW|FE_UNDERFLOW);
#endif

#if BGQ_QPX
	const master_args *args = argptr;
	const int j_max = args->j_max;
	const int k_max = args->k_max;
	const benchfunc_t benchfunc = args->benchfunc;
	const bgq_hmflags opts = args->opts;
	const bool nocom = opts & hm_nocom;
	const bool nooverlap = opts & hm_nooverlap;
	const bool nokamul = opts & hm_nokamul;
	const bool noprefetchlist = opts & hm_noprefetchlist;
	const bool noprefetchstream = opts & hm_noprefetchstream;
	const bool noprefetchexplicit = opts & hm_noprefetchexplicit;
	const bool noweylsend = opts & hm_noweylsend;
	const bool nobody = opts & hm_nobody;
	const bool nosurface = opts & hm_nosurface;
	const bool experimental = opts & hm_experimental;
	const bgq_hmflags implicitprefetch = opts & (hm_prefetchimplicitdisable | hm_prefetchimplicitoptimistic | hm_prefetchimplicitconfirmed);

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

	//if (g_proc_id==0) {
	//	L1P_StreamPolicy_t poli;
	//	L1P_CHECK(L1P_GetStreamPolicy(&poli));
	//	printf("Prefetch policy: %d\n",poli);
	//}

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

#endif
}

static void benchmark_free_worker(void *argptr, size_t tid, size_t threads) {
	mypapi_free();
}

static void HoppingMatrix_switch(bool isOdd, spinor *l, spinor *k, bgq_hmflags hmflags) {
	bool nocom = hmflags & hm_nocom;
	if (nocom) {
		Hopping_Matrix_nocom(isOdd, l, k);
	} else {
		Hopping_Matrix(isOdd, l, k);
	}
}

static double runcheck(bgq_hmflags hmflags, size_t k_max) {
	const size_t k = 0;
	// To ensure that zero is used in case of nocomm
	bgq_spinorfield_setup(&g_bgq_spinorfields[k], true, false, false, false, true, false);
	bgq_spinorfield_setup(&g_bgq_spinorfields[k + k_max], false, false, false, false, true, false);
	// Flow:
	// [k]isOdd -> [k+k_max]isEven -> [k]isOdd

	bgq_spinorfield_transfer(true, &g_bgq_spinorfields[k], g_spinor_field[k]);
#ifndef BGQ_COORDCHECK
	double compare_transfer = bgq_spinorfield_compare(true, &g_bgq_spinorfields[k], g_spinor_field[k], true);
	assert(compare_transfer == 0);
	// Must be exact copy
#endif

	bgq_HoppingMatrix(false, &g_bgq_spinorfields[k + k_max], &g_bgq_spinorfields[k], hmflags);
	HoppingMatrix_switch(false, g_spinor_field[k + k_max], g_spinor_field[k], hmflags);
#ifndef BGQ_COORDCHECK
	double compare_even = bgq_spinorfield_compare(false, &g_bgq_spinorfields[k + k_max], g_spinor_field[k + k_max], true);
	assert(compare_even < 0.01);
#endif

#ifndef BGQ_COORDCHECK
	bgq_spinorfield_transfer(false, &g_bgq_spinorfields[k + k_max], g_spinor_field[k + k_max]); 	// We don't want to accumulate errors
	compare_transfer = bgq_spinorfield_compare(false, &g_bgq_spinorfields[k + k_max], g_spinor_field[k + k_max], true);
	assert(compare_transfer == 0);
	// Must be exact copy
#endif

	bgq_HoppingMatrix(true, &g_bgq_spinorfields[k], &g_bgq_spinorfields[k + k_max], hmflags);
	HoppingMatrix_switch(true, g_spinor_field[k], g_spinor_field[k + k_max], hmflags);
#ifndef BGQ_COORDCHECK
	double compare_odd = bgq_spinorfield_compare(true, &g_bgq_spinorfields[k], g_spinor_field[k], true);
	assert(compare_odd < 0.01);
#else
	bgq_spinorfield_setup(&g_bgq_spinorfields[k], true, false, false, true, false, false); // Wait for data transmission
#endif

#ifndef BGQ_COORDCHECK
	return max(compare_even, compare_odd);
#else
	return 0;
#endif
}


typedef struct {
	int set;
	mypapi_counters result;
} mypapi_work_t;

static void mypapi_start_worker(void *arg_untyped, size_t tid, size_t threads) {
	mypapi_work_t *arg = arg_untyped;
	mypapi_start(arg->set);
}

static void mypapi_stop_worker(void *arg_untyped, size_t tid, size_t threads) {
	mypapi_work_t *arg = arg_untyped;
	arg->result = mypapi_stop();
}


static void donothing(void *arg, size_t tid, size_t threads) {
#if BGQ_QPX
	DelayTimeBase(1600*100);
#endif
}

static int benchmark_master(void *argptr) {
	master_args * const args = argptr;
	const int j_max = args->j_max;
	const int k_max = args->k_max;
	const benchfunc_t benchfunc = args->benchfunc;
	const bgq_hmflags opts = args->opts;
	const bool nocom = opts & hm_nocom;
	//const bool nooverlap = opts & hm_nooverlap;
	//const bool nokamul = opts & hm_nokamul;
	//const bool noprefetchlist = opts & hm_noprefetchlist;
	//const bool noprefetchstream = opts & hm_noprefetchstream;
	//const bool noprefetchexplicit = opts & hm_noprefetchexplicit;
	//const bool noweylsend = opts & hm_noweylsend;
	//const bool nobody = opts & hm_nobody;
	//const bool nosurface = opts & hm_nosurface;
	//const bool experimental = opts & hm_experimental;
	bool floatprecision = opts & hm_floatprecision;
	//const bgq_hmflags implicitprefetch = opts & (hm_prefetchimplicitdisable | hm_prefetchimplicitoptimistic | hm_prefetchimplicitconfirmed);
	bool withcheck = opts & hm_withcheck;

	// Setup thread options (prefetch setting, performance counters, etc.)
	bgq_master_call(&benchmark_setup_worker, argptr);


	double err = 0;
	if (withcheck) {
		err = runcheck(opts, k_max);
	}


	// Give us a fresh environment
	for (ucoord k = 0; k <= 2*k_max; k+=1) {
		if (g_bgq_spinorfields[k].isInitialized) {
			if (g_bgq_spinorfields[k].sec_collapsed_double)
				memset(g_bgq_spinorfields[k].sec_collapsed_double, 0, PHYSICAL_VOLUME * sizeof(*g_bgq_spinorfields[k].sec_collapsed_double));
		}
	}
	if (nocom) {
		for (ucoord d = 0; d < PHYSICAL_LD; d+=1) {
			memset(g_bgq_sec_recv_double[d], 0, bgq_section_size(bgq_direction2section(d, false)));
		}
	}
	random_spinor_field(g_spinor_field[0], VOLUME / 2, 0);
	bgq_spinorfield_transfer(true, &g_bgq_spinorfields[0], g_spinor_field[0]);


	uint64_t sumotime = 0;
	for (int i = 0; i < 20; i += 1) {
		uint64_t start_time = bgq_wcycles();
		donothing(NULL, 0, 0);
		uint64_t mid_time = bgq_wcycles();
		bgq_master_call(&donothing, NULL);
		uint64_t stop_time = bgq_wcycles();

		uint64_t time = (stop_time - mid_time) - (mid_time- start_time);
		sumotime += time;
	}
	double avgovhtime = (double)sumotime / 20.0;

	static mypapi_work_t mypapi_arg;
	double localsumtime = 0;
	double localsumsqtime = 0;
	uint64_t localsumcycles=0;
	uint64_t localsumflop = 0;
	mypapi_counters counters;
	counters.init = false;
	int iterations = 1; // Warmup phase
	iterations += j_max;
	if (iterations < 1 + MYPAPI_SETS)
		iterations = 1 + MYPAPI_SETS;

	for (int i = 0; i < iterations; i += 1) {
		//master_print("Starting iteration %d of %d\n", j+1, iterations);
		bool isWarmup = (i == 0);
		int j = i - 1;
		bool isPapi = !isWarmup && (i >= iterations - MYPAPI_SETS);
		int papiSet = i - (iterations - MYPAPI_SETS);
		bool isJMax = (0 <= j) && (j < j_max);

		double start_time;
		uint64_t start_cycles;
		uint64_t start_flop;
		if (isJMax) {
			start_flop = flopaccumulator;
		}
		if (isPapi) {
			bgq_master_sync();
			mypapi_arg.set = papiSet;
			bgq_master_call(mypapi_start_worker, &mypapi_arg);
		}
		if (isJMax) {
			start_time = MPI_Wtime();
			start_cycles = bgq_wcycles();
		}

		{
			// The main benchmark
			for (int k = 0; k < k_max; k += 1) {
				// Note that flops computation assumes that readWeyllayout is used
				benchfunc(opts, k, k_max);
			}

			bgq_master_sync(); // Wait for all threads to finish, to get worst thread timing
		}

		double end_time;
		uint64_t end_cycles;
		if (isJMax) {
			end_cycles = bgq_wcycles();
			end_time = MPI_Wtime();
		}
		if (isPapi) {
			bgq_master_call(mypapi_stop_worker, &mypapi_arg);
			counters = mypapi_merge_counters(&counters, &mypapi_arg.result);
		}

		if (isJMax) {
			double duration = end_time - start_time;
			localsumtime += duration;
			localsumsqtime += sqr(duration);
			localsumcycles += (end_cycles - start_cycles);
			localsumflop += (flopaccumulator - start_flop);
		}
	}

	bgq_master_call(&benchmark_free_worker, argptr);

	ucoord its = j_max;
	double localavgtime = localsumtime / its;
	double localavgsqtime = sqr(localavgtime);
	double localrmstime = sqrt((localsumsqtime / its) - localavgsqtime);
	double localcycles = (double)localsumcycles / (double)its;
	double localavgflop = (double)localsumflop / (double)its;

	double localtime[] = { localavgtime, localavgsqtime, localrmstime };
	double sumreduce[3] = { -1, -1, -1 };
	MPI_Allreduce(&localtime, &sumreduce, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	double sumtime = sumreduce[0];
	double sumsqtime = sumreduce[1];
	double sumrmstime = sumreduce[2];

	double avgtime = sumtime / g_nproc;
	double avglocalrms = sumrmstime / g_nproc;
	double rmstime = sqrt((sumsqtime / g_nproc) - sqr(avgtime));

	// Assume k_max lattice site updates, even+odd sites
	ucoord sites = LOCAL_LT * LOCAL_LX * LOCAL_LY * LOCAL_LZ;
	ucoord sites_body = PHYSICAL_BODY * PHYSICAL_LK * PHYSICAL_LP;
	ucoord sites_surface = PHYSICAL_SURFACE * PHYSICAL_LK * PHYSICAL_LP;
	assert(sites == sites_body+sites_surface);
	assert(sites == VOLUME);

	ucoord lup_body = k_max * sites_body;
	ucoord lup_surface = k_max * sites_surface;
	ucoord lup = lup_body + lup_surface;
	assert(lup == k_max * sites);

	benchstat *result = &args->result;
	result->avgtime = avgtime;
	result->localrmstime = avglocalrms;
	result->globalrmstime = rmstime;
	result->totcycles = localcycles;
	result->localavgflop = localavgflop;

	result->sites_surface = sites_surface;
	result->sites_body = sites_body;
	result->sites = sites;
	result->lup_surface = lup_surface;
	result->lup_body = lup_body;
	result->lup = lup;
	//result->flops = flops / avgtime;
	result->error = err;
	result->counters = counters;
	result->opts = opts;
	result->avgovhtime = avgovhtime;
	return EXIT_SUCCESS;
}

static benchstat runbench(benchfunc_t benchfunc, bgq_hmflags opts, int k_max, int j_max, int ompthreads) {
	omp_set_num_threads(ompthreads);
	master_args args = {
	        .j_max = j_max,
	        .k_max = k_max,
	        .benchfunc = benchfunc,
	        .opts = opts
	};
	int retcode = bgq_parallel(&benchmark_master, &args);
	assert(retcode == EXIT_SUCCESS);
	if (retcode != EXIT_SUCCESS) {
		exit(retcode);
	}

	return args.result;
}

static void print_repeat(const char * const str, const int count) {
	if (g_proc_id == 0) {
		for (int i = 0; i < count; i += 1) {
			printf("%s", str);
		}
	}
}

static bool kamuls[] = { false, true };
static char *kamuls_desc[] = { "dslash", "kamul" };

static bool sloppinesses[] = { false, true };
static char *sloppinesses_desc[] = { "double", "float" };

static int omp_threads[] = { 1, 2, 4, 8, 16, 32, 33, 48, 56, 64 };
static char *omp_threads_desc[] = { "1","2", "4", "8", "16", "32", "33", "48", "56", "64" };

#define DEFOPTS (hm_noprefetchlist | hm_nokamul)

static bgq_hmflags flags[] = {
        (DEFOPTS | hm_withcheck) & ~hm_nokamul,
        (DEFOPTS | hm_floatprecision | hm_withcheck) & ~hm_nokamul,
        DEFOPTS,
		DEFOPTS | hm_nospi,
		DEFOPTS | hm_nooverlap,
		DEFOPTS | hm_nocom,
		DEFOPTS | hm_nocom |             hm_nodistribute | hm_nodatamove,
		DEFOPTS | hm_nocom | hm_nobody |                   hm_nodatamove,
		DEFOPTS | hm_nocom | hm_nobody | hm_nodistribute                ,
		DEFOPTS | hm_nocom                               | hm_nodatamove,
		DEFOPTS | hm_nocom                               | hm_nodatamove | hm_floatprecision,
		DEFOPTS |            hm_nobody | hm_nodistribute | hm_nodatamove,
		DEFOPTS | hm_nospi | hm_nobody | hm_nodistribute | hm_nodatamove,
		DEFOPTS | hm_nocom | hm_nobody | hm_nodistribute | hm_nodatamove,
		DEFOPTS | hm_noprefetchstream | hm_prefetchimplicitdisable,
		DEFOPTS                       | hm_prefetchimplicitconfirmed,
		DEFOPTS | hm_noprefetchstream | hm_prefetchimplicitconfirmed,
		DEFOPTS | hm_noprefetchstream | hm_prefetchimplicitoptimistic
    };
static char* flags_desc[] = {
		"kamul dbl",
		"kamul sgl",
		"nokamul",
		"MPI",
		"+nooverlap",
		"+nocomm",
		"bodyonly",
		"distonly",
		"dmovonly",
		"volonly dbl",
		"volonly sgl",
		"SPI only",
		"MPI only",
		"idle",
		"pf disable",
		"pf stream",
		"pf confirmed",
		"pf optimistic"
	};

#define CELLWIDTH 15
#define SCELLWIDTH TOSTRING(CELLWIDTH)

static void print_stats(benchstat *stats) {
#if PAPI
	int threads = omp_get_num_threads();

	for (mypapi_interpretations j = 0; j < __pi_COUNT; j+=1) {
		printf("%10s|", "");
		char *desc = NULL;

		for (int i3 = 0; i3 < lengthof(flags); i3 += 1) {
			char str[80];
			str[0] = '\0';
			benchstat *stat = &stats[i3];
			bgq_hmflags opts = stat->opts;

			double avgtime = stat->avgtime;
			uint64_t lup = stat->lup;
			uint64_t flop = compute_flop(opts, stat->lup_body, stat->lup_surface);
			double flops = (double)flop/stat->avgtime;
			double localrms = stat->localrmstime / stat->avgtime;
			double globalrms = stat->globalrmstime / stat->avgtime;
			ucoord sites = stat->sites;

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

#if 0
			double nNecessaryInstr = 0;
			if (!(opts & hm_nobody))
			nNecessaryInstr += bodySites * (240/*QFMA*/+ 180/*QMUL+QADD*/+ 180/*LD+ST*/)/2;
			if (!(opts & hm_noweylsend))
			nNecessaryInstr += haloSites * (2*3*2/*QMUL+QADD*/+ 4*3/*LD*/+ 2*3/*ST*/)/2;
			if (!(opts & hm_nosurface))
			nNecessaryInstr += surfaceSites * (240/*QFMA*/+ 180/*QMUL+QADD*/+ 180/*LD+ST*/- 2*3*2/*QMUL+QADD*/- 4*3/*LD*/+ 2*3/*LD*/)/2;
			if (!(opts & hm_nokamul))
			nNecessaryInstr += sites * (8*2*3*1/*QFMA*/+ 8*2*3*1/*QMUL*/)/2;
#endif

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
			case pi_msecs:
				desc = "Iteration time";
				snprintf(str, sizeof(str), "%.3f mSecs",stat->avgtime/MILLI);
				break;
			case pi_cycpersite:
				desc = "per site update";
				snprintf(str, sizeof(str), "%.1f cyc", stat->totcycles / lup);
				break;
			case pi_instrpersite:
				desc = "instr per update";
				snprintf(str, sizeof(str), "%.1f", nInstructions / lup);
				break;
			case pi_fxupersite:
				desc = "FU instr per update";
				snprintf(str, sizeof(str), "%.1f", nAXUInstr / lup);
				break;
			case pi_flops:
				desc = "MFlop/s";
				snprintf(str, sizeof(str), "%.0f MFlop/s", flops/MEGA);
				break;
			case pi_flopsref:
				desc = "Speed";
				snprintf(str, sizeof(str), "%.0f MFlop/s", stat->localavgflop / (avgtime * MEGA));
				break;
			case pi_floppersite:
				desc = "Flop per site";
				snprintf(str, sizeof(str), "%.1f Flop", stat->localavgflop / sites);
				break;
			case pi_localrms:
				desc = "Thread RMS";
				snprintf(str, sizeof(str), "%.1f %%", 100.0*localrms);
				break;
			case pi_globalrms:
				desc = "Node RMS";
				snprintf(str, sizeof(str), "%.1f %%", 100.0*globalrms);
				break;
			case pi_avgovhtime:
				desc = "Threading overhead";
				snprintf(str, sizeof(str), "%.1f cyc", stat->avgovhtime);
				break;
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
				snprintf(str, sizeof(str), "%.3f cpi", nCycles / nInstructions);
				break;
			case pi_corecpi:
				desc = "Cycles per instruction (Core)";
				snprintf(str, sizeof(str), "%.3f cpi", nCoreCycles / nInstructions);
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
				snprintf(str, sizeof(str), "%.2f %%", 100 * nL1Hits / nCachableLoads);
				break;
			case pi_l1phitrate:
				desc = "L1P hit rate";
				snprintf(str, sizeof(str), "%.2f %%", 100 * nL1PHits / nL1PAccesses);
				break;
				//case pi_overhead:
				//desc = "Instr overhead";
				//snprintf(str, sizeof(str), "%.2f %%", 100 * (nInstructions - nNecessaryInstr) / nInstructions);
				//break;
			case pi_hitinl1p:
				desc = "Loads that hit in L1P";
				snprintf(str, sizeof(str), "%f %%" ,  100 * nL1PHits / nCachableLoads);
				break;
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


static void exec_table(benchfunc_t benchmark, bgq_hmflags additional_opts, int j_max, int k_max) {
//static void exec_table(bool sloppiness, hm_func_double hm_double, hm_func_float hm_float, bgq_hmflags additional_opts) {
	benchstat excerpt;

	if (g_proc_id == 0)
		printf("%10s|", "");
	for (int i3 = 0; i3 < lengthof(flags); i3 += 1) {
		if (g_proc_id == 0)
			printf("%-"SCELLWIDTH"s|", flags_desc[i3]);
	}
	if (g_proc_id == 0)
		printf("\n");
	print_repeat("-", 10 + 1 + (CELLWIDTH + 1) * lengthof(flags));
	if (g_proc_id == 0)
		printf("\n");
	for (int i2 = 0; i2 < lengthof(omp_threads); i2 += 1) {
		int threads = omp_threads[i2];

		if (g_proc_id == 0)
			printf("%-10s|", omp_threads_desc[i2]);

		benchstat stats[lengthof(flags)];
		for (int i3 = 0; i3 < lengthof(flags); i3 += 1) {
			bgq_hmflags hmflags = flags[i3];
			hmflags = hmflags | additional_opts;

			benchstat result = runbench(benchmark, hmflags, k_max, j_max, threads);
			stats[i3] = result;

			if (threads == 64 && i3 == 2) {
				excerpt = result;
			}

			char str[80] = { 0 };
			if (result.avgtime == 0)
				snprintf(str, sizeof(str), "~ %s", (result.error > 0.001) ? "X" : "");
			else
				snprintf(str, sizeof(str), "%.2f mlup/s%s", (double) result.lup / (result.avgtime * MEGA), (result.error > 0.001) ? "X" : "");
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

		print_repeat("-", 10 + 1 + (CELLWIDTH + 1) * lengthof(flags));
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

static void benchmark_hopmat(bgq_hmflags flags, int k, int k_max) {
	bgq_HoppingMatrix(false, &g_bgq_spinorfields[k + k_max], &g_bgq_spinorfields[k], flags);
	bgq_HoppingMatrix(true, &g_bgq_spinorfields[k], &g_bgq_spinorfields[k + k_max], flags);
}

static void benchmark_hopmatkernel(bgq_hmflags flags, int k, int k_max) {
	bool floppyprecision = flags & hm_floatprecision;

	bgq_spinorfield_layout layout = bgq_spinorfield_prepareRead(&g_bgq_spinorfields[k], true, true, !floppyprecision, floppyprecision, false);
	bgq_spinorfield_prepareWrite(&g_bgq_spinorfields[k+k_max], false, floppyprecision ? ly_weyl_float : ly_weyl_double);

	bgq_master_sync();
	static bgq_HoppingMatrix_workload work;
	work.isOdd_src = true;
	work.isOdd_dst = false;
	work.targetfield = &g_bgq_spinorfields[k+k_max];
	work.spinorfield = &g_bgq_spinorfields[k];
	work.ic_begin = 0;
	work.ic_end = PHYSICAL_VOLUME;
	work.noprefetchstream = flags & hm_noprefetchstream;
	bgq_HoppingMatrix_work(&work, flags & hm_nokamul, layout);
}

typedef struct {
	size_t k_max;
	bgq_hmflags opts;
	bool doSave;
} checkargs_t;

static int check_hopmat(void *arg_untyped) {
	checkargs_t *args = arg_untyped;
	ucoord k_max = args->k_max;
	ucoord k = 0;
	bgq_hmflags hmflags = args->opts;
	bool doSave = args->doSave;

	bgq_initbgqref();

	bgq_spinorfield_transfer(true, &g_bgq_spinorfields[k], g_spinor_field[k]);
	double compare_transfer = bgq_spinorfield_compare(true, &g_bgq_spinorfields[k], g_spinor_field[k], false);
#ifndef BGQ_COORDCHECK
	assert(compare_transfer == 0);
	// Must be exact copy
#endif

	for (ucoord z = 0; z < LOCAL_LZ ; z += 1) {
		for (ucoord y = 0; y < LOCAL_LY ; y += 1) {
			for (ucoord x = 0; x < LOCAL_LX ; x += 1) {
				for (ucoord t = 0; t < LOCAL_LT ; t += 1) {
					if (bgq_local2isOdd(t, x, y, z) != g_bgq_spinorfields[k].isOdd)
						continue;

					bgq_spinor ref = bgq_legacy_getspinor(g_spinor_field[k], t, x, y, z);
					bgq_spinor bgq = bgq_spinorfield_getspinor(&g_bgq_spinorfields[k], t, x, y, z);

					bgq_setdesc(BGQREF_SOURCE, "BGQREF_SOURCE");
					bgq_setrefvalue(t, x, y, z, BGQREF_SOURCE, ref.v[0].c[0]);
					bgq_setbgqvalue(t, x, y, z, BGQREF_SOURCE, bgq.v[0].c[0]);
				}
			}
		}
	}

	bgq_HoppingMatrix(false, &g_bgq_spinorfields[k + k_max], &g_bgq_spinorfields[k], hmflags);
	HoppingMatrix_switch(false, g_spinor_field[k + k_max], g_spinor_field[k], hmflags);

	for (ucoord z = 0; z < LOCAL_LZ ; z += 1) {
		for (ucoord y = 0; y < LOCAL_LY ; y += 1) {
			for (ucoord x = 0; x < LOCAL_LX ; x += 1) {
				for (ucoord t = 0; t < LOCAL_LT ; t += 1) {
					if (bgq_local2isOdd(t, x, y, z) != g_bgq_spinorfields[k + k_max].isOdd)
						continue;

					bgq_spinor ref = bgq_legacy_getspinor(g_spinor_field[k + k_max], t, x, y, z);
					bgq_spinor bgq = bgq_spinorfield_getspinor(&g_bgq_spinorfields[k + k_max], t, x, y, z);

					bgq_setdesc(BGQREF_RESULT, "BGQREF_RESULT");
					//assert(ref.v[1].c[0]!=-2);
					bgq_setrefvalue(t, x, y, z, BGQREF_RESULT, ref.v[0].c[0]);
					//assert(bgq.v[1].c[0]!=-2);
					bgq_setbgqvalue(t, x, y, z, BGQREF_RESULT, bgq.v[0].c[0]);
				}
			}
		}
	}

	if (doSave)
		bgq_savebgqref();
	double compare_even = bgq_spinorfield_compare(false, &g_bgq_spinorfields[k + k_max], g_spinor_field[k + k_max], false);
#ifndef BGQ_COORDCHECK
	assert(compare_even < 0.01);
#endif


	bgq_spinorfield_transfer(false, &g_bgq_spinorfields[k+k_max], g_spinor_field[k+k_max]);
	compare_transfer = bgq_spinorfield_compare(false, &g_bgq_spinorfields[k+k_max], g_spinor_field[k+k_max], false);
#ifndef BGQ_COORDCHECK
	assert(compare_transfer == 0);
#endif

	bgq_HoppingMatrix(true, &g_bgq_spinorfields[k], &g_bgq_spinorfields[k+k_max], hmflags);
	HoppingMatrix_switch(true, g_spinor_field[k], g_spinor_field[k+k_max], hmflags);

	double compare_odd = bgq_spinorfield_compare(true, &g_bgq_spinorfields[k], g_spinor_field[k], false);
#ifndef BGQ_COORDCHECK
	assert(compare_odd < 0.01);
#endif


	master_print("Comparison to reference version: even=%e odd=%e max difference\n", compare_even, compare_odd);
	return EXIT_SUCCESS;
}


static void exec_bench(int j_max, int k_max) {
	if (g_proc_id==0)
		bgq_qpx_unittest();

	bgq_indices_init();
	bgq_comm_mpi_init();
	bgq_comm_spi_init();
	bgq_initbgqref();

	bgq_spinorfields_init(2 * k_max + 1, 0);

	bgq_gaugefield_init();
	bgq_gaugefield_transferfrom(g_gauge_field);

	bool done = false;
#pragma omp parallel for schedule(static) ordered firstprivate(done)
	for (int j = 0; j < omp_get_max_threads(); j += 1)
	        {
		if (g_proc_id == 0 && !done) {
#pragma omp ordered
			{
				//printf("MK Here is omp_id %2d, cpu_id %2d running on SMT-Thread %d of Core %2d\n", omp_get_thread_num(), Kernel_ProcessorID(), Kernel_ProcessorThreadID(), Kernel_ProcessorCoreID());
				done = true;
			}
		}
	}

	uint64_t ws = 0;
	uint64_t indices = 0;
	uint64_t forcomm = 0;
	ws += PHYSICAL_VOLUME * sizeof(bgq_weylsite_double); // input spinorfield
	indices += PHYSICAL_VOLUME * sizeof(bgq_weylsite_double*); // sendptr
	indices += (2*LOCAL_HALO_T/PHYSICAL_LP+2*PHYSICAL_HALO_X+2*PHYSICAL_HALO_Y+2*PHYSICAL_HALO_Z)*sizeof(bgq_weylsite_double*); // consptr
	ws += PHYSICAL_VOLUME * sizeof(bgq_gaugesite); // gauge field
	forcomm += bgq_weyl_section_offset(sec_comm_end) - bgq_weyl_section_offset(sec_comm); // for communication

	ws += forcomm;
	ws += indices;
	uint64_t ws_write = ws + PHYSICAL_VOLUME * sizeof(bgq_weylsite_double); // target spinor
	master_print("Working set size: %.1fMB (%.1fMB incl target) (%.1fMB index, %.1fMB commbuf)\n", (double)ws/(1024.0*1024.0), (double)ws_write/(1024.0*1024.0),indices/MEBI,forcomm/MEBI);

	master_print("VOLUME=%d PHYSICAL_VOLUME=%zu PHYSICAL_BODY=%zu PHYSICAL_SURFACE=%zu\n", VOLUME, PHYSICAL_VOLUME, PHYSICAL_BODY, PHYSICAL_SURFACE);

	for (int k = 0; k < k_max; k += 1) {
		bgq_spinorfield_setup(&g_bgq_spinorfields[k],       true,  false, true, false, true, false);
		bgq_spinorfield_setup(&g_bgq_spinorfields[k+k_max], false, false, true, false, true, false);
	}

	for (int k = 0; k < k_max; k += 1) {
		/*initialize the pseudo-fermion fields*/
		random_spinor_field(g_spinor_field[k], VOLUME / 2, 0);
		bgq_spinorfield_transfer(true, &g_bgq_spinorfields[k], g_spinor_field[k]);
	}

	master_print("Double: ");
	checkargs_t checkargs_double = {
	        .k_max = k_max,
	        .opts = 0,
	        .doSave = false
	};
	//bgq_parallel(&check_hopmat, &checkargs_double);

	master_print("Float: ");
	checkargs_t checkargs_float = {
	        .k_max = k_max,
	        .opts = hm_floatprecision,
	        .doSave = true
	};
	bgq_parallel(&check_hopmat, &checkargs_float);

	master_print("Benchmark: hopmatkernel\n");
	exec_table(&benchmark_hopmatkernel, 0, j_max, k_max);
	print_repeat("\n", 2);

	master_print("Benchmark: HoppingMatrix\n");
	exec_table(&benchmark_hopmat, 0, j_max, k_max);
	print_repeat("\n", 2);

	master_print("Benchmarking done\n");
}


int main(int argc, char *argv[]) {
	int j, j_max = 5, k, k_max = 1;
#ifdef HAVE_LIBLEMON
	paramsXlfInfo *xlfInfo;
#endif
	int status = 0;

	static double t1, t2, dt, sdt, dts, qdt, sqdt;
	double antioptaway = 0.0;

#ifdef MPI
	static double dt2;

	DUM_DERI = 6;
	DUM_SOLVER = DUM_DERI+2;
	DUM_MATRIX = DUM_SOLVER+6;
	NO_OF_SPINORFIELDS = DUM_MATRIX+2;


//#  ifdef OMP
	int mpi_thread_provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_thread_provided);
//#  else
//  MPI_Init(&argc, &argv);
//#  endif
	MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);

#else
	g_proc_id = 0;
#endif

	g_rgi_C1 = 1.;
	bool legacy_spi = true;
	char *input_filename = "benchmark.input";
	int c = 0;
	while ((c = getopt(argc, argv, "vh?f:")) != -1) {
		switch (c) {
		case 'f':
			input_filename = malloc(strlen(optarg) + 1);
			strcpy(input_filename, optarg);
			break;
		case 'v':
			verbose = 1;
			break;
		case 'h':
			legacy_spi = false;
			break;
		case '?':
			default:
			//usage();
			exit(0);
		}
	}

	/* Read the input file */
	if ((status = read_input(input_filename)) != 0) {
		fprintf(stderr, "Could not find input file: %s\nAborting...\n", input_filename);
		exit(-1);
	}

#ifdef OMP
	if(omp_num_threads > 0)
	{
		omp_set_num_threads(omp_num_threads);
	}
	else {
		if( g_proc_id == 0 )
		printf("# No value provided for OmpNumThreads, using default (OMP_NUM_THREADS=%d)\n", omp_get_max_threads());

		//omp_num_threads = 1;
		//omp_set_num_threads(omp_num_threads);
	}

	init_omp_accumulators(omp_num_threads);
#endif

	tmlqcd_mpi_init(argc, argv);

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
	j = init_moment_field(VOLUME, VOLUMEPLUSRAND + g_dbw2rand);
	if (j != 0) {
		fprintf(stderr, "Not enough memory for moment fields! Aborting...\n");
		exit(0);
	}

	if (g_proc_id == 0) {
		fprintf(stdout, "# The number of processes is %d \n", g_nproc);
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
	if (legacy_spi) {
		master_print("MK calling init_dirac_halfspinor()\n");
	j = init_dirac_halfspinor();
	master_print("MK exit init_dirac_halfspinor()\n");
	if ( j!= 0) {
		fprintf(stderr, "Not enough memory for halfspinor fields! Aborting...\n");
		exit(0);
	}}
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
#if 0
	for (int t = 0; t <= 64; t+=1) {
		omp_set_num_threads(t);
#pragma omp parallel
		{
			size_t tid = omp_get_thread_num();
			for (int k = 0; k < 100; k+=1)
			if (tid==0) {
				ompbar();
			} else {
				ompbar2();
			}
		}
	}
#endif

	assert(even_odd_flag);
	exec_bench(j_max, k_max);
	return 0;
	/* END MK */

	while (sdt < 30.) {
#ifdef MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		t1 = gettime();
		antioptaway = 0.0;
		for (j = 0; j < j_max; j++) {
			for (k = 0; k < k_max; k++) {
				Hopping_Matrix(0, g_spinor_field[k + k_max], g_spinor_field[k]);
				Hopping_Matrix(1, g_spinor_field[k], g_spinor_field[k + k_max]);
				antioptaway += creal(g_spinor_field[2 * k_max][0].s0.c0);
			}
		}
		t2 = gettime();
		dt = t2 - t1;
#ifdef MPI
		MPI_Allreduce (&dt, &sdt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
		sdt = dt;
#endif
		qdt = dt * dt;
#ifdef MPI
		MPI_Allreduce (&qdt, &sqdt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
		printf("# The following result is just to make sure that the calculation is not optimized away: %e\n", antioptaway);
		printf("# Total compute time %e sec, variance of the time %e sec. (%d iterations).\n", sdt, sqdt, j_max);
		printf("# Communication switched on:\n# (%d Mflops [%d bit arithmetic])\n", (int) (1608.0f / sdt), (int) sizeof(spinor) / 3);
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
	return (0);
}
