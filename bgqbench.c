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
#include "bgq/bgq_gaugefield.h"
#include "bgq/bgq_spinorfield.h"
#include "bgq/bgq_HoppingMatrix.h"
#include "bgq/bgq_qpx.h"
#include "bgq/bgq_utils.h"
#include "bgq/bgq_dispatch.h"
#include <omp.h>
#include "bgq/mypapi.h"

#ifdef XLC
#include <l1p/pprefetch.h>
#include <l1p/sprefetch.h>
#endif

typedef bgq_spinorfield bgq_spinorfield_double;
typedef bgq_spinorfield bgq_spinorfield_float;

typedef bgq_gaugefield bgq_gaugefield_double;
typedef bgq_gaugefield bgq_gaugefield_float;

typedef bgq_gaugesite bgq_gaugesite_double;
typedef bgq_gaugesite bgq_gaugesite_float;

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

	double lups;
	double flops;

	double error;

	mypapi_counters counters;
	bgq_hmflags opts;
} benchstat;






typedef struct {
	int j_max;
	int k_max;
	benchfunc_t benchfunc;
	bgq_hmflags opts;
	benchstat result;
} master_args;

static void benchmark_setup_worker(void *argptr, size_t tid, size_t threads) {
#ifdef BGQ
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
	bgq_spinorfield_setup(&g_bgq_spinorfields[k], true, false, false, false, true);

	bgq_spinorfield_setup(&g_bgq_spinorfields[k + k_max], false, false, false, false, true);
#ifndef BGQ_COORDCHECK
	memset(g_bgq_spinorfields[k].sec_weyl, 0, bgq_weyl_section_offset(sec_end));
	memset(g_bgq_spinorfields[k + k_max].sec_weyl, 0, bgq_weyl_section_offset(sec_end));
#endif


	bgq_spinorfield_transfer(true, &g_bgq_spinorfields[k], g_spinor_field[k]);
#ifndef BGQ_COORDCHECK
	double compare_transfer = bgq_spinorfield_compare(false, &g_bgq_spinorfields[k], g_spinor_field[k], false);
	assert(compare_transfer == 0); // Must be exact copy
#endif

	bgq_HoppingMatrix(false, &g_bgq_spinorfields[k + k_max], &g_bgq_spinorfields[k], hmflags);
	HoppingMatrix_switch(false, g_spinor_field[k + k_max], g_spinor_field[k], hmflags);
#ifndef BGQ_COORDCHECK
	double compare_even = bgq_spinorfield_compare(false, &g_bgq_spinorfields[k + k_max], g_spinor_field[k + k_max], true);
	assert(compare_even < 0.01);
#endif


#ifndef BGQ_COORDCHECK
	bgq_spinorfield_transfer(false, &g_bgq_spinorfields[k + k_max], g_spinor_field[k + k_max]); 	// We don't want to accumulate errors
	compare_transfer = bgq_spinorfield_compare(false, &g_bgq_spinorfields[k + k_max], g_spinor_field[k + k_max], false);
	assert(compare_transfer == 0); // Must be exact copy
#endif

	bgq_HoppingMatrix(true, &g_bgq_spinorfields[k+2*k_max], &g_bgq_spinorfields[k + k_max], hmflags);
	HoppingMatrix_switch(true, g_spinor_field[k+2*k_max], g_spinor_field[k+k_max], hmflags);
#ifndef BGQ_COORDCHECK
	double compare_odd = bgq_spinorfield_compare(true, &g_bgq_spinorfields[k + 2*k_max], g_spinor_field[k + 2*k_max], true);
	assert(compare_odd < 0.01);
#endif


#ifndef BGQ_COORDCHECK
	return max(compare_even, compare_odd);
#else
	return 0;
#endif
}


static int benchmark_master(void *argptr) {
	master_args * const args = argptr;
	const int j_max = args->j_max;
	const int k_max = args->k_max;
	const benchfunc_t benchfunc = args->benchfunc;
	const bgq_hmflags opts = args->opts;
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
		//const bgq_hmflags implicitprefetch = opts & (hm_prefetchimplicitdisable | hm_prefetchimplicitoptimistic | hm_prefetchimplicitconfirmed);
bool withcheck = opts & hm_withcheck;

		// Setup thread options (prefetch setting etc.)
		bgq_master_call(&benchmark_setup_worker, argptr);


		double err = 0;
		if (withcheck) {
			err = runcheck(opts, k_max);
		}

		double localsumtime = 0;
		double localsumsqtime = 0;
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
			bool isJMax =(0 <= j) && (j < j_max);
			//bool isFirstJMax = (j==0);
			//bool isLastJMax = (j==j_max-1);

			double start_time;
			if (isJMax) {
				start_time = MPI_Wtime();
			}
			if (isPapi) {
				mypapi_start(papiSet);
			}

			{
				// The main benchmark
				for (int k = 0; k < k_max; k += 1) {
					benchfunc(opts, k, k_max);
				}
				bgq_master_sync(); // Wait for all threads to finish, to get worst thread timing
			}

			if (isPapi) {
				mypapi_counters curcounters = mypapi_stop();
				counters = mypapi_merge_counters(&counters, &curcounters);
			}
			if (isJMax) {
				double end_time = MPI_Wtime();
				double duration = end_time - start_time;
				localsumtime += duration;
				localsumsqtime += sqr(duration);
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
		int lups_body= PHYSICAL_BODY*PHYSICAL_LP*PHYSICAL_LK;
		int lups_surface=PHYSICAL_SURFACE*PHYSICAL_LP*PHYSICAL_LK;
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


	benchstat *result = &args->result;
	result->avgtime = avgtime;
	result->localrmstime = avglocalrms;
	result->globalrmstime = rmstime;
	result->lups = lups;
	result->flops = flops / avgtime;
	result->error = err;
	result->counters = counters;
	result->opts = opts;
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

static bool sloppinessess[] = { false, true };
static char *sloppinessess_desc[] = { "double", "float" };

static int omp_threads[] = { 1, 2, 4, 8, 16, 32, 33, 48, 56, 64 };
static char *omp_threads_desc[] = { "1", "2", "4", "8", "16", "32", "33", "48", "56", "64" };


static bgq_hmflags flags[] = {
		hm_withcheck
#if 0
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
#endif
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



#define CELLWIDTH 15
#define SCELLWIDTH TOSTRING(CELLWIDTH)

static void print_stats(benchstat *stats, int bodySites, int surfaceSites, int haloSites) {
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
			print_stats(stats, 0, 0, 0);
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
	bgq_HoppingMatrix(false, &g_bgq_spinorfields[k+k_max], &g_bgq_spinorfields[k], flags);
	bgq_HoppingMatrix(true, &g_bgq_spinorfields[k], &g_bgq_spinorfields[k+k_max], flags);
}


static void exec_bench(int j_max, int k_max) {
	bgq_indices_init();
	bgq_spinorfields_init(2*k_max+1, 0);

	bgq_gaugefield_init();
	bgq_gaugefield_transferfrom(g_gauge_field);

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
		bgq_spinorfield_transfer(true, &g_bgq_spinorfields[k], g_spinor_field[k]);
	}

	exec_table(&benchmark_hopmat , 0, j_max, k_max);
}


static void exec_bench_all() {
	print_repeat("\n", 2);

	for (int i0 = 0; i0 < lengthof(kamuls); i0 += 1) {
		//bool kamul = kamuls[i0];
		if (g_proc_id == 0) printf("Benchmark: %s\n", kamuls_desc[i0]);

		for (int i1 = 0; i1 < lengthof(sloppinessess); i1 += 1) {
			//bool sloppiness = sloppinessess[i1];

			if (g_proc_id == 0) printf("Benchmarking precision: %s\n", sloppinessess_desc[i1]);
			//exec_table(sloppiness, bgq_HoppingMatrix_double, bgq_HoppingMatrix_float, !kamul * hm_nokamul);
		}

		if (g_proc_id == 0) printf("\n");
	}

	if (g_proc_id == 0) printf("Benchmark done\n");
}

int main(int argc,char *argv[])
{
  int j,j_max=5,k,k_max = 1;
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



    while(sdt < 30.) {
#ifdef MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      t1 = gettime();
      antioptaway=0.0;
      for (j=0;j<j_max;j++) {
        for (k=0;k<k_max;k++) {
          Hopping_Matrix(0, g_spinor_field[k+k_max], g_spinor_field[k]);
          Hopping_Matrix(1, g_spinor_field[k], g_spinor_field[k+k_max]);
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
