/*
 * bgq_dispatch.c
 *
 *  Created on: Oct 15, 2012
 *      Author: meinersbur
 */

#define BGQ_DISPATCH_C_
#include "bgq_dispatch.h"

#include "bgq_qpx.h"

#ifdef BGQ
#include <l2/barrier.h>
#include <wu/wait.h>
#include <upci/upc_atomic.h>
#include <hwi/include/bqc/A2inlines.h>
#include <time.h> // nanosleep() system call
#endif

#include <omp.h>
#include <assert.h>
#include <stdbool.h>
#include <stdint.h>

//static L2_Barrier_t barrier;




static volatile bgq_worker_func g_bgq_dispatch_func;
static void * volatile g_bgq_dispatch_arg;
static volatile bool g_bgq_dispatch_terminate;
static volatile bool g_bgq_dispatch_sync;


int bgq_parallel(bgq_master_func master_func, void *master_arg) {
	assert(!omp_in_parallel() && "This starts the parallel section, do not call it within one");
	g_bgq_dispatch_func = NULL;
	g_bgq_dispatch_arg = NULL;
	g_bgq_dispatch_terminate = false;
	g_bgq_dispatch_sync = false;
	int master_result = 0;
	// We use OpenMP only to start the threads
	// Overhead of using OpenMP is too large
#pragma omp parallel
	{
		size_t tid = omp_get_thread_num();

		// Start workers, master thread will return
		bgq_worker();

		// Start control program in master
		if (tid == 0) {
			master_result = master_func(master_arg);

			// After program finishes, set flag so worker threads can terminate
			g_bgq_dispatch_func = NULL;
			g_bgq_dispatch_arg = NULL;
			g_bgq_dispatch_terminate = true;
			g_bgq_dispatch_sync = false;
			mbar();
		}

		// Wakeup workers to terminate
		bgq_worker();
	}
	return master_result;
}


void bgq_worker() {
	//assert(omp_in_parallel() && "Call this inside #pragma omp parallel");
	size_t threads = omp_get_num_threads();
	size_t tid = omp_get_thread_num();
	//assert((tid != 0) && "This function is for non-master threads only");


	if (tid == 0) {
		if (g_bgq_dispatch_func) {
			// There is some work to do, join the other workers
		} else {
			// Return so the master thread can execute the program
			return;
		}
	} else {
		// All others wait for commands from the master thread
	}


	while (true) {
		// Wait until every thread did its work
		// This doesn't need to be a barrier, waiting for submission of some work from the master is ok too
		// TODO: Hope OpenMP has a good implementation without busy waiting; if not, do some own work

#ifdef BGQ
		uint64_t ppc32 = mfspr(SPRN_PPR32);
		ThreadPriority_Low(); // Lower thread priority, so if busy waiting is used, do not impact other threads on core
#endif
#pragma omp barrier
#ifdef BGQ
		mtspr(SPRN_PPR32, ppc32); // Restore original priority
#endif
		// All threads should be here at the same time, including the master thread, which has issued some work, namely, calling a function

		if (g_bgq_dispatch_terminate) {
			// Exit program, or at least, the parallel section
			return;
		}

		if (g_bgq_dispatch_sync) {
			// This was meant just for synchronization between the threads, which already has been done
		} else {
			assert(g_bgq_dispatch_func);
			void *arg = g_bgq_dispatch_arg;
			g_bgq_dispatch_func(arg, tid, threads); //TODO: Shuffle tid to load-balance work?
		}

		if (tid==0) {
			// Let master thread continue the program
			break;
		}
		// All others, wait for the next command
	}
}


void bgq_master_call(bgq_worker_func func, void *arg) {
	assert(func);
	g_bgq_dispatch_func = func;
	g_bgq_dispatch_arg = arg;
	g_bgq_dispatch_terminate = false;
	g_bgq_dispatch_sync = false;
	mbar();

	bgq_worker();

	g_bgq_dispatch_func = NULL;
	g_bgq_dispatch_arg = NULL;
	g_bgq_dispatch_terminate = false;
}


void bgq_master_sync() {
	g_bgq_dispatch_func = NULL;
	g_bgq_dispatch_arg = NULL;
	g_bgq_dispatch_sync = true;
	g_bgq_dispatch_terminate = false;
	mbar();

	bgq_worker();

	g_bgq_dispatch_sync = false;
}


