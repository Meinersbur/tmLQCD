/*
 * bgq_dispatch.c
 *
 *  Created on: Oct 15, 2012
 *      Author: meinersbur
 */

#define BGQ_DISPATCH_C_
#include "bgq_dispatch.h"

#include "bgq_qpx.h"

#if BGQ_QPX
#include <l2/barrier.h>
//#include <wu/wait.h>
#include <upci/upc_atomic.h>
//#include <hwi/include/bqc/A2_inlines.h>
//#include <time.h> // nanosleep() system call
#endif

#include "../global.h"

#include <omp.h>
#include <assert.h>
#include <stdbool.h>
#include <stdint.h>

//static L2_Barrier_t barrier;


static char space[64*25+1]= {' '};


static volatile bgq_worker_func g_bgq_dispatch_func;
static void * volatile g_bgq_dispatch_arg;
static volatile bool g_bgq_dispatch_terminate;
static volatile bool g_bgq_dispatch_sync; // obsolete
//static volatile size_t g_bgq_dispatch_seq; // obsolete

static bool g_bgq_dispatch_pendingsync; // Accessed by master only


#if BGQ_QPX
static bool g_bgq_dispatch_barrier_initialized = false;
static L2_Barrier_t g_bgq_dispatch_barrier = L2_BARRIER_INITIALIZER;
#endif

static inline void bgq_thread_barrier() {
#if BGQ_QPX
	uint64_t savpri = Set_ThreadPriority_Low(); // Lower thread priority, so if busy waiting is used, do not impact other threads on core
	L2_Barrier(&g_bgq_dispatch_barrier, g_bgq_dispatch_threads);
	Restore_ThreadPriority(savpri);
#else
#pragma omp barrier
#endif
}


int bgq_parallel(bgq_master_func master_func, void *master_arg) {
	for (int i = 0; i < 64*25; i+=1)
		space[64*25] = ' ';
	space[64*25] = '\0';
	assert(!omp_in_parallel() && "This starts the parallel section, do not call it within one");
	g_bgq_dispatch_func = NULL;
	g_bgq_dispatch_arg = NULL;
	g_bgq_dispatch_terminate = false;
	g_bgq_dispatch_sync = false;
	//g_bgq_dispatch_seq = 0;
#if BGQ_QPX
	if (!g_bgq_dispatch_barrier_initialized) {
		Kernel_L2AtomicsAllocate(&g_bgq_dispatch_barrier, sizeof(g_bgq_dispatch_barrier));
		g_bgq_dispatch_barrier_initialized = true;
	}
#endif
	g_bgq_dispatch_threads = omp_get_max_threads();
	omp_num_threads = 1/*omp_get_num_threads()*/; // For legacy linalg (it depends on whether nested parallelism is enabled)

	int master_result = 0;
	// We use OpenMP only to start the threads
	// Overhead of using OpenMP is too large
#pragma omp parallel
	{
		size_t tid = omp_get_thread_num(); // Or

		// Start workers
		if (tid != 0) {
			g_bgq_dispatch_pendingsync = true;
			bgq_worker();
		}

		// Start control program in master
		if (tid == 0) {
			master_result = master_func(master_arg);

			bgq_master_sync();
			// After program finishes, set flag so worker threads can terminate
			g_bgq_dispatch_func = NULL;
			g_bgq_dispatch_arg = NULL;
			g_bgq_dispatch_terminate = true;
			g_bgq_dispatch_sync = false;
#if BGQ_QPX
			mbar();
#else
#pragma omp flush
#endif

			// Wakeup workers to terminate
			bgq_worker();
		}

		//printf("%*sEXIT: tid=%u\n", (int)tid*25, "",(int)tid);
	}
	g_bgq_dispatch_threads = 0;
	omp_num_threads = omp_get_max_threads();
	return master_result;
}


//static size_t count = 0;
//#pragma omp threadvar(count)

void bgq_worker() {
	//assert(omp_in_parallel() && "Call this inside #pragma omp parallel");
	assert(g_bgq_dispatch_threads == omp_get_num_threads());
	size_t threads = g_bgq_dispatch_threads;
	size_t tid = omp_get_thread_num(); // Or Kernel_ProcessorID()


	//assert((tid != 0) && "This function is for non-master threads only");
//size_t count = 0;
	while (true) {
		// Wait until every thread did its work
		// This doesn't need to be a barrier, waiting for submission of some work from the master is ok too
		// TODO: Hope OpenMP has a good implementation without busy waiting; if not, do some own work

		if (tid!=0) {
			// Guarantee that work is finished
			//TODO: can we implement this without this second barrier?
			// Required to ensure consistency of g_bgq_dispatch_sync, g_bgq_dispatch_terminate, g_bgq_dispatch_func
			bgq_thread_barrier(); // the sync barrier
		}


		// Master thread may write shared variables before this barrier
		bgq_thread_barrier();
#if BGQ_QPX
		mbar();
#else
#pragma omp flush
#endif
		// Worker threads read common variables after this barrier
		if (tid==0) {
			assert(!g_bgq_dispatch_pendingsync);
			g_bgq_dispatch_pendingsync = true;
		}

		//count += 1;
		// All threads should be here at the same time, including the master thread, which has issued some work, namely, calling a function

		if (g_bgq_dispatch_sync) {
			// This was meant just for synchronization between the threads, which already has been done
			//printf("%*sSYNC: tid=%u seq=%u\n", (int)tid*20, "", (int)tid, (int)g_bgq_dispatch_seq);
		} else if (g_bgq_dispatch_terminate) {
			// Exit program, or at least, the parallel section
			//printf("%*sTERM: tid=%u seq=%u\n", (int)tid*20, "",(int)tid, (int)g_bgq_dispatch_seq);
			return;
		} else {
			//printf("%*sCALL: tid=%u seq=%u\n", (int)tid*20, "",(int)tid, (int)g_bgq_dispatch_seq);
			assert(g_bgq_dispatch_func);
			void *arg = g_bgq_dispatch_arg;
			g_bgq_dispatch_func(arg, tid, threads); //TODO: Shuffle tid to load-balance work?
		}

		if (tid==0) {
			// Let master thread continue the program
			// Hint: master thread must call bgq_thread_barrier() sometime to release the workers from the following barrier
			return;
		}
		// All others, wait for the next command
	}
}


void bgq_master_call(bgq_worker_func func, void *arg) {
	assert(omp_get_thread_num()==0);
	assert(func);

	bgq_master_sync();

	//g_bgq_dispatch_seq += 1;
	//printf("MASTER CALL seq=%d--------------------------------------------------------\n", (int)g_bgq_dispatch_seq);
	g_bgq_dispatch_func = func;
	g_bgq_dispatch_arg = arg;
	g_bgq_dispatch_terminate = false;
	g_bgq_dispatch_sync = false;
#if BGQ_QPX
	mbar();
#else
#pragma omp flush
#endif

	bgq_worker();
}


void bgq_master_sync() {
	assert(omp_get_thread_num()==0);

	if (g_bgq_dispatch_pendingsync) {
		bgq_thread_barrier();
		//fflush(stdout);
		//printf("MASTER SYNC seq=%d--------------------------------------------------------\n", (int)g_bgq_dispatch_seq);
		//fflush(stdout);
		g_bgq_dispatch_pendingsync = false;
		return;
	} else {
		// Threads already at sync barrier
		return;
	}
}


