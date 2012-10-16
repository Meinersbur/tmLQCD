/*
 * bgq_dispatch.c
 *
 *  Created on: Oct 15, 2012
 *      Author: meinersbur
 */

#define BGQ_DISPATCH_C_
#include "bgq_dispatch.h"

#include <l2/barrier.h>
#include <wu/wait.h>
#include <upci/upc_atomic.h>
#include <hwi/include/bqc/A2inlines.h>
#include <time.h> // nanosleep() system call
#include <omp.h>

//static L2_Barrier_t barrier;




volatile bgq_dispatch_func g_bgq_dispatch_func;
volatile *void g_bgq_dispatch_arg;

void bgq_worker() {
	assert(omp_in_parallel() && "Call this inside #pragma omp parallel");
	size_t threads = omp_get_num_threads();
	size_t tid = omp_get_thread_num();
	//assert((tid != 0) && "This function is for non-master threads only");

#pragma omp master
{
	// Let the master thread do its work
	return;
}
// All others wait for commands from the master thread

	while (true) {

		// Wait until every thread did its work
		// This doesn't need to be a barrier, waiting for submission of some work from the master is ok too
		// TODO: Hope OpenMP has a good implementation without busy waiting; if not, do some own work
		uint64_t ppc32 = mfspr(SPRN_PPR32);
		ThreadPriority_Low(); // Lower thread priority, so if busy waiting is used, do not impact other threads on core
#pragma omp barrier (workerloop)
		mtspr(SPRN_PPR32, ppc32); // Restore original priority

		// All threads should be here at the same time, including the master thread, which has issued some work, namely, calling a function

	}
}

void bgq_master_call(bgq_dispatch_func func, void *arg) {

}

void bgq_master_nomorework() {

}


