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

//static L2_Barrier_t barrier;


void bgq_dispatch_HoppingMatrix_phase3() {

}



void bgq_dispatch_worklerloop() {
	while (true) {
		//WU_ArmWithAddress
	}
}

