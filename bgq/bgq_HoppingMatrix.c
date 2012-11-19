/*
 * bgq_HoppingMatrix.c
 *
 *  Created on: Oct 13, 2012
 *      Author: meinersbur
 */

#define BGQ_HOPPING_MATRIX_C_
#include "bgq_HoppingMatrix.h"

#include "bgq_field.h"
#include "bgq_spinorfield.h"
#include "bgq_qpx.h"
#include "bgq_dispatch.h"
#include "bgq_comm.h"
#include "bgq_workers.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define PRECISION double

static void bgq_HoppingMatrix_nokamul_worker_readFulllayout_double(void *arg, size_t tid, size_t threads) {
	//bgq_HoppingMatrix_worker(arg,tid,threads,false,true);

	bool const kamul = false;
	bool const readFulllayout = true;

#define PRECISION double
#define BGQ_HOPPINGMATRIXWORKER_INC_ 1
#include "bgq_HoppingMatrixWorker.inc.c"
#undef PRECISION
}
static void bgq_HoppingMatrix_nokamul_worker_readFulllayout_float(void *arg, size_t tid, size_t threads) {
	//bgq_HoppingMatrix_worker(arg,tid,threads,false,true);

	bool const kamul = false;
	bool const readFulllayout = true;

#define PRECISION float
#define BGQ_HOPPINGMATRIXWORKER_INC_ 1
#include "bgq_HoppingMatrixWorker.inc.c"
#undef PRECISION
}
static void bgq_HoppingMatrix_kamul_worker_readFulllayout_double(void *arg, size_t tid, size_t threads) {
	//bgq_HoppingMatrix_worker(arg,tid,threads,true,true);

	bool const kamul = true;
	bool const readFulllayout = true;

#define PRECISION double
#define BGQ_HOPPINGMATRIXWORKER_INC_ 1
#include "bgq_HoppingMatrixWorker.inc.c"
#undef PRECISION
}
static void bgq_HoppingMatrix_kamul_worker_readFulllayout_float(void *arg, size_t tid, size_t threads) {
	//bgq_HoppingMatrix_worker(arg,tid,threads,true,true);

	bool const kamul = true;
	bool const readFulllayout = true;

#define PRECISION float
#define BGQ_HOPPINGMATRIXWORKER_INC_ 1
#include "bgq_HoppingMatrixWorker.inc.c"
#undef PRECISION
}
static void bgq_HoppingMatrix_nokamul_worker_readWeyllayout_double(void *arg, size_t tid, size_t threads) {
	//bgq_HoppingMatrix_worker(arg,tid,threads,false,false);

	bool const kamul = false;
	bool const readFulllayout = false;

#define PRECISION double
#define BGQ_HOPPINGMATRIXWORKER_INC_ 1
#include "bgq_HoppingMatrixWorker.inc.c"
#undef PRECISION
}
static void bgq_HoppingMatrix_nokamul_worker_readWeyllayout_float(void *arg, size_t tid, size_t threads) {
	//bgq_HoppingMatrix_worker(arg,tid,threads,false,false);

	bool const kamul = false;
	bool const readFulllayout = false;

#define PRECISION float
#define BGQ_HOPPINGMATRIXWORKER_INC_ 1
#include "bgq_HoppingMatrixWorker.inc.c"
#undef PRECISION
}
static void bgq_HoppingMatrix_kamul_worker_readWeyllayout_double(void *arg, size_t tid, size_t threads) {
	//bgq_HoppingMatrix_worker(arg,tid,threads,true,false);

	bool const kamul = true;
	bool const readFulllayout = false;

#define PRECISION double
#define BGQ_HOPPINGMATRIXWORKER_INC_ 1
#include "bgq_HoppingMatrixWorker.inc.c"
#undef PRECISION
}
static void bgq_HoppingMatrix_kamul_worker_readWeyllayout_float(void *arg, size_t tid, size_t threads) {
	//bgq_HoppingMatrix_worker(arg,tid,threads,true,false);

	bool const kamul = true;
	bool const readFulllayout = false;

#define PRECISION float
#define BGQ_HOPPINGMATRIXWORKER_INC_ 1
#include "bgq_HoppingMatrixWorker.inc.c"
#undef PRECISION
}


void bgq_HoppingMatrix_work(bgq_HoppingMatrix_workload *work,  bool nokamul, bgq_spinorfield_layout layout) {
	assert(work);
	size_t sites = work->ic_end - work->ic_begin;
	uint64_t old = flopaccumulator;

	bgq_worker_func func;
	if (nokamul) {
		switch (layout) {
		case ly_full_double:
			func = &bgq_HoppingMatrix_nokamul_worker_readFulllayout_double;
			break;
		case ly_full_float:
			func = &bgq_HoppingMatrix_nokamul_worker_readFulllayout_float;
			break;
		case ly_weyl_double:
			func = &bgq_HoppingMatrix_nokamul_worker_readWeyllayout_double;
			break;
		case ly_weyl_float:
			func = &bgq_HoppingMatrix_nokamul_worker_readWeyllayout_float;
			break;
		default:
			UNREACHABLE
		}
	} else {
		switch (layout) {
		case ly_full_double:
			func = &bgq_HoppingMatrix_kamul_worker_readFulllayout_double;
			break;
		case ly_full_float:
			func = &bgq_HoppingMatrix_kamul_worker_readFulllayout_float;
			break;
		case ly_weyl_double:
			func = &bgq_HoppingMatrix_kamul_worker_readWeyllayout_double;
			break;
		case ly_weyl_float:
			func = &bgq_HoppingMatrix_kamul_worker_readWeyllayout_float;
			break;
		default:
			UNREACHABLE
		}
	}
	bgq_master_call(func, work);
	uint64_t flopPerSite = (
			/*accum spinor*/	7/*dirs*/ * (4 * 3)/*cmplx per spinor*/ * 2/*flops accum*/ +
			/*su3 mul*/			8/*dirs*/ * 2/*su3vec per weyl*/ * (6*9 + 2*3)/*flop per su3 mv-mul*/
	);
	if (!nokamul)
		flopPerSite += /*kamul*/			8/*dirs*/ * (2 * 3)/*cmplx per weyl*/ * 6/*flops cmplx mul*/;
	if (layout & ly_weyl)
		flopPerSite += /*weyl reduce*/     8/*dirs*/ * 2/*weyl per dir*/ * (2 * 3)/*cmplx per weyl*/ * 2/*flops*/;
	flopaccumulator += sites * PHYSICAL_LK * flopPerSite;
#if 0
	if (readFulllayout) {
		if (nokamul) {
			if (inputfield->isFulllayoutSloppy)
				bgq_master_call(&bgq_HoppingMatrix_nokamul_worker_readFulllayout_float, work);
			 else
				 bgq_master_call(&bgq_HoppingMatrix_nokamul_worker_readFulllayout_double, work);
			flopaccumulator += sites * PHYSICAL_LK * (
				/*weyl reduce*/     8/*dirs*/ * 2/*weyl per dir*/ * (2 * 3)/*cmplx per weyl*/ * 2/*flops*/ +
				/*su3 mul*/			8/*dirs*/ * 2/*su3vec per weyl*/ * (6*9 + 2*3)/*flop per su3 mv-mul*/
			);

		} else {
			if (inputfield->isFulllayoutSloppy)
				bgq_master_call(&bgq_HoppingMatrix_kamul_worker_readFulllayout_float, work);
			else
				bgq_master_call(&bgq_HoppingMatrix_kamul_worker_readFulllayout_double, work);
			flopaccumulator += sites * PHYSICAL_LK * (
				/*weyl reduce*/     8/*dirs*/ * 2/*weyl per dir*/ * (2 * 3)/*cmplx per weyl*/ * 2/*flops*/ +
				/*su3 mul*/			8/*dirs*/ * 2/*su3vec per weyl*/ * (6*9 + 2*3)/*flop per su3 mv-mul*/ +
				/*kamul*/			8/*dirs*/ * (2 * 3)/*cmplx per weyl*/ * 6/*flops cmplx mul*/
			);
		}
	} else {
		// readWeyl
		if (nokamul) {
			if (inputfield->isWeyllayoutSloppy)
				bgq_master_call(&bgq_HoppingMatrix_nokamul_worker_readWeyllayout_float, work);
			else
				bgq_master_call(&bgq_HoppingMatrix_nokamul_worker_readWeyllayout_double, work);
			flopaccumulator += sites * PHYSICAL_LK * (
				/*accum spinor*/	7/*dirs*/ * (4 * 3)/*cmplx per spinor*/ * 2/*flops accum*/ +
				/*weyl reduce*/     8/*dirs*/ * 2/*weyl per dir*/ * (2 * 3)/*cmplx per weyl*/ * 2/*flops*/ +
				/*su3 mul*/			8/*dirs*/ * 2/*su3vec per weyl*/ * (6*9 + 2*3)/*flop per su3 mv-mul*/
			);
		} else {
			if (inputfield->isWeyllayoutSloppy)
				bgq_master_call(&bgq_HoppingMatrix_kamul_worker_readWeyllayout_float, work);
			else
				bgq_master_call(&bgq_HoppingMatrix_kamul_worker_readWeyllayout_double, work);
			flopaccumulator += sites * PHYSICAL_LK * (
				/*accum spinor*/	7/*dirs*/ * (4 * 3)/*cmplx per spinor*/ * 2/*flops accum*/ +
				/*weyl reduce*/     8/*dirs*/ * 2/*weyl per dir*/ * (2 * 3)/*cmplx per weyl*/ * 2/*flops*/ +
				/*su3 mul*/			8/*dirs*/ * 2/*su3vec per weyl*/ * (6*9 + 2*3)/*flop per su3 mv-mul*/ +
				/*kamul*/			8/*dirs*/ * (2 * 3)/*cmplx per weyl*/ * 6/*flops cmplx mul*/
			);
		}
	}
#endif
	//master_print("nokamul=%d readFulllayout=%d sites=%zu flopaccum=%llu diff=%llu\n", nokamul, readFulllayout, sites, flopaccumulator, flopaccumulator-old);
}





void bgq_HoppingMatrix(bool isOdd, bgq_weylfield_controlblock *targetfield, bgq_weylfield_controlblock *spinorfield, bgq_hmflags opts) {
	assert(targetfield); // to be initialized

	assert(spinorfield);
	assert(spinorfield != targetfield);
	assert(spinorfield->isInitialized);
	assert(spinorfield->isOdd == !isOdd);
	bool readFullspinor;
	if (spinorfield->hasFullspinorData) {
		readFullspinor = true;
	} else if (spinorfield->hasWeylfieldData) {
		readFullspinor = false;
	} else {
		assert(!"Input does not contain data");
		UNREACHABLE
	}

	bool nocomm = opts & hm_nocom;
	bool nooverlap = opts & hm_nooverlap;
	bool nokamul = opts & hm_nokamul;
	bool nodistribute = opts & hm_nodistribute;
	bool nodatamove = opts & hm_nodatamove;
	bool nobody = opts & hm_nobody;
	bool nospi = opts & hm_nospi;
	bool noprefetchstream = opts & hm_noprefetchstream;
	bool floatprecision = opts & hm_floatprecision;

	bgq_spinorfield_setup(targetfield, isOdd, false, false, false, true, floatprecision);
	bgq_spinorfield_setup(spinorfield, !isOdd, readFullspinor, false, !readFullspinor, false, false);
	assert(targetfield->isOdd == isOdd);

	bgq_spinorfield_layout layout = bgq_spinorfield_bestLayout(spinorfield);

	// 0. Expect data from other neighbor node
	if (!nocomm) {
		bgq_comm_recv(nospi);
	}


	// 1. Distribute
	bgq_master_sync();
	static bgq_HoppingMatrix_workload work_surface;
	if (PHYSICAL_SURFACE > 0) {
		work_surface.isOdd_src = !isOdd;
		work_surface.isOdd_dst = isOdd;
		work_surface.targetfield = targetfield;
		work_surface.spinorfield = spinorfield;
		work_surface.ic_begin = bgq_surface2collapsed(0);
		work_surface.ic_end = bgq_surface2collapsed(PHYSICAL_SURFACE-1)+1;
		work_surface.noprefetchstream = noprefetchstream;
	}

	static bgq_HoppingMatrix_workload work_body;
	if (PHYSICAL_BODY > 0) {
		work_body.isOdd_src = !isOdd;
		work_body.isOdd_dst = isOdd;
		work_body.targetfield = targetfield;
		work_body.spinorfield = spinorfield;
		work_body.ic_begin = bgq_body2collapsed(0);
		work_body.ic_end = bgq_body2collapsed(PHYSICAL_BODY-1)+1;
		work_body.noprefetchstream = noprefetchstream;
	}

	// Compute surface and put data into the send buffers
	if ((PHYSICAL_SURFACE > 0) && !nodistribute) {
		bgq_HoppingMatrix_work(&work_surface, nokamul, layout);
	}


	if (BGQ_UNVECTORIZE && COMM_T && !nodatamove) {
		bgq_master_sync();
		static bgq_unvectorize_workload work_unvectorize;
		work_unvectorize.isOdd = isOdd;
		work_unvectorize.field = targetfield;
		if (layout & ly_sloppy)
			bgq_master_call(&bgq_HoppingMatrix_unvectorize_float, &work_unvectorize);
		else
			bgq_master_call(&bgq_HoppingMatrix_unvectorize_double, &work_unvectorize);
	}



// 2. Start communication
	if ((PHYSICAL_SURFACE > 0) && !nocomm) {
		bgq_master_sync(); // Wait for threads to finish surface before sending it
		//TODO: ensure there are no other communications pending
		bgq_comm_send(nospi);
		targetfield->waitingForRecv = true;
	}
	if ((PHYSICAL_SURFACE > 0) && !nodatamove) {
		targetfield->pendingDatamove = true;
	}
	targetfield->hmflags = opts;

	if (nooverlap) {
		// Do not wait until data is required, but do it here
		bgq_spinorfield_setup(targetfield, isOdd, false, false, true, false, false);
	}


// 3. Compute the body
	if ((PHYSICAL_BODY > 0) && !nobody) {
		bgq_HoppingMatrix_work(&work_body, nokamul, layout);
		if (!COMM_T && !nodatamove) {
			// Copy the data from HALO_T into the required locations
			bgq_master_sync();
			static bgq_work_datamove work_datamovet;
			work_datamovet.spinorfield = targetfield;
			work_datamovet.opts = opts;
			bgq_master_call(&bgq_HoppingMatrix_datamovet_worker, &work_datamovet);
		}
	}


// 4. Wait for the communication to finish
	/* Defer to bgq_spinorfield_setup as soon as the data is actually required */


// 5. Move received to correct location
	/* Done in bgq_spinorfield_setup whoever is using the field next*/


// 6. Compute the surface
	/* Done by funcs calling readWeyllayout */
}




