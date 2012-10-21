/*
 * bgq_HoppingMatrix.c
 *
 *  Created on: Oct 13, 2012
 *      Author: meinersbur
 */

#define BGQ_HOPPING_MATRIX_C_
#include "bgq_HoppingMatrix.h"

#include "bgq_field.h"
#include "bgq_qpx.h"
#include "bgq_dispatch.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define PRECISION double








static inline void bgq_HoppingMatrix_kernel_raw(bgq_weyl_ptr_t *targetptrs, bgq_weylsite *spinorsite, bgq_gaugesite *gaugesite) {
	bgq_su3_spinor_decl(result);
	//bgq_su3_spinor_zero(result);

	bgq_su3_weyl_prefetch(&spinorsite->tdown1); //TODO: currently has one cacheline more
	bgq_su3_matrix_prefetch(&gaugesite->su3[TDOWN]);

	// T+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tup);
		bgq_su3_mdecl(gauge_tup);

		bgq_su3_weyl_load_combine_double(weyl_tup, &spinorsite->tup1, &spinorsite->tup2);
		bgq_su3_matrix_load_double(gauge_tup, &gaugesite->su3[TUP]);
		bgq_su3_weyl_mvmul(weyl_tup, gauge_tup, weyl_tup);

		bgq_su3_expand_weyl_tup(result, weyl_tup);
	}

	bgq_su3_weyl_prefetch(&spinorsite->xup);
	bgq_su3_matrix_prefetch(&gaugesite->su3[XUP]);

	// T- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tdown);
		bgq_su3_mdecl(gauge_tdown);

		bgq_su3_weyl_load_combine_double(weyl_tdown, &spinorsite->tdown1, &spinorsite->tdown2);
		bgq_su3_matrix_load_double(gauge_tdown, &gaugesite->su3[TDOWN]);
		bgq_su3_weyl_mvmul(weyl_tdown, gauge_tdown, weyl_tdown);

		bgq_su3_accum_weyl_tdown(result, weyl_tdown);
	}

	bgq_su3_weyl_prefetch(&spinorsite->xdown);
	bgq_su3_matrix_prefetch(&gaugesite->su3[XDOWN]);

	// X+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xup);
		bgq_su3_mdecl(gauge_xup);

		bgq_su3_weyl_load_double(weyl_xup, &spinorsite->xup);
		bgq_su3_matrix_load_double(gauge_xup, &gaugesite->su3[XUP]);
		bgq_su3_weyl_mvmul(weyl_xup, gauge_xup, weyl_xup);

		bgq_su3_accum_weyl_xup(result, weyl_xup);
	}

	bgq_su3_weyl_prefetch(&spinorsite->yup);
	bgq_su3_matrix_prefetch(&gaugesite->su3[YUP]);

	// X- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xdown);
		bgq_su3_mdecl(gauge_xdown);

		bgq_su3_weyl_load_double(weyl_xdown, &spinorsite->xdown);
		bgq_su3_matrix_load_double(gauge_xdown, &gaugesite->su3[XDOWN]);
		bgq_su3_weyl_mvinvmul(weyl_xdown, gauge_xdown, weyl_xdown);

		bgq_su3_accum_weyl_xdown(result, weyl_xdown);
	}

	bgq_su3_weyl_prefetch(&spinorsite->ydown);
	bgq_su3_matrix_prefetch(&gaugesite->su3[YDOWN]);

	// Y+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_yup);
		bgq_su3_mdecl(gauge_yup);

		bgq_su3_weyl_load_double(weyl_yup, &spinorsite->yup);
		bgq_su3_matrix_load_double(gauge_yup, &gaugesite->su3[YUP]);
		bgq_su3_weyl_mvmul(weyl_yup, gauge_yup, weyl_yup);

		bgq_su3_accum_weyl_yup(result, weyl_yup);
	}

	bgq_su3_weyl_prefetch(&spinorsite->zup);
	bgq_su3_matrix_prefetch(&gaugesite->su3[ZUP]);

	// Y- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_ydown);
		bgq_su3_mdecl(gauge_ydown);

		bgq_su3_weyl_load_double(weyl_ydown, &spinorsite->ydown);
		bgq_su3_matrix_load_double(gauge_ydown, &gaugesite->su3[YDOWN]);
		bgq_su3_weyl_mvinvmul(weyl_ydown, gauge_ydown, weyl_ydown);

		bgq_su3_accum_weyl_ydown(result, weyl_ydown);
	}

	bgq_su3_weyl_prefetch(&spinorsite->zdown);
	bgq_su3_matrix_prefetch(&gaugesite->su3[ZDOWN]);

	// Z+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zup);
		bgq_su3_mdecl(gauge_zup);

		bgq_su3_weyl_load_double(weyl_zup, &spinorsite->zup);
		bgq_su3_matrix_load_double(gauge_zup, &gaugesite->su3[ZUP]);
		bgq_su3_weyl_mvmul(weyl_zup, gauge_zup, weyl_zup);

		bgq_su3_accum_weyl_zup(result, weyl_zup);
	}

	bgq_su3_weyl_prefetch(&(spinorsite + 1)->tup1); // next iteration; TODO: currently t-direction is 1 cacheline more of data
	bgq_su3_matrix_prefetch(&(gaugesite + 1)->su3[TUP]);

	// Z- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zdown);
		bgq_su3_mdecl(gauge_zdown);

		bgq_su3_weyl_load_double(weyl_zdown, &spinorsite->ydown);
		bgq_su3_matrix_load_double(gauge_zdown, &gaugesite->su3[ZDOWN]);
		bgq_su3_weyl_mvinvmul(weyl_zdown, gauge_zdown, weyl_zdown);

		bgq_su3_accum_weyl_zdown(result, weyl_zdown);
	}

	// Store the result
	bgq_spinorfield_weyl_store_raw(targetptrs, bgq_su3_spinor_vars(result));
}


typedef struct {
	bool isOdd;
	bgq_weylfield_controlblock *targetfield;
	bgq_weylfield_controlblock *spinorfield;
} bgq_HoppingMatrix_workload;


static void bgq_HoppingMatrix_worker_datamove(void *argptr, size_t tid, size_t threads) {
	//TODO: This code is relatively irregular kernel; should be possible to straighten it
	const bgq_HoppingMatrix_workload *args = argptr;
	const bgq_weylfield_controlblock *spinorfield = args->spinorfield;

	const size_t workload = PHYSICAL_HALO;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min(workload, begin+threadload);
	for (size_t i = begin; i<end; ) {
		WORKLOAD_DECL(i,workload);

		if (COMM_T && WORKLOAD_SPLIT(2*LOCAL_HALO_T)) {
			//TODO: Iteration workload imbalance
			// One DIM_T iteration loads/stores just half data as the other dimensions

			bgq_weyl_nonvec **destptrs;
			bgq_weyl_nonvec *weylbase;
			if (WORKLOAD_CHUNK(2)) {
				destptrs = spinorfield->destptrFromRecvTup;
				weylbase = spinorfield->sec_recv_tup;
			} else {
				destptrs = spinorfield->destptrFromRecvTdown;
				weylbase = spinorfield->sec_recv_tdown;
			}

			const size_t beginj = WORKLOAD_PARAM(LOCAL_HALO_T);
			const size_t endj = min(LOCAL_HALO_T, threadload);
			for (size_t j = beginj; j < endj; j+=1) {
				bgq_weyl_nonvec *weylsrc = &weylbase[j]; //TODO: Check strength reduction
				bgq_weyl_nonvec *destptr = destptrs[j];//TODO: Also prefetch destptrs

				bgq_su3_weyl_prefetch_double(&weylbase[j+1]); //TODO: prefetch just _nonvec
				bgq_su3_weyl_left_move_double(destptr, weylsrc);
			}
			i += (endj - beginj);
		} else {
			bgq_weyl_vec **destptrs;
			bgq_weyl_vec *weylbase;
			size_t count;

			if (COMM_X && WORKLOAD_SPLIT(2*LOCAL_HALO_X)) {
				if (WORKLOAD_CHUNK(2)) {
					destptrs = spinorfield->destptrFromRecvXup;
					weylbase = spinorfield->sec_recv_xup;
				} else {
					destptrs = spinorfield->destptrFromRecvXdown;
					weylbase = spinorfield->sec_recv_xdown;
				}
				count = LOCAL_HALO_X;
			} else if (COMM_Y && WORKLOAD_SPLIT(2*LOCAL_HALO_Y)) {
				if (WORKLOAD_CHUNK(2)) {
					destptrs = spinorfield->destptrFromRecvYup;
					weylbase = spinorfield->sec_recv_yup;
				} else {
					destptrs = spinorfield->destptrFromRecvYdown;
					weylbase = spinorfield->sec_recv_ydown;
				}
				count = LOCAL_HALO_Y;
			} else if(COMM_Z && WORKLOAD_SPLIT(2*LOCAL_HALO_Z)) {
				if (WORKLOAD_CHUNK(2)) {
					destptrs = spinorfield->destptrFromRecvZup;
					weylbase = spinorfield->sec_recv_zup;
				} else {
					destptrs = spinorfield->destptrFromRecvZdown;
					weylbase = spinorfield->sec_recv_zdown;
				}
				count = LOCAL_HALO_Z;
			} else {
				assert(!"Impossible to get here");
				exit(1);
			}

			const size_t beginj = WORKLOAD_PARAM(count);
			const size_t endj = min(count, threadload);
			for (size_t j = beginj; j < endj; j+=1) {
				bgq_weyl_vec *weylsrc = &weylbase[j]; //TODO: Check strength reduction
				bgq_weyl_vec *destptr = destptrs[j];

				bgq_su3_weyl_prefetch_double(&weylbase[j+1]);//TODO: Also prefetch destptrs
				bgq_su3_weyl_decl(weyl);
				bgq_su3_weyl_load_double(weyl, weylsrc);
				bgq_su3_weyl_store_double(destptr, weyl);
			}
			i += (endj - beginj);
		}

		WORKLOAD_CHECK
	}
}

static void bgq_HoppingMatrix_worker_surface(void *argptr, size_t tid, size_t threads) {
	const bgq_HoppingMatrix_workload *args = argptr;
	const bool isOdd = args->isOdd;
	const bgq_weylfield_controlblock *targetfield = args->targetfield;
	const bgq_weylfield_controlblock *spinorfield = args->spinorfield;

	const size_t workload = PHYSICAL_SURFACE;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min(workload, begin+threadload);
	for (size_t is = begin; is<end; is+=1) {
		//TODO: Check strength reduction
		bgq_weylsite *spinorsite = &spinorfield->sec_surface[is];
		bgq_gaugesite *gaugesite = &g_bgq_gaugefield_fromSurface[isOdd][is];
		bgq_weyl_ptr_t *destptrs = &targetfield->destptrFromSurface[is];

		bgq_HoppingMatrix_kernel_raw(destptrs, spinorsite, gaugesite); //TODO: Check if inlined
	}
}


static void bgq_HoppingMatrix_worker_body(void *argptr, size_t tid, size_t threads) {
	const bgq_HoppingMatrix_workload *args = argptr;
	const bool isOdd = args->isOdd;
	const bgq_weylfield_controlblock *targetfield = args->targetfield;
	const bgq_weylfield_controlblock *spinorfield = args->spinorfield;

	const size_t workload = PHYSICAL_BODY;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min(workload, begin+threadload);
	for (size_t ib = begin; ib<end; ib+=1) {
		//TODO: Check strength reduction
		bgq_weylsite *spinorsite = &spinorfield->sec_body[ib];
		bgq_gaugesite *gaugesite = &g_bgq_gaugefield_fromBody[isOdd][ib];
		bgq_weyl_ptr_t *destptrs = &targetfield->destptrFromBody[ib];

		bgq_HoppingMatrix_kernel_raw(destptrs, spinorsite, gaugesite); //TODO: Check if inlined
	}
}


void bgq_HoppingMatrix(bool isOdd, bgq_weylfield_controlblock *targetfield, bgq_weylfield_controlblock *spinorfield, bgq_hmflags opts) {
	assert(targetfield);
	bgq_spinorfield_setup(targetfield, isOdd, false, true,false,false);
	bgq_spinorfield_setup(spinorfield, !isOdd, false, false,true,false);
	assert(targetfield->isOdd == isOdd);

	assert(spinorfield);
	assert(spinorfield->isInitinialized);
	assert(spinorfield->isSloppy == false);
	assert(spinorfield->isOdd == !isOdd);
	assert(spinorfield->hasWeylfieldData);

// 1. Distribute
// Required before calling this function or done in bgq_spinorfield_setup

// 2. Start communication
	/* not yet implemented */

// 3. Compute the body
	bgq_HoppingMatrix_workload work = {
		.isOdd = isOdd,
		.targetfield = targetfield,
		.spinorfield = spinorfield
	};
	bgq_master_call(&bgq_HoppingMatrix_worker_body, &work);

// 4. Wait for the communication to finish
	/* Not yet implemented */

// 5. Move received to correct location
	bgq_master_call(&bgq_HoppingMatrix_worker_datamove, &work);

// 6. Compute the surface
	bgq_master_call(&bgq_HoppingMatrix_worker_surface, &work);
}




