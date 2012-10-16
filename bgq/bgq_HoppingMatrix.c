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

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>


inline void bgq_spinorfield_weyl_store_fromHalfvolume(bgq_weylfield_controlblock *targetfield, bool isOdd, size_t ih, bgq_su3_spinor_params(spinor)) {
	//bgq_spinorfield_reset(targetfield, isOdd, true, false);
	assert(targetfield->isInitinialized);
	assert(targetfield->isOdd == isOdd);
	assert(targetfield->hasWeylfieldData == true);

	bgq_weyl_ptr_t *weylptrs = targetfield->destptrFromHalfvolume[ih]; // TODO: Check that compiler does strength reduction after inline, otherwise do manually
	//TODO: probably compiler will li an offset for every destptrFromHalfvolume, can do better using addi

	// T+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tup);
		bgq_su3_reduce_weyl_tup(weyl_tup, spinor);
		bgq_su3_weyl_store_left_double(weylptrs->pd[P_TUP1], weyl_tup);
		bgq_su3_weyl_store_right_double(weylptrs->pd[P_TUP2], weyl_tup);
	}

	// T- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tdown);
		bgq_su3_reduce_weyl_tdown(weyl_tdown, spinor);
		bgq_su3_weyl_store_left_double(weylptrs->pd[P_TDOWN1], weyl_tup);
		bgq_su3_weyl_store_right_double(weylptrs->pd[P_TDOWN2], weyl_tup);
	}

	// X+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xup);
		bgq_su3_reduce_weyl_xup(weyl_xup, spinor);
		bgq_su3_weyl_store_double(weylptrs->pd[P_XUP], weyl_xup);
	}

	// X- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xdown);
		bgq_su3_reduce_weyl_xdown(weyl_xdown, spinor);
		bgq_su3_weyl_store_double(weylptrs->pd[P_XDOWN], weyl_xdown);
	}

	// Y+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_yup);
		bgq_su3_reduce_weyl_yup(weyl_yup, spinor);
		bgq_su3_weyl_store_double(weylptrs->pd[P_YUP], weyl_yup);
	}

	// Y- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_ydown);
		bgq_su3_reduce_weyl_ydown(weyl_ydown, spinor);
		bgq_su3_weyl_store_double(weylptrs->pd[P_YDOWN], weyl_ydown);
	}

	// Z+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zup);
		bgq_su3_reduce_weyl_zup(weyl_zup, spinor);
		bgq_su3_weyl_store_double(weylptrs->pd[P_ZUP], weyl_zup);
	}

	// Z- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zdown);
		bgq_su3_reduce_weyl_zdown(weyl_zdown, spinor);
		bgq_su3_weyl_store_double(weylptrs->pd[P_ZDOWN], weyl_zdown);
	}
}


inline void bgq_HoppingMatrix_kernel_surface() {
	bgq_su3_spinor_decl(result);
			bgq_su3_spinor_zero(result);

			// T+ //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_tup);
				bgq_su3_mdecl(gauge_tup);

				bgq_weyl_nonvec *weyladdr_tup1 = surfacedata;
				surfacedata += sizeof(bgq_weyl_nonvec);
				bgq_weyl_nonvec *weyladdr_tup2= surfacedata;
				surfacedata += sizeof(bgq_weyl_nonvec);
				bgq_su3_weyl_load_combine_double(weyl_tup, weyladdr_tup1, weyladdr_tup2);

				bgq_gaugesite *gaugesite_tup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, TUP, t1, t2, true,false);
				bgq_su3_matrix_load(gauge_tup, gaugesite_tup);
				//bgq_su3_matrix_flush(gaugesite_tup);
				bgq_su3_weyl_mvmul(weyl_tup,gauge_tup,weyl_tup);

				bgq_su3_expand_weyl_tup(result,weyl_tup);
			}

			bgq_su3_spinor_prefetch(addr + 2*4*3*2);
			bgq_su3_matrix_prefetch(gaugeaddr + 2*3*3*2);

			// T- //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_tdown);
				bgq_su3_mdecl(gauge_tdown);

				bgq_weyl_nonvec *weyladdr_tdown1 = surfacedata;
				surfacedata += sizeof(bgq_weyl_nonvec);
				bgq_weyl_nonvec *weyladdr_tdown2= surfacedata;
				surfacedata += sizeof(bgq_weyl_nonvec);
				bgq_su3_weyl_load_combine_double(weyl_tdown, weyladdr_tdown1, weyladdr_tdown2);

				bgq_gaugesite *gaugesite_tdown = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, TUP, t1-1, t2-1, true,false);
				bgq_su3_matrix_load(gauge_tdown, gaugesite_tdown);
				//bgq_su3_matrix_flush(gaugesite_tdown);
				bgq_su3_weyl_mvinvmul(weyl_tdown, gauge_tdown, weyl_tdown);

				bgq_su3_accum_weyl_tdown(result, weyl_tdown);
			}

			bgq_su3_spinor_prefetch(addr + 2*4*3*2);
			bgq_su3_matrix_prefetch(gaugeaddr + 2*3*3*2);

			// X+ //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_spinor_decl(spinor_xup);
				bgq_su3_weyl_decl(weyl_xup);
				bgq_su3_mdecl(gauge_xup);

				bgq_spinorsite *spinorsite_xup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x+1, y, z, t1, t2, true,false);
				bgq_su3_spinor_load(spinor_xup, spinorsite_xup);
				//bgq_su3_spinor_flush(spinorsite_xup);
				//bgq_su3_spinor_invalidate(spinorsite_xup);
				bgq_su3_reduce_weyl_xup(weyl_xup, spinor_xup);

				bgq_gaugesite *gaugesite_xup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, XUP, t1, t2, true,false);
				bgq_su3_matrix_load(gauge_xup, gaugesite_xup);
				//bgq_su3_matrix_flush(gaugesite_xup);
				//bgq_su3_matrix_invalidate(gaugesite_xup);
				bgq_su3_weyl_mvmul(weyl_xup, gauge_xup, weyl_xup);

				bgq_su3_accum_weyl_xup(result, weyl_xup);
			}

			bgq_su3_spinor_prefetch(addr + 2*4*3*2);
			bgq_su3_matrix_prefetch(gaugeaddr + 2*3*3*2);

			// X- //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_spinor_decl(spinor_xdown);
				bgq_su3_weyl_decl(weyl_xdown);
				bgq_su3_mdecl(gauge_xdown);

				bgq_spinorsite *spinorsite_xdown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x-1, y, z, t1, t2, true,false);
				bgq_su3_spinor_load(spinor_xdown, spinorsite_xdown);
				//bgq_su3_spinor_flush(spinorsite_xdown);
				//bgq_su3_spinor_invalidate(spinorsite_xdown);
				bgq_su3_reduce_weyl_xdown(weyl_xdown, spinor_xdown);

				bgq_gaugesite *gaugesite_xdown = BGQ_GAUGESITE(gaugefield, isOdd, tv, x-1, y, z, XUP, t1, t2, true,false);
				bgq_su3_matrix_load(gauge_xdown, gaugesite_xdown);
				//bgq_su3_matrix_flush(gaugesite_xdown);
				//bgq_su3_matrix_invalidate(gaugesite_xdown);
				bgq_su3_weyl_mvinvmul(weyl_xdown, gauge_xdown, weyl_xdown);

				bgq_su3_accum_weyl_xdown(result, weyl_xdown);
			}

			bgq_su3_spinor_prefetch(addr + 2*4*3*2);
			bgq_su3_matrix_prefetch(gaugeaddr + 2*3*3*2);

			// Y+ //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_spinor_decl(spinor_yup);
				bgq_su3_weyl_decl(weyl_yup);
				bgq_su3_mdecl(gauge_yup);

				bgq_spinorsite *spinorsite_yup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y+1, z, t1, t2, true,false);
				bgq_su3_spinor_load(spinor_yup, spinorsite_yup);
				//bgq_su3_spinor_flush(spinorsite_yup);
				//bgq_su3_spinor_invalidate(spinorsite_yup);
				bgq_su3_reduce_weyl_yup(weyl_yup, spinor_yup);

				bgq_gaugesite *gaugesite_yup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, YUP, t1, t2, true,false);
				bgq_su3_matrix_load(gauge_yup, gaugesite_yup);
				//bgq_su3_matrix_flush(gaugesite_yup);
				//bgq_su3_matrix_invalidate(gaugesite_yup);
				bgq_su3_weyl_mvmul(weyl_yup, gauge_yup, weyl_yup);

				bgq_su3_accum_weyl_yup(result, weyl_yup);
			}

			bgq_su3_spinor_prefetch(addr + 2*4*3*2);
			bgq_su3_matrix_prefetch(gaugeaddr + 2*3*3*2);

			// Y- //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_spinor_decl(spinor_ydown);
				bgq_su3_weyl_decl(weyl_ydown);
				bgq_su3_mdecl(gauge_ydown);

				bgq_spinorsite *spinorsite_ydown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y-1, z, t1, t2, true,false);
				bgq_su3_spinor_load(spinor_ydown, spinorsite_ydown);
				//bgq_su3_spinor_flush(spinorsite_ydown);
				//bgq_su3_spinor_invalidate(spinorsite_ydown);
				bgq_su3_reduce_weyl_ydown(weyl_ydown, spinor_ydown);

				bgq_gaugesite *gaugesite_ydown = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y-1, z, YUP, t1, t2, true,false);
				bgq_su3_matrix_load(gauge_ydown, gaugesite_ydown);
				//bgq_su3_matrix_flush(gaugesite_ydown);
				//bgq_su3_matrix_invalidate(gaugesite_ydown);
				bgq_su3_weyl_mvinvmul(weyl_ydown, gauge_ydown, weyl_ydown);

				bgq_su3_accum_weyl_ydown(result, weyl_ydown);
			}

			bgq_su3_spinor_prefetch(addr + 2*4*3*2);
			bgq_su3_matrix_prefetch(gaugeaddr + 2*3*3*2);

			// Z+ //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_spinor_decl(spinor_zup);
				bgq_su3_weyl_decl(weyl_zup);
				bgq_su3_mdecl(gauge_zup);

				bgq_spinorsite *spinorsite_zup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z+1, t1, t2, true,false);
				bgq_su3_spinor_load(spinor_zup, spinorsite_zup);
				//bgq_su3_spinor_flush(spinorsite_zup);
				//bgq_su3_spinor_invalidate(spinorsite_zup);
				bgq_su3_reduce_weyl_zup(weyl_zup, spinor_zup);

				bgq_gaugesite *gaugesite_zup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, ZUP, t1, t2, true,false);
				bgq_su3_matrix_load(gauge_zup, gaugesite_zup);
				//bgq_su3_matrix_flush(gaugesite_zup);
				//bgq_su3_matrix_invalidate(gaugesite_zup);
				bgq_su3_weyl_mvmul(weyl_zup, gauge_zup, weyl_zup);

				bgq_su3_accum_weyl_zup(result, weyl_zup);
			}

			bgq_su3_spinor_prefetch(addr + 2*4*3*2);
			bgq_su3_matrix_prefetch(gaugeaddr + 2*3*3*2);

			// Z- //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_spinor_decl(spinor_zdown);
				bgq_su3_weyl_decl(weyl_zdown);
				bgq_su3_mdecl(gauge_zdown);

				bgq_spinorsite *spinorsite_zdown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z-1, t1, t2, true, false);
				bgq_su3_spinor_load(spinor_zdown, spinorsite_zdown);
				//bgq_su3_spinor_flush(spinorsite_zdown);
				//bgq_su3_spinor_invalidate(spinorsite_zdown);
				bgq_su3_reduce_weyl_zdown(weyl_zdown, spinor_zdown);

				bgq_gaugesite *gaugesite_zdown = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z-1, ZUP, t1, t2, true, false);
				bgq_su3_matrix_load(gauge_zdown, gaugesite_zdown);
				//bgq_su3_matrix_flush(gaugesite_zdown);
				//bgq_su3_matrix_invalidate(gaugesite_zdown);
				bgq_su3_weyl_mvinvmul(weyl_zdown, gauge_zdown, weyl_zdown);

				bgq_su3_accum_weyl_zdown(result, weyl_zdown);
			}

			// Store the result
			//for (int i = 0; i < 4; i+=1) {
			bgq_spinorsite *targetsite = BGQ_SPINORSITE(targetfield, isOdd, tv, x, y, z, t1,t2, false, true);
			bgq_su3_spinor_zeroload(targetsite); /* no impact on performance */
			bgq_su3_spinor_store(targetsite, result);
			//bgq_su3_spinor_flush(targetsite); /* reduces performance */
			//}
}


void bgq_HoppingMatrix(bool isOdd, bgq_weylfield_controlblock *targetfield, bgq_weylfield_controlblock *spinorfield, bgq_gaugesite, int tid, int threads) {
	bgq_weylfield weyldata;

// 1. Distribute
	{
	size_t workload = PHYSICAL_VOLUME;
	size_t threadload = (PHYSICAL_VOLUME+threads-1)/threads;
	size_t begin = tid*threadload;
	size_t end = min(workload, (tid+1)*threadload);
	uint32_t *weylindex = bgq_distr_indices[begin*8];

	for (int xyztv = begin; xyztv < end; xyztv+=1) {
		bgq_su3_spinor_decl(spinor);
		// Get the spinor from somewhere
		bgq_HoppingMatrix_weylfield();
	}
	}



// 2. Start communication

// 3. Compute the body
	{
		size_t workload = PHYSICAL_BODY;
		size_t threadload = (workload+threads-1)/threads;
		size_t begin = tid*threadload;
		size_t end = min(workload,begin+threadload);
		COMPLEX_PRECISION *surfacedata = (char*)weyldata + bgq_weyl_section_offset(sec_surface);
		for (int i = begin; i < end; i+=1) {
			bgq_HoppingMatrix_kernel_surface();
		}
	}

// 4. Wait for communication to finish

// 5. Compute the surface
}



// Simple interface
void bgq_HoppingMatrix(bgq_weylfield_controlblock *targetfield, bgq_weylfield_controlblock *spinorfield) {
	assert(spinorfield->hasWeylfieldData);
	assert(!spinorfield->isSloppy);
	bool isOdd = !spinorfield->idOdd;

// 0. Configure the target field
	bgq_spinorfield_reset(targetfield, isOdd, true, false);

// 1. Start SPI communication

// 2. Compute body
}


