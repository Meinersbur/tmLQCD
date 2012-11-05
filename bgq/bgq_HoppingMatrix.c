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

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define PRECISION double






static inline void bgq_HoppingMatrix_kernel_raw(bgq_weyl_ptr_t *targetptrs, bgq_weylsite *spinorsite, bgq_gaugesite *gaugesite) {
	bgq_su3_spinor_decl(result);
	//bgq_su3_spinor_zero(result);

	bgq_su3_weyl_prefetch(&spinorsite->d[TDOWN]); //TODO: currently has one cacheline more
	bgq_su3_matrix_prefetch(&gaugesite->su3[TDOWN]);

	// T+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tup);
		bgq_su3_mdecl(gauge_tup);

		bgq_su3_weyl_load_double(weyl_tup, &spinorsite->d[TUP]);
		bgq_su3_matrix_load_double(gauge_tup, &gaugesite->su3[TUP]);
		bgq_su3_weyl_mvmul(weyl_tup, gauge_tup, weyl_tup);

		bgq_su3_expand_weyl_tup(result, weyl_tup);
	}

	bgq_su3_weyl_prefetch(&spinorsite->d[XUP]);
	bgq_su3_matrix_prefetch(&gaugesite->su3[XUP]);

	// T- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tdown);
		bgq_su3_mdecl(gauge_tdown);

		bgq_su3_weyl_load_double(weyl_tdown, &spinorsite->d[TDOWN]);
		bgq_su3_matrix_load_double(gauge_tdown, &gaugesite->su3[TDOWN]);
		bgq_su3_weyl_mvmul(weyl_tdown, gauge_tdown, weyl_tdown);

		bgq_su3_accum_weyl_tdown(result, weyl_tdown);
	}

	bgq_su3_weyl_prefetch(&spinorsite->d[XDOWN]);
	bgq_su3_matrix_prefetch(&gaugesite->su3[XDOWN]);

	// X+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xup);
		bgq_su3_mdecl(gauge_xup);

		bgq_su3_weyl_load_double(weyl_xup, &spinorsite->d[XUP]);
		bgq_su3_matrix_load_double(gauge_xup, &gaugesite->su3[XUP]);
		bgq_su3_weyl_mvmul(weyl_xup, gauge_xup, weyl_xup);

		bgq_su3_accum_weyl_xup(result, weyl_xup);
	}

	bgq_su3_weyl_prefetch(&spinorsite->d[YUP]);
	bgq_su3_matrix_prefetch(&gaugesite->su3[YUP]);

	// X- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xdown);
		bgq_su3_mdecl(gauge_xdown);

		bgq_su3_weyl_load_double(weyl_xdown, &spinorsite->d[XDOWN]);
		bgq_su3_matrix_load_double(gauge_xdown, &gaugesite->su3[XDOWN]);
		bgq_su3_weyl_mvinvmul(weyl_xdown, gauge_xdown, weyl_xdown);

		bgq_su3_accum_weyl_xdown(result, weyl_xdown);
	}

	bgq_su3_weyl_prefetch(&spinorsite->d[YDOWN]);
	bgq_su3_matrix_prefetch(&gaugesite->su3[YDOWN]);

	// Y+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_yup);
		bgq_su3_mdecl(gauge_yup);

		bgq_su3_weyl_load_double(weyl_yup, &spinorsite->d[YUP]);
		bgq_su3_matrix_load_double(gauge_yup, &gaugesite->su3[YUP]);
		bgq_su3_weyl_mvmul(weyl_yup, gauge_yup, weyl_yup);

		bgq_su3_accum_weyl_yup(result, weyl_yup);
	}

	bgq_su3_weyl_prefetch(&spinorsite->d[ZUP]);
	bgq_su3_matrix_prefetch(&gaugesite->su3[ZUP]);

	// Y- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_ydown);
		bgq_su3_mdecl(gauge_ydown);

		bgq_su3_weyl_load_double(weyl_ydown, &spinorsite->d[YDOWN]);
		bgq_su3_matrix_load_double(gauge_ydown, &gaugesite->su3[YDOWN]);
		bgq_su3_weyl_mvinvmul(weyl_ydown, gauge_ydown, weyl_ydown);

		bgq_su3_accum_weyl_ydown(result, weyl_ydown);
	}

	bgq_su3_weyl_prefetch(&spinorsite->d[ZDOWN]);
	bgq_su3_matrix_prefetch(&gaugesite->su3[ZDOWN]);

	// Z+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zup);
		bgq_su3_mdecl(gauge_zup);

		bgq_su3_weyl_load_double(weyl_zup, &spinorsite->d[ZUP]);
		bgq_su3_matrix_load_double(gauge_zup, &gaugesite->su3[ZUP]);
		bgq_su3_weyl_mvmul(weyl_zup, gauge_zup, weyl_zup);

		bgq_su3_accum_weyl_zup(result, weyl_zup);
	}

	bgq_su3_weyl_prefetch(&(spinorsite + 1)->d[TUP]); // next iteration; TODO: currently t-direction is 1 cacheline more of data
	bgq_su3_matrix_prefetch(&(gaugesite + 1)->su3[TUP]);

	// Z- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zdown);
		bgq_su3_mdecl(gauge_zdown);

		bgq_su3_weyl_load_double(weyl_zdown, &spinorsite->d[ZDOWN]);
		bgq_su3_matrix_load_double(gauge_zdown, &gaugesite->su3[ZDOWN]);
		bgq_su3_weyl_mvinvmul(weyl_zdown, gauge_zdown, weyl_zdown);

		bgq_su3_accum_weyl_zdown(result, weyl_zdown);
	}

	// Store the result
	bgq_HoppingMatrix_storeWeyllayout_raw(targetptrs, bgq_su3_spinor_vars(result));
}


typedef struct {
	bool isOdd_src;
	bool isOdd_dst;
	bgq_weylfield_controlblock *targetfield;
	bgq_weylfield_controlblock *spinorfield;
	ucoord ic_begin;
	ucoord ic_end;
	bool noprefetchstream;
} bgq_HoppingMatrix_workload;



static void bgq_HoppingMatrix_worker_surface(void *argptr, size_t tid, size_t threads) {
	const bgq_HoppingMatrix_workload *args = argptr;
	const bool isOdd = args->isOdd_src;
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
	const bool isOdd = args->isOdd_src;
	const bgq_weylfield_controlblock *targetfield = args->targetfield;
	const bgq_weylfield_controlblock *spinorfield = args->spinorfield;

	const size_t workload = PHYSICAL_SURFACE;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min(workload, begin+threadload);
	for (size_t ib = begin; ib<end; ib+=1) {
		//TODO: Check strength reduction
		bgq_weylsite *spinorsite = &spinorfield->sec_body[ib];
		bgq_gaugesite *gaugesite = &g_bgq_gaugefield_fromBody[isOdd][ib];
		bgq_weyl_ptr_t *destptrs = &targetfield->destptrFromBody[ib];
		for (size_t d = 0; d < PHYSICAL_LD; d+=1){
			assert((bgq_weyl_vec*)targetfield->sec_weyl <= destptrs->d[d] && destptrs->d[d] < (bgq_weyl_vec*)targetfield->sec_end);
		}

		bgq_HoppingMatrix_kernel_raw(destptrs, spinorsite, gaugesite); //TODO: Check if inlined
	}
}



static inline void bgq_HoppingMatrix_worker_readFulllayout(void * restrict arg, size_t tid, size_t threads, bool kamul, bool readFulllayout) {
	bgq_HoppingMatrix_workload *work = arg;
	bool isOdd = work->isOdd_src;
	bgq_weylfield_controlblock *spinorfield = work->spinorfield;
	bgq_weylfield_controlblock *targetfield = work->targetfield;
	ucoord ic_begin = work->ic_begin;
	ucoord ic_end = work->ic_end;
	bool noprefetchstream = work->noprefetchstream;

	bgq_vector4double_decl(qka0);
	bgq_complxval_splat(qka0,ka0);
	bgq_vector4double_decl(qka1);
	bgq_complxval_splat(qka1,ka1);
	bgq_vector4double_decl(qka2);
	bgq_complxval_splat(qka2,ka2);
	bgq_vector4double_decl(qka3);
	bgq_complxval_splat(qka3,ka3);

	const size_t workload = ic_end - ic_begin;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = ic_begin + tid*threadload;
	const size_t end = ic_end + min_sizet(workload, begin+threadload);

	if (!noprefetchstream) {
		bgq_prefetch_forward(&g_bgq_gaugefield_fromCollapsed[isOdd][begin]);
		if (readFulllayout) {
			bgq_prefetch_forward(&spinorfield->sec_fullspinor[begin]);
		} else {
			bgq_prefetch_forward(&spinorfield->sec_collapsed[begin]);
		}
		bgq_prefetch_forward(&targetfield->sendptr[begin]);
	}

	for (ucoord ic = begin; ic<end; ic+=1) {
		//TODO: Check optaway
		ucoord ih = bgq_collapsed2halfvolume(isOdd,ic);
		ucoord t1 = bgq_halfvolume2t1(isOdd,ih);
		ucoord t2 = bgq_halfvolume2t2(isOdd,ih);
		ucoord x = bgq_halfvolume2x(ih);
		ucoord y = bgq_halfvolume2y(ih);
		ucoord z = bgq_halfvolume2z(ih);

		bgq_gaugesite *gaugesite = &g_bgq_gaugefield_fromCollapsed[isOdd][ic];
		bgq_weyl_ptr_t *destptrs = &targetfield->sendptr[ic];

		//TODO: prefetching
		//TODO: Check inlining
		bgq_su3_spinor_decl(spinor);
		if (readFulllayout) {
			bgq_spinorsite *spinorsite = &spinorfield->sec_fullspinor[ic];
			bgq_su3_spinor_prefetch_double(&spinorfield->sec_fullspinor[ic+1]); // TODO: This prefetch is too early
			bgq_HoppingMatrix_loadFulllayout(spinor, spinorsite, t1, t2, x, y, z);
		} else {
			bgq_weylsite *weylsite = &spinorfield->sec_collapsed[ic];
			bgq_HoppingMatrix_loadWeyllayout(spinor, weylsite, t1, t2, x, y, z);
		}
		bgq_HoppingMatrix_compute_storeWeyllayout_alldir(destptrs, gaugesite, spinor, t1, t2, x, y, z, qka0,qka1,qka2,qka3,kamul);
	}
}

static void bgq_HoppingMatrix_nokamul_worker_readFulllayout(void *arg, size_t tid, size_t threads) {
	bgq_HoppingMatrix_worker_readFulllayout(arg,tid,threads,false,true);
}

static void bgq_HoppingMatrix_kamul_worker_readFulllayout(void *arg, size_t tid, size_t threads) {
	bgq_HoppingMatrix_worker_readFulllayout(arg,tid,threads,true,true);
}
static void bgq_HoppingMatrix_nokamul_worker_readWeyllayout(void *arg, size_t tid, size_t threads) {
	bgq_HoppingMatrix_worker_readFulllayout(arg,tid,threads,false,false);
}
static void bgq_HoppingMatrix_kamul_worker_readWeyllayout(void *arg, size_t tid, size_t threads) {
	bgq_HoppingMatrix_worker_readFulllayout(arg,tid,threads,true,false);
}


static void bgq_HoppingMatrix_worker_surface_precomm_readFulllayout(void *arg, size_t tid, size_t threads) {
	bgq_HoppingMatrix_workload *work = arg;
	bool isOdd = work->isOdd_src;
	bgq_weylfield_controlblock *spinorfield = work->spinorfield;
	bgq_weylfield_controlblock *targetfield = work->targetfield;
	const bool kamul = false;

	bgq_vector4double_decl(qka0);
	bgq_complxval_splat(qka0,ka0);
	bgq_vector4double_decl(qka1);
	bgq_complxval_splat(qka1,ka1);
	bgq_vector4double_decl(qka2);
	bgq_complxval_splat(qka2,ka2);
	bgq_vector4double_decl(qka3);
	bgq_complxval_splat(qka3,ka3);

	const size_t workload = PHYSICAL_SURFACE;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);
	for (ucoord is = begin; is<end; is+=1) {//TODO: Check removal
		ucoord ih = bgq_surface2halfvolume(isOdd, is);
		ucoord ic = bgq_surface2collapsed(is);
		ucoord t1 = bgq_halfvolume2t1(isOdd,ih);
		ucoord t2 = bgq_halfvolume2t2(isOdd,ih);
		ucoord x = bgq_halfvolume2x(ih);
		ucoord y = bgq_halfvolume2y(ih);
		ucoord z = bgq_halfvolume2z(ih);

		//TODO: Check strength reduction
		bgq_spinorsite *spinorsite = &spinorfield->sec_fullspinor_surface[is];
		bgq_gaugesite *gaugesite = &g_bgq_gaugefield_fromCollapsed[isOdd][is];
		bgq_weyl_ptr_t *destptrs = &targetfield->sendptr[ic];

		//TODO: prefetching
		//TODO: Check inlining
		bgq_su3_spinor_decl(spinor);
		bgq_HoppingMatrix_loadFulllayout(spinor, spinorsite, t1, t2, x, y, z);

		//bgq_HoppingMatrix_compute_storeWeyllayout(destptrs, gaugesite, spinor, t1, t2, x, y, z, false);
		bgq_HoppingMatrix_compute_storeWeyllayout_alldir(destptrs, gaugesite, spinor, t1, t2, x, y, z, qka0,qka1,qka2,qka3,kamul);
//#define BGQ_COMPUTEWEYL_INC_
//#include "bgq_ComputeWeyl.inc.c"
//#undef BGQ_COMPUTEWEYL_INC_
	}
}


static void bgq_HoppingMatrix_kamul_worker_surface_precomm_readFulllayout(void *arg, size_t tid, size_t threads) {
	bgq_HoppingMatrix_workload *work = arg;
	bool isOdd = work->isOdd_src;
	bgq_weylfield_controlblock *spinorfield = work->spinorfield;
	bgq_weylfield_controlblock *targetfield = work->targetfield;
	const bool kamul = true;

	bgq_vector4double_decl(qka0);
	bgq_complxval_splat(qka0,ka0);
	bgq_vector4double_decl(qka1);
	bgq_complxval_splat(qka1,ka1);
	bgq_vector4double_decl(qka2);
	bgq_complxval_splat(qka2,ka2);
	bgq_vector4double_decl(qka3);
	bgq_complxval_splat(qka3,ka3);

	const size_t workload = PHYSICAL_SURFACE;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);
	for (ucoord is = begin; is<end; is+=1) {//TODO: Check removal
		ucoord ih = bgq_surface2halfvolume(isOdd, is);
		ucoord ic = bgq_surface2collapsed(is);
		ucoord t1 = bgq_halfvolume2t1(isOdd,ih);
		ucoord t2 = bgq_halfvolume2t2(isOdd,ih);
		ucoord x = bgq_halfvolume2x(ih);
		ucoord y = bgq_halfvolume2y(ih);
		ucoord z = bgq_halfvolume2z(ih);

		//TODO: Check strength reduction
		bgq_spinorsite *spinorsite = &spinorfield->sec_fullspinor_surface[is];
		bgq_gaugesite *gaugesite = &g_bgq_gaugefield_fromCollapsed[isOdd][is];
		bgq_weyl_ptr_t *destptrs = &targetfield->sendptr[ic];

		//TODO: prefetching
		//TODO: Check inlining
		bgq_su3_spinor_decl(spinor);
		bgq_HoppingMatrix_loadFulllayout(spinor, spinorsite, t1, t2, x, y, z);
		//bgq_HoppingMatrix_compute_storeWeyllayout(destptrs, gaugesite, spinor, t1, t2, x, y, z, true);
		bgq_HoppingMatrix_compute_storeWeyllayout_alldir(destptrs, gaugesite, spinor, t1, t2, x, y, z, qka0,qka1,qka2,qka3,kamul);
	}
}


static void bgq_HoppingMatrix_worker_surface_precomm_readWeyllayout(void *arg, size_t tid, size_t threads) {
	bgq_HoppingMatrix_workload *work = arg;
	bool isOdd = work->isOdd_src;
	bgq_weylfield_controlblock *spinorfield = work->spinorfield;
	bgq_weylfield_controlblock *targetfield = work->targetfield;
	const bool kamul = false;

	bgq_vector4double_decl(qka0);
	bgq_complxval_splat(qka0,ka0);
	bgq_vector4double_decl(qka1);
	bgq_complxval_splat(qka1,ka1);
	bgq_vector4double_decl(qka2);
	bgq_complxval_splat(qka2,ka2);
	bgq_vector4double_decl(qka3);
	bgq_complxval_splat(qka3,ka3);

	const size_t workload = PHYSICAL_SURFACE;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);
	for (size_t is = begin; is<end; is+=1) {
		//TODO: Check optaway
		size_t ih = bgq_surface2halfvolume(isOdd,is);
		ucoord ic = bgq_surface2collapsed(is);
		size_t t1 = bgq_halfvolume2t1(isOdd,ih);
		size_t t2 = bgq_halfvolume2t2(isOdd,ih);
		size_t x = bgq_halfvolume2x(ih);
		size_t y = bgq_halfvolume2y(ih);
		size_t z = bgq_halfvolume2z(ih);

		//TODO: Check strength reduction
		bgq_weylsite *weylsite = &spinorfield->sec_collapsed[ic];
		bgq_gaugesite *gaugesite = &g_bgq_gaugefield_fromCollapsed[isOdd][ic];
		bgq_weyl_ptr_t *destptrs = &targetfield->sendptr[ic];

		//TODO: prefetching
		//TODO: Check inlining
		bgq_su3_spinor_decl(spinor);
		bgq_HoppingMatrix_loadWeyllayout(spinor, weylsite, t1, t2, x, y, z);
		//bgq_HoppingMatrix_compute_storeWeyllayout(destptrs, gaugesite, spinor, t1, t2, x, y, z, false);
		bgq_HoppingMatrix_compute_storeWeyllayout_alldir(destptrs, gaugesite, spinor, t1, t2, x, y, z, qka0,qka1,qka2,qka3,kamul);
	}
}

static void bgq_HoppingMatrix_kamul_worker_surface_precomm_readWeyllayout(void *arg, size_t tid, size_t threads) {
	bgq_HoppingMatrix_workload *work = arg;
	bool isOdd = work->isOdd_src;
	bgq_weylfield_controlblock *spinorfield = work->spinorfield;
	bgq_weylfield_controlblock *targetfield = work->targetfield;
	const bool kamul = true;

	bgq_vector4double_decl(qka0);
	bgq_complxval_splat(qka0,ka0);
	bgq_vector4double_decl(qka1);
	bgq_complxval_splat(qka1,ka1);
	bgq_vector4double_decl(qka2);
	bgq_complxval_splat(qka2,ka2);
	bgq_vector4double_decl(qka3);
	bgq_complxval_splat(qka3,ka3);

	const size_t workload = PHYSICAL_SURFACE;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);
	for (size_t is = begin; is<end; is+=1) {
		//TODO: Check optaway
		size_t ih = bgq_surface2halfvolume(isOdd,is);
		ucoord ic = bgq_surface2collapsed(is);
		size_t t1 = bgq_halfvolume2t1(isOdd,ih);
		size_t t2 = bgq_halfvolume2t2(isOdd,ih);
		size_t x = bgq_halfvolume2x(ih);
		size_t y = bgq_halfvolume2y(ih);
		size_t z = bgq_halfvolume2z(ih);

		//TODO: Check strength reduction
		bgq_weylsite *weylsite = &spinorfield->sec_collapsed[ic];
		bgq_gaugesite *gaugesite = &g_bgq_gaugefield_fromCollapsed[isOdd][ic];
		bgq_weyl_ptr_t *destptrs = &targetfield->sendptr[ic];

		//TODO: prefetching
		//TODO: Check inlining
		bgq_su3_spinor_decl(spinor);
		bgq_HoppingMatrix_loadWeyllayout(spinor, weylsite, t1, t2, x, y, z);
		//bgq_HoppingMatrix_compute_storeWeyllayout(destptrs, gaugesite, spinor, t1, t2, x, y, z, true);
		bgq_HoppingMatrix_compute_storeWeyllayout_alldir(destptrs, gaugesite, spinor, t1, t2, x, y, z, qka0,qka1,qka2,qka3,kamul);
	}
}



static void bgq_HoppingMatrix_worker_body_readFulllayout(void *arg, size_t tid, size_t threads) {
	bgq_HoppingMatrix_workload *work = arg;
	bool isOdd = work->isOdd_src;
	bgq_weylfield_controlblock *spinorfield = work->spinorfield;
	bgq_weylfield_controlblock *targetfield = work->targetfield;
	const bool kamul = false;

	bgq_vector4double_decl(qka0);
	bgq_complxval_splat(qka0,ka0);
	bgq_vector4double_decl(qka1);
	bgq_complxval_splat(qka1,ka1);
	bgq_vector4double_decl(qka2);
	bgq_complxval_splat(qka2,ka2);
	bgq_vector4double_decl(qka3);
	bgq_complxval_splat(qka3,ka3);

	const size_t workload = PHYSICAL_BODY;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);
	for (size_t ib = begin; ib<end; ib+=1) {
		//TODO: Check optaway
		size_t ih = bgq_body2halfvolume(isOdd, ib);
		ucoord ic = bgq_body2collapsed(ib);
		size_t t1 = bgq_halfvolume2t1(isOdd,ih);
		size_t t2 = bgq_halfvolume2t2(isOdd,ih);
		size_t x = bgq_halfvolume2x(ih);
		size_t y = bgq_halfvolume2y(ih);
		size_t z = bgq_halfvolume2z(ih);

		//TODO: Check strength reduction
		bgq_spinorsite *spinorsite = &spinorfield->sec_fullspinor[ic];
		bgq_gaugesite *gaugesite = &g_bgq_gaugefield_fromCollapsed[isOdd][ic];
		bgq_weyl_ptr_t *destptrs = &targetfield->sendptr[ic];

		//TODO: prefetching
		//TODO: Check inlining
		bgq_su3_spinor_decl(spinor);
		bgq_HoppingMatrix_loadFulllayout(spinor, spinorsite, t1, t2, x, y, z);
		//bgq_HoppingMatrix_compute_storeWeyllayout(destptrs, gaugesite, spinor, t1, t2, x, y, z, false);
		bgq_HoppingMatrix_compute_storeWeyllayout_alldir(destptrs, gaugesite, spinor, t1, t2, x, y, z, qka0,qka1,qka2,qka3,kamul);
	}
}




static void bgq_HoppingMatrix_kamul_worker_body_readFulllayout(void *arg, size_t tid, size_t threads) {
	bgq_HoppingMatrix_workload *work = arg;
	bool isOdd = work->isOdd_src;
	bgq_weylfield_controlblock *spinorfield = work->spinorfield;
	bgq_weylfield_controlblock *targetfield = work->targetfield;
	const bool kamul = true;

	bgq_vector4double_decl(qka0);
	bgq_complxval_splat(qka0,ka0);
	bgq_vector4double_decl(qka1);
	bgq_complxval_splat(qka1,ka1);
	bgq_vector4double_decl(qka2);
	bgq_complxval_splat(qka2,ka2);
	bgq_vector4double_decl(qka3);
	bgq_complxval_splat(qka3,ka3);

	const size_t workload = PHYSICAL_BODY;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);
	for (size_t ib = begin; ib<end; ib+=1) {
		//TODO: Check optaway
		size_t ih = bgq_body2halfvolume(isOdd, ib);
		ucoord ic = bgq_body2collapsed(ib);
		size_t t1 = bgq_halfvolume2t1(isOdd,ih);
		size_t t2 = bgq_halfvolume2t2(isOdd,ih);
		size_t x = bgq_halfvolume2x(ih);
		size_t y = bgq_halfvolume2y(ih);
		size_t z = bgq_halfvolume2z(ih);

		//TODO: Check strength reduction
		bgq_spinorsite *spinorsite = &spinorfield->sec_fullspinor[ic];
		bgq_gaugesite *gaugesite = &g_bgq_gaugefield_fromCollapsed[isOdd][ic];
		bgq_weyl_ptr_t *destptrs = &targetfield->sendptr[ic];

		//TODO: prefetching
		//TODO: Check inlining
		bgq_su3_spinor_decl(spinor);
		bgq_HoppingMatrix_loadFulllayout(spinor, spinorsite, t1, t2, x, y, z);
		//bgq_HoppingMatrix_compute_storeWeyllayout(destptrs, gaugesite, spinor, t1, t2, x, y, z, true);
		bgq_HoppingMatrix_compute_storeWeyllayout_alldir(destptrs, gaugesite, spinor, t1, t2, x, y, z, qka0,qka1,qka2,qka3,kamul);
	}
}

static void bgq_HoppingMatrix_worker_body_readWeyllayout(void *arg, size_t tid, size_t threads) {
	bgq_HoppingMatrix_workload *work = arg;
	bool isOdd = work->isOdd_src;
	bgq_weylfield_controlblock *spinorfield = work->spinorfield;
	bgq_weylfield_controlblock *targetfield = work->targetfield;
	const bool kamul = false;

	bgq_vector4double_decl(qka0);
	bgq_complxval_splat(qka0,ka0);
	bgq_vector4double_decl(qka1);
	bgq_complxval_splat(qka1,ka1);
	bgq_vector4double_decl(qka2);
	bgq_complxval_splat(qka2,ka2);
	bgq_vector4double_decl(qka3);
	bgq_complxval_splat(qka3,ka3);

	const size_t workload = PHYSICAL_BODY;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);
	for (size_t ib = begin; ib<end; ib+=1) {
		//TODO: Check optaway
		size_t ih = bgq_body2halfvolume(isOdd,ib);
		ucoord ic = bgq_body2collapsed(ib);
		size_t t1 = bgq_halfvolume2t1(isOdd,ih);
		size_t t2 = bgq_halfvolume2t2(isOdd,ih);
		size_t x = bgq_halfvolume2x(ih);
		size_t y = bgq_halfvolume2y(ih);
		size_t z = bgq_halfvolume2z(ih);

		//TODO: Check strength reduction
		bgq_weylsite *weylsite = &spinorfield->sec_body[ib];
		bgq_gaugesite *gaugesite = &g_bgq_gaugefield_fromCollapsed[isOdd][ic];
		bgq_weyl_ptr_t *destptrs = &targetfield->sendptr[ic];

		//TODO: prefetching
		//TODO: Check inlining
		bgq_su3_spinor_decl(spinor);
		bgq_HoppingMatrix_loadWeyllayout(spinor, weylsite, t1, t2, x, y, z);
		//bgq_HoppingMatrix_compute_storeWeyllayout(destptrs, gaugesite, spinor, t1, t2, x, y, z, false);
		bgq_HoppingMatrix_compute_storeWeyllayout_alldir(destptrs, gaugesite, spinor, t1, t2, x, y, z, qka0,qka1,qka2,qka3,kamul);
	}
}

static void bgq_HoppingMatrix_kamul_worker_body_readWeyllayout(void *arg, size_t tid, size_t threads) {
	bgq_HoppingMatrix_workload *work = arg;
	bool isOdd = work->isOdd_src;
	bgq_weylfield_controlblock *spinorfield = work->spinorfield;
	bgq_weylfield_controlblock *targetfield = work->targetfield;
	const bool kamul = true;

	bgq_vector4double_decl(qka0);
	bgq_complxval_splat(qka0,ka0);
	bgq_vector4double_decl(qka1);
	bgq_complxval_splat(qka1,ka1);
	bgq_vector4double_decl(qka2);
	bgq_complxval_splat(qka2,ka2);
	bgq_vector4double_decl(qka3);
	bgq_complxval_splat(qka3,ka3);

	const size_t workload = PHYSICAL_BODY;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);
	for (size_t ib = begin; ib<end; ib+=1) {
		//TODO: Check optaway
		size_t ih = bgq_body2halfvolume(isOdd,ib);
		ucoord ic = bgq_body2collapsed(ib);
		size_t t1 = bgq_halfvolume2t1(isOdd,ih);
		size_t t2 = bgq_halfvolume2t2(isOdd,ih);
		size_t x = bgq_halfvolume2x(ih);
		size_t y = bgq_halfvolume2y(ih);
		size_t z = bgq_halfvolume2z(ih);

		//TODO: Check strength reduction
		bgq_weylsite *weylsite = &spinorfield->sec_body[ib];
		bgq_gaugesite *gaugesite = &g_bgq_gaugefield_fromCollapsed[isOdd][ic];
		bgq_weyl_ptr_t *destptrs = &targetfield->sendptr[ic];

		//TODO: prefetching
		//TODO: Check inlining
		bgq_su3_spinor_decl(spinor);
		bgq_HoppingMatrix_loadWeyllayout(spinor, weylsite, t1, t2, x, y, z);
		//bgq_HoppingMatrix_compute_storeWeyllayout(destptrs, gaugesite, spinor, t1, t2, x, y, z, true);
		bgq_HoppingMatrix_compute_storeWeyllayout_alldir(destptrs, gaugesite, spinor, t1, t2, x, y, z, qka0,qka1,qka2,qka3,kamul);
	}
}


void bgq_HoppingMatrix(bool isOdd, bgq_weylfield_controlblock *targetfield, bgq_weylfield_controlblock *spinorfield, bgq_hmflags opts) {
	assert(targetfield);

	assert(spinorfield);
	assert(spinorfield->isInitinialized);
	assert(spinorfield->isSloppy == false);
	assert(spinorfield->isOdd == isOdd);
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

	bgq_spinorfield_setup(targetfield, isOdd, false, false, false, true);
	bgq_spinorfield_setup(spinorfield, !isOdd, readFullspinor, false, !readFullspinor, false);
	assert(targetfield->isOdd == isOdd);


	// 0. Expect data from other neighbor node
	if (!nocomm) {
		bgq_comm_recv(nospi);
	}


	// 1. Distribute
	bgq_master_sync();
	static bgq_HoppingMatrix_workload work_surface;
	work_surface.isOdd_src = !isOdd;
	work_surface.isOdd_dst = isOdd;
	work_surface.targetfield = targetfield;
	work_surface.spinorfield = spinorfield;
	work_surface.ic_begin=bgq_surface2collapsed(0);
	work_surface.ic_end=bgq_surface2collapsed(PHYSICAL_SURFACE);

	static bgq_HoppingMatrix_workload work_body ;
	work_body.isOdd_src = !isOdd;
	work_body.isOdd_dst = isOdd;
	work_body.targetfield = targetfield;
	work_body.spinorfield = spinorfield;
	work_body.ic_begin=bgq_body2collapsed(0);
	work_body.ic_end=bgq_body2collapsed(PHYSICAL_BODY);

	// Compute surface and put data into the send buffers
	if (!nodistribute) {
		if (readFullspinor) {
			if (nokamul)
				bgq_master_call(&bgq_HoppingMatrix_nokamul_worker_readFulllayout, &work_surface);
			else
				bgq_master_call(&bgq_HoppingMatrix_kamul_worker_readFulllayout, &work_surface);
		} else {
			// readWeyl
			if (nokamul)
				bgq_master_call(&bgq_HoppingMatrix_nokamul_worker_readWeyllayout, &work_surface);
			else
				bgq_master_call(&bgq_HoppingMatrix_kamul_worker_readWeyllayout, &work_surface);
		}
	}


// 2. Start communication
	if (!nocomm) {
		//TODO: ensure there are no other communications pending
		bgq_comm_send(nospi);
		targetfield->waitingForRecv = true;
	}
	if (!nodatamove) {
		targetfield->pendingDatamove = true;
	}
	targetfield->hmflags = opts;

	if (nooverlap) {
		// Do not wait until data is required, but do it here
		bgq_spinorfield_setup(targetfield, isOdd, false, false, true, false);
	}


// 3. Compute the body
	if (!nobody) {
		if (readFullspinor) {
			if (nokamul)
				bgq_master_call(&bgq_HoppingMatrix_nokamul_worker_readFulllayout, &work_body);
			else
				bgq_master_call(&bgq_HoppingMatrix_kamul_worker_readFulllayout, &work_body);
		} else {
			if (nokamul)
				bgq_master_call(&bgq_HoppingMatrix_nokamul_worker_readWeyllayout, &work_body);
			else
				bgq_master_call(&bgq_HoppingMatrix_kamul_worker_readWeyllayout, &work_body);
		}
	}


// 4. Wait for the communication to finish
	/* Defer to bgq_spinorfield_setup as soon as the data is actually required */


// 5. Move received to correct location
	/* Done in bgq_spinorfield_setup whoever is using the field next*/


// 6. Compute the surface
	/* Done by funcs calling readWeyllayout */
}




