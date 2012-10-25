/*
 * bgq_HoppingMatrix.h
 *
 *  Created on: Oct 13, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_HOPPINGMATRIX_H_
#define BGQ_HOPPINGMATRIX_H_

#include "bgq_field.h"
#include "bgq_qpx.h"
#include "bgq_utils.h"

#include <stdbool.h>

#ifndef BGQ_HOPPING_MATRIX_C_
#define EXTERN_INLINE EXTERN_INLINE_DECLARATION
#define EXTERN_FIELD extern
#define EXTERN_INIT(val)
#else
#define EXTERN_INLINE EXTERN_INLINE_DEFINITION
#define EXTERN_FIELD
#define EXTERN_INIT(val) = (val)
#endif


typedef enum {
	hm_nocom = 1 << 0,
	hm_nooverlap = 1 << 1,
	hm_nokamul = 1 << 2,
	hm_fixedoddness = 1 << 3,

	hm_noprefetchexplicit = 1 << 4,
	hm_noprefetchlist = 1 << 5,
	hm_noprefetchstream = 1 << 6,

	hm_noweylsend = 1 << 7,
	hm_nobody = 1 << 8,
	hm_nosurface = 1 << 9,

	hm_l1pnonstoprecord = 1 << 10,
	hm_experimental = 1 << 11,

	hm_prefetchimplicitdisable = 1 << 12,
	hm_prefetchimplicitoptimistic = 2 << 12,
	hm_prefetchimplicitconfirmed = 3 << 12
} bgq_hmflags;

void bgq_HoppingMatrix(bool isOdd, bgq_weylfield_controlblock *targetfield, bgq_weylfield_controlblock *spinorfield, bgq_hmflags opts);


EXTERN_INLINE void bgq_HoppingMatrix_loadFulllayout_raw(bgq_su3_spinor_params(*target), bgq_spinorsite *addr) {
	bgq_su3_spinor_load_double(*target,addr);
}


EXTERN_INLINE void bgq_HoppingMatrix_loadWeyllayout_raw(bgq_su3_spinor_params(*target), bgq_weylsite *spinorsite) {
	bgq_su3_spinor_decl(result);

	// T+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tup);
		bgq_su3_weyl_load_double(weyl_tup, &spinorsite->d[TUP]);
		bgq_su3_expand_weyl_tup(result, weyl_tup);
	}

	// T- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tdown);
		bgq_su3_weyl_load_double(weyl_tdown, &spinorsite->d[TDOWN]);
		bgq_su3_accum_weyl_tdown(result, weyl_tdown);
	}

	// X+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xup);
		bgq_su3_weyl_load_double(weyl_xup, &spinorsite->d[XUP]);
		bgq_su3_accum_weyl_xup(result, weyl_xup);
	}

	// X- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xdown);
		bgq_su3_weyl_load_double(weyl_xdown, &spinorsite->d[XDOWN]);
		bgq_su3_accum_weyl_xdown(result, weyl_xdown);
	}

	// Y+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_yup);
		bgq_su3_weyl_load_double(weyl_yup, &spinorsite->d[YUP]);
		bgq_su3_accum_weyl_yup(result, weyl_yup);
	}

	// Y- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_ydown);
		bgq_su3_weyl_load_double(weyl_ydown, &spinorsite->d[YDOWN]);
		bgq_su3_accum_weyl_ydown(result, weyl_ydown);
	}

	// Z+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zup);
		bgq_su3_weyl_load_double(weyl_zup, &spinorsite->d[ZUP]);
		bgq_su3_accum_weyl_zup(result, weyl_zup);
	}

	// Z- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zdown);
		bgq_su3_weyl_load_double(weyl_zdown, &spinorsite->d[ZDOWN]);
		bgq_su3_accum_weyl_zdown(result, weyl_zdown);
	}

	bgq_su3_spinor_mov(*target, result);
}


EXTERN_INLINE void bgq_HoppingMatrix_storeFulllayout_raw(bgq_spinorsite *targetptr, bgq_su3_spinor_params(spinor)) {
	//TODO: Callers should care for prefetching
	bgq_su3_spinor_store_double(targetptr, spinor);
}


EXTERN_INLINE void bgq_HoppingMatrix_storeWeyllayout_raw(bgq_weyl_ptr_t *targetptrs, bgq_su3_spinor_params(spinor)) {
	//TODO: prefetch targetptrs

	// T+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tup);
		bgq_su3_reduce_weyl_tup(weyl_tup, spinor);
		bgq_su3_weyl_store_double(targetptrs->d[TUP], weyl_tup);
	}

	// T- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tdown);
		bgq_su3_reduce_weyl_tdown(weyl_tdown, spinor);
		bgq_su3_weyl_store_double(targetptrs->d[TDOWN], weyl_tdown);
	}

	// X+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xup);
		bgq_su3_reduce_weyl_xup(weyl_xup, spinor);
		bgq_su3_weyl_store_double(targetptrs->d[XUP], weyl_xup);
	}

	// X- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xdown);
		bgq_su3_reduce_weyl_xdown(weyl_xdown, spinor);
		bgq_su3_weyl_store_double(targetptrs->d[XDOWN], weyl_xdown);
	}

	// Y+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_yup);
		bgq_su3_reduce_weyl_yup(weyl_yup, spinor);
		bgq_su3_weyl_store_double(targetptrs->d[YUP], weyl_yup);
	}

	// Y- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_ydown);
		bgq_su3_reduce_weyl_ydown(weyl_ydown, spinor);
		bgq_su3_weyl_store_double(targetptrs->d[YDOWN], weyl_ydown);
	}

	// Z+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zup);
		bgq_su3_reduce_weyl_zup(weyl_zup, spinor);
		bgq_su3_weyl_store_double(targetptrs->d[ZUP], weyl_zup);
	}

	// Z- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zdown);
		bgq_su3_reduce_weyl_zdown(weyl_zdown, spinor);
		bgq_su3_weyl_store_double(targetptrs->d[ZDOWN], weyl_zdown);
	}
}


EXTERN_INLINE void bgq_HoppingMatrix_compute_storeWeyllayout_raw(bgq_weyl_ptr_t *targetptrs, bgq_gaugesite *gaugesite, bgq_su3_spinor_params(spinor)) {
	//TODO: prefetch targetptrs

	// T+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tup);
		bgq_su3_reduce_weyl_tup(weyl_tup, spinor);

		bgq_su3_mdecl(gauge_tup);
		bgq_su3_matrix_load_double(gauge_tup, &gaugesite->su3[TUP]);
		bgq_su3_weyl_mvmul(weyl_tup, gauge_tup, weyl_tup);

		bgq_su3_weyl_store_double(targetptrs->d[TUP], weyl_tup);
	}

	// T- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tdown);
		bgq_su3_reduce_weyl_tdown(weyl_tdown, spinor);

		bgq_su3_mdecl(gauge_tdown);
		bgq_su3_matrix_load_double(gauge_tdown, &gaugesite->su3[TDOWN]);
		bgq_su3_weyl_mvmul(weyl_tdown, gauge_tdown, weyl_tdown);

		bgq_su3_weyl_store_double(targetptrs->d[TDOWN], weyl_tdown);
	}

	// X+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xup);
		bgq_su3_reduce_weyl_xup(weyl_xup, spinor);

		bgq_su3_mdecl(gauge_xup);
		bgq_su3_matrix_load_double(gauge_xup, &gaugesite->su3[XUP]);
		bgq_su3_weyl_mvmul(weyl_xup, gauge_xup, weyl_xup);

		bgq_su3_weyl_store_double(targetptrs->d[XUP], weyl_xup);
	}

	// X- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xdown);
		bgq_su3_reduce_weyl_xdown(weyl_xdown, spinor);

		bgq_su3_mdecl(gauge_xdown);
		bgq_su3_matrix_load_double(gauge_xdown, &gaugesite->su3[XDOWN]);
		bgq_su3_weyl_mvmul(weyl_xdown, gauge_xdown, weyl_xdown);

		bgq_su3_weyl_store_double(targetptrs->d[XDOWN], weyl_xdown);
	}

	// Y+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_yup);
		bgq_su3_reduce_weyl_yup(weyl_yup, spinor);

		bgq_su3_mdecl(gauge_yup);
		bgq_su3_matrix_load_double(gauge_yup, &gaugesite->su3[YUP]);
		bgq_su3_weyl_mvmul(weyl_yup, gauge_yup, weyl_yup);

		bgq_su3_weyl_store_double(targetptrs->d[YUP], weyl_yup);
	}

	// Y- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_ydown);
		bgq_su3_reduce_weyl_ydown(weyl_ydown, spinor);

		bgq_su3_mdecl(gauge_ydown);
		bgq_su3_matrix_load_double(gauge_ydown, &gaugesite->su3[YDOWN]);
		bgq_su3_weyl_mvmul(weyl_ydown, gauge_ydown, weyl_ydown);

		bgq_su3_weyl_store_double(targetptrs->d[YDOWN], weyl_ydown);
	}

	// Z+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zup);
		bgq_su3_reduce_weyl_zup(weyl_zup, spinor);

		bgq_su3_mdecl(gauge_zup);
		bgq_su3_matrix_load_double(gauge_zup, &gaugesite->su3[ZUP]);
		bgq_su3_weyl_mvmul(weyl_zup, gauge_zup, weyl_zup);

		bgq_su3_weyl_store_double(targetptrs->d[ZUP], weyl_zup);
	}

	// Z- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zdown);
		bgq_su3_reduce_weyl_zdown(weyl_zdown, spinor);

		bgq_su3_mdecl(gauge_zdown);
		bgq_su3_matrix_load_double(gauge_zdown, &gaugesite->su3[ZDOWN]);
		bgq_su3_weyl_mvmul(weyl_zdown, gauge_zdown, weyl_zdown);

		bgq_su3_weyl_store_double(targetptrs->d[ZDOWN], weyl_zdown);
	}
}

static void bgq_spinorfield_weyl_store_fromHalfvolume(bgq_weylfield_controlblock *targetfield, bool isOdd, size_t ih, bgq_su3_spinor_params(spinor)) {
	//bgq_spinorfield_reset(targetfield, isOdd, true, false);
	assert(targetfield->isInitinialized);
	assert(targetfield->isOdd == isOdd);

	bgq_weyl_ptr_t *weylptrs = &targetfield->destptrFromHalfvolume[ih]; // TODO: Check that compiler does strength reduction after inline, otherwise do manually
	for (size_t d = 0; d < PHYSICAL_LD; d+=1){
		assert((bgq_weyl_vec*)targetfield->sec_weyl <= weylptrs->d[d] && weylptrs->d[d] < (bgq_weyl_vec*)targetfield->sec_end);
	}
	//TODO: probably compiler will li an offset for every destptrFromHalfvolume, can do better using addi

	bgq_HoppingMatrix_storeWeyllayout_raw(weylptrs, bgq_su3_spinor_vars(spinor));
}






#undef EXTERN_INLINE
#undef EXTERN_FIELD
#undef EXTERN_INIT

#endif /* BGQ_HOPPINGMATRIX_H_ */
