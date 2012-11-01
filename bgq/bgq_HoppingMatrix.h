/*
 * bgq_HoppingMatrix.h
 *
 *  Created on: Oct 13, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_HOPPINGMATRIX_H_
#define BGQ_HOPPINGMATRIX_H_

#include "bgq_field.h"
#include "bgq_spinorfield.h"
#include "bgq_gaugefield.h"
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
	hm_prefetchimplicitconfirmed = 3 << 12,

	hm_withcheck = 1 << 14
} bgq_hmflags;

void bgq_HoppingMatrix(bool isOdd, bgq_weylfield_controlblock *targetfield, bgq_weylfield_controlblock *spinorfield, bgq_hmflags opts);


#define bgq_HoppingMatrix_loadFulllayout(spinor,addr,t_left,t_right,x,y,z) bgq_HoppingMatrix_loadFulllayout_raw(bgq_su3_spinor_vars(&spinor),addr,t_left,t_right,x,y,z)
EXTERN_INLINE void bgq_HoppingMatrix_loadFulllayout_raw(bgq_su3_spinor_params(*target), bgq_spinorsite *addr, size_t t1, size_t t2, size_t x, size_t y, size_t z) {
	bgq_su3_spinor_load_double(*target, addr);
	bgq_spinorqpx_expect(*target, t1, t2, x, y, z);
}


#define bgq_HoppingMatrix_loadWeyllayout(target, spinorsite, t_left, t_right, x,y,z) bgq_HoppingMatrix_loadWeyllayout_raw(bgq_su3_spinor_vars(&target), spinorsite, t_left, t_right, x,y,z)
EXTERN_INLINE void bgq_HoppingMatrix_loadWeyllayout_raw(bgq_su3_spinor_params(*target), bgq_weylsite *spinorsite, ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z) {
	bgq_su3_spinor_decl(result);

	// T+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tup);
		bgq_su3_weyl_load_double(weyl_tup, &spinorsite->d[TUP]);
		bgq_weylqpx_expect(weyl_tup, t1, t2, x, y, z, TUP, false);
		bgq_su3_expand_weyl_tup(result, weyl_tup);
	}

	// T- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tdown);
		bgq_su3_weyl_load_double(weyl_tdown, &spinorsite->d[TDOWN]);
		bgq_weylqpx_expect(weyl_tdown, t1, t2, x, y, z, TDOWN, false);
		bgq_su3_accum_weyl_tdown(result, weyl_tdown);
	}

	// X+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xup);
		bgq_su3_weyl_load_double(weyl_xup, &spinorsite->d[XUP]);
		bgq_weylqpx_expect(weyl_xup, t1, t2, x, y, z, XUP,false);
		bgq_su3_accum_weyl_xup(result, weyl_xup);
	}

	// X- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xdown);
		bgq_su3_weyl_load_double(weyl_xdown, &spinorsite->d[XDOWN]);
		bgq_weylqpx_expect(weyl_xdown, t1, t2, x, y, z, XDOWN, false);
		bgq_su3_accum_weyl_xdown(result, weyl_xdown);
	}

	// Y+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_yup);
		bgq_su3_weyl_load_double(weyl_yup, &spinorsite->d[YUP]);
		bgq_weylqpx_expect(weyl_yup, t1, t2, x, y, z, YUP,false);
		bgq_su3_accum_weyl_yup(result, weyl_yup);
	}

	// Y- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_ydown);
		bgq_su3_weyl_load_double(weyl_ydown, &spinorsite->d[YDOWN]);
		bgq_weylqpx_expect(weyl_ydown, t1, t2, x, y, z, YDOWN,false);
		bgq_su3_accum_weyl_ydown(result, weyl_ydown);
	}

	// Z+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zup);
		bgq_su3_weyl_load_double(weyl_zup, &spinorsite->d[ZUP]);
		bgq_weylqpx_expect(weyl_zup, t1, t2, x, y, z, ZUP,false);
		bgq_su3_accum_weyl_zup(result, weyl_zup);
	}

	// Z- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zdown);
		bgq_su3_weyl_load_double(weyl_zdown, &spinorsite->d[ZDOWN]);
		bgq_weylqpx_expect(weyl_zdown, t1, t2, x, y, z, ZDOWN,false);
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

#define bgq_HoppingMatrix_compute_storeWeyllayout(targetptrs,gaugesite,spinor,t1,t2,x,y,z) bgq_HoppingMatrix_compute_storeWeyllayout_raw(targetptrs,gaugesite,bgq_su3_spinor_vars(spinor),t1,t2,x,y,z)
EXTERN_INLINE void bgq_HoppingMatrix_compute_storeWeyllayout_raw(bgq_weyl_ptr_t *targetptrs, bgq_gaugesite *gaugesite, bgq_su3_spinor_params(spinor), ucoord t1,ucoord t2,ucoord x,ucoord y,ucoord z) {
	//TODO: prefetch targetptrs
	//bgq_spinorqpx_expect(spinor, t1, t2, x, y, z);

	// T+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tup);
		bgq_su3_reduce_weyl_tup(weyl_tup, spinor);

		bgq_su3_mdecl(gauge_tup);
		bgq_su3_matrix_load_double(gauge_tup, &gaugesite->su3[TUP]);
		bgq_gaugeqpx_expect(gauge_tup, t1, t2, x, y, z, TUP, true);
		bgq_su3_weyl_mvmul(weyl_tup, gauge_tup, weyl_tup);

		//bgq_weylvec_expect(*targetptrs->d[TUP], t1, t2, x, y, z, TUP, true);
		bgq_su3_weyl_store_double(targetptrs->d[TUP], weyl_tup);
		bgq_weylvec_written(targetptrs->d[TUP], t1, t2, x, y, z, TUP, true);
	}

	// T- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tdown);
		bgq_su3_reduce_weyl_tdown(weyl_tdown, spinor);

		bgq_su3_mdecl(gauge_tdown);
		bgq_su3_matrix_load_double(gauge_tdown, &gaugesite->su3[TDOWN]);
		bgq_gaugeqpx_expect(gauge_tdown, t1, t2, x, y, z, TDOWN, true);
		bgq_su3_weyl_mvmul(weyl_tdown, gauge_tdown, weyl_tdown);

		//bgq_weylvec_expect(*targetptrs->d[TDOWN], t1, t2, x, y, z, TDOWN, true);
		bgq_su3_weyl_store_double(targetptrs->d[TDOWN], weyl_tdown);
		bgq_weylvec_written(targetptrs->d[TDOWN], t1, t2, x, y, z, TDOWN, true);
	}

	// X+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xup);
		bgq_su3_reduce_weyl_xup(weyl_xup, spinor);

		bgq_su3_mdecl(gauge_xup);
		bgq_su3_matrix_load_double(gauge_xup, &gaugesite->su3[XUP]);
		bgq_gaugeqpx_expect(gauge_xup, t1, t2, x, y, z, XUP, true);
		bgq_su3_weyl_mvmul(weyl_xup, gauge_xup, weyl_xup);

		//bgq_weylvec_expect(*targetptrs->d[XUP], t1, t2, x, y, z, XUP, true);
		bgq_su3_weyl_store_double(targetptrs->d[XUP], weyl_xup);
		bgq_weylvec_written(targetptrs->d[XUP], t1, t2, x,y,z,XUP, true);
	}

	// X- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xdown);
		bgq_su3_reduce_weyl_xdown(weyl_xdown, spinor);

		bgq_su3_mdecl(gauge_xdown);
		bgq_su3_matrix_load_double(gauge_xdown, &gaugesite->su3[XDOWN]);
		bgq_gaugeqpx_expect(gauge_xdown, t1, t2, x, y, z, XDOWN, true);
		bgq_su3_weyl_mvmul(weyl_xdown, gauge_xdown, weyl_xdown);

		//bgq_weylvec_expect(*targetptrs->d[XDOWN], t1, t2, x, y, z, XDOWN, true);
		bgq_su3_weyl_store_double(targetptrs->d[XDOWN], weyl_xdown);
		bgq_weylvec_written(targetptrs->d[XDOWN], t1, t2, x,y,z,XDOWN, true);
	}

	// Y+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_yup);
		bgq_su3_reduce_weyl_yup(weyl_yup, spinor);

		bgq_su3_mdecl(gauge_yup);
		bgq_su3_matrix_load_double(gauge_yup, &gaugesite->su3[YUP]);
		bgq_gaugeqpx_expect(gauge_yup, t1, t2, x, y, z, YUP, true);
		bgq_su3_weyl_mvmul(weyl_yup, gauge_yup, weyl_yup);

		//bgq_weylvec_expect(*targetptrs->d[YUP], t1, t2, x, y, z, YUP, true);
		bgq_su3_weyl_store_double(targetptrs->d[YUP], weyl_yup);
		bgq_weylvec_written(targetptrs->d[YUP], t1, t2, x,y,z,YUP, true);
	}

	// Y- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_ydown);
		bgq_su3_reduce_weyl_ydown(weyl_ydown, spinor);

		bgq_su3_mdecl(gauge_ydown);
		bgq_su3_matrix_load_double(gauge_ydown, &gaugesite->su3[YDOWN]);
		bgq_gaugeqpx_expect(gauge_ydown, t1, t2, x, y, z, YDOWN, true);
		bgq_su3_weyl_mvmul(weyl_ydown, gauge_ydown, weyl_ydown);

		//bgq_weylvec_expect(*targetptrs->d[YDOWN], t1, t2, x, y, z, YDOWN, true);
		bgq_su3_weyl_store_double(targetptrs->d[YDOWN], weyl_ydown);
		bgq_weylvec_written(targetptrs->d[YDOWN], t1, t2, x,y,z,YDOWN, true);
	}

	// Z+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zup);
		bgq_su3_reduce_weyl_zup(weyl_zup, spinor);

		bgq_su3_mdecl(gauge_zup);
		bgq_su3_matrix_load_double(gauge_zup, &gaugesite->su3[ZUP]);
		bgq_gaugeqpx_expect(gauge_zup, t1, t2, x, y, z, ZUP, true);
		bgq_su3_weyl_mvmul(weyl_zup, gauge_zup, weyl_zup);

		//bgq_weylvec_expect(*targetptrs->d[ZUP], t1, t2, x, y, z, ZUP, true);
		bgq_su3_weyl_store_double(targetptrs->d[ZUP], weyl_zup);
		bgq_weylvec_written(targetptrs->d[ZUP], t1, t2, x,y,z,ZUP, true);
	}

	// Z- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zdown);
		bgq_su3_reduce_weyl_zdown(weyl_zdown, spinor);

		bgq_su3_mdecl(gauge_zdown);
		bgq_su3_matrix_load_double(gauge_zdown, &gaugesite->su3[ZDOWN]);
		bgq_gaugeqpx_expect(gauge_zdown, t1, t2, x, y, z, ZDOWN, true);
		bgq_su3_weyl_mvmul(weyl_zdown, gauge_zdown, weyl_zdown);

		//bgq_weylvec_expect(*targetptrs->d[ZDOWN], t1, t2, x, y, z, ZDOWN, true);
		bgq_su3_weyl_store_double(targetptrs->d[ZDOWN], weyl_zdown);
		bgq_weylvec_written(targetptrs->d[ZDOWN], t1, t2, x,y,z,ZDOWN, true);
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
