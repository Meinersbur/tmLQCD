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

#include "../boundary.h"

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


#define bgq_HoppingMatrix_loadWeyllayout(target, spinorsite, t1, t2, x,y,z) bgq_HoppingMatrix_loadWeyllayout_raw(bgq_su3_spinor_vars(&target), spinorsite, t1, t2, x, y, z)
EXTERN_INLINE void bgq_HoppingMatrix_loadWeyllayout_raw(bgq_su3_spinor_params(*target), bgq_weylsite *spinorsite, ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z) {
	bgq_su3_spinor_decl(result);

	// T+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tup);
		bgq_su3_weyl_load_double(weyl_tup, &spinorsite->d[TUP]);
		bgq_weylqpx_expect(weyl_tup, t1, t2, x, y, z, TUP, false);
				bgq_setdesc(BGQREF_TUP_RECV, "BGQREF_TUP_RECV");
				bgq_setbgqvalue(t1, x, y, z, BGQREF_TUP_RECV, bgq_cmplxval1(weyl_tup_v1_c0));
				bgq_setbgqvalue(t2, x, y, z, BGQREF_TUP_RECV, bgq_cmplxval2(weyl_tup_v1_c0));
		bgq_su3_expand_weyl_tup(result, weyl_tup);
	}

	// T- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tdown);
		bgq_su3_weyl_load_double(weyl_tdown, &spinorsite->d[TDOWN]);
		bgq_weylqpx_expect(weyl_tdown, t1, t2, x, y, z, TDOWN, false);
				bgq_setdesc(BGQREF_TDOWN_RECV, "BGQREF_TDOWN_RECV");
				bgq_setbgqvalue(t1, x, y, z, BGQREF_TDOWN_RECV, bgq_cmplxval1(weyl_tdown_v1_c0));
				bgq_setbgqvalue(t2, x, y, z, BGQREF_TDOWN_RECV, bgq_cmplxval2(weyl_tdown_v1_c0));
		bgq_su3_accum_weyl_tdown(result, weyl_tdown);
				bgq_setdesc(BGQREF_TDOWN_ACCUM, "BGQREF_TDOWN_ACCUM");
				bgq_setbgqvalue(t1, x, y, z, BGQREF_TDOWN_ACCUM, bgq_cmplxval1(result_v1_c0));
				bgq_setbgqvalue(t2, x, y, z, BGQREF_TDOWN_ACCUM, bgq_cmplxval2(result_v1_c0));
	}

	// X+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xup);
		bgq_su3_weyl_load_double(weyl_xup, &spinorsite->d[XUP]);
		bgq_weylqpx_expect(weyl_xup, t1, t2, x, y, z, XUP, false);
			bgq_setdesc(BGQREF_XUP_RECV, "BGQREF_XUP_RECV");
			bgq_setbgqvalue(t1, x, y, z, BGQREF_XUP_RECV, bgq_cmplxval1(weyl_xup_v1_c0));
			bgq_setbgqvalue(t2, x, y, z, BGQREF_XUP_RECV, bgq_cmplxval2(weyl_xup_v1_c0));
			ucoord index = bgq_offset2index(bgq_fieldpointer2offset(&spinorsite->d[XUP]));
			assert(bgq_cmplxval1(weyl_xup_v1_c0)!=0);
			assert(bgq_cmplxval2(weyl_xup_v1_c0)!=0);
		bgq_su3_accum_weyl_xup(result, weyl_xup);
			bgq_setdesc(BGQREF_XUP_ACCUM, "BGQREF_XUP_ACCUM");
			bgq_setbgqvalue(t1, x, y, z, BGQREF_XUP_ACCUM, bgq_cmplxval1(result_v1_c0));
			bgq_setbgqvalue(t2, x, y, z, BGQREF_XUP_ACCUM, bgq_cmplxval2(result_v1_c0));
	}

	// X- //////////////////////////////////////////////////////////////////////////
	{
	    if (t1==0 && x==0 && y==0 && z==0) {
	    	int a = 0;
	    }


		bgq_su3_weyl_decl(weyl_xdown);
		bgq_su3_weyl_load_double(weyl_xdown, &spinorsite->d[XDOWN]);
		bgq_weylqpx_expect(weyl_xdown, t1, t2, x, y, z, XDOWN, false);
			bgq_setdesc(BGQREF_XDOWN_RECV, "BGQREF_XDOWN_RECV");
			bgq_setbgqvalue(t1, x, y, z, BGQREF_XDOWN_RECV, bgq_cmplxval1(weyl_xdown_v1_c0));
			bgq_setbgqvalue(t2, x, y, z, BGQREF_XDOWN_RECV, bgq_cmplxval2(weyl_xdown_v1_c0));
		bgq_su3_accum_weyl_xdown(result, weyl_xdown);
			bgq_setdesc(BGQREF_XDOWN_ACCUM, "BGQREF_XDOWN_ACCUM");
			bgq_setbgqvalue(t1, x, y, z, BGQREF_XDOWN_ACCUM, bgq_cmplxval1(result_v1_c0));
			bgq_setbgqvalue(t2, x, y, z, BGQREF_XDOWN_ACCUM, bgq_cmplxval2(result_v1_c0));
	}

	// Y+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_yup);
		bgq_su3_weyl_load_double(weyl_yup, &spinorsite->d[YUP]);
		bgq_weylqpx_expect(weyl_yup, t1, t2, x, y, z, YUP,false);
			bgq_setdesc(BGQREF_YUP_RECV, "BGQREF_YUP_RECV");
			bgq_setbgqvalue(t1, x, y, z, BGQREF_YUP_RECV, bgq_cmplxval1(weyl_yup_v1_c0));
			bgq_setbgqvalue(t2, x, y, z, BGQREF_YUP_RECV, bgq_cmplxval2(weyl_yup_v1_c0));
		bgq_su3_accum_weyl_yup(result, weyl_yup);
			bgq_setdesc(BGQREF_YUP_ACCUM, "BGQREF_YUP_ACCUM");
			bgq_setbgqvalue(t1, x, y, z, BGQREF_YUP_ACCUM, bgq_cmplxval1(result_v1_c0));
			bgq_setbgqvalue(t2, x, y, z, BGQREF_YUP_ACCUM, bgq_cmplxval2(result_v1_c0));
	}

	// Y- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_ydown);
		bgq_su3_weyl_load_double(weyl_ydown, &spinorsite->d[YDOWN]);
		bgq_weylqpx_expect(weyl_ydown, t1, t2, x, y, z, YDOWN, false);
		bgq_su3_accum_weyl_ydown(result, weyl_ydown);
			bgq_setdesc(BGQREF_YDOWN_ACCUM, "BGQREF_YDOWN_ACCUM");
			bgq_setbgqvalue(t1, x, y, z, BGQREF_YDOWN_ACCUM, bgq_cmplxval1(result_v1_c0));
			bgq_setbgqvalue(t2, x, y, z, BGQREF_YDOWN_ACCUM, bgq_cmplxval2(result_v1_c0));
	}

	// Z+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zup);
		bgq_su3_weyl_load_double(weyl_zup, &spinorsite->d[ZUP]);
		bgq_weylqpx_expect(weyl_zup, t1, t2, x, y, z, ZUP, false);
			bgq_setdesc(BGQREF_ZUP_RECV, "BGQREF_ZUP_RECV");
			bgq_setbgqvalue(t1, x, y, z, BGQREF_ZUP_RECV, bgq_cmplxval1(weyl_zup_v1_c0));
			bgq_setbgqvalue(t2, x, y, z, BGQREF_ZUP_RECV, bgq_cmplxval2(weyl_zup_v1_c0));
		bgq_su3_accum_weyl_zup(result, weyl_zup);
			bgq_setdesc(BGQREF_ZUP_ACCUM, "BGQREF_ZUP_ACCUM");
			bgq_setbgqvalue(t1, x, y, z, BGQREF_ZUP_ACCUM, bgq_cmplxval1(result_v1_c0));
			bgq_setbgqvalue(t2, x, y, z, BGQREF_ZUP_ACCUM, bgq_cmplxval2(result_v1_c0));
	}

	// Z- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zdown);
		bgq_su3_weyl_load_double(weyl_zdown, &spinorsite->d[ZDOWN]);
		bgq_weylqpx_expect(weyl_zdown, t1, t2, x, y, z, ZDOWN,false);
		bgq_su3_accum_weyl_zdown(result, weyl_zdown);
	}

	bgq_setdesc(BGQREF_ACCUM, "BGQREF_ACCUM");
	bgq_setbgqvalue(t1, x, y, z, BGQREF_ACCUM, bgq_cmplxval1(result_v1_c0));
	bgq_setbgqvalue(t2, x, y, z, BGQREF_ACCUM, bgq_cmplxval2(result_v1_c0));
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


#define bgq_HoppingMatrix_compute_storeWeyllayout(targetptrs,gaugesite,spinor,t1,t2,x,y,z,kamul) bgq_HoppingMatrix_compute_storeWeyllayout_raw(targetptrs,gaugesite,bgq_su3_spinor_vars(spinor),t1,t2,x,y,z,kamul)
EXTERN_INLINE void bgq_HoppingMatrix_compute_storeWeyllayout_raw(bgq_weyl_ptr_t *targetptrs, bgq_gaugesite *gaugesite, bgq_su3_spinor_params(spinor), ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z, bool kamul) {
	//TODO: prefetch targetptrs
	//bgq_spinorqpx_expect(spinor, t1, t2, x, y, z);

    if (t1==0 && x==3 && y==0 && z==0) {
    	int a = 0;
    }

	// T+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tup);
		bgq_su3_reduce_weyl_tdown(weyl_tup, spinor);

		bgq_su3_mdecl(gauge_tup);
		bgq_su3_matrix_load_double(gauge_tup, &gaugesite->su3[TUP]);
		bgq_gaugeqpx_expect(gauge_tup, t1, t2, x, y, z, TUP, true);
				bgq_setdesc(BGQREF_TDOWN_GAUGE, "BGQREF_TDOWN_GAUGE");
				bgq_setbgqvalue_src(t1, x, y, z, TUP, BGQREF_TDOWN_GAUGE, bgq_cmplxval1(gauge_tup_c00));
				bgq_setbgqvalue_src(t2, x, y, z, TUP, BGQREF_TDOWN_GAUGE, bgq_cmplxval2(gauge_tup_c00));
		bgq_su3_weyl_mvinvmul(weyl_tup, gauge_tup, weyl_tup);
		if (kamul) { //TODO: Check conditional constant
			//TODO: Check code motion
			bgq_vector4double_decl(qka0);
			bgq_complxval_splat(qka0,ka0);
			bgq_su3_weyl_cmul(weyl_tup, qka0, weyl_tup);
		}
				bgq_setdesc(BGQREF_TDOWN_KAMUL,"BGQREF_TDOWN_KAMUL");
				bgq_setbgqvalue_src(t1, x, y, z, TUP, BGQREF_TDOWN_KAMUL, bgq_cmplxval1(weyl_tup_v1_c0));
				bgq_setbgqvalue_src(t2, x, y, z, TUP, BGQREF_TDOWN_KAMUL, bgq_cmplxval2(weyl_tup_v1_c0));

		//bgq_weylvec_expect(*targetptrs->d[TUP], t1, t2, x, y, z, TUP, true);
		bgq_su3_weyl_store_double(targetptrs->d[TUP], weyl_tup);
		bgq_weylvec_written(targetptrs->d[TUP], t1, t2, x, y, z, TUP, true);
	}

	// T- /////////////////////////////////////////////////////////////////////////
	{
		bgq_setdesc(BGQREF_TUP_SOURCE,"BGQREF_TUP_SOURCE");
				bgq_setbgqvalue_src(t1, x, y, z, TDOWN, BGQREF_TUP_SOURCE, bgq_cmplxval1(spinor_v1_c0));
				bgq_setbgqvalue_src(t2, x, y, z, TDOWN, BGQREF_TUP_SOURCE, bgq_cmplxval2(spinor_v1_c0));

		bgq_su3_weyl_decl(weyl_tdown);
		bgq_su3_reduce_weyl_tup(weyl_tdown, spinor);
				bgq_setdesc(BGQREF_TUP, "BGQREF_TUP");
				bgq_setbgqvalue_src(t1, x, y, z, TDOWN, BGQREF_TUP, bgq_cmplxval1(weyl_tdown_v1_c0));
				bgq_setbgqvalue_src(t2, x, y, z, TDOWN, BGQREF_TUP, bgq_cmplxval2(weyl_tdown_v1_c0));

		bgq_su3_mdecl(gauge_tdown);
		bgq_su3_matrix_load_double(gauge_tdown, &gaugesite->su3[TDOWN]);
		bgq_gaugeqpx_expect(gauge_tdown, t1, t2, x, y, z, TDOWN, true);
				bgq_setdesc(BGQREF_TUP_GAUGE, "BGQREF_TUP_GAUGE");
				bgq_setbgqvalue_src(t1, x, y, z, TDOWN, BGQREF_TUP_GAUGE, bgq_cmplxval1(gauge_tdown_c00));
				bgq_setbgqvalue_src(t2, x, y, z, TDOWN, BGQREF_TUP_GAUGE, bgq_cmplxval2(gauge_tdown_c00));
		bgq_su3_weyl_mvmul(weyl_tdown, gauge_tdown, weyl_tdown);
				bgq_setdesc(BGQREF_TUP_WEYL,"BGQREF_TUP_WEYL");
				bgq_setbgqvalue_src(t1, x, y, z, TDOWN, BGQREF_TUP_WEYL, bgq_cmplxval1(weyl_tdown_v1_c0));
				bgq_setbgqvalue_src(t2, x, y, z, TDOWN, BGQREF_TUP_WEYL, bgq_cmplxval2(weyl_tdown_v1_c0));
		if (kamul) {
			bgq_vector4double_decl(qka0);
			bgq_complxval_splat(qka0,ka0);
			bgq_su3_weyl_cmul(weyl_tdown, qka0, weyl_tdown);
		}
			bgq_setdesc(BGQREF_TUP_KAMUL,"BGQREF_TUP_KAMUL");
			bgq_setbgqvalue_src(t1, x, y, z, TDOWN, BGQREF_TUP_KAMUL, bgq_cmplxval1(weyl_tdown_v1_c0));
			bgq_setbgqvalue_src(t2, x, y, z, TDOWN, BGQREF_TUP_KAMUL, bgq_cmplxval2(weyl_tdown_v1_c0));

		//bgq_weylvec_expect(*targetptrs->d[TDOWN], t1, t2, x, y, z, TDOWN, true);
		bgq_su3_weyl_store_double(targetptrs->d[TDOWN], weyl_tdown);
		bgq_weylvec_written(targetptrs->d[TDOWN], t1, t2, x, y, z, TDOWN, true);
	}

	// X+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xup);
		bgq_su3_reduce_weyl_xdown(weyl_xup, spinor);
				bgq_setdesc(BGQREF_XDOWN, "BGQREF_XDOWN");
				bgq_setbgqvalue_src(t1, x, y, z, XUP, BGQREF_XDOWN, bgq_cmplxval1(weyl_xup_v1_c0));
				bgq_setbgqvalue_src(t2, x, y, z, XUP, BGQREF_XDOWN, bgq_cmplxval2(weyl_xup_v1_c0));

		bgq_su3_mdecl(gauge_xup);
		bgq_su3_matrix_load_double(gauge_xup, &gaugesite->su3[XUP]);
		bgq_gaugeqpx_expect(gauge_xup, t1, t2, x, y, z, XUP, true);
				bgq_setdesc(BGQREF_XDOWN_GAUGE, "BGQREF_XDOWN_GAUGE");
				bgq_setbgqvalue_src(t1, x, y, z, XUP, BGQREF_XDOWN_GAUGE, bgq_cmplxval1(gauge_xup_c02));
				bgq_setbgqvalue_src(t2, x, y, z, XUP, BGQREF_XDOWN_GAUGE, bgq_cmplxval2(gauge_xup_c02));
		bgq_su3_weyl_mvinvmul(weyl_xup, gauge_xup, weyl_xup);
				bgq_setdesc(BGQREF_XDOWN_WEYL,"BGQREF_XDOWN_WEYL");
				bgq_setbgqvalue_src(t1, x, y, z, XUP, BGQREF_XDOWN_WEYL, bgq_cmplxval1(weyl_xup_v1_c0));
				bgq_setbgqvalue_src(t2, x, y, z, XUP, BGQREF_XDOWN_WEYL, bgq_cmplxval2(weyl_xup_v1_c0));
		if (kamul) {
			bgq_vector4double_decl(qka1);
			bgq_complxval_splat(qka1, ka1);
			bgq_su3_weyl_cmul(weyl_xup, qka1, weyl_xup);
		}
				bgq_setdesc(BGQREF_XDOWN_KAMUL,"BGQREF_XDOWN_KAMUL");
				bgq_setbgqvalue_src(t1, x, y, z, XUP, BGQREF_XDOWN_KAMUL, bgq_cmplxval1(weyl_xup_v1_c0));
				bgq_setbgqvalue_src(t2, x, y, z, XUP, BGQREF_XDOWN_KAMUL, bgq_cmplxval2(weyl_xup_v1_c0));

		//bgq_weylvec_expect(*targetptrs->d[XUP], t1, t2, x, y, z, XUP, true);
		bgq_su3_weyl_store_double(targetptrs->d[XUP], weyl_xup);
		bgq_weylvec_written(targetptrs->d[XUP], t1, t2, x,y,z,XUP, true);
	}

	// X- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xdown);
		bgq_su3_reduce_weyl_xup(weyl_xdown, spinor);

		bgq_su3_mdecl(gauge_xdown);
		bgq_su3_matrix_load_double(gauge_xdown, &gaugesite->su3[XDOWN]);
		bgq_gaugeqpx_expect(gauge_xdown, t1, t2, x, y, z, XDOWN, true);
		bgq_su3_weyl_mvmul(weyl_xdown, gauge_xdown, weyl_xdown);
		if (kamul) {
			bgq_vector4double_decl(qka1);
			bgq_complxval_splat(qka1, ka1);
			bgq_su3_weyl_cmul(weyl_xdown, qka1, weyl_xdown);
		}
				bgq_setdesc(BGQREF_XUP_KAMUL,"BGQREF_XUP_KAMUL");
				bgq_setbgqvalue_src(t1, x, y, z, XDOWN, BGQREF_XUP_KAMUL, bgq_cmplxval1(weyl_xdown_v1_c0));
				bgq_setbgqvalue_src(t2, x, y, z, XDOWN, BGQREF_XUP_KAMUL, bgq_cmplxval2(weyl_xdown_v1_c0));

		//bgq_weylvec_expect(*targetptrs->d[XDOWN], t1, t2, x, y, z, XDOWN, true);
		size_t offset = bgq_fieldpointer2offset(targetptrs->d[XDOWN]);
		if (offset == 280192) {
			int c = 0;
		}
		offset = bgq_fieldpointer2offset(targetptrs->d[XDOWN]);
		ucoord index = bgq_offset2index(offset);
		if (index == 50) {
			int a = 0;
		}
		if (bgq_cmplxval1(weyl_xdown_v1_c0)==0) {
			int b = 0;
		}
		bgq_su3_weyl_store_double(targetptrs->d[XDOWN], weyl_xdown);
		bgq_weylvec_written(targetptrs->d[XDOWN], t1, t2, x,y,z,XDOWN, true);
	}

	// Y+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_yup);
		bgq_su3_reduce_weyl_ydown(weyl_yup, spinor);

		bgq_su3_mdecl(gauge_yup);
		bgq_su3_matrix_load_double(gauge_yup, &gaugesite->su3[YUP]);
		bgq_gaugeqpx_expect(gauge_yup, t1, t2, x, y, z, YUP, true);
		bgq_su3_weyl_mvinvmul(weyl_yup, gauge_yup, weyl_yup);
		if (kamul) {
			bgq_vector4double_decl(qka2);
			bgq_complxval_splat(qka2, ka2);
			bgq_su3_weyl_cmul(weyl_yup, qka2, weyl_yup);
		}

		//bgq_weylvec_expect(*targetptrs->d[YUP], t1, t2, x, y, z, YUP, true);
		bgq_su3_weyl_store_double(targetptrs->d[YUP], weyl_yup);
		bgq_weylvec_written(targetptrs->d[YUP], t1, t2, x,y,z,YUP, true);
	}

	// Y- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_ydown);
		bgq_su3_reduce_weyl_yup(weyl_ydown, spinor);

		bgq_su3_mdecl(gauge_ydown);
		bgq_su3_matrix_load_double(gauge_ydown, &gaugesite->su3[YDOWN]);
		bgq_gaugeqpx_expect(gauge_ydown, t1, t2, x, y, z, YDOWN, true);
		bgq_su3_weyl_mvmul(weyl_ydown, gauge_ydown, weyl_ydown);
		if (kamul) {
			bgq_vector4double_decl(qka2);
			bgq_complxval_splat(qka2, ka2);
			bgq_su3_weyl_cmul(weyl_ydown, qka2, weyl_ydown);
		}

		//bgq_weylvec_expect(*targetptrs->d[YDOWN], t1, t2, x, y, z, YDOWN, true);
		bgq_su3_weyl_store_double(targetptrs->d[YDOWN], weyl_ydown);
		bgq_weylvec_written(targetptrs->d[YDOWN], t1, t2, x,y,z,YDOWN, true);
	}

	// Z+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zup);
		bgq_su3_reduce_weyl_zdown(weyl_zup, spinor);

		bgq_su3_mdecl(gauge_zup);
		bgq_su3_matrix_load_double(gauge_zup, &gaugesite->su3[ZUP]);
		bgq_gaugeqpx_expect(gauge_zup, t1, t2, x, y, z, ZUP, true);
		bgq_su3_weyl_mvinvmul(weyl_zup, gauge_zup, weyl_zup);
		if (kamul) {
			bgq_vector4double_decl(qka3);
			bgq_complxval_splat(qka3, ka3);
			bgq_su3_weyl_cmul(weyl_zup, qka3, weyl_zup);
		}

		//bgq_weylvec_expect(*targetptrs->d[ZUP], t1, t2, x, y, z, ZUP, true);
		bgq_su3_weyl_store_double(targetptrs->d[ZUP], weyl_zup);
		bgq_weylvec_written(targetptrs->d[ZUP], t1, t2, x,y,z,ZUP, true);
	}

	// Z- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zdown);
		bgq_su3_reduce_weyl_zup(weyl_zdown, spinor);
			bgq_setdesc(BGQREF_ZUP,"BGQREF_ZUP");
			bgq_setbgqvalue_src(t1, x, y, z, ZDOWN, BGQREF_ZUP, bgq_cmplxval1(weyl_zdown_v1_c0));
			bgq_setbgqvalue_src(t2, x, y, z, ZDOWN, BGQREF_ZUP, bgq_cmplxval2(weyl_zdown_v1_c0));

		bgq_su3_mdecl(gauge_zdown);
		bgq_su3_matrix_load_double(gauge_zdown, &gaugesite->su3[ZDOWN]);
		bgq_gaugeqpx_expect(gauge_zdown, t1, t2, x, y, z, ZDOWN, true);
				bgq_setdesc(BGQREF_ZUP_GAUGE, "BGQREF_ZUP_GAUGE");
				bgq_setbgqvalue_src(t1, x, y, z, ZDOWN, BGQREF_ZUP_GAUGE, bgq_cmplxval1(gauge_zdown_c00));
				bgq_setbgqvalue_src(t2, x, y, z, ZDOWN, BGQREF_ZUP_GAUGE, bgq_cmplxval2(gauge_zdown_c00));
		bgq_su3_weyl_mvmul(weyl_zdown, gauge_zdown, weyl_zdown);
		if (kamul) {
			bgq_vector4double_decl(qka3);
			bgq_complxval_splat(qka3, ka3);
			bgq_su3_weyl_cmul(weyl_zdown, qka3, weyl_zdown);
		}
		bgq_setdesc(BGQREF_ZUP_KAMUL,"BGQREF_ZUP_KAMUL");
		bgq_setbgqvalue_src(t1, x, y, z, ZDOWN, BGQREF_ZUP_KAMUL, bgq_cmplxval1(weyl_zdown_v1_c0));
		bgq_setbgqvalue_src(t2, x, y, z, ZDOWN, BGQREF_ZUP_KAMUL, bgq_cmplxval2(weyl_zdown_v1_c0));

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
