/*
 * bgq_readWeyllayout.inc.c
 *
 *  Created on: Nov 13, 2012
 *      Author: meinersbur
 */



#ifndef BGQ_READWEYLLAYOUT_INC_
#include "bgq_qpx.h"
#include "bgq_spinorfield.h"

void bgq_readWeyllayout(bgq_su3_spinor_params(/*out*/spinor), bgq_weylsite *weylsite, ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z)
#endif
{

	// T+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tup);
		bgq_su3_weyl_load_double(weyl_tup, &weylsite->d[TUP]);
		//bgq_su3_weyl_valgen(weyl_tup);
		bgq_weylqpx_expect(weyl_tup, t1, t2, x, y, z, TUP, false);
				bgq_setdesc(BGQREF_TUP_RECV, "BGQREF_TUP_RECV");
				bgq_setbgqvalue(t1, x, y, z, BGQREF_TUP_RECV, bgq_cmplxval1(weyl_tup_v1_c0));
				bgq_setbgqvalue(t2, x, y, z, BGQREF_TUP_RECV, bgq_cmplxval2(weyl_tup_v1_c0));
		bgq_su3_expand_weyl_tup(spinor, weyl_tup);
	}

	//bgq_su3_weyl_prefetch_double(&spinorsite->d[XUP]);

	// T- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tdown);
		bgq_su3_weyl_load_double(weyl_tdown, &weylsite->d[TDOWN]);
		//bgq_su3_weyl_valgen(weyl_tdown);
		bgq_weylqpx_expect(weyl_tdown, t1, t2, x, y, z, TDOWN, false);
				bgq_setdesc(BGQREF_TDOWN_RECV, "BGQREF_TDOWN_RECV");
				bgq_setbgqvalue(t1, x, y, z, BGQREF_TDOWN_RECV, bgq_cmplxval1(weyl_tdown_v1_c0));
				bgq_setbgqvalue(t2, x, y, z, BGQREF_TDOWN_RECV, bgq_cmplxval2(weyl_tdown_v1_c0));
		bgq_su3_accum_weyl_tdown(spinor, weyl_tdown);
				bgq_setdesc(BGQREF_TDOWN_ACCUM, "BGQREF_TDOWN_ACCUM");
				bgq_setbgqvalue(t1, x, y, z, BGQREF_TDOWN_ACCUM, bgq_cmplxval1(spinor_v1_c0));
				bgq_setbgqvalue(t2, x, y, z, BGQREF_TDOWN_ACCUM, bgq_cmplxval2(spinor_v1_c0));
	}

	//bgq_su3_weyl_prefetch_double(&spinorsite->d[XDOWN]);

	// X+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xup);
		bgq_su3_weyl_load_double(weyl_xup, &weylsite->d[XUP]);
		//bgq_su3_weyl_valgen(weyl_xup);
		bgq_weylqpx_expect(weyl_xup, t1, t2, x, y, z, XUP, false);
			bgq_setdesc(BGQREF_XUP_RECV, "BGQREF_XUP_RECV");
			bgq_setbgqvalue(t1, x, y, z, BGQREF_XUP_RECV, bgq_cmplxval1(weyl_xup_v1_c0));
			bgq_setbgqvalue(t2, x, y, z, BGQREF_XUP_RECV, bgq_cmplxval2(weyl_xup_v1_c0));
		bgq_su3_accum_weyl_xup(spinor, weyl_xup);
			bgq_setdesc(BGQREF_XUP_ACCUM, "BGQREF_XUP_ACCUM");
			bgq_setbgqvalue(t1, x, y, z, BGQREF_XUP_ACCUM, bgq_cmplxval1(spinor_v1_c0));
			bgq_setbgqvalue(t2, x, y, z, BGQREF_XUP_ACCUM, bgq_cmplxval2(spinor_v1_c0));
	}

	//bgq_su3_weyl_prefetch_double(&spinorsite->d[YUP]);

	// X- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xdown);
		bgq_su3_weyl_load_double(weyl_xdown, &weylsite->d[XDOWN]);
		//bgq_su3_weyl_valgen(weyl_xdown);
		bgq_weylqpx_expect(weyl_xdown, t1, t2, x, y, z, XDOWN, false);
			bgq_setdesc(BGQREF_XDOWN_RECV, "BGQREF_XDOWN_RECV");
			bgq_setbgqvalue(t1, x, y, z, BGQREF_XDOWN_RECV, bgq_cmplxval1(weyl_xdown_v1_c0));
			bgq_setbgqvalue(t2, x, y, z, BGQREF_XDOWN_RECV, bgq_cmplxval2(weyl_xdown_v1_c0));
		bgq_su3_accum_weyl_xdown(spinor, weyl_xdown);
			bgq_setdesc(BGQREF_XDOWN_ACCUM, "BGQREF_XDOWN_ACCUM");
			bgq_setbgqvalue(t1, x, y, z, BGQREF_XDOWN_ACCUM, bgq_cmplxval1(spinor_v1_c0));
			bgq_setbgqvalue(t2, x, y, z, BGQREF_XDOWN_ACCUM, bgq_cmplxval2(spinor_v1_c0));
	}

	//bgq_su3_weyl_prefetch_double(&spinorsite->d[YDOWN]);

	// Y+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_yup);
		bgq_su3_weyl_load_double(weyl_yup, &weylsite->d[YUP]);
		//bgq_su3_weyl_valgen(weyl_yup);
		bgq_weylqpx_expect(weyl_yup, t1, t2, x, y, z, YUP,false);
			bgq_setdesc(BGQREF_YUP_RECV, "BGQREF_YUP_RECV");
			bgq_setbgqvalue(t1, x, y, z, BGQREF_YUP_RECV, bgq_cmplxval1(weyl_yup_v1_c0));
			bgq_setbgqvalue(t2, x, y, z, BGQREF_YUP_RECV, bgq_cmplxval2(weyl_yup_v1_c0));
		bgq_su3_accum_weyl_yup(spinor, weyl_yup);
			bgq_setdesc(BGQREF_YUP_ACCUM, "BGQREF_YUP_ACCUM");
			bgq_setbgqvalue(t1, x, y, z, BGQREF_YUP_ACCUM, bgq_cmplxval1(spinor_v1_c0));
			bgq_setbgqvalue(t2, x, y, z, BGQREF_YUP_ACCUM, bgq_cmplxval2(spinor_v1_c0));
	}

	//bgq_su3_weyl_prefetch_double(&spinorsite->d[ZUP]);

	// Y- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_ydown);
		bgq_su3_weyl_load_double(weyl_ydown, &weylsite->d[YDOWN]);
		//bgq_su3_weyl_valgen(weyl_ydown);
		assert(bgq_cmplxval1(weyl_ydown_v1_c0)!=0);
		bgq_weylqpx_expect(weyl_ydown, t1, t2, x, y, z, YDOWN, false);
		bgq_su3_accum_weyl_ydown(spinor, weyl_ydown);
			bgq_setdesc(BGQREF_YDOWN_ACCUM, "BGQREF_YDOWN_ACCUM");
			bgq_setbgqvalue(t1, x, y, z, BGQREF_YDOWN_ACCUM, bgq_cmplxval1(spinor_v1_c0));
			bgq_setbgqvalue(t2, x, y, z, BGQREF_YDOWN_ACCUM, bgq_cmplxval2(spinor_v1_c0));
	}

	//bgq_su3_weyl_prefetch_double(&spinorsite->d[ZDOWN]);

	// Z+ //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zup);
		bgq_su3_weyl_load_double(weyl_zup, &weylsite->d[ZUP]);
		//bgq_su3_weyl_valgen(weyl_zup);
		assert(bgq_cmplxval1(weyl_zup_v1_c0)!=0);
		bgq_weylqpx_expect(weyl_zup, t1, t2, x, y, z, ZUP, false);
			bgq_setdesc(BGQREF_ZUP_RECV, "BGQREF_ZUP_RECV");
			bgq_setbgqvalue(t1, x, y, z, BGQREF_ZUP_RECV, bgq_cmplxval1(weyl_zup_v1_c0));
			bgq_setbgqvalue(t2, x, y, z, BGQREF_ZUP_RECV, bgq_cmplxval2(weyl_zup_v1_c0));
		bgq_su3_accum_weyl_zup(spinor, weyl_zup);
			bgq_setdesc(BGQREF_ZUP_ACCUM, "BGQREF_ZUP_ACCUM");
			bgq_setbgqvalue(t1, x, y, z, BGQREF_ZUP_ACCUM, bgq_cmplxval1(spinor_v1_c0));
			bgq_setbgqvalue(t2, x, y, z, BGQREF_ZUP_ACCUM, bgq_cmplxval2(spinor_v1_c0));
	}

	//bgq_su3_weyl_prefetch_double(&spinorsite->d[TUP]);

	// Z- //////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zdown);
		bgq_su3_weyl_load_double(weyl_zdown, &weylsite->d[ZDOWN]);
		//bgq_su3_weyl_valgen(weyl_zdown);
		bgq_weylqpx_expect(weyl_zdown, t1, t2, x, y, z, ZDOWN,false);
		bgq_su3_accum_weyl_zdown(spinor, weyl_zdown);
	}


			bgq_setdesc(BGQREF_ACCUM, "BGQREF_ACCUM");
			bgq_setbgqvalue(t1, x, y, z, BGQREF_ACCUM, bgq_cmplxval1(spinor_v1_c0));
			bgq_setbgqvalue(t2, x, y, z, BGQREF_ACCUM, bgq_cmplxval2(spinor_v1_c0));
	//bgq_su3_spinor_mov(*target, result);
}

#undef BGQ_READWEYLLAYOUT_INC_
