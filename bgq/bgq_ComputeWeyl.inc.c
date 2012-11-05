/*
 * IBM xlc refuses to inline functions that big, therefore we have to include the function using the preprocessor
 */

#ifndef BGQ_COMPUTEWEYL_INC_
#include "bgq_qpx.h"

#include <stdbool.h>
void bgq_HoppingMatrix_compute_storeWeyllayout_raw(bgq_weyl_ptr_t *targetptrs, bgq_gaugesite *gaugesite, bgq_su3_spinor_params(spinor), ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z, bool kamul)
#endif
{

	//TODO: prefetch targetptrs
	//bgq_spinorqpx_expect(spinor, t1, t2, x, y, z);

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
			bgq_su3_weyl_cmul(weyl_xdown, qka1, weyl_xdown);
		}
				bgq_setdesc(BGQREF_XUP_KAMUL,"BGQREF_XUP_KAMUL");
				bgq_setbgqvalue_src(t1, x, y, z, XDOWN, BGQREF_XUP_KAMUL, bgq_cmplxval1(weyl_xdown_v1_c0));
				bgq_setbgqvalue_src(t2, x, y, z, XDOWN, BGQREF_XUP_KAMUL, bgq_cmplxval2(weyl_xdown_v1_c0));

		//bgq_weylvec_expect(*targetptrs->d[XDOWN], t1, t2, x, y, z, XDOWN, true);
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
