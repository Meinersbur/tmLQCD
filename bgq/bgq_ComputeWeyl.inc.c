/*
 * IBM xlc refuses to inline functions that big, therefore we have to include the function using the preprocessor
 */


#ifndef BGQ_COMPUTEWEYL_INC_
#include "bgq_qpx.h"
#include "bgq_spinorfield.h"
#include "bgq_gaugefield.h"

#include <stdbool.h>

#define PRECISION double

void bgq_HoppingMatrix_compute_storeWeyllayout_raw(bgq_weyl_ptr_t *targetptrs, bgq_gaugesite *gaugesite, bgq_su3_spinor_params(spinor), bgq_params(qka0), bgq_params(qka1), bgq_params(qka2), bgq_params(qka3),ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z, bool kamul)
#endif

#ifndef BGQ_COMPUTEWEYL_INSERTPREFETCH
#define BGQ_COMPUTEWEYL_INSERTPREFETCH
#endif

{
	asm volatile ("");
	bgq_prefetch(&targetptrs->d[TUP]);
	bgq_su3_matrixnext_prefetch(gaugesite);

	// T+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_tup);
		bgq_su3_reduce_weyl_tdown(weyl_tup, spinor);

		bgq_su3_mdecl(gauge_tup);
		bgq_qvlfuxa(gauge_tup_c00, gaugesite, 32);
		bgq_qvlfuxa(gauge_tup_c01, gaugesite, 32);
		bgq_qvlfuxa(gauge_tup_c02, gaugesite, 32);
		bgq_qvlfuxa(gauge_tup_c10, gaugesite, 32);
		bgq_qvlfuxa(gauge_tup_c11, gaugesite, 32);
		bgq_qvlfuxa(gauge_tup_c12, gaugesite, 32);
		bgq_qvlfuxa(gauge_tup_c20, gaugesite, 32);
		bgq_qvlfuxa(gauge_tup_c21, gaugesite, 32);
		bgq_qvlfuxa(gauge_tup_c22, gaugesite, 32);
				bgq_gaugeqpx_expect(gauge_tup, t1, t2, x, y, z, TUP, true);
				bgq_setdesc(BGQREF_TDOWN_GAUGE, "BGQREF_TDOWN_GAUGE");
				bgq_setbgqvalue_src(t1, x, y, z, TUP, BGQREF_TDOWN_GAUGE, bgq_cmplxval1(gauge_tup_c00));
				bgq_setbgqvalue_src(t2, x, y, z, TUP, BGQREF_TDOWN_GAUGE, bgq_cmplxval2(gauge_tup_c00));

		bgq_su3_weyl_mvinvmul(weyl_tup, gauge_tup, weyl_tup);
		if (kamul) {
			bgq_su3_weyl_cmul(weyl_tup, qka0, weyl_tup);
		}
				bgq_setdesc(BGQREF_TDOWN_KAMUL,"BGQREF_TDOWN_KAMUL");
				bgq_setbgqvalue_src(t1, x, y, z, TUP, BGQREF_TDOWN_KAMUL, bgq_cmplxval1(weyl_tup_v1_c0));
				bgq_setbgqvalue_src(t2, x, y, z, TUP, BGQREF_TDOWN_KAMUL, bgq_cmplxval2(weyl_tup_v1_c0));

		bgq_su3_weyl_zeroload(targetptrs->d[TUP]);
		bgq_su3_weyl_store(targetptrs->d[TUP], weyl_tup);
		bgq_weylvec_written(targetptrs->d[TUP], t1, t2, x, y, z, TUP, true);
	}

	asm volatile ("");
	//bgq_prefetch(&targetptrs->d[TDOWN]);
	bgq_su3_matrixnext_prefetch(gaugesite);

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
		bgq_qvlfuxa(gauge_tdown_c00, gaugesite, 32);
		bgq_qvlfuxa(gauge_tdown_c01, gaugesite, 32);
		bgq_qvlfuxa(gauge_tdown_c02, gaugesite, 32);
		bgq_qvlfuxa(gauge_tdown_c10, gaugesite, 32);
		bgq_qvlfuxa(gauge_tdown_c11, gaugesite, 32);
		bgq_qvlfuxa(gauge_tdown_c12, gaugesite, 32);
		bgq_qvlfuxa(gauge_tdown_c20, gaugesite, 32);
		bgq_qvlfuxa(gauge_tdown_c21, gaugesite, 32);
		bgq_qvlfuxa(gauge_tdown_c22, gaugesite, 32);
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

		bgq_su3_weyl_zeroload(targetptrs->d[TDOWN]);
		bgq_su3_weyl_store(targetptrs->d[TDOWN], weyl_tdown);
		bgq_weylvec_written(targetptrs->d[TDOWN], t1, t2, x, y, z, TDOWN, true);
	}

	asm volatile ("");
	//bgq_prefetch(&targetptrs->d[XUP]);
	bgq_su3_matrixnext_prefetch(gaugesite);

	// X+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xup);
		bgq_su3_reduce_weyl_xdown(weyl_xup, spinor);

		bgq_su3_mdecl(gauge_xup);
		bgq_qvlfuxa(gauge_xup_c00, gaugesite, 32);
		bgq_qvlfuxa(gauge_xup_c01, gaugesite, 32);
		bgq_qvlfuxa(gauge_xup_c02, gaugesite, 32);
		bgq_qvlfuxa(gauge_xup_c10, gaugesite, 32);
		bgq_qvlfuxa(gauge_xup_c11, gaugesite, 32);
		bgq_qvlfuxa(gauge_xup_c12, gaugesite, 32);
		bgq_qvlfuxa(gauge_xup_c20, gaugesite, 32);
		bgq_qvlfuxa(gauge_xup_c21, gaugesite, 32);
		bgq_qvlfuxa(gauge_xup_c22, gaugesite, 32);

		bgq_su3_weyl_mvinvmul(weyl_xup, gauge_xup, weyl_xup);
		if (kamul) {
			bgq_su3_weyl_cmul(weyl_xup, qka1, weyl_xup);
		}

		bgq_su3_weyl_zeroload(targetptrs->d[XUP]);
		bgq_su3_weyl_store(targetptrs->d[XUP], weyl_xup);
		bgq_weylvec_written(targetptrs->d[XUP], t1, t2, x,y,z,XUP, true);
	}

	asm volatile ("");
	//bgq_prefetch(&targetptrs->d[XDOWN]);
	bgq_su3_matrixnext_prefetch(gaugesite);

	// X- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xdown);
		bgq_su3_reduce_weyl_xup(weyl_xdown, spinor);

		bgq_su3_mdecl(gauge_xdown);
		bgq_qvlfuxa(gauge_xdown_c00, gaugesite, 32);
		bgq_qvlfuxa(gauge_xdown_c01, gaugesite, 32);
		bgq_qvlfuxa(gauge_xdown_c02, gaugesite, 32);
		bgq_qvlfuxa(gauge_xdown_c10, gaugesite, 32);
		bgq_qvlfuxa(gauge_xdown_c11, gaugesite, 32);
		bgq_qvlfuxa(gauge_xdown_c12, gaugesite, 32);
		bgq_qvlfuxa(gauge_xdown_c20, gaugesite, 32);
		bgq_qvlfuxa(gauge_xdown_c21, gaugesite, 32);
		bgq_qvlfuxa(gauge_xdown_c22, gaugesite, 32);

		bgq_su3_weyl_mvmul(weyl_xdown, gauge_xdown, weyl_xdown);
		if (kamul) {
			bgq_su3_weyl_cmul(weyl_xdown, qka1, weyl_xdown);
		}

		bgq_su3_weyl_zeroload(targetptrs->d[XDOWN]);
		bgq_su3_weyl_store(targetptrs->d[XDOWN], weyl_xdown);
		bgq_weylvec_written(targetptrs->d[XDOWN], t1, t2, x,y,z,XDOWN, true);
	}

	asm volatile ("");
	//bgq_prefetch(&targetptrs->d[YUP]);
	bgq_su3_matrixnext_prefetch(gaugesite);

	// Y+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_yup);
		bgq_su3_reduce_weyl_ydown(weyl_yup, spinor);

		bgq_su3_mdecl(gauge_yup);
		bgq_qvlfuxa(gauge_yup_c00, gaugesite, 32);
		bgq_qvlfuxa(gauge_yup_c01, gaugesite, 32);
		bgq_qvlfuxa(gauge_yup_c02, gaugesite, 32);
		bgq_qvlfuxa(gauge_yup_c10, gaugesite, 32);
		bgq_qvlfuxa(gauge_yup_c11, gaugesite, 32);
		bgq_qvlfuxa(gauge_yup_c12, gaugesite, 32);
		bgq_qvlfuxa(gauge_yup_c20, gaugesite, 32);
		bgq_qvlfuxa(gauge_yup_c21, gaugesite, 32);
		bgq_qvlfuxa(gauge_yup_c22, gaugesite, 32);

		bgq_su3_weyl_mvinvmul(weyl_yup, gauge_yup, weyl_yup);
		if (kamul) {
			bgq_su3_weyl_cmul(weyl_yup, qka2, weyl_yup);
		}

		bgq_su3_weyl_zeroload(targetptrs->d[YUP]);
		bgq_su3_weyl_store(targetptrs->d[YUP], weyl_yup);
		bgq_weylvec_written(targetptrs->d[YUP], t1, t2, x,y,z,YUP, true);
	}

	asm volatile ("");
	//bgq_prefetch(&targetptrs->d[YDOWN]);
	bgq_su3_matrixnext_prefetch(gaugesite);

	// Y- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_ydown);
		bgq_su3_reduce_weyl_yup(weyl_ydown, spinor);

		bgq_su3_mdecl(gauge_ydown);
		bgq_qvlfuxa(gauge_ydown_c00, gaugesite, 32);
		bgq_qvlfuxa(gauge_ydown_c01, gaugesite, 32);
		bgq_qvlfuxa(gauge_ydown_c02, gaugesite, 32);
		bgq_qvlfuxa(gauge_ydown_c10, gaugesite, 32);
		bgq_qvlfuxa(gauge_ydown_c11, gaugesite, 32);
		bgq_qvlfuxa(gauge_ydown_c12, gaugesite, 32);
		bgq_qvlfuxa(gauge_ydown_c20, gaugesite, 32);
		bgq_qvlfuxa(gauge_ydown_c21, gaugesite, 32);
		bgq_qvlfuxa(gauge_ydown_c22, gaugesite, 32);

		bgq_su3_weyl_mvmul(weyl_ydown, gauge_ydown, weyl_ydown);
		if (kamul) {
			bgq_su3_weyl_cmul(weyl_ydown, qka2, weyl_ydown);
		}

		bgq_su3_weyl_zeroload(targetptrs->d[YDOWN]);
		bgq_su3_weyl_store(targetptrs->d[YDOWN], weyl_ydown);
		bgq_weylvec_written(targetptrs->d[YDOWN], t1, t2, x,y,z,YDOWN, true);
	}

	asm volatile ("");
	//bgq_prefetch(&targetptrs->d[ZUP]);
	bgq_su3_matrixnext_prefetch(gaugesite);

	// Z+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zup);
		bgq_su3_reduce_weyl_zdown(weyl_zup, spinor);

		bgq_su3_mdecl(gauge_zup);
		bgq_qvlfuxa(gauge_zup_c00, gaugesite, 32);
		bgq_qvlfuxa(gauge_zup_c01, gaugesite, 32);
		bgq_qvlfuxa(gauge_zup_c02, gaugesite, 32);
		bgq_qvlfuxa(gauge_zup_c10, gaugesite, 32);
		bgq_qvlfuxa(gauge_zup_c11, gaugesite, 32);
		bgq_qvlfuxa(gauge_zup_c12, gaugesite, 32);
		bgq_qvlfuxa(gauge_zup_c20, gaugesite, 32);
		bgq_qvlfuxa(gauge_zup_c21, gaugesite, 32);
		bgq_qvlfuxa(gauge_zup_c22, gaugesite, 32);

		bgq_su3_weyl_mvinvmul(weyl_zup, gauge_zup, weyl_zup);
		if (kamul) {
			bgq_su3_weyl_cmul(weyl_zup, qka3, weyl_zup);
		}

		bgq_su3_weyl_zeroload(targetptrs->d[ZUP]);
		bgq_su3_weyl_store(targetptrs->d[ZUP], weyl_zup);
		bgq_weylvec_written(targetptrs->d[ZUP], t1, t2, x,y,z,ZUP, true);
	}

	asm volatile ("");
	//bgq_prefetch(&targetptrs->d[ZDOWN]);
	BGQ_COMPUTEWEYL_INSERTPREFETCH

	// Z- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zdown);
		bgq_su3_reduce_weyl_zup(weyl_zdown, spinor);

		bgq_su3_mdecl(gauge_zdown);
		bgq_qvlfuxa(gauge_zdown_c00, gaugesite, 32);
		bgq_qvlfuxa(gauge_zdown_c01, gaugesite, 32);
		bgq_qvlfuxa(gauge_zdown_c02, gaugesite, 32);
		bgq_qvlfuxa(gauge_zdown_c10, gaugesite, 32);
		bgq_qvlfuxa(gauge_zdown_c11, gaugesite, 32);
		bgq_qvlfuxa(gauge_zdown_c12, gaugesite, 32);
		bgq_qvlfuxa(gauge_zdown_c20, gaugesite, 32);
		bgq_qvlfuxa(gauge_zdown_c21, gaugesite, 32);
		bgq_qvlfuxa(gauge_zdown_c22, gaugesite, 32);

		bgq_su3_weyl_mvmul(weyl_zdown, gauge_zdown, weyl_zdown);
		if (kamul) {
			bgq_su3_weyl_cmul(weyl_zdown, qka3, weyl_zdown);
		}

		bgq_su3_weyl_zeroload(targetptrs->d[ZDOWN]);
		bgq_su3_weyl_store(targetptrs->d[ZDOWN], weyl_zdown);
		bgq_weylvec_written(targetptrs->d[ZDOWN], t1, t2, x,y,z,ZDOWN, true);
	}
	asm volatile ("");
}

#undef BGQ_COMPUTEWEYL_INSERTPREFETCH
#undef BGQ_COMPUTEWEYL_INC_
