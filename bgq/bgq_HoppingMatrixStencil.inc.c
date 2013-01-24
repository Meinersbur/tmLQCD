/*
 * bgq_HoppingMatrixStencil.inc.c
 *
 *  Created on: Jan 23, 2013
 *      Author: meinersbur
 */

#ifndef BGQ_HOPPINGMATRIXSTENCIL_INC_
#include "bgq_qpx.h"
#include "bgq_spinorfield.h"
#include "bgq_gaugefield.h"

#include <stdbool.h>

#define PRECISION double
#define GETSPINORFIELDPTR(DIRECTION, direction) (pspinor_ ## direction)

void bgq_HoppingMatrix_stencil_raw(bgq_spinor_vec *pspinor_tup, bgq_spinor_vec *pspinor_tdown, bgq_spinor_vec *pspinor_xup, bgq_spinor_vec *pspinor_xdown, bgq_spinor_vec *pspinor_yup, bgq_spinor_vec *pspinor_ydown, bgq_spinor_vec *pspinor_zup, bgq_spinor_vec *pspinor_zdown, bgq_gaugesite *gaugesite, bgq_spinor_vec *target, bgq_params(qka0), bgq_params(qka1), bgq_params(qka2), bgq_params(qka3), ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z, bool kamul)
#endif
{
#ifndef NDEBUG
#endif

	bgq_su3_spinor_decl(result);

	// T+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_spinor_decl(spinor_tup);
		bgq_spinor_vec *spinorptr_tup = GETSPINORFIELDPTR(TUP, tup);
		bgq_su3_spinor_load(spinor_tup, spinorptr_tup);

		bgq_su3_weyl_decl(weyl_tup);
		bgq_su3_reduce_weyl_tup(weyl_tup, spinor_tup);

		bgq_su3_mdecl(gauge_tup);
		bgq_qvlfduxa(gauge_tup_c00, gaugesite, 32);
		bgq_qvlfduxa(gauge_tup_c01, gaugesite, 32);
		bgq_qvlfduxa(gauge_tup_c02, gaugesite, 32);
		bgq_qvlfduxa(gauge_tup_c10, gaugesite, 32);
		bgq_qvlfduxa(gauge_tup_c11, gaugesite, 32);
		bgq_qvlfduxa(gauge_tup_c12, gaugesite, 32);
		bgq_qvlfduxa(gauge_tup_c20, gaugesite, 32);
		bgq_qvlfduxa(gauge_tup_c21, gaugesite, 32);
		bgq_qvlfduxa(gauge_tup_c22, gaugesite, 32);

		bgq_su3_weyl_mvmul(weyl_tup, gauge_tup, weyl_tup);

		if (kamul) {
			bgq_su3_weyl_cmul(weyl_tup, qka0, weyl_tup);
		}

		bgq_su3_expand_weyl_tup(result, weyl_tup);
	}


	// T- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_spinor_decl(spinor_tdown);
		bgq_su3_spinor_load(spinor_tdown, (GETSPINORFIELDPTR(TDOWN, tdown)));

		bgq_su3_weyl_decl(weyl_tdown);
		bgq_su3_reduce_weyl_tdown(weyl_tdown, spinor_tdown);

		bgq_su3_mdecl(gauge_tdown);
		bgq_qvlfduxa(gauge_tdown_c00, gaugesite, 32);
		bgq_qvlfduxa(gauge_tdown_c01, gaugesite, 32);
		bgq_qvlfduxa(gauge_tdown_c02, gaugesite, 32);
		bgq_qvlfduxa(gauge_tdown_c10, gaugesite, 32);
		bgq_qvlfduxa(gauge_tdown_c11, gaugesite, 32);
		bgq_qvlfduxa(gauge_tdown_c12, gaugesite, 32);
		bgq_qvlfduxa(gauge_tdown_c20, gaugesite, 32);
		bgq_qvlfduxa(gauge_tdown_c21, gaugesite, 32);
		bgq_qvlfduxa(gauge_tdown_c22, gaugesite, 32);

		bgq_su3_weyl_mvmul(weyl_tdown, gauge_tdown, weyl_tdown);

		if (kamul) {
			bgq_su3_weyl_cmul(weyl_tdown, qka0, weyl_tdown);
		}

		bgq_su3_accum_weyl_tdown(result, weyl_tdown);
	}


	// X+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_spinor_decl(spinor_xup);
		bgq_su3_spinor_load(spinor_xup, (GETSPINORFIELDPTR(XUP, xup)));

		bgq_su3_weyl_decl(weyl_xup);
		bgq_su3_reduce_weyl_xup(weyl_xup, spinor_xup);

		bgq_su3_mdecl(gauge_xup);
		bgq_qvlfduxa(gauge_xup_c00, gaugesite, 32);
		bgq_qvlfduxa(gauge_xup_c01, gaugesite, 32);
		bgq_qvlfduxa(gauge_xup_c02, gaugesite, 32);
		bgq_qvlfduxa(gauge_xup_c10, gaugesite, 32);
		bgq_qvlfduxa(gauge_xup_c11, gaugesite, 32);
		bgq_qvlfduxa(gauge_xup_c12, gaugesite, 32);
		bgq_qvlfduxa(gauge_xup_c20, gaugesite, 32);
		bgq_qvlfduxa(gauge_xup_c21, gaugesite, 32);
		bgq_qvlfduxa(gauge_xup_c22, gaugesite, 32);

		bgq_su3_weyl_mvmul(weyl_xup, gauge_xup, weyl_xup);

		if (kamul) {
			bgq_su3_weyl_cmul(weyl_xup, qka1, weyl_xup);
		}

		bgq_su3_accum_weyl_xup(result, weyl_xup);
	}


	// X- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_spinor_decl(spinor_xdown);
		bgq_su3_spinor_load(spinor_xdown, (GETSPINORFIELDPTR(XDOWN, xdown)));

		bgq_su3_weyl_decl(weyl_xdown);
		bgq_su3_reduce_weyl_xdown(weyl_xdown, spinor_xdown);

		bgq_su3_mdecl(gauge_xdown);
		bgq_qvlfduxa(gauge_xdown_c00, gaugesite, 32);
		bgq_qvlfduxa(gauge_xdown_c01, gaugesite, 32);
		bgq_qvlfduxa(gauge_xdown_c02, gaugesite, 32);
		bgq_qvlfduxa(gauge_xdown_c10, gaugesite, 32);
		bgq_qvlfduxa(gauge_xdown_c11, gaugesite, 32);
		bgq_qvlfduxa(gauge_xdown_c12, gaugesite, 32);
		bgq_qvlfduxa(gauge_xdown_c20, gaugesite, 32);
		bgq_qvlfduxa(gauge_xdown_c21, gaugesite, 32);
		bgq_qvlfduxa(gauge_xdown_c22, gaugesite, 32);

		bgq_su3_weyl_mvmul(weyl_xdown, gauge_xdown, weyl_xdown);

		if (kamul) {
			bgq_su3_weyl_cmul(weyl_xdown, qka1, weyl_xdown);
		}

		bgq_su3_accum_weyl_xdown(result, weyl_xdown);
	}


	// Y+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_spinor_decl(spinor_yup);
		bgq_su3_spinor_load(spinor_yup, (GETSPINORFIELDPTR(YUP, yup)));

		bgq_su3_weyl_decl(weyl_yup);
		bgq_su3_reduce_weyl_yup(weyl_yup, spinor_yup);

		bgq_su3_mdecl(gauge_yup);
		bgq_qvlfduxa(gauge_yup_c00, gaugesite, 32);
		bgq_qvlfduxa(gauge_yup_c01, gaugesite, 32);
		bgq_qvlfduxa(gauge_yup_c02, gaugesite, 32);
		bgq_qvlfduxa(gauge_yup_c10, gaugesite, 32);
		bgq_qvlfduxa(gauge_yup_c11, gaugesite, 32);
		bgq_qvlfduxa(gauge_yup_c12, gaugesite, 32);
		bgq_qvlfduxa(gauge_yup_c20, gaugesite, 32);
		bgq_qvlfduxa(gauge_yup_c21, gaugesite, 32);
		bgq_qvlfduxa(gauge_yup_c22, gaugesite, 32);

		bgq_su3_weyl_mvmul(weyl_yup, gauge_yup, weyl_yup);

		if (kamul) {
			bgq_su3_weyl_cmul(weyl_yup, qka2, weyl_yup);
		}

		bgq_su3_accum_weyl_yup(result, weyl_yup);
	}


	// Y- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_spinor_decl(spinor_ydown);
		bgq_su3_spinor_load(spinor_ydown, (GETSPINORFIELDPTR(YDOWN, ydown)));

		bgq_su3_weyl_decl(weyl_ydown);
		bgq_su3_reduce_weyl_ydown(weyl_ydown, spinor_ydown);

		bgq_su3_mdecl(gauge_ydown);
		bgq_qvlfduxa(gauge_ydown_c00, gaugesite, 32);
		bgq_qvlfduxa(gauge_ydown_c01, gaugesite, 32);
		bgq_qvlfduxa(gauge_ydown_c02, gaugesite, 32);
		bgq_qvlfduxa(gauge_ydown_c10, gaugesite, 32);
		bgq_qvlfduxa(gauge_ydown_c11, gaugesite, 32);
		bgq_qvlfduxa(gauge_ydown_c12, gaugesite, 32);
		bgq_qvlfduxa(gauge_ydown_c20, gaugesite, 32);
		bgq_qvlfduxa(gauge_ydown_c21, gaugesite, 32);
		bgq_qvlfduxa(gauge_ydown_c22, gaugesite, 32);

		bgq_su3_weyl_mvmul(weyl_ydown, gauge_ydown, weyl_ydown);

		if (kamul) {
			bgq_su3_weyl_cmul(weyl_ydown, qka2, weyl_ydown);
		}

		bgq_su3_accum_weyl_ydown(result, weyl_ydown);
	}


	// Z+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_spinor_decl(spinor_zup);
		bgq_su3_spinor_load(spinor_zup, (GETSPINORFIELDPTR(ZUP, zup)));

		bgq_su3_weyl_decl(weyl_zup);
		bgq_su3_reduce_weyl_zup(weyl_zup, spinor_zup);

		bgq_su3_mdecl(gauge_zup);
		bgq_qvlfduxa(gauge_zup_c00, gaugesite, 32);
		bgq_qvlfduxa(gauge_zup_c01, gaugesite, 32);
		bgq_qvlfduxa(gauge_zup_c02, gaugesite, 32);
		bgq_qvlfduxa(gauge_zup_c10, gaugesite, 32);
		bgq_qvlfduxa(gauge_zup_c11, gaugesite, 32);
		bgq_qvlfduxa(gauge_zup_c12, gaugesite, 32);
		bgq_qvlfduxa(gauge_zup_c20, gaugesite, 32);
		bgq_qvlfduxa(gauge_zup_c21, gaugesite, 32);
		bgq_qvlfduxa(gauge_zup_c22, gaugesite, 32);

		bgq_su3_weyl_mvmul(weyl_zup, gauge_zup, weyl_zup);

		if (kamul) {
			bgq_su3_weyl_cmul(weyl_zup, qka1, weyl_zup);
		}

		bgq_su3_accum_weyl_zup(result, weyl_zup);
	}


	// Z- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_spinor_decl(spinor_zdown);
		bgq_su3_spinor_load(spinor_zdown, (GETSPINORFIELDPTR(ZDOWN, zdown)));

		bgq_su3_weyl_decl(weyl_zdown);
		bgq_su3_reduce_weyl_zdown(weyl_zdown, spinor_zdown);

		bgq_su3_mdecl(gauge_zdown);
		bgq_qvlfduxa(gauge_zdown_c00, gaugesite, 32);
		bgq_qvlfduxa(gauge_zdown_c01, gaugesite, 32);
		bgq_qvlfduxa(gauge_zdown_c02, gaugesite, 32);
		bgq_qvlfduxa(gauge_zdown_c10, gaugesite, 32);
		bgq_qvlfduxa(gauge_zdown_c11, gaugesite, 32);
		bgq_qvlfduxa(gauge_zdown_c12, gaugesite, 32);
		bgq_qvlfduxa(gauge_zdown_c20, gaugesite, 32);
		bgq_qvlfduxa(gauge_zdown_c21, gaugesite, 32);
		bgq_qvlfduxa(gauge_zdown_c22, gaugesite, 32);

		bgq_su3_weyl_mvmul(weyl_zdown, gauge_zdown, weyl_zdown);

		if (kamul) {
			bgq_su3_weyl_cmul(weyl_zdown, qka1, weyl_zdown);
		}

		bgq_su3_accum_weyl_zdown(result, weyl_zdown);
	}


	bgq_su3_spinor_store(target, result);

}
#undef BGQ_HOPPINGMATRIXSTENCIL_INC_
