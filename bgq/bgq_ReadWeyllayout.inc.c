/*
 * bgq_readWeyllayout.inc.c
 *
 *  Created on: Nov 13, 2012
 *      Author: meinersbur
 */



#ifndef BGQ_READWEYLLAYOUT_INC_
#include "bgq_qpx.h"
#include "bgq_spinorfield.h"

#define PRECISION double

void bgq_readWeyllayout(bgq_su3_spinor_params(/*out*/spinor), bgq_weylsite *weylsite, ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z)
#endif

#ifndef BGQ_READWEYLLAYOUT_INSERTPREFETCH
#define BGQ_READWEYLLAYOUT_INSERTPREFETCH
#endif

{

	bgq_su3_weylnext_prefetch(weylsite);

	// TUP
	{
		bgq_su3_weyl_decl(weylnext_tup);
		bgq_qvlfuxa(weylnext_tup_v0_c0, weylsite, 32);
		bgq_qvlfuxa(weylnext_tup_v0_c1, weylsite, 32);
		bgq_qvlfuxa(weylnext_tup_v0_c2, weylsite, 32);
		bgq_qvlfuxa(weylnext_tup_v1_c0, weylsite, 32);
		bgq_qvlfuxa(weylnext_tup_v1_c1, weylsite, 32);
		bgq_qvlfuxa(weylnext_tup_v1_c2, weylsite, 32);
		bgq_su3_expand_weyl_tup(spinor, weylnext_tup);
	}

	bgq_su3_weylnext_prefetch(weylsite);

	// TDOWN
	{
		bgq_su3_weyl_decl(weylnext_tdown);
		bgq_qvlfuxa(weylnext_tdown_v0_c0, weylsite, 32);
		bgq_qvlfuxa(weylnext_tdown_v0_c1, weylsite, 32);
		bgq_qvlfuxa(weylnext_tdown_v0_c2, weylsite, 32);
		bgq_qvlfuxa(weylnext_tdown_v1_c0, weylsite, 32);
		bgq_qvlfuxa(weylnext_tdown_v1_c1, weylsite, 32);
		bgq_qvlfuxa(weylnext_tdown_v1_c2, weylsite, 32);
		bgq_su3_accum_weyl_tdown(spinor, weylnext_tdown);
	}

	bgq_su3_weylnext_prefetch(weylsite);

	// X+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weylnext_xup);
		bgq_qvlfuxa(weylnext_xup_v0_c0, weylsite, 32);
		bgq_qvlfuxa(weylnext_xup_v0_c1, weylsite, 32);
		bgq_qvlfuxa(weylnext_xup_v0_c2, weylsite, 32);
		bgq_qvlfuxa(weylnext_xup_v1_c0, weylsite, 32);
		bgq_qvlfuxa(weylnext_xup_v1_c1, weylsite, 32);
		bgq_qvlfuxa(weylnext_xup_v1_c2, weylsite, 32);
		bgq_su3_accum_weyl_xup(spinor, weylnext_xup);
	}

	bgq_su3_weylnext_prefetch(weylsite);

	// X- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weylnext_xdown);
		bgq_qvlfuxa(weylnext_xdown_v0_c0, weylsite, 32);
		bgq_qvlfuxa(weylnext_xdown_v0_c1, weylsite, 32);
		bgq_qvlfuxa(weylnext_xdown_v0_c2, weylsite, 32);
		bgq_qvlfuxa(weylnext_xdown_v1_c0, weylsite, 32);
		bgq_qvlfuxa(weylnext_xdown_v1_c1, weylsite, 32);
		bgq_qvlfuxa(weylnext_xdown_v1_c2, weylsite, 32);
		bgq_su3_accum_weyl_xdown(spinor, weylnext_xdown);
	}

	bgq_su3_weylnext_prefetch(weylsite);

	// Y+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weylnext_yup);
		bgq_qvlfuxa(weylnext_yup_v0_c0, weylsite, 32);
		bgq_qvlfuxa(weylnext_yup_v0_c1, weylsite, 32);
		bgq_qvlfuxa(weylnext_yup_v0_c2, weylsite, 32);
		bgq_qvlfuxa(weylnext_yup_v1_c0, weylsite, 32);
		bgq_qvlfuxa(weylnext_yup_v1_c1, weylsite, 32);
		bgq_qvlfuxa(weylnext_yup_v1_c2, weylsite, 32);
		bgq_su3_accum_weyl_yup(spinor, weylnext_yup);
	}

	bgq_su3_weylnext_prefetch(weylsite);

	// Y- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weylnext_ydown);
		bgq_qvlfuxa(weylnext_ydown_v0_c0, weylsite, 32);
		bgq_qvlfuxa(weylnext_ydown_v0_c1, weylsite, 32);
		bgq_qvlfuxa(weylnext_ydown_v0_c2, weylsite, 32);
		bgq_qvlfuxa(weylnext_ydown_v1_c0, weylsite, 32);
		bgq_qvlfuxa(weylnext_ydown_v1_c1, weylsite, 32);
		bgq_qvlfuxa(weylnext_ydown_v1_c2, weylsite, 32);
		bgq_su3_accum_weyl_ydown(spinor, weylnext_ydown);
	}

	bgq_su3_weylnext_prefetch(weylsite);

	// Z+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weylnext_zup);
		bgq_qvlfuxa(weylnext_zup_v0_c0, weylsite, 32);
		bgq_qvlfuxa(weylnext_zup_v0_c1, weylsite, 32);
		bgq_qvlfuxa(weylnext_zup_v0_c2, weylsite, 32);
		bgq_qvlfuxa(weylnext_zup_v1_c0, weylsite, 32);
		bgq_qvlfuxa(weylnext_zup_v1_c1, weylsite, 32);
		bgq_qvlfuxa(weylnext_zup_v1_c2, weylsite, 32);
		bgq_su3_accum_weyl_zup(spinor, weylnext_zup);
	}

	BGQ_READWEYLLAYOUT_INSERTPREFETCH

	// Z- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weylnext_zdown);
		bgq_qvlfuxa(weylnext_zdown_v0_c0, weylsite, 32);
		bgq_qvlfuxa(weylnext_zdown_v0_c1, weylsite, 32);
		bgq_qvlfuxa(weylnext_zdown_v0_c2, weylsite, 32);
		bgq_qvlfuxa(weylnext_zdown_v1_c0, weylsite, 32);
		bgq_qvlfuxa(weylnext_zdown_v1_c1, weylsite, 32);
		bgq_qvlfuxa(weylnext_zdown_v1_c2, weylsite, 32);
		bgq_su3_accum_weyl_zdown(spinor, weylnext_zdown);
	}

}

#undef BGQ_READWEYLLAYOUT_INSERTPREFETCH
#undef BGQ_READWEYLLAYOUT_INC_
