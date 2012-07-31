
// BGQ_HM_TLINE_FLUSHLINE
// BGQ_HM_TLINE_RAGGEDLINE

#ifndef BGQ_HM_TUP_TLINEINDENT
#error Must define the line indention (0, 1 or -1 for runtime-conditional)
#define BGQ_HM_TUP_TLINEINDENT -1
#endif

#ifndef BGQ_HM_TUP_RIGHTWRAPAROUND
#define BGQ_HM_TUP_RIGHTWRAPAROUND 0
#endif

#ifndef BGQ_HM_TUP_COMPUTE
#define BGQ_HM_TUP_COMPUTE 1
#endif

#ifndef BGQ_HM_TUP_ACCUMULATE
#define BGQ_HM_TUP_ACCUMULATE 0
#endif

#ifndef BGQ_HM_TUP_WRITECARRYSPINOR
#define BGQ_HM_TUP_WRITECARRYSPINOR 0
#endif

#ifndef BGQ_HM_DIR_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void HoppingMatrix_site_tup(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int x, int y, int z, int tv, int k) {
	bgq_su3_spinor_decl(result);
#endif

	{
#if (BGQ_HM_TUP_TLINEINDENT==-1) || (BGQ_HM_TUP_TLINEINDENT==1)
#if BGQ_HM_TUP_RIGHTWRAPAROUND
		const int tv_right = 0;
#else
		const int tv_right = tv+1;
#endif
#endif

		// Load the input spinor
		bgq_su3_spinor_decl(spinor_tup);
		// # = stencil site to update (either even or odd sites, depending on what isOdd says)
		// _ = neighbor sites to read (therefore !isOdd)
		// (  2  ) = vector site with tv-number
		// top line = vector sites of isOdd sites
		// bottom line = vector sites of !isOdd sites

#if BGQ_HM_TUP_TLINEINDENT==-1
		if ((x+y+z)%2 == isOdd) {
#endif
#if (BGQ_HM_TUP_TLINEINDENT==-1) || (BGQ_HM_TUP_TLINEINDENT==0)
		// (  0  ) (  1  ) (  2  )
		// |# _ # _ # _ # _ # _ # _|
		//   (  0  ) (  1  ) (  2  )
		// T_UP = tv
		// T_DOWN = half tv, half tv-1
		assert((x+y+z)%2 == isOdd);

		bgq_spinorsite_double *spinorsite_tup = BGQ_SPINORSITE(spinorfield, !isOdd, x, y, z, tv);
		bgq_su3_spinor_double_load(spinor_tup, spinorsite_tup);
#endif
#if BGQ_HM_TUP_TLINEINDENT==-1
		} else {
#endif
#if (BGQ_HM_TUP_TLINEINDENT==-1) || (BGQ_HM_TUP_TLINEINDENT==1)
		//   (  0  ) (  1  ) (  2  )
		// |_ # _ # _ # _ # _ # _ #|
		// (  0  ) (  1  ) (  2  )
		// T_UP = half tv, half tv+1
		// T_DOWN = tv
		assert((x+y+z)%2 == !isOdd);

		bgq_su3_spinor_decl_rightonly(spinor_tup_mid);
		bgq_spinorsite_double *spinorsite_tup_mid = BGQ_SPINORSITE(spinorfield, !isOdd, x, y, z, tv);
		bgq_su3_spinor_double_load_right_torightonly(spinor_tup_mid, spinorsite_tup_mid);

		bgq_su3_spinor_decl_leftonly(spinor_tup_right);
		bgq_spinorsite_double *spinorsite_tup_right = BGQ_SPINORSITE(spinorfield, !isOdd, x, y, z, tv_right);
		bgq_su3_spinor_double_load_left_toleftonly(spinor_tup_right, spinorsite_tup_right);

		bgq_su3_spinor_merge(spinor_tup, spinor_tup_mid, spinor_tup_right);
#endif
#if BGQ_HM_TUP_TLINEINDENT==-1
		}
#endif

		// Compute its halfspinor
		bgq_su3_weyl_decl(weyl_tup);
		bgq_su3_vpiadd(weyl_tup_v0, spinor_tup_v0, spinor_tup_v2);
		bgq_su3_vpisub(weyl_tup_v1, spinor_tup_v1, spinor_tup_v3);

#if BGQ_HM_TUP_COMPUTE
		bgq_su3_mdecl(gauge_tup);

		bgq_gaugesite_double *gaugesite_tup;
#if BGQ_HM_TUP_TLINEINDENT==-1
		if ((x+y+z)%2 == isOdd) {
#endif
#if (BGQ_HM_TUP_TLINEINDENT==-1) || (BGQ_HM_TUP_TLINEINDENT==0)
		gaugesite_tup = BGQ_GAUGESITE(gaugefield, isOdd, x, y, z, tv, T_UP);
#endif
#if BGQ_HM_TUP_TLINEINDENT==-1
		} else {
#endif
#if (BGQ_HM_TUP_TLINEINDENT==-1) || (BGQ_HM_TUP_TLINEINDENT==1)
		gaugesite_tup = BGQ_GAUGESITE(gaugefield, isOdd, x, y, z, tv_right, T_RAGGED_UP);
#endif
#if BGQ_HM_TUP_TLINEINDENT==-1
		}
#endif
		bgq_su3_matrix_double_load(gauge_tup, gaugesite_tup);

		bgq_su3_mvmul(weyl_tup_v0, gauge_tup, weyl_tup_v0);
		bgq_su3_mvmul(weyl_tup_v1, gauge_tup, weyl_tup_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(weyl_tup_v0, qka3, weyl_tup_v0);
		bgq_su3_cvmul(weyl_tup_v1, qka3, weyl_tup_v1);
#endif
#endif

#if BGQ_HM_TUP_ACCUMULATE
		bgq_su3_vadd(result_v0, result_v0, weyl_tup_v0);
		bgq_su3_vadd(result_v1, result_v1, weyl_tup_v1);
		bgq_su3_vpisub(result_v2, result_v2, weyl_tup_v0);
		bgq_su3_vpiadd(result_v3, result_v3, weyl_tup_v1);
#endif

#if BGQ_HM_TUP_WRITECARRYSPINOR
		bgq_su3_spinor_mov(spinor_tcarry, spinor_tup);
#endif
	}

#ifndef BGQ_HM_DIR_NOFUNC
}
#endif

#undef BGQ_HM_TUP_TLINEINDENT
#undef BGQ_HM_TUP_RIGHTWRAPAROUND
#undef BGQ_HM_TUP_COMPUTE
#undef BGQ_HM_TUP_ACCUMULATE
#undef BGQ_HM_TUP_WRITECARRYSPINOR
