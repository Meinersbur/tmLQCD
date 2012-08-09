
#ifndef BGQ_HM_TDOWN_TLINEINDENT
#define BGQ_HM_TDOWN_TLINEINDENT -1
#endif

#ifndef BGQ_HM_TDOWN_COMPUTE
#define BGQ_HM_TDOWN_COMPUTE 0
#endif

#ifndef BGQ_HM_TDOWN_ACCUMULATE
#define BGQ_HM_TDOWN_ACCUMULATE 0
#endif


#ifndef BGQ_HM_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void bgq_HoppingMatrix_tdown(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int tv, int x, int y, int z, int t1, int t2) {
	bgq_su3_spinor_decl(result);
#endif
	{

		// Load the input spinor
	bgq_su3_spinor_decl(spinor_tdown);
	#if (BGQ_HM_TDOWN_TLINEINDENT==-1)
		if (((tv+x+y+z+isOdd)&1) == 0) {
	#endif
	#if (BGQ_HM_TDOWN_TLINEINDENT==-1) || (BGQ_HM_TDOWN_TLINEINDENT==0)
			bgq_spinorsite_double *spinorsite_tdown_left = BGQ_SPINORSITE_RIGHT(spinorfield, !isOdd, tv-1, x, y, z, t1-3,t1-1, false);
			bgq_su3_spinor_decl_rightonly(spinor_tdown_left);
			bgq_su3_spinor_double_load_right_torightonly(spinor_tdown_left, spinorsite_tdown_left);

			bgq_spinorsite_double *spinorsite_tdown_mid = BGQ_SPINORSITE_LEFT(spinorfield, !isOdd, tv, x, y, z, t2-1,t2+1, false);
			bgq_su3_spinor_decl_leftonly(spinor_tdown_mid);
			bgq_su3_spinor_double_load_left_toleftonly(spinor_tdown_mid, spinorsite_tdown_mid);

			bgq_su3_spinor_merge(spinor_tdown, spinor_tdown_left, spinor_tdown_mid);
	#endif
	#if (BGQ_HM_TDOWN_TLINEINDENT==-1)
		} else {
	#endif
	#if (BGQ_HM_TDOWN_TLINEINDENT==-1) || (BGQ_HM_TDOWN_TLINEINDENT==1)
			bgq_spinorsite_double *spinorsite_tdown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z, t1-1, t2-1,false);
			bgq_su3_spinor_double_load(spinor_tdown, spinorsite_tdown);
	#endif
	#if (BGQ_HM_TDOWN_TLINEINDENT==-1)
		}
	#endif


		// Compute its halfspinor
		bgq_su3_weyl_decl(weyl_tdown);
		bgq_su3_vsub(weyl_tdown_v0, spinor_tdown_v0, spinor_tdown_v2);
		bgq_su3_vsub(weyl_tdown_v1, spinor_tdown_v1, spinor_tdown_v3);


#if BGQ_HM_TDOWN_COMPUTE
		bgq_su3_mdecl(gauge_tdown);
		bgq_gaugesite_double *gaugesite_tdown = BGQ_GAUGESITE(gaugefield, !isOdd, tv-1, x, y, z, TUP_SHIFT, t1-1, t2-1, true,false);
		bgq_su3_matrix_double_load(gauge_tdown, gaugesite_tdown);

		bgq_su3_mvinvmul(weyl_tdown_v0, gauge_tdown, weyl_tdown_v0);
		bgq_su3_mvinvmul(weyl_tdown_v1, gauge_tdown, weyl_tdown_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(weyl_tdown_v0, qka0, weyl_tdown_v0);
		bgq_su3_cvmul(weyl_tdown_v1, qka0, weyl_tdown_v1);
#endif
#endif


#if BGQ_HM_TDOWN_ACCUMULATE
		bgq_su3_vadd(result_v0, result_v0, weyl_tdown_v0);
		bgq_su3_vadd(result_v1, result_v1, weyl_tdown_v1);
		bgq_su3_vsub(result_v2, result_v2, weyl_tdown_v0);
		bgq_su3_vsub(result_v3, result_v3, weyl_tdown_v1);
#endif


	}
#ifndef BGQ_HM_NOFUNC
}
#endif


#undef BGQ_HM_TDOWN_TLINEINDENT
#undef BGQ_HM_TDOWN_COMPUTE
#undef BGQ_HM_TDOWN_ACCUMULATE

