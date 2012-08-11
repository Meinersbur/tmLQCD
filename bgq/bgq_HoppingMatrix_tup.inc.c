
#ifndef BGQ_HM_TUP_FIRST
#define BGQ_HM_TUP_FIRST 0
#endif

#ifndef BGQ_HM_TUP_TLINEINDENT
#define BGQ_HM_TUP_TLINEINDENT -1
#endif

#ifndef BGQ_HM_TUP_COMPUTE
#define BGQ_HM_TUP_COMPUTE 0
#endif

#ifndef BGQ_HM_TUP_ACCUMULATE
#define BGQ_HM_TUP_ACCUMULATE 0
#endif


#ifndef BGQ_HM_DIR_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void HoppingMatrix_tup(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int tv, int x, int y, int z, int t1, int t2) {
	bgq_su3_spinor_decl(result);
#else
{
#endif


	// Load the input spinor
	bgq_su3_spinor_decl(spinor_tup);
#if (BGQ_HM_TUP_TLINEINDENT==-1)
	if (((tv+x+y+z+isOdd)&1) == 0) {
#endif
#if (BGQ_HM_TUP_TLINEINDENT==-1) || (BGQ_HM_TUP_TLINEINDENT==0)
		bgq_spinorsite_double *spinorsite_tup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z, t1+1, t2+1, true,false);
		bgq_su3_spinor_double_load(spinor_tup, spinorsite_tup);
#endif
#if (BGQ_HM_TUP_TLINEINDENT==-1)
	} else {
#endif
#if (BGQ_HM_TUP_TLINEINDENT==-1) || (BGQ_HM_TUP_TLINEINDENT==1)
		bgq_spinorsite_double *spinorsite_tup_mid = BGQ_SPINORSITE_LEFT(spinorfield, !isOdd, tv, x, y, z, t1-1,t1+1, true,false);
		bgq_su3_spinor_decl_rightonly(spinor_tup_mid);
		bgq_su3_spinor_double_load_right_torightonly(spinor_tup_mid, spinorsite_tup_mid);

		bgq_spinorsite_double *spinorsite_tup_right = BGQ_SPINORSITE_RIGHT(spinorfield, !isOdd, tv+1, x, y, z, t2+1,t2+3, true,false);
		bgq_su3_spinor_decl_leftonly(spinor_tup_right);
		bgq_su3_spinor_double_load_left_toleftonly(spinor_tup_right, spinorsite_tup_right);

		bgq_su3_spinor_merge(spinor_tup, spinor_tup_mid, spinor_tup_right);
#endif
#if (BGQ_HM_TUP_TLINEINDENT==-1)
	}
#endif


	// Compute its halfspinor
	bgq_su3_weyl_decl(weyl_tup);
	bgq_su3_vadd(weyl_tup_v0, spinor_tup_v0, spinor_tup_v2);
	bgq_su3_vadd(weyl_tup_v1, spinor_tup_v1, spinor_tup_v3);


#if BGQ_HM_TUP_COMPUTE
	bgq_gaugesite_double *gaugesite_tup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, TUP_SHIFT, t1, t2, true,false);
	bgq_su3_mdecl(gauge_tup);
	bgq_su3_matrix_double_load(gauge_tup, gaugesite_tup);

	bgq_su3_mvmul(weyl_tup_v0, gauge_tup, weyl_tup_v0);
	bgq_su3_mvmul(weyl_tup_v1, gauge_tup, weyl_tup_v1);

#ifndef BGQ_HM_NOKAMUL
	bgq_su3_cvmul(weyl_tup_v0, qka0, weyl_tup_v0);
	bgq_su3_cvmul(weyl_tup_v1, qka0, weyl_tup_v1);
#endif
#endif


#if BGQ_HM_TUP_ACCUMULATE
#if BGQ_HM_TUP_FIRST
		bgq_su3_vmov(result_v0, weyl_tup_v0);
		bgq_su3_vmov(result_v1, weyl_tup_v1);
		bgq_su3_vmov(result_v2, weyl_tup_v0);
		bgq_su3_vmov(result_v3, weyl_tup_v1);
#else
		bgq_su3_vadd(result_v0, result_v0, weyl_tup_v0);
		bgq_su3_vadd(result_v1, result_v1, weyl_tup_v1);
		bgq_su3_vadd(result_v2, result_v2, weyl_tup_v0);
		bgq_su3_vadd(result_v3, result_v3, weyl_tup_v1);
#endif
#endif


}


#undef BGQ_HM_TUP_FIRST
#undef BGQ_HM_TUP_TLINEINDENT
#undef BGQ_HM_TUP_COMPUTE
#undef BGQ_HM_TUP_ACCUMULATE

