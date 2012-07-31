
#ifndef BGQ_HM_ZUP_WEYLREAD
#define BGQ_HM_ZUP_WEYLREAD 0
#endif

#ifndef BGQ_HM_ZUP_COMPUTE
#define BGQ_HM_ZUP_COMPUTE 0
#endif

#ifndef BGQ_HM_ZUP_WEYL_SEND
#define BGQ_HM_ZUP_WEYL_SEND 0
#endif

#ifndef BGQ_HM_ZUP_ACCUMULATE
#define BGQ_HM_ZUP_ACCUMULATE 0
#endif


#ifndef BGQ_HM_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void HoppingMatrix_site(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefieldeo_double gaugefield, bool isOdd, int x, int y, int z, int tv, int k) {
	bgq_su3_spinor_decl(result);
#endif

	{
		bgq_su3_weyl_decl(weyl_zup);
#if BGQ_HM_ZUP_WEYLREAD==-1
		if (z==PHYSICAL_LZ-1) {
#endif
#if (BGQ_HM_ZUP_WEYLREAD==-1) || (BGQ_HM_ZUP_WEYLREAD==1)
		bgq_weylsite_double *weylsite_zup = BGQ_WEYLSITE_Z(weylxchange_yup_recv_double, !isOdd, x, y, z+1, tv);
		bgq_su3_weyl_double_load(weyl_zup, weylsite_zup);
#endif
#if BGQ_HM_ZUP_WEYLREAD==-1
		} else {
#endif
#if (BGQ_HM_ZUP_WEYLREAD==-1) || (BGQ_HM_ZUP_WEYLREAD==0)
		// Load the input spinor
		bgq_su3_spinor_decl(spinor_zup);
		bgq_spinorsite_double *spinorsite_zup = BGQ_SPINORSITE(spinorfield, !isOdd, x, y, z+1, tv);
		bgq_su3_spinor_double_load(spinor_zup, spinorsite_zup);

		// Compute its halfspinor
		bgq_su3_vadd(weyl_zup_v0, spinor_zup_v0, spinor_zup_v3);
		bgq_su3_vsub(weyl_zup_v1, spinor_zup_v1, spinor_zup_v2);
#endif
#if BGQ_HM_ZUP_WEYLREAD==-1
		}
#endif

#if BGQ_HM_ZUP_ACCUMULATE
		bgq_gaugesite_double *gaugesite_zup = BGQ_GAUGESITE(gaugefield, isOdd, x, y, z, tv, Z_UP);
		bgq_su3_mdecl(gauge_zup);
		bgq_su3_matrix_double_load(gauge_zup, gaugesite_zup);

		bgq_su3_mvmul(weyl_zup_v0, gauge_zup, weyl_zup_v0);
		bgq_su3_mvmul(weyl_zup_v1, gauge_zup, weyl_zup_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(weyl_zup_v0, qka2, weyl_zup_v0);
		bgq_su3_cvmul(weyl_zup_v1, qka2, weyl_zup_v1);
#endif
#endif

#if BGQ_HM_ZUP_WEYL_SEND
		// Store the halfspinor to be transfered to the neighbor node
		bgq_su3_weyl_double_store(weylxchange_zup_send_double, weyl_zup);
#endif

#if BGQ_HM_ZUP_ACCUMULATE
		bgq_su3_vadd(result_v0, result_v0, weyl_zup_v0);
		bgq_su3_vadd(result_v1, result_v1, weyl_zup_v1);
		bgq_su3_vsub(result_v2, result_v2, weyl_zup_v1);
		bgq_su3_vadd(result_v3, result_v3, weyl_zup_v0);
#endif
	}

#ifndef BGQ_HM_NOFUNC
}
#endif

#undef BGQ_HM_ZUP_WEYLREAD
#undef BGQ_HM_ZUP_COMPUTE
#undef BGQ_HM_ZUP_WEYL_SEND
#undef BGQ_HM_ZUP_ACCUMULATE
