#ifndef BGQ_HM_TDOWN_WEYLREAD
#define BGQ_HM_TDOWN_WEYLREAD 0
#endif

#ifndef BGQ_HM_TDOWN_COMPUTE
#define BGQ_HM_TDOWN_COMPUTE 0
#endif

#ifndef BGQ_HM_TDOWN_WEYL_SEND
#define BGQ_HM_TDOWN_WEYL_SEND 0
#endif

#ifndef BGQ_HM_TDOWN_ACCUMULATE
#define BGQ_HM_TDOWN_ACCUMULATE 0
#endif

#ifndef BGQ_HM_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void bgq_HoppingMatrix_site_tdown(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int t, int x, int y, int zv, int z1, int z2) {
	bgq_su3_spinor_decl(result);
#endif
	{

		bgq_su3_weyl_decl(weyl_tdown);
#if BGQ_HM_TDOWN_WEYLREAD==-1
		if (x==0) {
#endif
#if (BGQ_HM_TDOWN_WEYLREAD==-1) || (BGQ_HM_TDOWN_WEYLREAD==1)
		bgq_weylsite_double *weylsite_xdown = BGQ_WEYLSITE_X(weylxchange_recv_double[T_DOWN], !isOdd, x-1, y, z, tv);
		bgq_su3_weyl_double_load(weyl_xdown, weylsite_xdown);
#endif
#if BGQ_HM_TDOWN_WEYLREAD==-1
	} else {
#endif
#if (BGQ_HM_TDOWN_WEYLREAD==-1) || (BGQ_HM_TDOWN_WEYLREAD==0)
		// Load the input spinor
		bgq_su3_spinor_decl(spinor_tdown);
		bgq_spinorsite_double *spinorsite_tdown = BGQ_SPINORSITE(spinorfield, !isOdd, t-1, x, y, zv, z1,z2);
		bgq_su3_spinor_double_load(spinor_tdown, spinorsite_tdown);

		// Compute its halfspinor
		bgq_su3_vsub(weyl_tdown_v0, spinor_tdown_v0, spinor_tdown_v2);
		bgq_su3_vsub(weyl_tdown_v1, spinor_tdown_v1, spinor_tdown_v3);
#endif
#if BGQ_HM_TDOWN_WEYLREAD==-1
	}
#endif

#if BGQ_HM_TDOWN_COMPUTE
		bgq_su3_mdecl(gauge_tdown);
		bgq_gaugesite_double *gaugesite_tdown = BGQ_GAUGESITE(gaugefield, !isOdd, t-1, x, y, zv, T_UP);
		bgq_su3_matrix_double_load(gauge_tdown, gaugesite_tdown);

		bgq_su3_mvinvmul(weyl_tdown_v0, gauge_tdown, weyl_tdown_v0);
		bgq_su3_mvinvmul(weyl_tdown_v1, gauge_tdown, weyl_tdown_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(weyl_tdown_v0, qka3, weyl_tdown_v0);
		bgq_su3_cvmul(weyl_tdown_v1, qka3, weyl_tdown_v1);
#endif
#endif


#if BGQ_HM_TDOWN_WEYL_SEND
		// Store the halfspinor to be transfered to the neighbor node
		bgq_su3_weyl_double_store(weylxchange_send_double[T_UP], weyl_tdown);
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


#undef BGQ_HM_TDOWN_WEYLREAD
#undef BGQ_HM_TDOWN_COMPUTE
#undef BGQ_HM_TDOWN_WEYL_SEND
#undef BGQ_HM_TDOWN_ACCUMULATE

