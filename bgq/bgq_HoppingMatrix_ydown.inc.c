
#ifndef BGQ_HM_YDOWN_WEYLREAD
#define BGQ_HM_YDOWN_WEYLREAD 0
#endif

#ifndef BGQ_HM_YDOWN_COMPUTE
#define BGQ_HM_YDOWN_COMPUTE 0
#endif

#ifndef BGQ_HM_YDOWN_WEYL_SEND
#define BGQ_HM_YDOWN_WEYL_SEND 0
#endif

#ifndef BGQ_HM_YDOWN_ACCUMULATE
#define BGQ_HM_YDOWN_ACCUMULATE 0
#endif


#ifndef BGQ_HM_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void HoppingMatrix_site(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int x, int y, int z, int tv, int k) {
	bgq_su3_spinor_decl(result);
#endif

	{
		bgq_su3_weyl_decl(weyl_ydown);
#if BGQ_HM_YDOWN_WEYLREAD
		bgq_weylsite_double *weylsite_ydown = BGQ_WEYLSITE_Y(weylxchange_xdown_recv_double, !isOdd, x, y-1, z, tv);
		bgq_su3_weyl_double_load(weyl_ydown, weylsite_ydown);
#else
		// Load the input spinor
		bgq_su3_spinor_decl(spinor_ydown);
		bgq_spinorsite_double *spinorsite_ydown = BGQ_SPINORSITE(spinorfield, !isOdd, x, y-1, z, tv);
		bgq_su3_spinor_double_load(spinor_ydown, spinorsite_ydown);

		// Compute its halfspinor

		bgq_su3_vpisub(weyl_ydown_v0, spinor_ydown_v0, spinor_ydown_v3);
		bgq_su3_vpisub(weyl_ydown_v1, spinor_ydown_v1, spinor_ydown_v2);
#endif

#if BGQ_HM_YDOWN_COMPUTE
		bgq_su3_mdecl(gauge_ydown);
		bgq_gaugesite_double *gaugesite_ydown = BGQ_GAUGESITE(gaugefield, !isOdd, x, y-1, z, tv, Y_UP);
		bgq_su3_matrix_double_load(gauge_ydown, gaugesite_ydown);

		bgq_su3_mvinvmul(weyl_ydown_v0, gauge_ydown, weyl_ydown_v0);
		bgq_su3_mvinvmul(weyl_ydown_v1, gauge_ydown, weyl_ydown_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(weyl_ydown_v0, qka1, weyl_ydown_v0);
		bgq_su3_cvmul(weyl_ydown_v1, qka1, weyl_ydown_v1);
#endif
#endif

#if BGQ_HM_YDOWN_WEYL_SEND
		// Store the halfspinor to be transfered to the neighbor node
		bgq_su3_weyl_double_store(weylxchange_ydown_send_double, weyl_ydown);
#endif

#if BGQ_HM_YDOWN_ACCUMULATE
		bgq_su3_vadd(result_v0, result_v0, weyl_ydown_v0);
		bgq_su3_vadd(result_v1, result_v1, weyl_ydown_v1);
		bgq_su3_vpiadd(result_v2, result_v2, weyl_ydown_v1);
		bgq_su3_vpiadd(result_v3, result_v3, weyl_ydown_v0);
#endif

	}

#ifndef BGQ_HM_NOFUNC
}
#endif

#undef BGQ_HM_YDOWN_WEYLREAD
#undef BGQ_HM_YDOWN_COMPUTE
#undef BGQ_HM_YDOWN_WEYL_SEND
#undef BGQ_HM_YDOWN_ACCUMULATE
