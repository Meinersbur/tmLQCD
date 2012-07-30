
#ifndef BGQ_HM_ZDOWN_WEYLREAD
#define BGQ_HM_ZDOWN_WEYLREAD 0
#endif

#ifndef BGQ_HM_ZDOWN_COMPUTE
#define BGQ_HM_ZDOWN_COMPUTE 0
#endif

#ifndef BGQ_HM_ZDOWN_WEYL_SEND
#define BGQ_HM_ZDOWN_WEYL_SEND 0
#endif

#ifndef BGQ_HM_ZDOWN_ACCUMULATE
#define BGQ_HM_ZDOWN_ACCUMULATE 0
#endif


#ifndef BGQ_HM_DIR_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void HoppingMatrix_site(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int x, int y, int z, int tv, int k) {
	bgq_su3_spinor_decl(result);
#endif

	{
		bgq_su3_weyl_decl(weyl_zdown);
#if BGQ_HM_ZDOWN_WEYLREAD
		bgq_weylsite_double *weylsite_zdown = BGQ_WEYLSITE_Z(weylxchange_xdown_recv_double, !isOdd, x, y, z-1, tv);
		bgq_su3_weyl_double_load(weyl_zdown, weylsite_zdown);
#else
		// Load the input spinor
		bgq_su3_spinor_decl(spinor_zdown);
		bgq_spinorsite_double *spinorsite_zdown = BGQ_SPINORSITE(spinorfield, !isOdd, x, y, z-1, tv);
		bgq_su3_spinor_double_load(spinor_zdown, spinorsite_zdown);

		// Compute its halfspinor
		bgq_su3_vsub(weyl_zdown_v0, spinor_zdown_v0, spinor_zdown_v3);
		bgq_su3_vadd(weyl_zdown_v1, spinor_zdown_v1, spinor_zdown_v2);
#endif

#if BGQ_HM_ZDOWN_COMPUTE
		bgq_su3_mdecl(gauge_zdown);
		bgq_gaugesite_double *gaugesite_zdown = BGQ_GAUGESITE(gaugefield, !isOdd, x, y, z-1, tv, Z_UP);
		bgq_su3_matrix_double_load(gauge_zdown, gaugesite_zdown);

		bgq_su3_mvinvmul(weyl_zdown_v0, gauge_zdown, weyl_zdown_v0);
		bgq_su3_mvinvmul(weyl_zdown_v1, gauge_zdown, weyl_zdown_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(weyl_zdown_v0, qka2, weyl_zdown_v0);
		bgq_su3_cvmul(weyl_zdown_v1, qka2, weyl_zdown_v1);
#endif
#endif

#if BGQ_HM_ZDOWN_WEYL_SEND
		// Store the halfspinor to be transfered to the neighbor node
		bgq_su3_weyl_double_store(weylxchange_zdown_send_double, weyl_zdown);
#endif

#if BGQ_HM_ZDOWN_ACCUMULATE
		bgq_su3_vadd(result_v0, result_v0, weyl_zdown_v0);
		bgq_su3_vadd(result_v1, result_v1, weyl_zdown_v1);
		bgq_su3_vadd(result_v2, result_v2, weyl_zdown_v1);
		bgq_su3_vsub(result_v3, result_v3, weyl_zdown_v0);
#endif
	}

#ifndef BGQ_HM_DIR_NOFUNC
}
#endif

#undef BGQ_HM_ZDOWN_WEYLREAD
#undef BGQ_HM_ZDOWN_COMPUTE
#undef BGQ_HM_ZDOWN_WEYL_SEND
#undef BGQ_HM_ZDOWN_ACCUMULATE
