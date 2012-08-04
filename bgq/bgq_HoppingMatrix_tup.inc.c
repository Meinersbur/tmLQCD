
#ifndef BGQ_HM_TUP_WEYLREAD
#define BGQ_HM_TUP_WEYLREAD 0
#endif

#ifndef BGQ_HM_TUP_COMPUTE
#define BGQ_HM_TUP_COMPUTE 0
#endif

#ifndef BGQ_HM_TUP_WEYL_SEND
#define BGQ_HM_TUP_WEYL_SEND 0
#endif

#ifndef BGQ_HM_TUP_ACCUMULATE
#define BGQ_HM_TUP_ACCUMULATE 0
#endif

#ifndef BGQ_HM_DIR_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void HoppingMatrix_site_tup(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int t, int x, int y, int z, int zv, int z1, int z2) {
	bgq_su3_spinor_decl(result);
#endif
	{


		bgq_su3_weyl_decl(weyl_tup);
#if BGQ_HM_TUP_WEYLREAD==-1
		if (t==PHYSICAL_LT-1) {
#endif
#if (BGQ_HM_TUP_WEYLREAD==-1) || (BGQ_HM_TUP_WEYLREAD==1)
			bgq_weylsite_double *weylsite_tup = BGQ_WEYLSITE_T(weylxchange_recv_double[T_UP], !isOdd, t+1, x, y, zv);
			bgq_su3_weyl_double_load(weyl_tup, weylsite_tup);
#endif
#if BGQ_HM_TUP_WEYLREAD==-1
		} else {
#endif
#if (BGQ_HM_TUP_WEYLREAD==-1) || (BGQ_HM_TUP_WEYLREAD==0)
			// Load the input spinor
			bgq_su3_spinor_decl(spinor_tup);
			bgq_spinorsite_double *spinorsite_tup = BGQ_SPINORSITE(spinorfield, !isOdd, t+1, x, y, zv, z1,z2);
			bgq_su3_spinor_double_load(spinor_tup, spinorsite_tup);

			// Compute its halfspinor
			bgq_su3_vadd(weyl_tup_v0, spinor_tup_v0, spinor_tup_v2);
			bgq_su3_vadd(weyl_tup_v1, spinor_tup_v1, spinor_tup_v3);
#endif
#if (BGQ_HM_TUP_WEYLREAD==-1)
		}
#endif


#if BGQ_HM_TUP_COMPUTE
		bgq_su3_mdecl(gauge_tup);
		bgq_gaugesite_double *gaugesite_tup = BGQ_GAUGESITE(gaugefield, isOdd, t, x, y, zv, T_UP);
		bgq_su3_matrix_double_load(gauge_tup, gaugesite_tup);

		bgq_su3_mvmul(weyl_tup_v0, gauge_tup, weyl_tup_v0);
		bgq_su3_mvmul(weyl_tup_v1, gauge_tup, weyl_tup_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(weyl_tup_v0, qka0, weyl_tup_v0);
		bgq_su3_cvmul(weyl_tup_v1, qka0, weyl_tup_v1);
#endif
#endif


#if BGQ_HM_TUP_WEYL_SEND
		// Store the halfspinor to be transfered to the neighbor node
		bgq_su3_weyl_double_store(weylxchange_send_double[T_DOWN], weyl_tup);
#endif


#if BGQ_HM_TUP_ACCUMULATE
		bgq_su3_vmov(result_v0, weyl_tup_v0);
		bgq_su3_vmov(result_v1, weyl_tup_v1);
		bgq_su3_vmov(result_v2, weyl_tup_v0);
		bgq_su3_vmov(result_v3, weyl_tup_v1);
#endif


	}
#ifndef BGQ_HM_DIR_NOFUNC
}
#endif


#undef BGQ_HM_TUP_WEYLREAD
#undef BGQ_HM_TUP_WEYL_SEND
#undef BGQ_HM_TUP_COMPUTE
#undef BGQ_HM_TUP_ACCUMULATE
