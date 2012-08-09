
#ifndef BGQ_HM_ZUP_RIGHTWRAPAROUND
#define BGQ_HM_ZUP_RIGHTWRAPAROUND -1 /* Runtime-conditional */
#endif

#ifndef BGQ_HM_ZUP_COMPUTE
#define BGQ_HM_ZUP_COMPUTE 0
#endif

#ifndef BGQ_HM_ZUP_ACCUMULATE
#define BGQ_HM_ZUP_ACCUMULATE 0
#endif


#ifndef BGQ_HM_DIR_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void bgq_HoppingMatrix_zup(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int tv, int x, int y, int z, int t1, int t2) {
	bgq_su3_spinor_decl(result);
#endif
	{


#if (BGQ_HM_ZUP_RIGHTWRAPAROUND==0)
		assert(z!=PHYSICAL_LZ-1);
		const int z_right = z+1;
#elif (BGQ_HM_ZUP_RIGHTWRAPAROUND==1)
		assert(z==PHYSICAL_LZ-1);
		const int z_right = 0;
#else
		const int z_right = mod(z+1, PHYSICAL_LZ);
#endif


		// Load the input spinor
		bgq_su3_spinor_decl(spinor_zup);
		bgq_spinorsite_double *spinorsite_zup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z_right, t1,t2, false);
		bgq_su3_spinor_double_load(spinor_zup, spinorsite_zup);


		// Compute its halfspinor
		bgq_su3_weyl_decl(weyl_zup);
		bgq_su3_vpiadd(weyl_zup_v0, spinor_zup_v0, spinor_zup_v2);
		bgq_su3_vpisub(weyl_zup_v1, spinor_zup_v1, spinor_zup_v3);


#if BGQ_HM_ZUP_COMPUTE
		bgq_gaugesite_double *gaugesite_zup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, ZUP, t1,t2, true,false);
		bgq_su3_mdecl(gauge_zup);
		bgq_su3_matrix_double_load(gauge_zup, gaugesite_zup);

		bgq_su3_mvmul(weyl_zup_v0, gauge_zup, weyl_zup_v0);
		bgq_su3_mvmul(weyl_zup_v1, gauge_zup, weyl_zup_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(weyl_zup_v0, qka3, weyl_zup_v0);
		bgq_su3_cvmul(weyl_zup_v1, qka3, weyl_zup_v1);
#endif
#endif


#if BGQ_HM_ZUP_ACCUMULATE
		bgq_su3_vadd(result_v0, result_v0, weyl_zup_v0);
		bgq_su3_vadd(result_v1, result_v1, weyl_zup_v1);
		bgq_su3_vpisub(result_v2, result_v2, weyl_zup_v0);
		bgq_su3_vpiadd(result_v3, result_v3, weyl_zup_v1);
#endif


	}
#ifndef BGQ_HM_DIR_NOFUNC
}
#endif


#undef BGQ_HM_ZUP_COMPUTE
#undef BGQ_HM_ZUP_ACCUMULATE
