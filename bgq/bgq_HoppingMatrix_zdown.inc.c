
#ifndef BGQ_HM_ZDOWN_LEFTWRAPAROUND
#define BGQ_HM_ZDOWN_LEFTWRAPAROUND -1 /* Runtime-conditional */
#endif

#ifndef BGQ_HM_ZDOWN_COMPUTE
#define BGQ_HM_ZDOWN_COMPUTE 0
#endif

#ifndef BGQ_HM_ZDOWN_ACCUMULATE
#define BGQ_HM_ZDOWN_ACCUMULATE 0
#endif


#ifndef BGQ_HM_DIR_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void bgq_HoppingMatrix_zdown(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int tv, int x, int y, int z, int t1, int t2) {
	bgq_su3_spinor_decl(result);
#endif
	{


#if (BGQ_HM_ZDOWN_LEFTWRAPAROUND==0)
		assert(z!=0);
		const int z_left = z-1;
#elif (BGQ_HM_ZDOWN_LEFTWRAPAROUND==1)
		assert(zv==0);
		const int z_left = PHYSICAL_LZ-1;
#elif (BGQ_HM_ZDOWN_LEFTWRAPAROUND==-1)
		const int z_left = mod(z-1, PHYSICAL_LZ);
#endif


		// Load the input spinor
		bgq_su3_spinor_decl(spinor_zdown);
		bgq_spinorsite_double *spinorsite_zdown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z_left, t1,t2, true,false);
		bgq_su3_spinor_double_load(spinor_zdown, spinorsite_zdown);


		// Compute its halfspinor
		bgq_su3_weyl_decl(weyl_zdown);
		bgq_su3_vpisub(weyl_zdown_v0, spinor_zdown_v0, spinor_zdown_v2);
		bgq_su3_vpiadd(weyl_zdown_v1, spinor_zdown_v1, spinor_zdown_v3);


#if BGQ_HM_ZDOWN_COMPUTE
		bgq_su3_mdecl(gauge_zdown);
		bgq_gaugesite_double *gaugesite_zdown = BGQ_GAUGESITE(gaugefield, !isOdd, tv, x, y, z-1, ZUP, t1,t2, true,false);
		bgq_su3_matrix_double_load(gauge_zdown, gaugesite_zdown);

		bgq_su3_mvinvmul(weyl_zdown_v0, gauge_zdown, weyl_zdown_v0);
		bgq_su3_mvinvmul(weyl_zdown_v1, gauge_zdown, weyl_zdown_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(weyl_zdown_v0, qka3, weyl_zdown_v0);
		bgq_su3_cvmul(weyl_zdown_v1, qka3, weyl_zdown_v1);
#endif
#endif


#if BGQ_HM_ZDOWN_ACCUMULATE
		bgq_su3_vadd(result_v0, result_v0, weyl_zdown_v0);
		bgq_su3_vadd(result_v1, result_v1, weyl_zdown_v1);
		bgq_su3_vpiadd(result_v2, result_v2, weyl_zdown_v0);
		bgq_su3_vpisub(result_v3, result_v3, weyl_zdown_v1);
#endif


	}
#ifndef BGQ_HM_DIR_NOFUNC
}
#endif


#undef BGQ_HM_ZDOWN_COMPUTE
#undef BGQ_HM_ZDOWN_ACCUMULATE
