#ifndef BGQ_HM_ZDOWN_ZLINEINDENT
#error Must define the line indention (0, 1 or -1 for runtime-conditional)
#define BGQ_HM_ZDOWN_ZLINEINDENT -1
#endif

#ifndef BGQ_HM_ZDOWN_LEFTWRAPAROUND
#define BGQ_HM_ZDOWN_LEFTWRAPAROUND -1 /* Runtime-conditional */
#endif

#ifndef BGQ_HM_ZDOWN_COMPUTE
#define BGQ_HM_ZDOWN_COMPUTE 0
#endif

#ifndef BGQ_HM_ZDOWN_ACCUMULATE
#define BGQ_HM_ZDOWN_ACCUMULATE 0
#endif

#ifndef BGQ_HM_ZDOWN_READCARRY
#define BGQ_HM_ZDOWN_READCARRY 0
#endif

#ifndef BGQ_HM_ZDOWN_WRITECARRY
#define BGQ_HM_ZDOWN_WRITECARRY 0
#endif

#ifndef BGQ_HM_DIR_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void bgq_HoppingMatrix_zdown(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int t, int x, int y, int zv) {
	bgq_su3_spinor_decl(result);
#endif
	{

#if (BGQ_HM_ZDOWN_LEFTWRAPAROUND==0)
		const int zv_left = zv-1;
#elif (BGQ_HM_ZDOWN_LEFTWRAPAROUND==1)
		assert(zv==0);
		const int zv_left = PHYSICAL_LZV-1;
#elif (BGQ_HM_ZDOWN_LEFTWRAPAROUND==-1)
		const int zv_left = (PHYSICAL_LZV + zv - 1) % PHYSICAL_LZV;
#endif

		bgq_su3_weyl_decl(weyl_zdown);
#if BGQ_HM_ZDOWN_READCARRY
		bgq_su3_weyl_mov(weyl_zdown, weyl_zcarry);
#else
		// Load the input spinor
		bgq_su3_spinor_decl(spinor_zdown);
#if (BGQ_HM_ZDOWN_ZLINEINDENT==-1)
		if ((t + x + y) % 2 == isOdd) {
#endif
#if (BGQ_HM_ZDOWN_ZLINEINDENT==-1) || (BGQ_HM_ZDOWN_ZLINEINDENT==0)
			assert((t+x+y)%2 == isOdd);

			bgq_su3_spinor_decl_rightonly(spinor_zdown_left);
			bgq_spinorsite_double *spinorsite_zdown_left = BGQ_SPINORSITE(spinorfield, !isOdd, t, x, y, zv_left, z2-1, z2+1);
			bgq_su3_spinor_double_load_right_torightonly(spinor_zdown_left, spinorsite_zdown_left);

			bgq_su3_spinor_decl_leftonly(spinor_zdown_mid);
			bgq_spinorsite_double *spinorsite_zdown_mid = BGQ_SPINORSITE(spinorfield, !isOdd, t, x, y, zv, z1-3, z1-1);
			bgq_su3_spinor_double_load_left_toleftonly(spinor_zdown_mid, spinorsite_zdown_mid);

			bgq_su3_spinor_merge(spinor_zdown, spinor_zdown_left, spinor_zdown_mid);
#endif
#if BGQ_HM_ZDOWN_ZLINEINDENT==-1
		} else {
#endif
#if (BGQ_HM_ZDOWN_ZLINEINDENT==-1) || (BGQ_HM_ZDOWN_ZLINEINDENT==1)
			assert((t+x+y)%2 == !isOdd);

			bgq_spinorsite_double *spinorsite_zdown = BGQ_SPINORSITE(spinorfield, !isOdd, t, x, y, zv, z1-1,z2-1);
			bgq_su3_spinor_double_load(spinor_zdown, spinorsite_zdown);
#endif
#if (BGQ_HM_ZDOWN_ZLINEINDENT==-1)
		}
#endif

#if BGQ_HM_ZDOWN_WRITECARRY
		// Compute the halfspinor for the zup of the next (downwards-)iteration
		//TODO: Confirm correct computation
		bgq_su3_vadd(weyl_zcarry_v0, spinor_zdown_v0, spinor_zdown_v3);
		bgq_su3_vsub(weyl_zcarry_v1, spinor_zdown_v1, spinor_zdown_v2);
#endif

#if BGQ_HM_ZDOWN_ACCUMULATE /* unused otherwise */
		// Compute its halfspinor
		bgq_su3_vsub(weyl_zdown_v0, spinor_zdown_v0, spinor_zdown_v3);
		bgq_su3_vadd(weyl_zdown_v1, spinor_zdown_v1, spinor_zdown_v2);
#endif
#endif

#if BGQ_HM_ZDOWN_COMPUTE
		bgq_su3_mdecl(gauge_zdown);
		bgq_gaugesite_double *gaugesite_zdown = BGQ_GAUGESITE(gaugefield, !isOdd, t, x, y, zv, Z_UP_SHIFT);
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
		bgq_su3_vadd(result_v2, result_v2, weyl_zdown_v1);
		bgq_su3_vsub(result_v3, result_v3, weyl_zdown_v0);
#endif

	}
#ifndef BGQ_HM_DIR_NOFUNC
}
#endif

#undef BGQ_HM_ZDOWN_ZLINEINDENT
#undef BGQ_HM_ZDOWN_COMPUTE
#undef BGQ_HM_ZDOWN_ACCUMULATE
#undef BGQ_HM_ZDOWN_WRITECARRY
#undef BGQ_HM_ZDOWN_READCARRY
