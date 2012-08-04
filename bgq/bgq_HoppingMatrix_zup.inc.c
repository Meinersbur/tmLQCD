
#ifndef BGQ_HM_ZUP_ZLINEINDENT
#error Must define the line indention (0, 1 or -1 for runtime-conditional)
#define BGQ_HM_ZUP_ZLINEINDENT -1
#endif

#ifndef BGQ_HM_ZUP_RIGHTWRAPAROUND
#define BGQ_HM_ZUP_RIGHTWRAPAROUND -1 /* Runtime-conditional */
#endif

#ifndef BGQ_HM_ZUP_COMPUTE
#define BGQ_HM_ZUP_COMPUTE 0
#endif

#ifndef BGQ_HM_ZUP_ACCUMULATE
#define BGQ_HM_ZUP_ACCUMULATE 0
#endif


#ifndef BGQ_HM_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void bgq_HoppingMatrix_zup(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int t, int x, int y, int z, int zv, int z1, int z2) {
	bgq_su3_spinor_decl(result);
#endif
	{

#if (BGQ_HM_ZUP_RIGHTWRAPAROUND==0)
		const int zv_right = zv+1;
#elif (BGQ_HM_ZUP_RIGHTWRAPAROUND==1)
		assert(zv==PHYSICAL_LZV-1);
		const int zv_right = 0;
#else
		const int zv_right = mod(zv + 1, PHYSICAL_LZV);
#endif


		// Load the input spinor
		bgq_su3_spinor_decl(spinor_zup);
		// # = stencil site to update (either even or odd sites, depending on what isOdd says)
		// _ = neighbor sites to read (therefore !isOdd)
		// (  2  ) = vector site with tv-number
		// top line = vector sites of isOdd sites
		// bottom line = vector sites of !isOdd sites
#if (BGQ_HM_ZUP_ZLINEINDENT==-1)
		if (((t + x + y)&1) == isOdd) {
#endif
#if (BGQ_HM_ZUP_ZLINEINDENT==-1) || (BGQ_HM_ZUP_ZLINEINDENT==0)
			// (  0  ) (  1  ) (  2  )
			// |# _ # _ # _ # _ # _ # _|
			//   (  0  ) (  1  ) (  2  )
			// T_UP = zv
			// T_DOWN = merge2(zv, zv-1)
			assert(((t+x+y)&1) == isOdd);

			bgq_spinorsite_double *spinorsite_zup = BGQ_SPINORSITE(spinorfield, !isOdd, t, x, y, zv, z1+1, z2+1);
			bgq_su3_spinor_double_load(spinor_zup, spinorsite_zup);
#endif
#if BGQ_HM_ZUP_ZLINEINDENT==-1
		} else {
#endif
#if (BGQ_HM_ZUP_ZLINEINDENT==-1) || (BGQ_HM_ZUP_ZLINEINDENT==1)
			//   (  0  ) (  1  ) (  2  )
			// |_ # _ # _ # _ # _ # _ #|
			// (  0  ) (  1  ) (  2  )
			// T_UP = half tv, half tv+1
			// T_DOWN = tv
			assert(((t+x+y)&1) == !isOdd);

			bgq_su3_spinor_decl_rightonly(spinor_zup_mid);
			bgq_spinorsite_double *spinorsite_zup_mid = BGQ_SPINORSITE(spinorfield, !isOdd, t, x, y, zv, z1-1,z1+1);
			bgq_su3_spinor_double_load_right_torightonly(spinor_zup_mid, spinorsite_zup_mid);

			bgq_su3_spinor_decl_leftonly(spinor_zup_right);
			bgq_spinorsite_double *spinorsite_zup_right = BGQ_SPINORSITE(spinorfield, !isOdd, t, x, y, zv_right, z2-1,z2+1);
			bgq_su3_spinor_double_load_left_toleftonly(spinor_zup_right, spinorsite_zup_right);

			bgq_su3_spinor_merge(spinor_zup, spinor_zup_mid, spinor_zup_right);
#endif
#if (BGQ_HM_ZUP_ZLINEINDENT==-1)
		}
#endif

		// Compute its halfspinor
		bgq_su3_weyl_decl(weyl_zup);
		bgq_su3_vadd(weyl_zup_v0, spinor_zup_v0, spinor_zup_v3);
		bgq_su3_vsub(weyl_zup_v1, spinor_zup_v1, spinor_zup_v2);


#if BGQ_HM_ZUP_COMPUTE
		bgq_gaugesite_double *gaugesite_zup = BGQ_GAUGESITE(gaugefield, isOdd, t, x, y, zv, Z_UP);
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
		bgq_su3_vsub(result_v2, result_v2, weyl_zup_v1);
		bgq_su3_vadd(result_v3, result_v3, weyl_zup_v0);
#endif

	}
#ifndef BGQ_HM_NOFUNC
}
#endif

#undef BGQ_HM_ZUP_ZLINEINDENT
#undef BGQ_HM_ZUP_COMPUTE
#undef BGQ_HM_ZUP_ACCUMULATE
