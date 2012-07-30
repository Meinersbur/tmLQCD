
// BGQ_HM_TLINE_FLUSHLINE
// BGQ_HM_TLINE_RAGGEDLINE

#if !defined(BGQ_HM_TLINE_FLUSHLINE)
#error Need to define flush- or ragged line
#endif

#ifndef BGQ_HM_TDOWN_LEFTWRAPAROUND
#define BGQ_HM_TDOWN_LEFTWRAPAROUND 0
#endif

#ifndef BGQ_HM_TDOWN_COMPUTE
#define BGQ_HM_TDOWN_COMPUTE 0
#endif

#ifndef BGQ_HM_TDOWN_ACCUMULATE
#define BGQ_HM_TDOWN_ACCUMULATE 0
#endif

#ifndef BGQ_HM_TDOWN_READCARRYSPINOR
#define BGQ_HM_TDOWN_READCARRYSPINOR 0
#endif

#ifndef BGQ_HM_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void HoppingMatrix_site(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefieldeo_double gaugefield, bool isOdd, int x, int y, int z, int tv, int k) {
	bgq_su3_spinor_decl(result);
#endif

	{
#if BGQ_HM_TLINE_FLUSHLINE
#if BGQ_HM_TDOWN_LEFTWRAPAROUND
		const int tv_left = PHYSICAL_LTV-1;
#else
		const int tv_left = tv-1;
#endif
#endif

		// Load the input spinor
		bgq_su3_spinor_decl(spinor_tdown);
#if BGQ_HM_TDOWN_READCARRYSPINOR
		bgq_su3_spinor_mov(spinor_tdown, spinor_tcarry);
#else
#if BGQ_HM_TLINE_FLUSHLINE
		// (  0  ) (  1  ) (  2  )
		// |# _ # _ # _ # _ # _ # _|
		//   (  0  ) (  1  ) (  2  )
		// T_UP = tv
		// T_DOWN = half tv, half tv-1
		assert((x+y+z)%2 == isOdd);

		bgq_spinorsite_double *spinorsite_tdown_left = BGQ_SPINORSITE(spinorfield, !isOdd, x, y, z, tv_left);
		bgq_su3_spinor_decl_rightonly(spinor_tdown_left);
		bgq_su3_spinor_double_load_right_torightonly(spinor_tdown_left, spinorsite_tdown_left);

		bgq_spinorsite_double *spinorsite_tdown_mid = BGQ_SPINORSITE(spinorfield, !isOdd, x, y, z, tv);
		bgq_su3_spinor_decl_leftonly(spinor_tdown_mid);
		bgq_su3_spinor_double_load_left_toleftonly(spinor_tdown_mid, spinorsite_tdown_mid);

		bgq_su3_spinor_merge(spinor_tdown, spinor_tdown_left/*toright*/, spinor_tdown_mid/*toleft*/);
#else
		//   (  0  ) (  1  ) (  2  )
		// |_ # _ # _ # _ # _ # _ #|
		// (  0  ) (  1  ) (  2  )
		// T_UP = half tv, half tv+1
		// T_DOWN = tv
		assert((x+y+z)%2 == !isOdd);

		bgq_spinorsite_double *spinorsite_tdown = BGQ_SPINORSITE(spinorfield, !isOdd, x, y, z, tv);
		bgq_su3_spinor_double_load(spinor_tdown, spinorsite_tdown);
#endif
#endif

		bgq_su3_weyl_decl(weyl_tdown);
		bgq_su3_vpisub(weyl_tdown_v0, spinor_tdown_v0, spinor_tdown_v2);
		bgq_su3_vpiadd(weyl_tdown_v1, spinor_tdown_v1, spinor_tdown_v3);

#if BGQ_HM_TDOWN_COMPUTE
		bgq_su3_mdecl(gauge_tdown);
#if BGQ_HM_TLINE_FLUSHLINE
		bgq_gaugesite_double *gaugesite_tdown = BGQ_GAUGESITE(gaugefield, isOdd, x, y, z, tv_left, T_RAGGED_UP);
#else
		bgq_gaugesite_double *gaugesite_tdown = BGQ_GAUGESITE(gaugefield, isOdd, x, y, z, tv, T_UP);
#endif
		bgq_su3_matrix_double_load(gauge_tdown, gaugesite_tdown);

		bgq_su3_mvinvmul(weyl_tdown_v0, gauge_tdown, weyl_tdown_v0);
		bgq_su3_mvinvmul(weyl_tdown_v1, gauge_tdown, weyl_tdown_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(weyl_tdown_v0, qka3, weyl_tdown_v0);
		bgq_su3_cvmul(weyl_tdown_v1, qka3, weyl_tdown_v1);
#endif
#endif

#if BGQ_HM_TDOWN_ACCUMULATE
		bgq_su3_vadd(result_v0, result_v0, weyl_tdown_v0);
		bgq_su3_vadd(result_v1, result_v1, weyl_tdown_v1);
		bgq_su3_vpiadd(result_v2, result_v2, weyl_tdown_v0);
		bgq_su3_vpisub(result_v3, result_v3, weyl_tdown_v1);
#endif
	}


#ifndef BGQ_HM_NOFUNC
}
#endif

#undef BGQ_HM_TDOWN_LEFTWRAPAROUND
#undef BGQ_HM_TDOWN_COMPUTE
#undef BGQ_HM_TDOWN_ACCUMULATE
#undef BGQ_HM_TDOWN_READCARRYSPINOR
