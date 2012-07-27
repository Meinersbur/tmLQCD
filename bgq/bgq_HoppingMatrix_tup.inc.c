/*
 * bgq_HoppingMatrix_xup.inc.c
 *
 *  Created on: Jul 27, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_HM_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void HoppingMatrix_site(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefieldeo_double gaugefield, bool isOdd, int x, int y, int z, int tv, int k) {
	bgq_su3_spinor_decl(result);
	#endif


	{

		bgq_su3_spinor_decl(spinor_tup);
#if BGQ_HM_T_ALTERNATING
#elif BGQ_HM_T_EVENLINE
#elif BGQ_HM_T_ODDLINE
#else
		// runtime condition

		// # = stencil site to update (either even or odd sites, depending on what isOdd says)
		// _ = neighbor sites to read (therefore !isOdd)
		// (  2  ) = vector site with tv-number
		// top line = vector sites of isOdd sites
		// bottom line = vector sites of !isOdd sites
		if ((x + y + z) % 2 == isOdd) {
			// (  0  ) (  1  ) (  2  )
			// |# _ # _ # _ # _ # _ # _|
			//   (  0  ) (  1  ) (  2  )
			// T_UP = tv
			// T_DOWN = half tv, half tv-1

			bgq_spinorfield_double *spinor_tup = bgq_spinorsite_double_physical_pointer(spinorfield, !isOdd, x, y, z, tv);
			bgq_su3_spinor_double_load(spinor_tup, spinor_tup);
		} else {
			//   (  0  ) (  1  ) (  2  )
			// |_ # _ # _ # _ # _ # _ #|
			// (  0  ) (  1  ) (  2  )
			// T_UP = half tv, half tv+1
			// T_DOWN = tv

			bgq_spinorfield_double *spinor_tup_mid = bgq_spinorsite_double_physical_pointer(spinorfield, !isOdd, x, y, z, tv);
			bgq_spinorfield_double *spinor_tup_right = bgq_spinorsite_double_physical_pointer(spinorfield, !isOdd, x, y, z, tv + 1);

			bgq_su3_spinor_decl(spinor_tup_mid);
			bgq_su3_spinor_decl(spinor_tup_right);
			bgq_su3_spinor_double_load_left(spinor_tup_right, spinor_tup_right);
			bgq_su3_spinor_double_load_right(spinor_tup_mid, spinor_tup_mid);
			bgq_su3_spinor_merge(spinor_tup, spinor_tup_mid, spinor_tup_right);
		}
#endif

		bgq_su3_weyl_decl(weyl_tup);
		bgq_su3_vpiadd(weyl_tup_v0, spinor_tup_v0, spinor_tup_v2);
		bgq_su3_vipsub(weyl_tup_v1, spinor_tup_v1, spinor_tup_v3);

		bgq_su3_mdecl(gauge_tup);

		bgq_gaugesite_double *gauge_tup = bgq_gaugesiteeo_double_physical_pointer(gaugefield, isOdd, x, y, z, tv, T_UP);
		bgq_su3_matrix_double_load(gauge_tup, gauge_tup);

		bgq_su3_weyl_decl(chi_tup);
		bgq_su3_mvmul(chi_tup_v0, gauge_tup, weyl_tup_v0);
		bgq_su3_mvmul(chi_tup_v1, gauge_tup, weyl_tup_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(chi_tup_v0, ka2, chi_tup_v0);
		bgq_su3_cvmul(chi_tup_v1, ka2, chi_tup_v1);
#endif

		bgq_su3_vadd(result_v0, result_v0, chi_tup_v0);
		bgq_su3_vadd(result_v1, result_v1, chi_tup_v1);
		bgq_su3_vpisub(result_v2, result_v2, chi_tup_v0);
		bgq_su3_vpiadd(result_v3, result_v3, chi_tup_v1);
	}

#ifndef BGQ_HM_NOFUNC
}
#endif

