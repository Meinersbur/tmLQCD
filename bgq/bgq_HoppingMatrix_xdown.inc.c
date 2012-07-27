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
		bgq_su3_spinor_decl(spinor_xdown);
		bgq_spinorfield_double *spinor_xdown = bgq_spinorsite_double_physical_pointer(spinorfield, !isOdd, x - 1, y, z, tv);
		bgq_su3_spinor_double_load(spinor_xdown, spinor_xdown);

		bgq_su3_weyl_decl(weyl_xdown);
		bgq_su3_vsub(weyl_xdown_v0, spinor_xdown_v0, spinor_xdown_v2);
		bgq_su3_vsub(weyl_xdown_v1, spinor_xdown_v1, spinor_xdown_v3);

		bgq_su3_mdecl(gauge_xdown);
		bgq_gaugesite_double *gauge_xdown = bgq_gaugesiteeo_double_physical_pointer(gaugefield, !isOdd, x - 1, y, z, tv, X_UP);
		bgq_su3_matrix_double_load(gauge_xdown, gauge_xdown);

		bgq_su3_weyl_decl(chi_xdown);
		bgq_su3_mvinvmul(chi_xdown_v0, gauge_xdown, weyl_xdown_v0);
		bgq_su3_mvinvmul(chi_xdown_v1, gauge_xdown, weyl_xdown_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(chi_xdown_v0, ka0, chi_xdown_v0);
		bgq_su3_cvmul(chi_xdown_v1, ka0, chi_xdown_v1);
#endif

		bgq_su3_vadd(result_v0, result_v0, chi_xdown_v0);
		bgq_su3_vadd(result_v1, result_v1, chi_xdown_v1);
		bgq_su3_vsub(result_v2, result_v2, chi_xdown_v0);
		bgq_su3_vsub(result_v3, result_v3, chi_xdown_v1);
	}

#ifndef BGQ_HM_NOFUNC
}
#endif

