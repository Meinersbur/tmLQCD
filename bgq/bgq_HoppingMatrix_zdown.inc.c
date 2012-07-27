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
		bgq_su3_spinor_decl(spinor_zdown);
		bgq_spinorfield_double *spinor_zdown = bgq_spinorsite_double_physical_pointer(spinorfield, !isOdd, x, y, z - 1, tv);
		bgq_su3_spinor_double_load(spinor_zdown, spinor_zdown);

		bgq_su3_weyl_decl(weyl_zdown);
		bgq_su3_vsub(weyl_zdown_v0, spinor_zdown_v0, spinor_zdown_v3);
		bgq_su3_vadd(weyl_zdown_v1, spinor_zdown_v1, spinor_zdown_v2);

		bgq_su3_mdecl(gauge_zdown);
		bgq_gaugesite_double *gauge_zdown = bgq_gaugesiteeo_double_physical_pointer(gaugefield, !isOdd, x, y, z - 1, tv, Z_UP);
		bgq_su3_matrix_double_load(gauge_zdown, gauge_zdown);

		bgq_su3_weyl_decl(chi_zdown);
		bgq_su3_mvinvmul(chi_zdown_v0, gauge_zdown, weyl_zdown_v0);
		bgq_su3_mvinvmul(chi_zdown_v1, gauge_zdown, weyl_zdown_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(chi_zdown_v0, ka1, chi_zdown_v0);
		bgq_su3_cvmul(chi_zdown_v1, ka1, chi_zdown_v1);
#endif

		bgq_su3_vadd(result_v0, result_v0, chi_zdown_v0);
		bgq_su3_vadd(result_v1, result_v1, chi_zdown_v1);
		bgq_su3_vadd(result_v2, result_v2, chi_zdown_v1);
		bgq_su3_vsub(result_v3, result_v3, chi_zdown_v0);
	}

#ifndef BGQ_HM_NOFUNC
}
#endif

