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
		bgq_su3_spinor_decl(spinor_ydown);
		bgq_spinorfield_double *spinor_ydown = bgq_spinorsite_double_physical_pointer(spinorfield, !isOdd, x, y - 1, z, tv);
		bgq_su3_spinor_double_load(spinor_ydown, spinor_ydown);

		bgq_su3_weyl_decl(weyl_ydown);
		bgq_su3_vpisub(weyl_ydown_v0, spinor_ydown_v0, spinor_ydown_v3);
		bgq_su3_vpisub(weyl_ydown_v1, spinor_ydown_v1, spinor_ydown_v2);

		bgq_su3_mdecl(gauge_ydown);
		bgq_gaugesite_double *gauge_ydown = bgq_gaugesiteeo_double_physical_pointer(gaugefield, !isOdd, x, y - 1, z, tv, Y_UP);
		bgq_su3_matrix_double_load(gauge_ydown, gauge_ydown);

		bgq_su3_weyl_decl(chi_ydown);
		bgq_su3_mvinvmul(chi_ydown_v0, gauge_ydown, weyl_ydown_v0);
		bgq_su3_mvinvmul(chi_ydown_v1, gauge_ydown, weyl_ydown_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(chi_ydown_v0, ka1, chi_ydown_v0);
		bgq_su3_cvmul(chi_ydown_v1, ka1, chi_ydown_v1);
#endif

		bgq_su3_vadd(result_v0, result_v0, chi_ydown_v0);
		bgq_su3_vadd(result_v1, result_v1, chi_ydown_v1);
		bgq_su3_vpiadd(result_v2, result_v2, chi_ydown_v1);
		bgq_su3_vpiadd(result_v3, result_v3, chi_ydown_v0);
	}

#ifndef BGQ_HM_NOFUNC
}
#endif

