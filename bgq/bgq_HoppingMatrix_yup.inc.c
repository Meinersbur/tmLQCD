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
		bgq_su3_spinor_decl(spinor_yup);
		bgq_spinorfield_double *spinor_yup = bgq_spinorsite_double_physical_pointer(spinorfield, !isOdd, x, y + 1, z, tv);
		bgq_su3_spinor_double_load(spinor_yup, spinor_yup);

		bgq_su3_weyl_decl(weyl_yup);
		bgq_su3_vpiadd(weyl_yup_v0, spinor_yup_v0, spinor_yup_v3);
		bgq_su3_vpiadd(weyl_yup_v1, spinor_yup_v1, spinor_yup_v2);

		bgq_su3_mdecl(gauge_yup);
		bgq_gaugesite_double *gauge_yup = bgq_gaugesiteeo_double_physical_pointer(gaugefield, isOdd, x, y, z, tv, Y_UP);
		bgq_su3_matrix_double_load(gauge_yup, gauge_yup);

		bgq_su3_weyl_decl(chi_yup);
		bgq_su3_mvmul(chi_yup_v0, gauge_yup, weyl_yup_v0);
		bgq_su3_mvmul(chi_yup_v1, gauge_yup, weyl_yup_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(chi_yup_v0, ka1, chi_yup_v0);
		bgq_su3_cvmul(chi_yup_v1, ka1, chi_yup_v1);
#endif

		bgq_su3_vadd(result_v0, result_v0, chi_yup_v0);
		bgq_su3_vadd(result_v1, result_v1, chi_yup_v1);
		bgq_su3_vpisub(result_v2, result_v2, chi_yup_v1);
		bgq_su3_vpisub(result_v3, result_v3, chi_yup_v0);
	}

#ifndef BGQ_HM_NOFUNC
}
#endif

