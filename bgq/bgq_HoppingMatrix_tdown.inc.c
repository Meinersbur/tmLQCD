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
		bgq_su3_spinor_decl(spinor_tdown);

#if BGQ_HM_T_ALTERNATING
#elif BGQ_HM_T_EVENLINE
#elif BGQ_HM_T_ODDLINE
#else
		if ((x + y + z) % 2 == isOdd) {
			bgq_spinorfield_double *spinor_tdown_left = bgq_spinorsite_double_physical_pointer(spinorfield, !isOdd, x, y, z, tv - 1);
			bgq_spinorfield_double *spinor_tdown_mid = bgq_spinorsite_double_physical_pointer(spinorfield, !isOdd, x, y, z, tv);

			bgq_su3_spinor_decl(spinor_tdown_left);
			bgq_su3_spinor_decl(spinor_tdown_mid);
			bgq_su3_spinor_double_load_right(spinor_tdown_left, spinor_tdown_left);
			bgq_su3_spinor_double_load_left(spinor_tdown_mid, spinor_tdown_mid);

			bgq_su3_spinor_merge(spinor_tdown, spinor_tdown_left, spinor_tdown_mid);
		} else {
			bgq_spinorfield_double *spinor_tdown = bgq_spinorsite_double_physical_pointer(spinorfield, !isOdd, x, y, z, tv);
			bgq_su3_spinor_double_load(spinor_tdown, spinor_tdown);
		}
#endif

		bgq_su3_weyl_decl(weyl_tdown);
		bgq_su3_vpisub(weyl_tdown_v0, spinor_tdown_v0, spinor_tdown_v2);
		bgq_su3_vpiadd(weyl_tdown_v1, spinor_tdown_v1, spinor_tdown_v3);

		bgq_su3_mdecl(gauge_tdown);
#if BGQ_HM_T_ALTERNATING
#elif BGQ_HM_T_EVENLINE
#elif BGQ_HM_T_ODDLINE
#else
		if ((x + y + z) % 2 == isOdd) {
			bgq_gaugesite_double *gauge_tdown_left = bgq_gaugesiteeo_double_physical_pointer(gaugefield, !isOdd, x, y, z, tv - 1, T_UP);
			bgq_gaugesite_double *gauge_tdown_mid = bgq_gaugesiteeo_double_physical_pointer(gaugefield, !isOdd, x, y, z, tv, T_UP);
			bgq_su3_mdecl(gauge_tdown_left);
			bgq_su3_mdecl(gauge_tdown_mid);

			bgq_su3_matrix_double_load_right(gauge_tdown_left, gauge_tdown_left);
			bgq_su3_matrix_double_load_left(gauge_tdown_mid, gauge_tdown_mid);

			bgq_su3_mmerge(gauge_tdown, gauge_tdown_left, gauge_tdown_mid);
		} else {
			bgq_gaugesite_double *gauge_tdown = bgq_gaugesiteeo_double_physical_pointer(gaugefield, !isOdd, x, y, z, tv, T_UP);
			bgq_su3_matrix_double_load(gauge_tdown, gauge_tdown);
		}
#endif
		bgq_su3_weyl_decl(chi_tdown);
		bgq_su3_mvinvmul(chi_tdown_v0, gauge_tdown, weyl_tdown_v0);
		bgq_su3_mvinvmul(chi_tdown_v1, gauge_tdown, weyl_tdown_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(chi_tdown_v0, ka1, chi_tdown_v0);
		bgq_su3_cvmul(chi_tdown_v1, ka1, chi_tdown_v1);
#endif

		bgq_su3_vadd(result_v0, result_v0, chi_tdown_v0);
		bgq_su3_vadd(result_v1, result_v1, chi_tdown_v1);
		bgq_su3_vpiadd(result_v2, result_v2, chi_tdown_v0);
		bgq_su3_vpisub(result_v3, result_v3, chi_tdown_v1);
	}


#ifndef BGQ_HM_NOFUNC
}
#endif

