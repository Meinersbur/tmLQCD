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
		// Load the input spinor
		bgq_su3_spinor_decl(spinor_xup);
		bgq_spinorfield_double *spinor_xup = bgq_spinorsite_double_physical_pointer(spinorfield, !isOdd, x + 1, y, z, tv);
		bgq_su3_spinor_double_load(spinor_xup, spinor_xup);

		// Compute its halfspinor
		bgq_su3_weyl_decl(weyl_xup);
		bgq_su3_vadd(weyl_xup_v0, spinor_xup_v0, spinor_xup_v2);
		bgq_su3_vadd(weyl_xup_v1, spinor_xup_v1, spinor_xup_v3);

#if BGQ_HM_BORDER
		// Store the halfspinor to be transfered to the neighbor node
		bgq_su3_weyldouble_store(weylxchange_xup_send_double, weyl_xup);
#else
		// Load the interaction matrix between the lattice sites
		bgq_su3_mdecl(gauge_xup);
		bgq_gaugesite_double *gauge_xup = bgq_gaugesiteeo_double_physical_pointer(gaugefield, isOdd, x, y, z, tv, X_UP);
		bgq_su3_matrix_double_load(gauge_xup, gauge_xup);

		// Multiply the halfspinor with the matrix
		bgq_su3_weyl_decl(chi_xup);
		bgq_su3_mvmul(chi_xup_v0, gauge_xup, weyl_xup_v0);
		bgq_su3_mvmul(chi_xup_v1, gauge_xup, weyl_xup_v1);

#ifndef BGQ_HM_NOKAMUL
		// Multiply with custom constant
		bgq_su3_cvmul(chi_xup_v0, ka0, chi_xup_v0);
		bgq_su3_cvmul(chi_xup_v1, ka0, chi_xup_v1);
#endif

		// Add up at the output lattice site
		bgq_su3_vmov(result_v0, chi_xup_v0);
		bgq_su3_vmov(result_v1, chi_xup_v1);
		bgq_su3_vmov(result_v2, chi_xup_v0);
		bgq_su3_vmov(result_v3, chi_xup_v1);
	}
#endif

#ifndef BGQ_HM_NOFUNC
}
#endif

