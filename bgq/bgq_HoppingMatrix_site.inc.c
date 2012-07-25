#ifndef BGQ_HM_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void HoppingMatrix_site(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefieldeo_double gaugefield, bool isOdd, int x, int y, int z, int tv, int k) {
#endif

	bgq_su3_spinor_decl(result);

// direction X_UP /////////////////////////////////////////////////////////////

	{
		bgq_su3_spinor_decl(spinor_xup);
		bgq_spinorfield_double *spinor_xup = bgq_spinorsite_double_physical_pointer(spinorfield, !isOdd, x + 1, y, z, tv);
		bgq_su3_spinor_double_load(spinor_xup, spinor_xup);

		bgq_su3_weyl_decl(weyl_xup);
		bgq_su3_vadd(weyl_xup_v0, spinor_xup_v0, spinor_xup_v2);
		bgq_su3_vadd(weyl_xup_v1, spinor_xup_v1, spinor_xup_v3);

		bgq_su3_mdecl(gauge_xup);
		bgq_gaugesite_double *gauge_xup = bgq_gaugesiteeo_double_physical_pointer(gaugefield, isOdd, x, y, z, tv, X_UP);
		bgq_su3_matrix_double_load(gauge_xup, gauge_xup);

		bgq_su3_weyl_decl(chi_xup);
		bgq_su3_mvmul(chi_xup_v0, gauge_xup, weyl_xup_v0);
		bgq_su3_mvmul(chi_xup_v1, gauge_xup, weyl_xup_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(chi_xup_v0, ka0, chi_xup_v0);
		bgq_su3_cvmul(chi_xup_v1, ka0, chi_xup_v1);
#endif

		bgq_su3_vmov(result_v0, chi_xup_v0);
		bgq_su3_vmov(result_v1, chi_xup_v1);
		bgq_su3_vmov(result_v2, chi_xup_v0);
		bgq_su3_vmov(result_v3, chi_xup_v1);
	}

// direction X_DOWN /////////////////////////////////////////////////////////////

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

// direction Y_UP /////////////////////////////////////////////////////////////

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

// direction Y_DOWN /////////////////////////////////////////////////////////////

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

// direction Z_UP /////////////////////////////////////////////////////////////

	{
		bgq_su3_spinor_decl(spinor_zup);
		bgq_spinorfield_double *spinor_zup = bgq_spinorsite_double_physical_pointer(spinorfield, !isOdd, x, y, z + 1, tv);
		bgq_su3_spinor_double_load(spinor_zup, spinor_zup);

		bgq_su3_weyl_decl(weyl_zup);
		bgq_su3_vadd(weyl_zup_v0, spinor_zup_v0, spinor_zup_v3);
		bgq_su3_vsub(weyl_zup_v1, spinor_zup_v1, spinor_zup_v2);

		bgq_su3_mdecl(gauge_zup);
		bgq_gaugesite_double *gauge_zup = bgq_gaugesiteeo_double_physical_pointer(gaugefield, isOdd, x, y, z, tv, Z_UP);
		bgq_su3_matrix_double_load(gauge_zup, gauge_zup);

		bgq_su3_weyl_decl(chi_zup);
		bgq_su3_mvmul(chi_zup_v0, gauge_yup, weyl_zup_v0);
		bgq_su3_mvmul(chi_zup_v1, gauge_yup, weyl_zup_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(chi_zup_v0, ka2, chi_zup_v0);
		bgq_su3_cvmul(chi_zup_v1, ka2, chi_zup_v1);
#endif

		bgq_su3_vadd(result_v0, result_v0, chi_zup_v0);
		bgq_su3_vadd(result_v1, result_v1, chi_zup_v1);
		bgq_su3_vsub(result_v2, result_v2, chi_zup_v1);
		bgq_su3_vadd(result_v3, result_v3, chi_zup_v0);
	}

	// direction Z_DOWN /////////////////////////////////////////////////////////////

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

	// direction T_UP /////////////////////////////////////////////////////////////

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

	// direction T_DOWN /////////////////////////////////////////////////////////////

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

///////////////////////////////////////////////////////////////////////////////
// Store the spinor

	bgq_spinorsite_double *spinor_target = bgq_spinorsite_double_physical_pointer(spinorfield, isOdd, x, y, z, tv);
	bgq_su3_spinor_double_store(spinor_target, result);

#ifndef BGQ_HM_NOFUNC
}
#endif

