/*
 * bgq_operators.inc.c
 *
 *  Created on: Sep 19, 2012
 *      Author: meinersbur
 */

#define BGQ_OPERATORS_INC_C_

#ifndef PRECISION
//#define BGQ_PRECISION 64
//#include "bgq_field.inc.h"
#include "bgq_field_double.h"
#endif

#include "bgq_field.inc.h"
#include "bgq.h"
#include <omp.h>
#include "bgq_field.h"
#include "../global.h"


void bgq_mul_one_pm_imu_sub_mul_gamma5(bgq_spinorfield l, bgq_spinorfield k, bgq_spinorfield j, bool isOdd, double sign) {
	assert(!omp_in_parallel());

	sign = (sign >= 0) ? 1 : -1;

	bgq_vector4double_decl(z);
	bgq_cconst(z, 1, sign * g_mu);

	bgq_vector4double_decl(w);
	bgq_cconst(w, 1, -sign * g_mu);

#pragma omp parallel
	{
#pragma omp for schedule(static) nowait
		for (int txy = 0; txy < VOLUME_ZLINES; txy += 1) {
			WORKLOAD_DECL(txy, VOLUME_SITES);
			const int tv = WORKLOAD_PARAM(PHYSICAL_LTV);
			const int x = WORKLOAD_PARAM(PHYSICAL_LX);
			const int y = WORKLOAD_PARAM(PHYSICAL_LY);
			WORKLOAD_CHECK

			for (int z = 0; z < PHYSICAL_LZ; z += 1) {
				const int t1 = ((isOdd + x + y + z) & 1) + tv * PHYSICAL_LP * PHYSICAL_LK;
				const int t2 = t1 + 2;

				bgq_spinorsite *site_r = BGQ_SPINORSITE(k, isOdd, tv, x, y, z, t1, t2, true, false);
				bgq_su3_spinor_decl(r);
				bgq_su3_spinor_load(r, site_r);

				bgq_spinorsite *site_s = BGQ_SPINORSITE(j, isOdd, tv, x, y, z, t1, t2, true, false);
				bgq_su3_spinor_decl(s);
				bgq_su3_spinor_load(s, site_s);

				bgq_su3_spinor_decl(phi);
				bgq_su3_cvmul(phi_v0, z, r_v0);
				bgq_su3_cvmul(phi_v1, z, r_v1);
				bgq_su3_cvmul(phi_v2, w, r_v2);
				bgq_su3_cvmul(phi_v3, w, r_v3);

				bgq_su3_spinor_decl(result);
				bgq_su3_vsub(result_v0, phi_v0, s_v0);
				bgq_su3_vsub(result_v1, phi_v1, s_v1);
				bgq_su3_vsub(result_v2, s_v2, phi_v2);
				bgq_su3_vsub(result_v3, s_v3, phi_v3);

				bgq_spinorsite *site_t = BGQ_SPINORSITE(l, isOdd, tv, x, y, z, t1, t2, true, false);
				bgq_su3_spinor_store(site_t, result);
			}
		}
	}
}


void bgq_mul_one_pm_imu_inv(bgq_spinorfield spinorfield, bool isOdd, double sign) {
	assert(!omp_in_parallel());

	double nrm = 1.0 / (1.0 + g_mu * g_mu);
	sign = (sign >= 0) ? 1 : -1;

	bgq_vector4double_decl(z);
	bgq_cconst(z, nrm, sign * nrm * g_mu);

	bgq_vector4double_decl(w);
	bgq_cconst(w, nrm, -sign * nrm * g_mu);

#pragma omp parallel
	{
#pragma omp for schedule(static) nowait
		for (int txy = 0; txy < VOLUME_ZLINES; txy += 1) {
			WORKLOAD_DECL(txy, VOLUME_SITES);
			const int tv = WORKLOAD_PARAM(PHYSICAL_LTV);
			const int x = WORKLOAD_PARAM(PHYSICAL_LX);
			const int y = WORKLOAD_PARAM(PHYSICAL_LY);
			WORKLOAD_CHECK

			for (int z = 0; z < PHYSICAL_LZ; z += 1) {
				const int t1 = ((isOdd + x + y + z) & 1) + tv * PHYSICAL_LP * PHYSICAL_LK;
				const int t2 = t1 + 2;

				bgq_spinorsite *spinorsite = BGQ_SPINORSITE(spinorfield, isOdd, tv, x, y, z, t1, t2, true, false);
				bgq_su3_spinor_decl(spinor);
				bgq_su3_spinor_load(spinor, spinorsite);

				bgq_su3_spinor_decl(result);
				bgq_su3_cvmul(result_v0, z, spinor_v0);
				bgq_su3_cvmul(result_v1, z, spinor_v1);
				bgq_su3_cvmul(result_v2, w, spinor_v2);
				bgq_su3_cvmul(result_v3, w, spinor_v3);

				bgq_su3_spinor_store(spinorsite, result);
			}
		}
	}
}


double bgq_scalar_prod(bgq_spinorfield spinorfield1, bgq_spinorfield spinorfield2, bool isOdd, bool parallel) {
	assert(!omp_in_parallel());

	bgq_vector4double_decl(shared_sum);
	bgq_zero(shared_sum);

#pragma omp parallel
	{
		bgq_vector4double_decl(thread_sum);
		bgq_zero(thread_sum);

#pragma omp for schedule(static) nowait
		for (int txy = 0; txy < VOLUME_ZLINES; txy+=1) {
			WORKLOAD_DECL(txy, VOLUME_SITES);
			const int tv = WORKLOAD_PARAM(PHYSICAL_LTV);
			const int x = WORKLOAD_PARAM(PHYSICAL_LX);
			const int y = WORKLOAD_PARAM(PHYSICAL_LY);
			WORKLOAD_CHECK

			for (int z = 0; z < PHYSICAL_LZ; z += 1) {
				const int t1 = ((isOdd + x + y + z) & 1) + tv * PHYSICAL_LP * PHYSICAL_LK;
				const int t2 = t1 + 2;

				bgq_spinorsite *spinorsite1 = BGQ_SPINORSITE(spinorfield1, isOdd, tv, x, y, z, t1, t2, true, false);
				bgq_su3_spinor_decl(spinor1);
				bgq_su3_spinor_load(spinor1, spinorsite1);

				bgq_spinorsite *spinorsite2 = BGQ_SPINORSITE(spinorfield2, isOdd, tv, x, y, z, t1, t2, true, false);
				bgq_su3_spinor_decl(spinor2);
				bgq_su3_spinor_load(spinor2, spinorsite2);

				// re = a.re * b.re + a.im * b.im + re
				// im = a.re * b.im - a.im * b.re + im
				bgq_xmadd(thread_sum, spinor1_v0_c0, spinor2_v0_c0, thread_sum);
				bgq_rxxcpnmadd(thread_sum, spinor1_v0_c0, spinor2_v0_c0, thread_sum);

				bgq_xmadd(thread_sum, spinor1_v0_c1, spinor2_v0_c1, thread_sum);
				bgq_rxxcpnmadd(thread_sum, spinor1_v0_c1, spinor2_v0_c1, thread_sum);

				bgq_xmadd(thread_sum, spinor1_v0_c2, spinor2_v0_c2, thread_sum);
				bgq_rxxcpnmadd(thread_sum, spinor1_v0_c2, spinor2_v0_c2, thread_sum);

				bgq_xmadd(thread_sum, spinor1_v1_c0, spinor2_v1_c0, thread_sum);
				bgq_rxxcpnmadd(thread_sum, spinor1_v1_c0, spinor2_v1_c0, thread_sum);

				bgq_xmadd(thread_sum, spinor1_v1_c1, spinor2_v1_c1, thread_sum);
				bgq_rxxcpnmadd(thread_sum, spinor1_v1_c1, spinor2_v1_c1, thread_sum);

				bgq_xmadd(thread_sum, spinor1_v1_c2, spinor2_v1_c2, thread_sum);
				bgq_rxxcpnmadd(thread_sum, spinor1_v1_c2, spinor2_v1_c2, thread_sum);

				bgq_xmadd(thread_sum, spinor1_v2_c0, spinor2_v2_c0, thread_sum);
				bgq_rxxcpnmadd(thread_sum, spinor1_v2_c0, spinor2_v2_c0, thread_sum);

				bgq_xmadd(thread_sum, spinor1_v2_c1, spinor2_v2_c1, thread_sum);
				bgq_rxxcpnmadd(thread_sum, spinor1_v2_c1, spinor2_v2_c1, thread_sum);

				bgq_xmadd(thread_sum, spinor1_v2_c2, spinor2_v2_c2, thread_sum);
				bgq_rxxcpnmadd(thread_sum, spinor1_v2_c2, spinor2_v2_c2, thread_sum);

				bgq_xmadd(thread_sum, spinor1_v3_c0, spinor2_v3_c0, thread_sum);
				bgq_rxxcpnmadd(thread_sum, spinor1_v3_c0, spinor2_v3_c0, thread_sum);

				bgq_xmadd(thread_sum, spinor1_v3_c1, spinor2_v3_c1, thread_sum);
				bgq_rxxcpnmadd(thread_sum, spinor1_v3_c1, spinor2_v3_c1, thread_sum);

				bgq_xmadd(thread_sum, spinor1_v3_c2, spinor2_v3_c2, thread_sum);
				bgq_rxxcpnmadd(thread_sum, spinor1_v3_c2, spinor2_v3_c2, thread_sum);
			}
		}

#ifdef XLC
#pragma tm_atomic
#else
#pragma omp critical
#endif
		{
			// not using Kahan summation at the moment
			bgq_add(shared_sum, thread_sum, shared_sum);
		}
	}

	complexdouble result = 0;
	complexdouble node_result = bgq_cmplxval1(shared_sum) + bgq_cmplxval2(shared_sum);
	if (parallel) {
		MPI_Allreduce(&node_result, &result, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD );
	} else {
		result = node_result;
	}

	return result;
}


complexdouble bgq_square_norm(bgq_spinorfield spinorfield, bool isOdd, bool parallel) {
	assert(!omp_in_parallel());

	bgq_vector4double_decl(shared_sum);
	bgq_zero(shared_sum);

#pragma omp parallel
	{
		bgq_vector4double_decl(thread_sum);
		bgq_zero(thread_sum);

#pragma omp for schedule(static) nowait
		for (int txy = 0; txy < VOLUME_ZLINES; txy+=1) {
			WORKLOAD_DECL(txy, VOLUME_SITES);
			const int tv = WORKLOAD_PARAM(PHYSICAL_LTV);
			const int x = WORKLOAD_PARAM(PHYSICAL_LX);
			const int y = WORKLOAD_PARAM(PHYSICAL_LY);
			WORKLOAD_CHECK

			for (int z = 0; z < PHYSICAL_LZ; z += 1) {
				const int t1 = ((isOdd + x + y + z) & 1) + tv * PHYSICAL_LP * PHYSICAL_LK;
				const int t2 = t1 + 2;

				bgq_spinorsite *spinorsite = BGQ_SPINORSITE(spinorfield, isOdd, tv, x, y, z, t1, t2, true, false);
				bgq_su3_spinor_decl(spinor);
				bgq_su3_spinor_load(spinor, spinorsite);

				//TODO: prefetch
				bgq_madd(thread_sum, spinor_v0_c0, spinor_v0_c0, thread_sum);
				bgq_madd(thread_sum, spinor_v0_c1, spinor_v0_c1, thread_sum);
				bgq_madd(thread_sum, spinor_v0_c2, spinor_v0_c2, thread_sum);
				bgq_madd(thread_sum, spinor_v1_c0, spinor_v1_c0, thread_sum);
				bgq_madd(thread_sum, spinor_v1_c1, spinor_v1_c1, thread_sum);
				bgq_madd(thread_sum, spinor_v1_c2, spinor_v1_c2, thread_sum);
				bgq_madd(thread_sum, spinor_v2_c0, spinor_v2_c0, thread_sum);
				bgq_madd(thread_sum, spinor_v2_c1, spinor_v2_c1, thread_sum);
				bgq_madd(thread_sum, spinor_v2_c2, spinor_v2_c2, thread_sum);
				bgq_madd(thread_sum, spinor_v3_c0, spinor_v3_c0, thread_sum);
				bgq_madd(thread_sum, spinor_v3_c1, spinor_v3_c1, thread_sum);
				bgq_madd(thread_sum, spinor_v3_c2, spinor_v3_c2, thread_sum);
			}
		}

#ifdef XLC
#pragma tm_atomic
#else
#pragma omp critical
#endif
		{
			// not using Kahan summation at the moment
			bgq_add(shared_sum, thread_sum, shared_sum);
		}
	}

	double result = 0;
	double node_result = bgq_elem0(shared_sum) + bgq_elem1(shared_sum) + bgq_elem2(shared_sum) + bgq_elem3(shared_sum);
	if (parallel) {
		MPI_Allreduce(&node_result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	} else {
		result = node_result;
	}

	return result;
}

