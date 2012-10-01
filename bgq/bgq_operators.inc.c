/*
 * bgq_operators.inc.c
 *
 *  Created on: Sep 19, 2012
 *      Author: meinersbur
 */

#define BGQ_OPERATORS_INC_C_

#ifndef BGQ_PRECISION
#define BGQ_PRECISION 64
#include "bgq_precisionselect.inc.c"
#endif

#include "bgq_operators.inc.h"
//#include "bgq_field.inc.h"
#include "bgq.h"
#include <omp.h>
#include "bgq_field.h"
#include "../global.h"
#include <omp.h>


void bgq_spinorfield_mul_weyl_complex(bgq_spinorfield targetfield, bgq_spinorfield spinorfield, bool isOdd, complexdouble z, complexdouble w) {
	assert(!omp_in_parallel());
	bgq_spinorfield_setOdd(targetfield, isOdd, true);
	bgq_spinorfield_setOdd(spinorfield, isOdd, false);

	bgq_vector4double_decl(qz);
	bgq_cconst(qz, creal(z), cimag(z));

	bgq_vector4double_decl(qw);
	bgq_cconst(qw, creal(w), cimag(w));

#pragma omp parallel for schedule(static) firstprivate(bgq_vars(qz), bgq_vars(qw))
	for (int txy = 0; txy < VOLUME_ZLINES; txy += 1) {
		WORKLOAD_DECL(txy, VOLUME_ZLINES);
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
			bgq_su3_cvmul(result_v0, qz, spinor_v0);
			bgq_su3_cvmul(result_v1, qz, spinor_v1);
			bgq_su3_cvmul(result_v2, qw, spinor_v2);
			bgq_su3_cvmul(result_v3, qw, spinor_v3);

			bgq_spinorsite *targetsite = BGQ_SPINORSITE(targetfield, isOdd, tv, x, y, z, t1, t2, true, false);
			bgq_su3_spinor_store(targetsite, result);
		}
	}
}






void bgq_mul_r(bgq_spinorfield targetfield, double c, bgq_spinorfield spinorfield, bool isOdd) {
	assert(!omp_in_parallel());
	//assert(targetfield != spinorfield);
	bgq_spinorfield_setOdd(targetfield, isOdd, true);
	bgq_spinorfield_setOdd(spinorfield, isOdd, false);

	bgq_vector4double_decl(qc);
	bgq_cconst(qc, c, c);

#pragma omp parallel for schedule(static) firstprivate(bgq_vars(qc))
	for (int txy = 0; txy < VOLUME_ZLINES; txy += 1) {
		WORKLOAD_DECL(txy, VOLUME_ZLINES);
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
			bgq_mul(result_v0_c0, qc, spinor_v0_c0);
			bgq_mul(result_v0_c1, qc, spinor_v0_c1);
			bgq_mul(result_v0_c2, qc, spinor_v0_c2);
			bgq_mul(result_v1_c0, qc, spinor_v1_c0);
			bgq_mul(result_v1_c1, qc, spinor_v1_c1);
			bgq_mul(result_v1_c2, qc, spinor_v1_c2);
			bgq_mul(result_v2_c0, qc, spinor_v2_c0);
			bgq_mul(result_v2_c1, qc, spinor_v2_c1);
			bgq_mul(result_v2_c2, qc, spinor_v2_c2);
			bgq_mul(result_v3_c0, qc, spinor_v3_c0);
			bgq_mul(result_v3_c1, qc, spinor_v3_c1);
			bgq_mul(result_v3_c2, qc, spinor_v3_c2);

			bgq_spinorsite *targetsite = BGQ_SPINORSITE(targetfield, isOdd, tv, x, y, z, t1, t2, true, false);
			bgq_su3_spinor_store(targetsite, result);
		}
	}
}


void bgq_assign_add_mul_r(bgq_spinorfield targetfield, bgq_spinorfield spinorfield, bool isOdd, double c) {
	assert(!omp_in_parallel());
	assert(targetfield != spinorfield);
	bgq_spinorfield_setOdd(targetfield, isOdd, false);
	bgq_spinorfield_setOdd(spinorfield, isOdd, false);

	bgq_vector4double_decl(qc);
	bgq_cconst(qc, c, c);

#pragma omp parallel for schedule(static) firstprivate(bgq_vars(qc))
	for (int txy = 0; txy < VOLUME_ZLINES; txy += 1) {
		WORKLOAD_DECL(txy, VOLUME_ZLINES);
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

			bgq_spinorsite *targetsite = BGQ_SPINORSITE(targetfield, isOdd, tv, x, y, z, t1, t2, true, false);
			bgq_su3_spinor_decl(targetspinor);
			bgq_su3_spinor_load(targetspinor, targetsite);

			bgq_su3_spinor_decl(result);
			bgq_madd(result_v0_c0, qc, spinor_v0_c0, targetspinor_v0_c0);
			bgq_madd(result_v0_c1, qc, spinor_v0_c1, targetspinor_v0_c1);
			bgq_madd(result_v0_c2, qc, spinor_v0_c2, targetspinor_v0_c2);
			bgq_madd(result_v1_c0, qc, spinor_v1_c0, targetspinor_v1_c0);
			bgq_madd(result_v1_c1, qc, spinor_v1_c1, targetspinor_v1_c1);
			bgq_madd(result_v1_c2, qc, spinor_v1_c2, targetspinor_v1_c2);
			bgq_madd(result_v2_c0, qc, spinor_v2_c0, targetspinor_v2_c0);
			bgq_madd(result_v2_c1, qc, spinor_v2_c1, targetspinor_v2_c1);
			bgq_madd(result_v2_c2, qc, spinor_v2_c2, targetspinor_v2_c2);
			bgq_madd(result_v3_c0, qc, spinor_v3_c0, targetspinor_v3_c0);
			bgq_madd(result_v3_c1, qc, spinor_v3_c1, targetspinor_v3_c1);
			bgq_madd(result_v3_c2, qc, spinor_v3_c2, targetspinor_v3_c2);

			bgq_su3_spinor_store(targetsite, result);
		}
	}
}


void bgq_assign(bgq_spinorfield targetfield, bgq_spinorfield spinorfield, bool isOdd) {
	bgq_spinorfield_setOdd(targetfield, isOdd, true);
	bgq_spinorfield_setOdd(spinorfield, isOdd, false);
	assert(targetfield != spinorfield);

#pragma omp parallel
	{
		// Parallel memcpy; could do some optimization like alignment
		int n = omp_get_num_threads();
		int i = omp_get_thread_num();

		size_t size = sizeof(bgq_spinorsite) * VOLUME_SITES;
		size_t chunksize = (size + n - 1/* round up */) / n;
		char *to_start = (char*)targetfield + chunksize*i;
		char *from_start = (char*)spinorfield + chunksize*i;

		if ((from_start + chunksize) > (char*)spinorfield) {
			chunksize = (char*)spinorfield + size - from_start;
		}
		assert(chunksize >= 0);
		memcpy(to_start, from_start, chunksize);
	}
}


void bgq_plus(bgq_spinorfield targetfield, bgq_spinorfield lhsfield, bgq_spinorfield rhsfield, bool isOdd) {
	assert(!omp_in_parallel());
	bgq_spinorfield_setOdd(targetfield, isOdd, true);
	bgq_spinorfield_setOdd(lhsfield, isOdd, false);
	bgq_spinorfield_setOdd(rhsfield, isOdd, false);

#pragma omp parallel for schedule(static)
	for (int txy = 0; txy < VOLUME_ZLINES; txy += 1) {
		WORKLOAD_DECL(txy, VOLUME_ZLINES);
		const int tv = WORKLOAD_PARAM(PHYSICAL_LTV);
		const int x = WORKLOAD_PARAM(PHYSICAL_LX);
		const int y = WORKLOAD_PARAM(PHYSICAL_LY);
		WORKLOAD_CHECK

		for (int z = 0; z < PHYSICAL_LZ; z += 1) {
			const int t1 = ((isOdd + x + y + z) & 1) + tv * PHYSICAL_LP * PHYSICAL_LK;
			const int t2 = t1 + 2;

			bgq_spinorsite *lhssite = BGQ_SPINORSITE(lhsfield, isOdd, tv, x, y, z, t1, t2, true, false);
			bgq_su3_spinor_decl(lhs);
			bgq_su3_spinor_load(lhs, lhssite);

			bgq_spinorsite *rhssite = BGQ_SPINORSITE(lhsfield, isOdd, tv, x, y, z, t1, t2, true, false);
			bgq_su3_spinor_decl(rhs);
			bgq_su3_spinor_load(rhs, rhssite);

			bgq_su3_spinor_decl(result);
			bgq_add(result_v0_c0, lhs_v0_c0, rhs_v0_c0);
			bgq_add(result_v0_c1, lhs_v0_c1, rhs_v0_c1);
			bgq_add(result_v0_c2, lhs_v0_c2, rhs_v0_c2);
			bgq_add(result_v1_c0, lhs_v1_c0, rhs_v1_c0);
			bgq_add(result_v1_c1, lhs_v1_c1, rhs_v1_c1);
			bgq_add(result_v1_c2, lhs_v1_c2, rhs_v1_c2);
			bgq_add(result_v2_c0, lhs_v2_c0, rhs_v2_c0);
			bgq_add(result_v2_c1, lhs_v2_c1, rhs_v2_c1);
			bgq_add(result_v2_c2, lhs_v2_c2, rhs_v2_c2);
			bgq_add(result_v3_c0, lhs_v3_c0, rhs_v3_c0);
			bgq_add(result_v3_c1, lhs_v3_c1, rhs_v3_c1);
			bgq_add(result_v3_c2, lhs_v3_c2, rhs_v3_c2);

			bgq_spinorsite *targetsite = BGQ_SPINORSITE(targetfield, isOdd, tv, x, y, z, t1, t2, true, false);
			bgq_su3_spinor_store(targetsite, result);
		}
	}
}


void bgq_diff(bgq_spinorfield targetfield, bgq_spinorfield lhsfield, bgq_spinorfield rhsfield, bool isOdd) {
	assert(!omp_in_parallel());
	bgq_spinorfield_setOdd(targetfield, isOdd, true);
	bgq_spinorfield_setOdd(lhsfield, isOdd, false);
	bgq_spinorfield_setOdd(rhsfield, isOdd, false);

#pragma omp parallel for schedule(static)
	for (int txy = 0; txy < VOLUME_ZLINES; txy += 1) {
		WORKLOAD_DECL(txy, VOLUME_ZLINES);
		const int tv = WORKLOAD_PARAM(PHYSICAL_LTV);
		const int x = WORKLOAD_PARAM(PHYSICAL_LX);
		const int y = WORKLOAD_PARAM(PHYSICAL_LY);
		WORKLOAD_CHECK

		for (int z = 0; z < PHYSICAL_LZ; z += 1) {
			const int t1 = ((isOdd + x + y + z) & 1) + tv * PHYSICAL_LP * PHYSICAL_LK;
			const int t2 = t1 + 2;

			bgq_spinorsite *lhssite = BGQ_SPINORSITE(lhsfield, isOdd, tv, x, y, z, t1, t2, true, false);
			bgq_su3_spinor_decl(lhs);
			bgq_su3_spinor_load(lhs, lhssite);

			bgq_spinorsite *rhssite = BGQ_SPINORSITE(lhsfield, isOdd, tv, x, y, z, t1, t2, true, false);
			bgq_su3_spinor_decl(rhs);
			bgq_su3_spinor_load(rhs, rhssite);

			bgq_su3_spinor_decl(result);
			bgq_sub(result_v0_c0, lhs_v0_c0, rhs_v0_c0);
			bgq_sub(result_v0_c1, lhs_v0_c1, rhs_v0_c1);
			bgq_sub(result_v0_c2, lhs_v0_c2, rhs_v0_c2);
			bgq_sub(result_v1_c0, lhs_v1_c0, rhs_v1_c0);
			bgq_sub(result_v1_c1, lhs_v1_c1, rhs_v1_c1);
			bgq_sub(result_v1_c2, lhs_v1_c2, rhs_v1_c2);
			bgq_sub(result_v2_c0, lhs_v2_c0, rhs_v2_c0);
			bgq_sub(result_v2_c1, lhs_v2_c1, rhs_v2_c1);
			bgq_sub(result_v2_c2, lhs_v2_c2, rhs_v2_c2);
			bgq_sub(result_v3_c0, lhs_v3_c0, rhs_v3_c0);
			bgq_sub(result_v3_c1, lhs_v3_c1, rhs_v3_c1);
			bgq_sub(result_v3_c2, lhs_v3_c2, rhs_v3_c2);

			bgq_spinorsite *targetsite = BGQ_SPINORSITE(targetfield, isOdd, tv, x, y, z, t1, t2, true, false);
			bgq_su3_spinor_store(targetsite, result);
		}
	}
}


void bgq_gamma5(bgq_spinorfield targetfield, bgq_spinorfield spinorfield, bool isOdd) {
	assert(!omp_in_parallel());
	bgq_spinorfield_setOdd(targetfield, isOdd, true);
	bgq_spinorfield_setOdd(spinorfield, isOdd, false);
#pragma omp parallel for schedule(static)
	for (int txy = 0; txy < VOLUME_ZLINES; txy += 1) {
		WORKLOAD_DECL(txy, VOLUME_ZLINES);
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
			bgq_su3_vmov(result_v0, spinor_v0);
			bgq_su3_vmov(result_v1, spinor_v1);
			bgq_neg(result_v2_c0, spinor_v2_c0);
			bgq_neg(result_v2_c1, spinor_v2_c1);
			bgq_neg(result_v2_c2, spinor_v2_c2);
			bgq_neg(result_v3_c0, spinor_v3_c0);
			bgq_neg(result_v3_c1, spinor_v3_c1);
			bgq_neg(result_v3_c2, spinor_v3_c2);

			bgq_spinorsite *targetsite = BGQ_SPINORSITE(targetfield, isOdd, tv, x, y, z, t1, t2, true, false);
			bgq_su3_spinor_store(targetsite, result);
		}
	}
}


void bgq_assign_mul_add_r(bgq_spinorfield targetfield, double c, bgq_spinorfield spinorfield, bool isOdd) {
	assert(!omp_in_parallel());
	bgq_spinorfield_setOdd(targetfield, isOdd, false);

	bgq_vector4double_decl(fact);
	bgq_cconst(fact, c, c);

#pragma omp parallel for schedule(static)
	for (int txy = 0; txy < VOLUME_ZLINES; txy += 1) {
		WORKLOAD_DECL(txy, VOLUME_ZLINES);
		const int tv = WORKLOAD_PARAM(PHYSICAL_LTV);
		const int x = WORKLOAD_PARAM(PHYSICAL_LX);
		const int y = WORKLOAD_PARAM(PHYSICAL_LY);
		WORKLOAD_CHECK

		for (int z = 0; z < PHYSICAL_LZ; z += 1) {
			const int t1 = ((isOdd + x + y + z) & 1) + tv * PHYSICAL_LP * PHYSICAL_LK;
			const int t2 = t1 + 2;

			bgq_spinorsite *targetsite = BGQ_SPINORSITE(targetfield, isOdd, tv, x, y, z, t1, t2, true, false);
			bgq_su3_spinor_decl(spinor);
			bgq_su3_spinor_load(spinor, targetsite);

			bgq_spinorsite *addsite = BGQ_SPINORSITE(spinorfield, isOdd, tv, x, y, z, t1, t2, true, false);
			bgq_su3_spinor_decl(addme);
			bgq_su3_spinor_load(addme, addsite);

			bgq_su3_spinor_decl(result);
			bgq_madd(result_v0_c0, fact, spinor_v0_c0, addme_v0_c0);
			bgq_madd(result_v0_c1, fact, spinor_v0_c1, addme_v0_c1);
			bgq_madd(result_v0_c2, fact, spinor_v0_c2, addme_v0_c2);
			bgq_madd(result_v1_c0, fact, spinor_v1_c0, addme_v1_c0);
			bgq_madd(result_v1_c1, fact, spinor_v1_c1, addme_v1_c1);
			bgq_madd(result_v1_c2, fact, spinor_v1_c2, addme_v1_c2);
			bgq_madd(result_v2_c0, fact, spinor_v2_c0, addme_v2_c0);
			bgq_madd(result_v2_c1, fact, spinor_v2_c1, addme_v2_c1);
			bgq_madd(result_v2_c2, fact, spinor_v2_c2, addme_v2_c2);
			bgq_madd(result_v3_c0, fact, spinor_v3_c0, addme_v3_c0);
			bgq_madd(result_v3_c1, fact, spinor_v3_c1, addme_v3_c1);
			bgq_madd(result_v3_c2, fact, spinor_v3_c2, addme_v3_c2);

			bgq_su3_spinor_store(targetsite, result);
		}
	}

	bgq_spinorfield_setOdd(targetfield, isOdd, true);
	bgq_spinorfield_setOdd(spinorfield, isOdd, false);
}


void bgq_mul_one_pm_imu_sub_mul_gamma5(bgq_spinorfield l, bgq_spinorfield k, bgq_spinorfield j, bool isOdd, double sign) {
	assert(!omp_in_parallel());
	bgq_spinorfield_setOdd(l, isOdd, true);
	bgq_spinorfield_setOdd(k, isOdd, false);
	bgq_spinorfield_setOdd(j, isOdd, false);

	sign = (sign >= 0) ? 1 : -1;

	bgq_vector4double_decl(qz);
	bgq_cconst(qz, 1, sign * g_mu);

	bgq_vector4double_decl(qw);
	bgq_cconst(qw, 1, -sign * g_mu);

#pragma omp parallel for schedule(static)
	for (int txy = 0; txy < VOLUME_ZLINES; txy += 1) {
		WORKLOAD_DECL(txy, VOLUME_ZLINES);
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
			bgq_su3_cvmul(phi_v0, qz, r_v0);
			bgq_su3_cvmul(phi_v1, qz, r_v1);
			bgq_su3_cvmul(phi_v2, qw, r_v2);
			bgq_su3_cvmul(phi_v3, qw, r_v3);

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


complexdouble bgq_scalar_prod(bgq_spinorfield spinorfield1, bgq_spinorfield spinorfield2, bool isOdd, bool parallel) {
	assert(!omp_in_parallel());
	bgq_spinorfield_setOdd(spinorfield1, isOdd, false);
	bgq_spinorfield_setOdd(spinorfield2, isOdd, false);

	bgq_vector4double_decl(shared_sum);
	bgq_zero(shared_sum);

#pragma omp parallel
	{
		bgq_vector4double_decl(thread_sum);
		bgq_zero(thread_sum);

#pragma omp for schedule(static) nowait
		for (int txy = 0; txy < VOLUME_ZLINES; txy+=1) {
			WORKLOAD_DECL(txy, VOLUME_ZLINES);
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
#pragma tm_atomic /* safe_mode */
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
		MPI_Allreduce(&node_result, &result, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
	} else {
		result = node_result;
	}

	return result;
}


double bgq_scalar_prod_r(bgq_spinorfield spinorfield1, bgq_spinorfield spinorfield2, bool isOdd, bool parallel) {
	assert(!omp_in_parallel());
	assert(spinorfield1 != spinorfield2); // Use bgq_scalar_prod is they are equal
	bgq_spinorfield_setOdd(spinorfield1, isOdd, false);
	bgq_spinorfield_setOdd(spinorfield2, isOdd, false);

	bgq_vector4double_decl(shared_sum);
	bgq_zero(shared_sum);

#pragma omp parallel
	{
		bgq_vector4double_decl(thread_sum);
		bgq_zero(thread_sum);

#pragma omp for schedule(static) nowait
		for (int txy = 0; txy < VOLUME_ZLINES; txy+=1) {
			WORKLOAD_DECL(txy, VOLUME_ZLINES);
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

				bgq_madd(thread_sum, spinor1_v0_c0, spinor2_v0_c0, thread_sum);
				bgq_madd(thread_sum, spinor1_v0_c1, spinor2_v0_c1, thread_sum);
				bgq_madd(thread_sum, spinor1_v0_c2, spinor2_v0_c2, thread_sum);
				bgq_madd(thread_sum, spinor1_v1_c0, spinor2_v1_c0, thread_sum);
				bgq_madd(thread_sum, spinor1_v1_c1, spinor2_v1_c1, thread_sum);
				bgq_madd(thread_sum, spinor1_v1_c2, spinor2_v1_c2, thread_sum);
				bgq_madd(thread_sum, spinor1_v2_c0, spinor2_v2_c0, thread_sum);
				bgq_madd(thread_sum, spinor1_v2_c1, spinor2_v2_c1, thread_sum);
				bgq_madd(thread_sum, spinor1_v2_c2, spinor2_v2_c2, thread_sum);
				bgq_madd(thread_sum, spinor1_v3_c0, spinor2_v3_c0, thread_sum);
				bgq_madd(thread_sum, spinor1_v3_c1, spinor2_v3_c1, thread_sum);
				bgq_madd(thread_sum, spinor1_v3_c2, spinor2_v3_c2, thread_sum);
			}
		}

#ifdef XLC
#pragma tm_atomic /* safe_mode */
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
		MPI_Allreduce(&node_result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	} else {
		result = node_result;
	}

	return result;
}


double bgq_square_norm(bgq_spinorfield spinorfield, bool isOdd, bool parallel) {
	assert(!omp_in_parallel());
	bgq_spinorfield_setOdd(spinorfield, isOdd, false);

	bgq_vector4double_decl(shared_sum);
	bgq_zero(shared_sum);

#pragma omp parallel
	{
		bgq_vector4double_decl(thread_sum);
		bgq_zero(thread_sum);

#pragma omp for schedule(static) nowait
		for (int txy = 0; txy < VOLUME_ZLINES; txy+=1) {
			WORKLOAD_DECL(txy, VOLUME_ZLINES);
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
#pragma tm_atomic /* safe_mode */
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
		MPI_Allreduce(&node_result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	} else {
		result = node_result;
	}

	return result;
}



