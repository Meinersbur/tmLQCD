/*
 * bgq_operators.c
 *
 *  Created on: Sep 23, 2012
 *      Author: meinersbur
 */

#include "bgq_field_double.h"


#if BGQ_REPLACE
void mul_one_pm_imu_sub_mul_gamma5(spinor * const l, spinor * const k, spinor * const j, const double _sign) {
	bgq_spinorfield spinorfield_l = bgq_translate_spinorfield(l);
	bgq_spinorfield spinorfield_k = bgq_translate_spinorfield(k);
	bgq_spinorfield spinorfield_j = bgq_translate_spinorfield(j);


	bool isOdd = bgq_spinorfield_isOdd(spinorfield_l);
	assert(isOdd == bgq_spinorfield_isOdd(spinorfield_k));
	assert(isOdd == bgq_spinorfield_isOdd(spinorfield_j));
	bgq_mul_one_pm_imu_sub_mul_gamma5(spinorfield_l, spinorfield_k, spinorfield_j, _sign);
}


void mul_one_pm_imu_inv(spinor * const l, const double _sign) {
	bgq_spinorfield spinorfield = bgq_translate_spinorfield(l);

	bool isOdd = bgq_spinorfield_isOdd(spinorfield);
	bgq_mul_one_pm_imu_inv(spinorfield, isOdd, _sign);
}


complex scalar_prod(spinor * const S, spinor * const R, const int N, const int parallel) {
	bgq_spinorfield spinorfield1 = bgq_translate_spinorfield(S);
	bgq_spinorfield spinorfield2 = bgq_translate_spinorfield(R);
	assert(N==VOLUME/2);

	bool isOdd = bgq_spinorfield_isOdd(spinorfield1);
	assert(isOdd == bgq_spinorfield_isOdd(spinorfield2));
	return bgq_scalar_prod_double(spinorfield1, spinorfield2, isOdd, parallel);
}


double square_norm(spinor * const P, const int N, const int parallel) {
	bgq_spinorfield spinorfield = bgq_translate_spinorfield(P);
	assert(N==VOLUME/2);

	bool isOdd = bgq_spinorfield_isOdd(spinorfield);
	bgq_square_norm_double(spinorfield, parallel);
}


#endif
