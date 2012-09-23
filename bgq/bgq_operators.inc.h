/*
 * bgq_operators.inc.h
 *
 *  Created on: Sep 23, 2012
 *      Author: meinersbur
 */

#include "bgq_field.inc.h"

#define bgq_mul_one_pm_imu_sub_mul_gamma5 NAME2(bgq_mul_one_pm_imu_sub_mul_gamma5,PRECISION)
void bgq_mul_one_pm_imu_sub_mul_gamma5(bgq_spinorfield l, bgq_spinorfield k, bgq_spinorfield j, bool isOdd, double sign);

#define bgq_mul_one_pm_imu_inv NAME2(bgq_mul_one_pm_imu_inv,PRECISION)
void bgq_mul_one_pm_imu_inv(bgq_spinorfield spinorfield, bool isOdd, double sign);

#define bgq_scalar_prod NAME2(bgq_scalar_prod,PRECISION)
complexdouble bgq_scalar_prod(bgq_spinorfield spinorfield1, bgq_spinorfield spinorfield2, bool isOdd, bool parallel);

#define bgq_square_norm NAME2(bgq_square_norm,PRECISION)
double bgq_square_norm(bgq_spinorfield spinorfield, bool isOdd, bool parallel);
