/*
 * bgq_operators.inc.h
 *
 *  Created on: Sep 23, 2012
 *      Author: meinersbur
 */

#include "bgq_field.inc.h"

// targetfield[].v[0,1] = z * spinorfield[].v[0,1]
// targetfield[].v[2,3] = w * spinorfield[].v[2,3]
#define bgq_spinorfield_mul_weyl_complex NAME2(bgq_spinorfield_mul_weyl_complex,PRECISION)
void bgq_spinorfield_mul_weyl_complex(bgq_spinorfield targetfield, bgq_spinorfield spinorfield, bool isOdd, complexdouble z, complexdouble w);

#define bgq_mul_r NAME2(bgq_mul_r,PRECISION)
void bgq_mul_r(bgq_spinorfield targetfield, double c, bgq_spinorfield spinorfield, bool isOdd);

//#define bgq_assign_mul_one_pm_imu NAME2(bgq_assign_mul_one_pm_imu,PRECISION)
//void bgq_assign_mul_one_pm_imu(bgq_spinorfield targetfield, bgq_spinorfield spinorfield, bool isOdd, double sign);

#define bgq_assign_add_mul_r NAME2(bgq_assign_add_mul_r,PRECISION)
void bgq_assign_add_mul_r(bgq_spinorfield targetfield, bgq_spinorfield spinorfield, bool isOdd, double c);

#define bgq_assign NAME2(bgq_assign,PRECISION)
void bgq_assign(bgq_spinorfield targetfield, bgq_spinorfield spinorfield, bool isOdd);

#define bgq_plus NAME2(bgq_plus,PRECISION)
void bgq_plus(bgq_spinorfield targetfield, bgq_spinorfield lhsfield, bgq_spinorfield rhsfield, bool isOdd);

#define bgq_diff NAME2(bgq_diff,PRECISION)
void bgq_diff(bgq_spinorfield targetfield, bgq_spinorfield lhsfield, bgq_spinorfield rhsfield, bool isOdd);

#define bgq_gamma5 NAME2(bgq_gamma5,PRECISION)
void bgq_gamma5(bgq_spinorfield targetfield, bgq_spinorfield spinorfield, bool isOdd);

#define bgq_assign_mul_add_r NAME2(bgq_assign_mul_add_r,PRECISION)
void bgq_assign_mul_add_r(bgq_spinorfield targerfield, double c, bgq_spinorfield spinorfield, bool isOdd);

//#define bgq_assign_mul_one_pm_imu_inv NAME2(bgq_assign_mul_one_pm_imu_inv,PRECISION)
//void bgq_assign_mul_one_pm_imu_inv(bgq_spinorfield targetfield, bgq_spinorfield spinorfield, bool isOdd, double sign);

#define bgq_mul_one_pm_imu_sub_mul_gamma5 NAME2(bgq_mul_one_pm_imu_sub_mul_gamma5,PRECISION)
void bgq_mul_one_pm_imu_sub_mul_gamma5(bgq_spinorfield l, bgq_spinorfield k, bgq_spinorfield j, bool isOdd, double sign);

//#define bgq_mul_one_pm_imu_inv NAME2(bgq_mul_one_pm_imu_inv,PRECISION)
//void bgq_mul_one_pm_imu_inv(bgq_spinorfield spinorfield, bool isOdd, double sign);

#define bgq_scalar_prod NAME2(bgq_scalar_prod,PRECISION)
complexdouble bgq_scalar_prod(bgq_spinorfield spinorfield1, bgq_spinorfield spinorfield2, bool isOdd, bool parallel);

#define bgq_scalar_prod_r NAME2(bgq_scalar_prod_r,PRECISION)
double bgq_scalar_prod_r(bgq_spinorfield spinorfield1, bgq_spinorfield spinorfield2, bool isOdd, bool parallel);

#define bgq_square_norm NAME2(bgq_square_norm,PRECISION)
double bgq_square_norm(bgq_spinorfield spinorfield, bool isOdd, bool parallel);
