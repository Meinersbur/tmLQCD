/*
 * bgq_stdoperators.c
 *
 *  Created on: Nov 21, 2012
 *      Author: meinersbur
 */

#include "bgq_stdoperators.h"

#include "bgq_field.h"
#include "bgq_qpx.h"



static inline void bgq_site_diff(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), ucoord ic) {
	bgq_su3_spinor_decl(result);

	bgq_su3_vsub(result_v0, spinor1_v0, spinor2_v0);
	bgq_su3_vsub(result_v1, spinor1_v1, spinor2_v1);
	bgq_su3_vsub(result_v2, spinor1_v2, spinor2_v2);
	bgq_su3_vsub(result_v3, spinor1_v3, spinor2_v3);

	bgq_su3_spinor_mov(*target, result);
}


static inline void bgq_site_set(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor), ucoord ic) {
	bgq_su3_spinor_decl(result);

	bgq_su3_spinor_mov(result, spinor);

	bgq_su3_spinor_mov(*target, result);
}


static inline void bgq_site_imul(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor), bgq_params(qz), bgq_params(qw), ucoord ic) {
	bgq_su3_spinor_decl(result);

	bgq_su3_cvmul(result_v0,qz,spinor_v0);
	bgq_su3_cvmul(result_v1,qz,spinor_v1);
	bgq_su3_cvmul(result_v2,qw,spinor_v2);
	bgq_su3_cvmul(result_v3,qw,spinor_v3);

	bgq_su3_spinor_mov(*target, result);
}



#define OPERATOR_NAME bgq_spinorfield_diff
#define OPERATOR_ARGFIELDS 2
#define OPERATOR_VECSITEFUNC bgq_site_diff
#include "bgq_operator.inc.c"


#define OPERATOR_NAME bgq_spinorfield_set
#define OPERATOR_ARGFIELDS 0
#define OPERATOR_VECSITEFUNC bgq_site_set
#define OPERATOR_EXTRAPARMS bgq_su3_spinor_params(spinor)
#define OPERATOR_EXTRAARGS bgq_su3_spinor_vars(spinor)
#include "bgq_operator.inc.c"


#define OPERATOR_NAME bgq_spinorfield_imul_raw
#define OPERATOR_ARGFIELDS 1
#define OPERATOR_VECSITEFUNC bgq_site_imul
#define OPERATOR_EXTRAPARMS bgq_params(qz),bgq_params(qw)
#define OPERATOR_EXTRAARGS bgq_vars(qz),bgq_vars(qw)
#include "bgq_operator.inc.c"

void bgq_spinorfield_imul_double(bgq_weylfield_controlblock *targetfield, bool isOdd, bgq_weylfield_controlblock *sourcefield, complex_double z, complex_double w) {
	bgq_vector4double_decl(qz);
	bgq_cconst(qz, creal(z), cimag(z));
	bgq_vector4double_decl(qw);
	bgq_cconst(qw, creal(w), cimag(w));
	bgq_spinorfield_imul_raw_double(targetfield, isOdd, sourcefield, bgq_vars(qz), bgq_vars(qw));
}

void bgq_spinorfield_imul_float(bgq_weylfield_controlblock *targetfield, bool isOdd, bgq_weylfield_controlblock *sourcefield, complex_double z, complex_double w) {
	bgq_vector4double_decl(qz);
	bgq_cconst(qz, creal(z), cimag(z));
	bgq_vector4double_decl(qw);
	bgq_cconst(qw, creal(w), cimag(w));
	bgq_spinorfield_imul_raw_float(targetfield, isOdd, sourcefield, bgq_vars(qz), bgq_vars(qw));
}


#if BGQ_REPLACE
void assign_mul_one_pm_imu_inv(spinor * const l, spinor * const k, const double _sign, const int N) {
	bgq_weylfield_controlblock *targetfield = bgq_translate_spinorfield(l);
	bgq_weylfield_controlblock *sourcefield = bgq_translate_spinorfield(k);
	bool isOdd = bgq_spinorfield_isOdd(sourcefield);

	 _Complex double z,w;
	  int ix;
	  double sign=-1.;
	  spinor *r, *s;
	  double nrm = 1./(1.+g_mu*g_mu);

	  if(_sign < 0.){
	    sign = 1.;
	  }

	  z = nrm + (sign * nrm * g_mu) * I;
	  w = conj(z);

	  bgq_spinorfield_imul_double(targetfield, isOdd, sourcefield, z, w);
}
#endif
