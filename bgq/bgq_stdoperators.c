/*
 * bgq_stdoperators.c
 *
 *  Created on: Nov 21, 2012
 *      Author: meinersbur
 */

#include "bgq_stdoperators.h"

#include "bgq_field.h"
#include "bgq_qpx.h"



static inline void bgq_site_sub(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), ucoord ic) {
	bgq_su3_spinor_decl(result);

	bgq_su3_vsub(result_v0, spinor1_v0, spinor2_v0);
	bgq_su3_vsub(result_v1, spinor1_v1, spinor2_v1);
	bgq_su3_vsub(result_v2, spinor1_v2, spinor2_v2);
	bgq_su3_vsub(result_v3, spinor1_v3, spinor2_v3);

	bgq_su3_spinor_mov(*target, result);
}


static inline void bgq_site_cmul_plain_sub(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), bgq_params(qc), ucoord ic) {
	bgq_su3_spinor_decl(result);
	bgq_su3_spinor_decl(cmul1);

	// Unfortunately, there is no fused-multiply-sub for this
	bgq_mul(cmul1_v0_c0, qc, spinor1_v0_c0);
	bgq_mul(cmul1_v0_c1, qc, spinor1_v0_c1);
	bgq_mul(cmul1_v0_c2, qc, spinor1_v0_c2);
	bgq_mul(cmul1_v1_c0, qc, spinor1_v1_c0);
	bgq_mul(cmul1_v1_c1, qc, spinor1_v1_c1);
	bgq_mul(cmul1_v1_c2, qc, spinor1_v1_c2);
	bgq_mul(cmul1_v2_c0, qc, spinor1_v2_c0);
	bgq_mul(cmul1_v2_c1, qc, spinor1_v2_c1);
	bgq_mul(cmul1_v2_c2, qc, spinor1_v2_c2);
	bgq_mul(cmul1_v3_c0, qc, spinor1_v3_c0);
	bgq_mul(cmul1_v3_c1, qc, spinor1_v3_c1);
	bgq_mul(cmul1_v3_c2, qc, spinor1_v3_c2);

	bgq_su3_vsub(result_v0, cmul1_v0, spinor2_v0);
	bgq_su3_vsub(result_v1, cmul1_v1, spinor2_v1);
	bgq_su3_vsub(result_v2, cmul1_v2, spinor2_v2);
	bgq_su3_vsub(result_v3, cmul1_v3, spinor2_v3);

	bgq_su3_spinor_mov(*target, result);
}


static inline void bgq_site_icjgmul_plain_sub(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), bgq_params(qc), bgq_params(qw), ucoord ic) {
	bgq_su3_spinor_decl(result);
	bgq_su3_spinor_decl(imul1);

	assert(bgq_cmplxval1(qc) == bgq_cmplxval2(qc));

	bgq_su3_cvmul(imul1_v0, qc, spinor1_v0);
	bgq_su3_cvmul(imul1_v1, qc, spinor1_v1);
	bgq_su3_vsub(result_v0, imul1_v0, spinor2_v0);
	bgq_su3_vsub(result_v1, imul1_v1, spinor2_v1);

#if 1
	assert(bgq_cmplxval1(qw) == bgq_cmplxval2(qw));
	assert(bgq_cmplxval1(qc) == -bgq_cmplxval1(qw));

// re = spinor2.re + qw.re * spinor1.re + qw.im * spinor1.im
// im = spinor2.im + qw.re * spinor1.im - qw.im * spinor1.re

	bgq_xmadd(result_v2_c0, qw, spinor1_v2_c0, spinor2_v2_c0);
	bgq_xxcpnmadd(result_v2_c0, spinor1_v2_c0, qw, result_v2_c0);
	bgq_xmadd(result_v2_c1, qw, spinor1_v2_c1, spinor2_v2_c1);
	bgq_xxcpnmadd(result_v2_c1, spinor1_v2_c1, qw, result_v2_c1);
	bgq_xmadd(result_v2_c2, qw, spinor1_v2_c2, spinor2_v2_c2);
	bgq_xxcpnmadd(result_v2_c2, spinor1_v2_c2, qw, result_v2_c2);

	bgq_xmadd(result_v3_c0, qw, spinor1_v3_c0, spinor2_v3_c0);
	bgq_xxcpnmadd(result_v3_c0, spinor1_v3_c0, qw, result_v3_c0);
	bgq_xmadd(result_v3_c1, qw, spinor1_v3_c1, spinor2_v3_c1);
	bgq_xxcpnmadd(result_v3_c1, spinor1_v3_c1, qw, result_v3_c1);
	bgq_xmadd(result_v3_c2, qw, spinor1_v3_c2, spinor2_v3_c2);
	bgq_xxcpnmadd(result_v3_c2, spinor1_v3_c2, qw, result_v3_c2);
#else
	bgq_su3_cjgvmul(imul1_v2, qc, spinor1_v2);
	bgq_su3_cjgvmul(imul1_v3, qc, spinor1_v3);
	bgq_su3_vsub(result_v2, spinor2_v2, imul1_v2);
	bgq_su3_vsub(result_v3, spinor2_v3, imul1_v3);
#endif

	bgq_su3_spinor_mov(*target, result);
}


static inline void bgq_site_set(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor), ucoord ic) {
	bgq_su3_spinor_decl(result);

	bgq_su3_spinor_mov(result, spinor);

	bgq_su3_spinor_mov(*target, result);
}


static inline void bgq_site_cmul(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor), bgq_params(qc), ucoord ic) {
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


static inline void bgq_site_cmul_plain_add(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), bgq_params(qc), ucoord ic) {
	bgq_su3_spinor_decl(result);

	bgq_madd(result_v0_c0, qc, spinor1_v0_c0, spinor2_v0_c0);
	bgq_madd(result_v0_c1, qc, spinor1_v0_c1, spinor2_v0_c1);
	bgq_madd(result_v0_c2, qc, spinor1_v0_c2, spinor2_v0_c2);
	bgq_madd(result_v1_c0, qc, spinor1_v1_c0, spinor2_v1_c0);
	bgq_madd(result_v1_c1, qc, spinor1_v1_c1, spinor2_v1_c1);
	bgq_madd(result_v1_c2, qc, spinor1_v1_c2, spinor2_v1_c2);
	bgq_madd(result_v2_c0, qc, spinor1_v2_c0, spinor2_v2_c0);
	bgq_madd(result_v2_c1, qc, spinor1_v2_c1, spinor2_v2_c1);
	bgq_madd(result_v2_c2, qc, spinor1_v2_c2, spinor2_v2_c2);
	bgq_madd(result_v3_c0, qc, spinor1_v3_c0, spinor2_v3_c0);
	bgq_madd(result_v3_c1, qc, spinor1_v3_c1, spinor2_v3_c1);
	bgq_madd(result_v3_c2, qc, spinor1_v3_c2, spinor2_v3_c2);

	bgq_su3_spinor_mov(*target, result);
}


static inline void bgq_site_gamma5(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor), ucoord ic) {
	bgq_su3_spinor_decl(result);

	bgq_mov(result_v0_c0, spinor_v0_c0);
	bgq_mov(result_v0_c1, spinor_v0_c1);
	bgq_mov(result_v0_c2, spinor_v0_c2);
	bgq_mov(result_v1_c0, spinor_v1_c0);
	bgq_mov(result_v1_c1, spinor_v1_c1);
	bgq_mov(result_v1_c2, spinor_v1_c2);
	bgq_neg(result_v2_c0, spinor_v2_c0);
	bgq_neg(result_v2_c1, spinor_v2_c1);
	bgq_neg(result_v2_c2, spinor_v2_c2);
	bgq_neg(result_v3_c0, spinor_v3_c0);
	bgq_neg(result_v3_c1, spinor_v3_c1);
	bgq_neg(result_v3_c2, spinor_v3_c2);

	bgq_su3_spinor_mov(*target, result);
}



#define OPERATOR_NAME bgq_spinorfield_sub
#define OPERATOR_ARGFIELDS 2
#define OPERATOR_VECSITEFUNC bgq_site_sub
#include "bgq_operator.inc.c"

#if BGQ_REPLACE
void diff(spinor * const Q, const spinor * const R, const spinor * const S, const int N) {
	assert(N == VOLUME/2);
	bgq_weylfield_controlblock *targetfield = bgq_translate_spinorfield(Q);
	bgq_weylfield_controlblock *sourcefield1 = bgq_translate_spinorfield(R);
	bgq_weylfield_controlblock *sourcefield2 = bgq_translate_spinorfield(S);

	bgq_spinorfield_sub_double(targetfield, tri_unknown, sourcefield1, sourcefield2);
}
#endif

#define OPERATOR_NAME bgq_spinorfield_set
#define OPERATOR_ARGFIELDS 0
#define OPERATOR_VECSITEFUNC bgq_site_set
#define OPERATOR_EXTRAPARMS bgq_su3_spinor_params(spinor)
#define OPERATOR_EXTRAARGS bgq_su3_spinor_vars(spinor)
#include "bgq_operator.inc.c"


#define OPERATOR_NAME bgq_spinorfield_copy_raw
#define OPERATOR_ARGFIELDS 1
#define OPERATOR_VECSITEFUNC bgq_site_set
#include "bgq_operator.inc.c"

void bgq_spinorfield_copy(bgq_weylfield_controlblock *target, bgq_spinorfield_layout targetLayout, bgq_weylfield_controlblock *source) {
	assert(target);
	assert(source);
	assert(source != target);

	switch (targetLayout) {
	case ly_full_double:
		if (source->has_fulllayout_double) {
			bgq_spinorfield_prepareRead(source, tri_unknown, false, true, false, false, false);
			bgq_spinorfield_prepareWrite(target, source->isOdd, ly_full_double, false);
			bgq_master_memcpy(target->sec_fullspinor_double, source->sec_fullspinor_double, PHYSICAL_VOLUME * sizeof(*source->sec_fullspinor_double));
		} else {
			bgq_spinorfield_copy_raw_double(target, source->isOdd, source);
		}
		break;
	case ly_full_float:
		if (source->has_fulllayout_float) {
			bgq_spinorfield_prepareRead(source, tri_unknown, false, false, true, false, false);
			bgq_spinorfield_prepareWrite(target, source->isOdd, ly_full_float, false);
			bgq_master_memcpy(target->sec_fullspinor_float, source->sec_fullspinor_float, PHYSICAL_VOLUME * sizeof(*source->sec_fullspinor_float));
		} else {
			bgq_spinorfield_copy_raw_float(target, source->isOdd, source);
		}
		break;
	case ly_legacy:
		bgq_spinorfield_prepareRead(source, false, false, false, false, false, true);
		bgq_spinorfield_prepareWrite(target, source->isOdd, ly_legacy, false);
		bgq_master_memcpy(target->legacy_field, source->legacy_field, VOLUME/2 * sizeof(*source->legacy_field));
		break;
	default:
		master_error(1, "Cannot copy to layout");
	}

	assert(bgq_spinorfield_prepareRead(target, source->isOdd, true, true, true, true, true) == targetLayout);
}

#if BGQ_REPLACE
void assign(spinor * const R, spinor * const S, const int N) {
	assert(N == VOLUME/2);
	bgq_weylfield_controlblock *target = bgq_translate_spinorfield(R);
	bgq_weylfield_controlblock *source = bgq_translate_spinorfield(S);

	bgq_spinorfield_copy(target, ly_full_double, source);
}
#endif


#define OPERATOR_NAME bgq_spinorfield_imul_raw
#define OPERATOR_ARGFIELDS 1
#define OPERATOR_VECSITEFUNC bgq_site_imul
#define OPERATOR_EXTRAPARMS bgq_params(qz),bgq_params(qw)
#define OPERATOR_EXTRAARGS bgq_vars(qz),bgq_vars(qw)
#include "bgq_operator.inc.c"

void bgq_spinorfield_imul_double(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield, complex_double z, complex_double w) {
	bgq_vector4double_decl(qz);
	bgq_cconst(qz, creal(z), cimag(z));
	bgq_vector4double_decl(qw);
	bgq_cconst(qw, creal(w), cimag(w));
	bgq_spinorfield_imul_raw_double(targetfield, isOdd, sourcefield, bgq_vars(qz), bgq_vars(qw));
}

void bgq_spinorfield_imul_float(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield, complex_double z, complex_double w) {
	bgq_vector4double_decl(qz);
	bgq_cconst(qz, creal(z), cimag(z));
	bgq_vector4double_decl(qw);
	bgq_cconst(qw, creal(w), cimag(w));
	bgq_spinorfield_imul_raw_float(targetfield, isOdd, sourcefield, bgq_vars(qz), bgq_vars(qw));
}



#define OPERATOR_NAME bgq_spinorfield_cmul_plain_add_raw
#define OPERATOR_ARGFIELDS 2
#define OPERATOR_VECSITEFUNC bgq_site_cmul_plain_add
#define OPERATOR_EXTRAPARMS bgq_params(qc)
#define OPERATOR_EXTRAARGS bgq_vars(qc)
#include "bgq_operator.inc.c"


void bgq_spinorfield_cmul_plain_add(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *spinorfield1, bgq_weylfield_controlblock *spinorfield2, double c) {
	bgq_vector4double_decl(qc);
	bgq_cconst(qc, c, c);
	bgq_spinorfield_cmul_plain_add_raw_double(targetfield, isOdd, spinorfield1, spinorfield2, bgq_vars(qc));
}



#if BGQ_REPLACE
void assign_mul_one_pm_imu_inv(spinor * const l, spinor * const k, const double _sign, const int N) {
	assert(N == VOLUME/2);
	bgq_weylfield_controlblock *targetfield = bgq_translate_spinorfield(l);
	bgq_weylfield_controlblock *sourcefield = bgq_translate_spinorfield(k);

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

	  bgq_spinorfield_imul_double(targetfield, tri_unknown, sourcefield, z, w);
}


void assign_mul_add_r(spinor * const R, const double c, const spinor * const S, const int N) {
	assert(N == VOLUME/2);
	bgq_weylfield_controlblock *targetfield = bgq_translate_spinorfield(R);
	bgq_weylfield_controlblock *sourcefield = bgq_translate_spinorfield(S);

	bgq_spinorfield_cmul_plain_add(targetfield, tri_unknown, targetfield, sourcefield, c);
}
#endif


#define OPERATOR_NAME bgq_spinorfield_gamma5
#define OPERATOR_ARGFIELDS 1
#define OPERATOR_VECSITEFUNC bgq_site_gamma5
#include "bgq_operator.inc.c"

#if BGQ_REPLACE
void gamma5(spinor * const l, spinor * const k, const int V) {
	assert(V == VOLUME/2);
	bgq_weylfield_controlblock *targetfield = bgq_translate_spinorfield(l);
	bgq_weylfield_controlblock *sourcefield = bgq_translate_spinorfield(k);

	bgq_spinorfield_gamma5_double(targetfield, tri_unknown, sourcefield);
}
#endif




#define OPERATOR_NAME bgq_spinorfield_cmul_raw
#define OPERATOR_ARGFIELDS 1
#define OPERATOR_VECSITEFUNC bgq_site_cmul
#define OPERATOR_EXTRAPARMS bgq_params(qc)
#define OPERATOR_EXTRAARGS bgq_vars(qc)
#include "bgq_operator.inc.c"

void bgq_spinorfield_cmul_double(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield, double c) {
	bgq_vector4double_decl(qc);
	bgq_cconst(qc, c, c);
	bgq_spinorfield_cmul_raw_double(targetfield, isOdd, sourcefield, bgq_vars(qc));
}

void bgq_spinorfield_cmul_float(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield, double c) {
	bgq_vector4double_decl(qc);
	bgq_cconst(qc, c, c);
	bgq_spinorfield_cmul_raw_float(targetfield, isOdd, sourcefield, bgq_vars(qc));
}



#define OPERATOR_NAME bgq_spinorfield_cmul_plain_sub_raw
#define OPERATOR_ARGFIELDS 2
#define OPERATOR_VECSITEFUNC bgq_site_cmul_plain_sub
#define OPERATOR_EXTRAPARMS bgq_params(qc)
#define OPERATOR_EXTRAARGS bgq_vars(qc)
#include "bgq_operator.inc.c"

void bgq_spinorfield_cmul_plain_sub_double(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield1, bgq_weylfield_controlblock *sourcefield2, double c) {
	bgq_vector4double_decl(qc);
	bgq_cconst(qc, c, c);
	bgq_spinorfield_cmul_plain_sub_raw_double(targetfield, isOdd, sourcefield1, sourcefield2,bgq_vars(qc));
}

void bgq_spinorfield_cmul_plain_sub_float(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield1, bgq_weylfield_controlblock *sourcefield2, double c) {
	bgq_vector4double_decl(qc);
	bgq_cconst(qc, c, c);
	bgq_spinorfield_cmul_plain_sub_raw_float(targetfield, isOdd, sourcefield1, sourcefield2, bgq_vars(qc));
}


#define OPERATOR_NAME bgq_spinorfield_icjgmul_plain_sub_raw
#define OPERATOR_ARGFIELDS 2
#define OPERATOR_VECSITEFUNC bgq_site_icjgmul_plain_sub
#define OPERATOR_EXTRAPARMS bgq_params(qc),bgq_params(qw)
#define OPERATOR_EXTRAARGS bgq_vars(qc),bgq_vars(qw)
#include "bgq_operator.inc.c"

void bgq_spinorfield_icjgmul_plain_sub_double(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield1, bgq_weylfield_controlblock *sourcefield2, complex_double c) {
	bgq_vector4double_decl(qc);
	bgq_cconst(qc, creal(c), cimag(c));
	bgq_vector4double_decl(qw);
	bgq_neg(qw, qc);
	bgq_spinorfield_icjgmul_plain_sub_raw_double(targetfield, isOdd, sourcefield1, sourcefield2, bgq_vars(qc), bgq_vars(qw));
}

void bgq_spinorfield_icjgmul_plain_sub_float(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield1, bgq_weylfield_controlblock *sourcefield2, complex_double c) {
	bgq_vector4double_decl(qc);
	bgq_cconst(qc, creal(c), cimag(c));
	bgq_vector4double_decl(qw);
	bgq_neg(qw, qc);
	bgq_spinorfield_icjgmul_plain_sub_raw_float(targetfield, isOdd, sourcefield1, sourcefield2, bgq_vars(qc), bgq_vars(qw));
}


