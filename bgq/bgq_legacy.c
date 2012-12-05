/*
 * bgq_legacy.c
 *
 *  Created on: Dec 5, 2012
 *      Author: meinersbur
 */

#if BGQ_REPLACE

#include "bgq_stdoperators.h"
#include "bgq_stdreductions.h"
#include "bgq_HoppingMatrix.h"


double assign_mul_add_r_and_square(spinor * const R, const double c, const spinor * const S, const int N, const int parallel) {
	assert(N == VOLUME/2);
	bgq_weylfield_controlblock *targetfield = bgq_translate_spinorfield(R);
	bgq_weylfield_controlblock *sourcefield = bgq_translate_spinorfield(S);

	bgq_spinorfield_rmul_plain_add_double(targetfield, tri_unknown, targetfield, sourcefield, c);
	if (parallel)
		return bgq_spinorfield_sqrnorm_global(tri_unknown, targetfield);
	else
		return bgq_spinorfield_sqrnorm_local(tri_unknown, targetfield);
}


_Complex double scalar_prod(spinor *S, spinor *R, int N, int parallel) {
	assert(N == LOCAL_VOLUME/2);
	bgq_weylfield_controlblock *field1 = bgq_translate_spinorfield(S);
	bgq_weylfield_controlblock *field2 = bgq_translate_spinorfield(R);

	if (parallel)
		return bgq_spinorfield_innerprod_global(tri_unknown, field1, field2);
	else
		return bgq_spinorfield_innerprod_local(tri_unknown, field1, field2);
}


double scalar_prod_r(spinor *S, const spinor * const R, int N, int parallel) {
	assert(N == LOCAL_VOLUME/2);
	bgq_weylfield_controlblock *field1 = bgq_translate_spinorfield(S);
	bgq_weylfield_controlblock *field2 = bgq_translate_spinorfield(R);

	if (parallel)
		return bgq_spinorfield_innerprod_r_global(tri_unknown, field1, field2);
	else
		return bgq_spinorfield_innerprod_r_local(tri_unknown, field1, field2);
}


double square_norm(spinor *P, int N, int parallel) {
	assert(N == LOCAL_VOLUME/2);
	bgq_weylfield_controlblock *field = bgq_translate_spinorfield(P);

	if (parallel)
		return bgq_spinorfield_sqrnorm_global(tri_unknown, field);
	else
		return bgq_spinorfield_sqrnorm_local(tri_unknown, field);
}


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

	bgq_spinorfield_rmul_plain_add_double(targetfield, tri_unknown, targetfield, sourcefield, c);
}

void assign_add_mul_r(spinor * const P, spinor * const Q, const double c, const int N) {
	assert(N == VOLUME/2);
	bgq_weylfield_controlblock *targetfield = bgq_translate_spinorfield(P);
	bgq_weylfield_controlblock *sourcefield = bgq_translate_spinorfield(Q);

	bgq_spinorfield_rmul_plain_add_double(targetfield, tri_unknown, sourcefield, targetfield, c);
}


void gamma5(spinor * const l, spinor * const k, const int V) {
	assert(V == VOLUME/2);
	bgq_weylfield_controlblock *targetfield = bgq_translate_spinorfield(l);
	bgq_weylfield_controlblock *sourcefield = bgq_translate_spinorfield(k);

	bgq_spinorfield_gamma5_double(targetfield, tri_unknown, sourcefield);
}


void assign(spinor * const R, spinor * const S, const int N) {
	assert(N == VOLUME/2);
	bgq_weylfield_controlblock *target = bgq_translate_spinorfield(R);
	bgq_weylfield_controlblock *source = bgq_translate_spinorfield(S);

	bgq_spinorfield_copy(target, ly_full_double, source);
}


void diff(spinor * const Q, const spinor * const R, const spinor * const S, const int N) {
	assert(N == VOLUME/2);
	bgq_weylfield_controlblock *targetfield = bgq_translate_spinorfield(Q);
	bgq_weylfield_controlblock *sourcefield1 = bgq_translate_spinorfield(R);
	bgq_weylfield_controlblock *sourcefield2 = bgq_translate_spinorfield(S);

	bgq_spinorfield_sub_double(targetfield, tri_unknown, sourcefield1, sourcefield2);
}


void zero_spinor_field(spinor * const k, const int N) {
	assert(N == VOLUME/2);
	bgq_weylfield_controlblock *field = bgq_translate_spinorfield(k);
	bgq_spinorfield_zero(field, tri_unknown);
}


void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k) {
	bgq_weylfield_controlblock *targetfield = bgq_translate_spinorfield(l);
	bgq_weylfield_controlblock *sourcefield = bgq_translate_spinorfield(k);

	bgq_HoppingMatrix(ieo, targetfield, sourcefield, 0);
}


void Hopping_Matrix_nocom(const int ieo, spinor * const l, spinor * const k) {
	bgq_weylfield_controlblock *targetfield = bgq_translate_spinorfield(l);
	bgq_weylfield_controlblock *sourcefield = bgq_translate_spinorfield(k);

	bgq_HoppingMatrix(ieo, targetfield, sourcefield, hm_nocom);
}


void tm_times_Hopping_Matrix(const int ieo, spinor * const l, spinor * const k, double complex const cfactor) {
	bgq_weylfield_controlblock *targetfield = bgq_translate_spinorfield(l);
	bgq_weylfield_controlblock *sourcefield = bgq_translate_spinorfield(k);

	bgq_HoppingMatrix(ieo, targetfield, sourcefield, 0);
	bgq_spinorfield_rmul_double(targetfield, ieo, targetfield, cfactor);
}


void tm_sub_Hopping_Matrix(const int ieo, spinor * const l, spinor * p, spinor * const k, complex double const cfactor) {
	bgq_weylfield_controlblock *targetfield = bgq_translate_spinorfield(l);
	bgq_weylfield_controlblock *sourcefield = bgq_translate_spinorfield(k);
	bgq_weylfield_controlblock *sourcefield_sub = bgq_translate_spinorfield(p);

	bgq_HoppingMatrix(ieo, targetfield, sourcefield, 0);
	bgq_spinorfield_icjgmul_plain_sub_double(targetfield, ieo, sourcefield_sub, targetfield, cfactor);
}

#endif

