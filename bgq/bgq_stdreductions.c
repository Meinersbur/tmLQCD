/*
 * bgq_stdreductions.c
 *
 *  Created on: Nov 10, 2012
 *      Author: meinersbur
 */

#include "bgq_qpx.h"
#include "bgq_field.h"


typedef struct {
	bgq_vector4double_decl(sum);
} bgq_vectorsum_t;



typedef struct {
	bgq_vector4double_decl(ks); // sum
	bgq_vector4double_decl(kc); // compensation or carried error (i.e. ks+kc yields better approximation than just ks)
} bgq_vectorkahan_t;

static inline bgq_vectorsum_t bgq_reduce_initzero() {
	bgq_vectorsum_t result;
	bgq_zero(result.sum);
	return result;
}

static inline bgq_vectorkahan_t bgq_reduce_initkahan() {
	bgq_vectorkahan_t result;
	bgq_zero(result.ks);
	bgq_zero(result.kc);
	return result;
}


#define bgq_kahan_add(prev, val) bgq_kahan_add_raw(prev, bgq_vars(val))
static inline bgq_vectorkahan_t bgq_kahan_add_raw(bgq_vectorkahan_t prev, bgq_params(val)) {
	bgq_vector4double_decl(ks);
	bgq_vector4double_decl(kc);
	bgq_mov(ks, prev.ks);
	bgq_mov(kc, prev.kc);

	bgq_vector4double_decl(tr); // val+compensation
	bgq_vector4double_decl(ts); // next sum
	bgq_vector4double_decl(tt); // val+next compensation

	bgq_add(tr, val, kc); // val + carried_error
	bgq_add(ts, tr, ks); // val + carried_error + sum - error
	bgq_sub(tt, ts, ks); // val + carried_error - error
	bgq_mov(ks, ts);
	bgq_sub(kc, tr, tt); // error

	bgq_vectorkahan_t result;
    bgq_mov(result.ks, ks);
    bgq_mov(result.kc, kc);
    return result;
}


static inline void bgq_reduce_prod(bgq_vectorsum_t *sum, bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), ucoord ic) {
	bgq_vector4double_decl(result);

	bgq_mov(result, sum->sum);

	// re = a.re * b.re + a.im * b.im + re
    // im = a.re * b.im - a.im * b.re + im
	bgq_xmadd(result, spinor1_v0_c0, spinor2_v0_c0, result);
    bgq_rxxcpnmadd(result, spinor1_v0_c0, spinor2_v0_c0, result);

    bgq_xmadd(result, spinor1_v0_c1, spinor2_v0_c1, result);
    bgq_rxxcpnmadd(result, spinor1_v0_c1, spinor2_v0_c1, result);

    bgq_xmadd(result, spinor1_v0_c2, spinor2_v0_c2, result);
    bgq_rxxcpnmadd(result, spinor1_v0_c2, spinor2_v0_c2, result);

    bgq_xmadd(result, spinor1_v1_c0, spinor2_v1_c0, result);
    bgq_rxxcpnmadd(result, spinor1_v1_c0, spinor2_v1_c0, result);

    bgq_xmadd(result, spinor1_v1_c1, spinor2_v1_c1, result);
    bgq_rxxcpnmadd(result, spinor1_v1_c1, spinor2_v1_c1, result);

    bgq_xmadd(result, spinor1_v1_c2, spinor2_v1_c2, result);
    bgq_rxxcpnmadd(result, spinor1_v1_c2, spinor2_v1_c2, result);

    bgq_xmadd(result, spinor1_v2_c0, spinor2_v2_c0, result);
    bgq_rxxcpnmadd(result, spinor1_v2_c0, spinor2_v2_c0, result);

    bgq_xmadd(result, spinor1_v2_c1, spinor2_v2_c1, result);
    bgq_rxxcpnmadd(result, spinor1_v2_c1, spinor2_v2_c1, result);

    bgq_xmadd(result, spinor1_v2_c2, spinor2_v2_c2, result);
    bgq_rxxcpnmadd(result, spinor1_v2_c2, spinor2_v2_c2, result);

    bgq_xmadd(result, spinor1_v3_c0, spinor2_v3_c0, result);
    bgq_rxxcpnmadd(result, spinor1_v3_c0, spinor2_v3_c0, result);

    bgq_xmadd(result, spinor1_v3_c1, spinor2_v3_c1, result);
    bgq_rxxcpnmadd(result, spinor1_v3_c1, spinor2_v3_c1, result);

    bgq_xmadd(result, spinor1_v3_c2, spinor2_v3_c2, result);
    bgq_rxxcpnmadd(result, spinor1_v3_c2, spinor2_v3_c2, result);

    bgq_mov(sum->sum, result);
}


static inline void bgq_reduce_prod_r(bgq_vectorsum_t *accumulator, bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), ucoord ic) {
	bgq_vector4double_decl(result);
	bgq_mov(result, accumulator->sum);

	bgq_madd(result, spinor1_v0_c0, spinor2_v0_c0, result);
	bgq_madd(result, spinor1_v0_c1, spinor2_v0_c1, result);
	bgq_madd(result, spinor1_v0_c2, spinor2_v0_c2, result);
	bgq_madd(result, spinor1_v1_c0, spinor2_v1_c0, result);
	bgq_madd(result, spinor1_v1_c1, spinor2_v1_c1, result);
	bgq_madd(result, spinor1_v1_c2, spinor2_v1_c2, result);
	bgq_madd(result, spinor1_v2_c0, spinor2_v2_c0, result);
	bgq_madd(result, spinor1_v2_c1, spinor2_v2_c1, result);
	bgq_madd(result, spinor1_v2_c2, spinor2_v2_c2, result);
	bgq_madd(result, spinor1_v3_c0, spinor2_v3_c0, result);
	bgq_madd(result, spinor1_v3_c1, spinor2_v3_c1, result);
	bgq_madd(result, spinor1_v3_c2, spinor2_v3_c2, result);

	bgq_mov(accumulator->sum, result);
}




static inline void bgq_reduce_norm(bgq_vectorkahan_t *accumulator, bgq_su3_spinor_params(spinor), ucoord ic) {
	bgq_vector4double_decl(siteresult);
	bgq_mul(siteresult, spinor_v0_c0, spinor_v0_c0);
	bgq_madd(siteresult, spinor_v0_c1, spinor_v0_c1, siteresult);
	bgq_madd(siteresult, spinor_v0_c2, spinor_v0_c2, siteresult);
	bgq_madd(siteresult, spinor_v1_c0, spinor_v1_c0, siteresult);
	bgq_madd(siteresult, spinor_v1_c1, spinor_v1_c1, siteresult);
	bgq_madd(siteresult, spinor_v1_c2, spinor_v1_c2, siteresult);
	bgq_madd(siteresult, spinor_v2_c0, spinor_v2_c0, siteresult);
	bgq_madd(siteresult, spinor_v2_c1, spinor_v2_c1, siteresult);
	bgq_madd(siteresult, spinor_v2_c2, spinor_v2_c2, siteresult);
	bgq_madd(siteresult, spinor_v3_c0, spinor_v3_c0, siteresult);
	bgq_madd(siteresult, spinor_v3_c1, spinor_v3_c1, siteresult);
	bgq_madd(siteresult, spinor_v3_c2, spinor_v3_c2, siteresult);

	*accumulator = bgq_kahan_add(*accumulator, siteresult);
}


static inline void bgq_combine_add(bgq_vectorsum_t *arg1, bgq_vectorsum_t arg2) {
	bgq_vector4double_decl(result);
	bgq_mov(result, arg1->sum);

	bgq_add(result, result, arg2.sum);

    bgq_mov(arg1->sum, result);
}


static inline void bgq_combine_kahan(bgq_vectorkahan_t *arg1, bgq_vectorkahan_t arg2) {
	bgq_vectorkahan_t result = *arg1;

#if 1
	result = bgq_kahan_add(result, arg2.kc);
	result = bgq_kahan_add(result, arg2.ks);
#else
	bgq_add(result.ks, arg1->ks, arg2.ks);
	bgq_add(result.kc, arg1->kc, arg2.kc);
#endif

	*arg1 = result;
}


#define REDUCTION_NAME bgq_spinorfield_innerprod_raw
#define REDUCTION_ARGFIELDS 2
#define REDUCTION_RETURNTYPE(...) bgq_vectorsum_t (__VA_ARGS__)
#define REDUCTION_VARINIT bgq_reduce_initzero
#define REDUCTION_SITEREDUCEFUNC bgq_reduce_prod
#define REDUCTION_COMBINEFUNC bgq_combine_add
#include "bgq_reduction.inc.c"

complex_double bgq_spinorfield_innerprod_local(bgq_weylfield_controlblock *field1, bgq_weylfield_controlblock *field2) {
	bgq_vectorsum_t rawresult = bgq_spinorfield_innerprod_raw(field1, field2);
	complex_double localresult = bgq_cmplxval1(rawresult.sum) + bgq_cmplxval2(rawresult.sum);
	return localresult;
}

complex_double bgq_spinorfield_innerprod_global(bgq_weylfield_controlblock *field1, bgq_weylfield_controlblock *field2) {
	complex_double localresult = bgq_spinorfield_innerprod_local(field1, field2);
	complex_double globalresult;
	MPI_Allreduce(&localresult, &globalresult, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
	return globalresult;
}


#define REDUCTION_NAME bgq_innerprod_r_raw
#define REDUCTION_ARGFIELDS 2
#define REDUCTION_RETURNTYPE(...) bgq_vectorsum_t (__VA_ARGS__)
#define REDUCTION_VARINIT bgq_reduce_initzero
#define REDUCTION_SITEREDUCEFUNC bgq_reduce_prod_r
#define REDUCTION_COMBINEFUNC bgq_combine_add
#include "bgq_reduction.inc.c"

double bgq_spinorfield_innerprod_r_local(bgq_weylfield_controlblock *field1, bgq_weylfield_controlblock *field2) {
	bgq_vectorsum_t rawresult = bgq_innerprod_r_raw(field1, field2);
	double localresult = bgq_elem0(rawresult.sum) +  bgq_elem1(rawresult.sum) +  bgq_elem2(rawresult.sum) +  bgq_elem3(rawresult.sum);
	return localresult;
}

double bgq_spinorfield_innerprod_r_global(bgq_weylfield_controlblock *field1, bgq_weylfield_controlblock *field2) {
	double localresult = bgq_spinorfield_innerprod_r_local(field1, field2);
	double globalresult;
	MPI_Allreduce(&localresult, &globalresult, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return globalresult;
}


#define REDUCTION_NAME bgq_spinorfield_sqrnorm_raw
#define REDUCTION_ARGFIELDS 1
#define REDUCTION_RETURNTYPE(...) bgq_vectorkahan_t (__VA_ARGS__)
#define REDUCTION_VARINIT bgq_reduce_initkahan
#define REDUCTION_SITEREDUCEFUNC bgq_reduce_norm
#define REDUCTION_COMBINEFUNC bgq_combine_kahan
#include "bgq_reduction.inc.c"

double bgq_spinorfield_sqrnorm_local(bgq_weylfield_controlblock *field) {
	bgq_vectorkahan_t rawresult = bgq_spinorfield_sqrnorm_raw(field);
	double localresult = bgq_elem0(rawresult.kc) + bgq_elem1(rawresult.kc) + bgq_elem2(rawresult.kc) + bgq_elem3(rawresult.kc); // This is probably useless
	localresult += bgq_elem0(rawresult.ks);
	localresult += bgq_elem1(rawresult.ks);
	localresult += bgq_elem2(rawresult.ks);
	localresult += bgq_elem3(rawresult.ks);
	return localresult;
}

double bgq_spinorfield_sqrnorm_global(bgq_weylfield_controlblock *field) {
	double localresult = bgq_spinorfield_sqrnorm_local(field);
	double globalresult;
	MPI_Allreduce(&localresult, &globalresult, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return globalresult;
}




#if BGQ_REPLACE
_Complex double scalar_prod(spinor *S, spinor *R, int N, int parallel) {
	assert(N == LOCAL_VOLUME/2);
	bgq_weylfield_controlblock *field1 = (bgq_weylfield_controlblock*)S;
	bgq_weylfield_controlblock *field2 = (bgq_weylfield_controlblock*)R;

	if (parallel)
		return bgq_spinorfield_innerprod_global(field1, field2);
	else
		return bgq_spinorfield_innerprod_local(field1, field2);
}


double scalar_prod_r(spinor *S, const spinor * const R, int N, int parallel) {
	assert(N == LOCAL_VOLUME/2);
	bgq_weylfield_controlblock *field1 = (bgq_weylfield_controlblock*)S;
	bgq_weylfield_controlblock *field2 = (bgq_weylfield_controlblock*)R;

	if (parallel)
		return bgq_spinorfield_innerprod_r_global(field1, field2);
	else
		return bgq_spinorfield_innerprod_r_local(field1, field2);
}


double square_norm(spinor *P, int N, int parallel) {
	assert(N == LOCAL_VOLUME/2);
	bgq_weylfield_controlblock *field = (bgq_weylfield_controlblock*)P;

	if (parallel)
		return bgq_spinorfield_sqrnorm_global(field);
	else
		return bgq_spinorfield_sqrnorm_local(field);
}

#endif

