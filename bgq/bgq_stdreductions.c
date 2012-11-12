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


static inline bgq_vectorsum_t bgq_reduce_prod(bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), ucoord tv, ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z) {
	bgq_vector4double_decl(result);

	// re = a.re * b.re + a.im * b.im + re
    // im = a.re * b.im - a.im * b.re + im
	bgq_xmul(result, spinor1_v0_c0, spinor2_v0_c0);
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

    bgq_vectorsum_t vresult;
    bgq_mov(vresult.sum,result);
	return vresult;
}


static inline bgq_vectorsum_t bgq_reduce_prod_r(bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), ucoord tv, ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z) {
	bgq_vector4double_decl(result);

	bgq_mul(result, spinor1_v0_c0, spinor2_v0_c0);
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

    bgq_vectorsum_t vresult;
    bgq_mov(vresult.sum,result);
	return vresult;
}




static inline bgq_vectorsum_t bgq_reduce_norm(bgq_su3_spinor_params(spinor), ucoord tv, ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z) {
	bgq_vector4double_decl(result);

	bgq_mul(result, spinor_v0_c0, spinor_v0_c0);
	bgq_madd(result, spinor_v0_c1, spinor_v0_c1, result);
	bgq_madd(result, spinor_v0_c2, spinor_v0_c2, result);
	bgq_madd(result, spinor_v1_c0, spinor_v1_c0, result);
	bgq_madd(result, spinor_v1_c1, spinor_v1_c1, result);
	bgq_madd(result, spinor_v1_c2, spinor_v1_c2, result);
	bgq_madd(result, spinor_v2_c0, spinor_v2_c0, result);
	bgq_madd(result, spinor_v2_c1, spinor_v2_c1, result);
	bgq_madd(result, spinor_v2_c2, spinor_v2_c2, result);
	bgq_madd(result, spinor_v3_c0, spinor_v3_c0, result);
	bgq_madd(result, spinor_v3_c1, spinor_v3_c1, result);
	bgq_madd(result, spinor_v3_c2, spinor_v3_c2, result);

    bgq_vectorsum_t vresult;
    bgq_mov(vresult.sum,result);
	return vresult;
}


static inline bgq_vectorsum_t bgq_combine_add(bgq_vectorsum_t arg1, bgq_vectorsum_t arg2) {
	bgq_vectorsum_t vresult;
	bgq_add(vresult.sum, arg1.sum, arg2.sum);
	return vresult;
}


#define REDUCTION_NAME bgq_innerprod
#define REDUCTION_ARGFIELDS 2
#define REDUCTION_RETURNTYPE(...) bgq_vectorsum_t (__VA_ARGS__)
#define REDUCTION_SITEREDUCEFUNC bgq_reduce_prod
#define REDUCTION_COMBINEFUNC bgq_combine_add
#include "bgq_reduction.inc.c"

#define REDUCTION_NAME bgq_innerprod_r
#define REDUCTION_ARGFIELDS 2
#define REDUCTION_RETURNTYPE(...) bgq_vectorsum_t (__VA_ARGS__)
#define REDUCTION_SITEREDUCEFUNC bgq_reduce_prod_r
#define REDUCTION_COMBINEFUNC bgq_combine_add
#include "bgq_reduction.inc.c"

#define REDUCTION_NAME bgq_norm
#define REDUCTION_ARGFIELDS 1
#define REDUCTION_RETURNTYPE(...) bgq_vectorsum_t (__VA_ARGS__)
#define REDUCTION_SITEREDUCEFUNC bgq_reduce_norm
#define REDUCTION_COMBINEFUNC bgq_combine_add
#include "bgq_reduction.inc.c"



#if BGQ_REPLACE
_Complex double scalar_prod(spinor *S, spinor *R, int N, int parallel) {
	assert(N == LOCAL_VOLUME/2);
	bgq_weylfield_controlblock *field1 = (bgq_weylfield_controlblock*)S;
	bgq_weylfield_controlblock *field2 = (bgq_weylfield_controlblock*)R;

	bgq_vectorsum_t vresult = bgq_innerprod(field1,field2);
	complexdouble localresult = bgq_cmplxval1(vresult.sum) + bgq_cmplxval2(vresult.sum);

	if (parallel) {
		complexdouble globalresult = 0;
		MPI_Allreduce(&localresult, &localresult, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD );
		return globalresult;
	} else {
		return localresult;
	}
}


 double scalar_prod_r(spinor *S, const spinor * const R, int N, int parallel) {
	assert(N == LOCAL_VOLUME/2);
	bgq_weylfield_controlblock *field1 = (bgq_weylfield_controlblock*)S;
	bgq_weylfield_controlblock *field2 = (bgq_weylfield_controlblock*)R;

	bgq_vectorsum_t vresult = bgq_innerprod_r(field1,field2);
	double localresult = bgq_elem0(vresult.sum) + bgq_elem1(vresult.sum) + bgq_elem2(vresult.sum) + bgq_elem3(vresult.sum);

	if (parallel) {
		double globalresult = 0;
		MPI_Allreduce(&localresult, &localresult, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
		return globalresult;
	} else {
		return localresult;
	}
}



double square_norm(spinor *P, int N, int parallel) {
	assert(N == LOCAL_VOLUME/2);
	bgq_weylfield_controlblock *field = (bgq_weylfield_controlblock*)P;

	bgq_vectorsum_t vresult = bgq_norm(field);
	double localresult = bgq_elem0(vresult.sum) + bgq_elem1(vresult.sum) + bgq_elem2(vresult.sum) + bgq_elem3(vresult.sum);

	if (parallel) {
		double globalresult = 0;
		MPI_Allreduce(&localresult, &localresult, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		return globalresult;
	} else {
		return localresult;
	}
}
#endif

