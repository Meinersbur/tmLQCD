/*
 * bgq.h
 *
 *  Created on: Jul 25, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_H_
#define BGQ_H_

#include <mpi.h>
#include "complex_c99.h"

#ifndef XLC
//typedef double vector4double[4];
typedef struct {
	double q[4];
} vector4double;

#define bgq_vector4double_decl(name) \
	double NAME2(name,q0); \
	double NAME2(name,q1); \
	double NAME2(name,q2); \
	double NAME2(name,q3)

#define bgq_vector4double_decl_leftonly(name) \
	double NAME2(name,q0); \
	double NAME2(name,q1)

#define bgq_vector4double_decl_rightonly(name) \
	double NAME2(name,q2); \
	double NAME2(name,q3)

#define bgq_lda(dst,offset,addr)                                 \
	assert( (((size_t)addr) + offset) % 32 == 0);                \
	NAME2(dst,q0) = ((vector4double*)(((char*)addr) + offset))->q[0]; \
	NAME2(dst,q1) = ((vector4double*)(((char*)addr) + offset))->q[1]; \
	NAME2(dst,q2) = ((vector4double*)(((char*)addr) + offset))->q[2]; \
	NAME2(dst,q3) = ((vector4double*)(((char*)addr) + offset))->q[3]
//TODO: setting the 5 least significant bits to zero

#define bgq_ld2a(dst,offset,addr) \
	assert( (((size_t)addr) + offset) % 16 == 0);                \
	dst##_q0 = ((vector4double*)(((char*)addr) + offset))->q[0]; \
	dst##_q1 = ((vector4double*)(((char*)addr) + offset))->q[1]; \
	dst##_q2 = dst##_q0;                                         \
	dst##_q3 = dst##_q1

#define bgq_ld2a_leftonly(dst,offset,addr) \
	assert( (((size_t)addr) + offset) % 16 == 0);                \
	dst##_q0 = ((vector4double*)(((char*)addr) + offset))->q[0]; \
	dst##_q1 = ((vector4double*)(((char*)addr) + offset))->q[1]

#define bgq_ld2a_rightonly(dst,offset,addr) \
	assert( (((size_t)addr) + offset) % 16 == 0);                \
	dst##_q2 = ((vector4double*)(((char*)addr) + offset))->q[0]; \
	dst##_q3 = ((vector4double*)(((char*)addr) + offset))->q[1]

#define bgq_sta(src,offset,addr) \
	assert( (((size_t)addr) + offset) % 32 == 0);                \
	((vector4double*)(((char*)addr) + offset))->q[0] = src##_q0; \
	((vector4double*)(((char*)addr) + offset))->q[1] = src##_q1; \
	((vector4double*)(((char*)addr) + offset))->q[2] = src##_q2; \
	((vector4double*)(((char*)addr) + offset))->q[3] = src##_q3; \

#define bgq_add(dst,lhs,rhs)        \
	dst##_q0 = lhs##_q0 + rhs##_q0; \
	dst##_q1 = lhs##_q1 + rhs##_q1; \
	dst##_q2 = lhs##_q2 + rhs##_q2; \
	dst##_q3 = lhs##_q3 + rhs##_q3

#define bgq_sub(dst,lhs,rhs)        \
	dst##_q0 = lhs##_q0 - rhs##_q0; \
	dst##_q1 = lhs##_q1 - rhs##_q1; \
	dst##_q2 = lhs##_q2 - rhs##_q2; \
	dst##_q3 = lhs##_q3 - rhs##_q3

#define bgq_xxnpmadd(dst,a,b,c)                     \
	{                                               \
		bgq_vector4double_decl(MAKENAME5(xxnpmadd,dst,a,b,v)); \
		MAKENAME6(xxnpmadd,dst,a,b,v,q0) = - (CONCAT(a,_q1) * CONCAT(b,_q1) - CONCAT(c,_q0)); \
		MAKENAME6(xxnpmadd,dst,a,b,v,q1) =    CONCAT(a,_q0) * CONCAT(b,_q1) + CONCAT(c,_q1) ; \
		MAKENAME6(xxnpmadd,dst,a,b,v,q2) = - (CONCAT(a,_q3) * CONCAT(b,_q3) - CONCAT(c,_q2)); \
		MAKENAME6(xxnpmadd,dst,a,b,v,q3) =    CONCAT(a,_q2) * CONCAT(b,_q3) + CONCAT(c,_q3) ; \
		bgq_mov(dst,MAKENAME5(xxnpmadd,dst,a,b,v));                                           \
	}

#define bgq_iadd(dst,lhs,rhs)         \
	{                                 \
		bgq_vector4double_decl(tmp);  \
		tmp_q0 = lhs##_q0 - rhs##_q1; \
		tmp_q1 = lhs##_q1 + rhs##_q0; \
		tmp_q2 = lhs##_q2 - rhs##_q3; \
		tmp_q3 = lhs##_q3 + rhs##_q2; \
		bgq_mov(dst,tmp);             \
	}

#define bgq_isub(dst,lhs,rhs)         \
	{                                 \
		bgq_vector4double_decl(tmp);  \
		tmp_q0 = lhs##_q0 + rhs##_q1; \
		tmp_q1 = lhs##_q1 - rhs##_q0; \
		tmp_q2 = lhs##_q2 + rhs##_q3; \
		tmp_q3 = lhs##_q3 - rhs##_q2; \
		bgq_mov(dst,tmp);             \
	}

#define bgq_merge2(dst, a23_to01, b01_to23) \
	{                                       \
		bgq_vector4double_decl(MAKENAME4(merge2,dst,a23_to01, b01_to23));      \
		MAKENAME5(merge2,dst,a23_to01,b01_to23,q0) = NAME2(a23_to01,q2);      \
		MAKENAME5(merge2,dst,a23_to01,b01_to23,q1) = NAME2(a23_to01,q3);      \
		MAKENAME5(merge2,dst,a23_to01,b01_to23,q2) = NAME2(b01_to23,q0);      \
		MAKENAME5(merge2,dst,a23_to01,b01_to23,q3) = NAME2(b01_to23,q1);      \
		bgq_mov(dst,MAKENAME4(merge2,dst,a23_to01,b01_to23));                   \
	}

#define bgq_xmul(dst,lhs,rhs)              \
	{                                      \
		bgq_vector4double_decl(MAKENAME4(xmul,dst,lhs,rhs));                \
		MAKENAME5(xmul,dst,lhs,rhs,q0) = CONCAT(lhs,_q0) * CONCAT(rhs,_q0); \
		MAKENAME5(xmul,dst,lhs,rhs,q1) = CONCAT(lhs,_q0) * CONCAT(rhs,_q1); \
		MAKENAME5(xmul,dst,lhs,rhs,q2) = CONCAT(lhs,_q2) * CONCAT(rhs,_q2); \
		MAKENAME5(xmul,dst,lhs,rhs,q3) = CONCAT(lhs,_q2) * CONCAT(rhs,_q3); \
		bgq_mov(dst,MAKENAME4(xmul,dst,lhs,rhs));                           \
	}

#define bgq_xmadd(dst,a,b,c)                      \
	{                                             \
		bgq_vector4double_decl(MAKENAME5(xmul,dst,a,b,c));                            \
		MAKENAME6(xmul,dst,a,b,c,q0) = CONCAT(a,_q0) * CONCAT(b,_q0) + CONCAT(c,_q0); \
		MAKENAME6(xmul,dst,a,b,c,q1) = CONCAT(a,_q0) * CONCAT(b,_q1) + CONCAT(c,_q1); \
		MAKENAME6(xmul,dst,a,b,c,q2) = CONCAT(a,_q2) * CONCAT(b,_q2) + CONCAT(c,_q2); \
		MAKENAME6(xmul,dst,a,b,c,q3) = CONCAT(a,_q2) * CONCAT(b,_q3) + CONCAT(c,_q3); \
		bgq_mov(dst,MAKENAME5(xmul,dst,a,b,c));                                       \
	}

#define bgq_mov(dst,src) \
	NAME2(dst,q0) = NAME2(src,q0); \
	NAME2(dst,q1) = NAME2(src,q1); \
	NAME2(dst,q2) = NAME2(src,q2); \
	NAME2(dst,q3) = NAME2(src,q3)

#define bgq_cconst(dst,re,im) \
	dst##_q0 = re;            \
	dst##_q1 = im;            \
	dst##_q2 = re;            \
	dst##_q3 = im

#define bgq_xxcpnmadd(dst,a,b,c)                      \
	{                                                 \
		bgq_vector4double_decl(MAKENAME5(xmul,dst,a,b,c));   \
		MAKENAME6(xmul,dst,a,b,c,q0) =    CONCAT(a,_q1) * CONCAT(b,_q1) + CONCAT(c,_q0) ; \
		MAKENAME6(xmul,dst,a,b,c,q1) = - (CONCAT(a,_q0) * CONCAT(b,_q1) - CONCAT(c,_q1)); \
		MAKENAME6(xmul,dst,a,b,c,q2) =    CONCAT(a,_q3) * CONCAT(b,_q3) + CONCAT(c,_q2) ; \
		MAKENAME6(xmul,dst,a,b,c,q3) = - (CONCAT(a,_q2) * CONCAT(b,_q3) - CONCAT(c,_q3)); \
		bgq_mov(dst,MAKENAME5(xmul,dst,a,b,c));                                           \
	}

#else

#define bgq_vector4double_decl(name) \
	vector4double name

#define bgq_vector4double_decl_leftonly(name) \
	vector4double name

#define bgq_vector4double_decl_rightonly(name) \
	vector4double name

#define bgq_lda(dst,offset,addr) \
	dst = vec_lda(offset,addr)

#define bgq_ld2a(dst,offset,addr) \
	dst = vec_ld2a(offset, addr)

#define bgq_ld2a_leftonly(dst,offset,addr) \
	bgq_ld2a(dst,offset,addr)

#define bgq_ld2a_rightonly(dst,offset,addr) \
	bgq_ld2a(dst,offset,addr)

#define bgq_sta(src,offset,addr) \
	vec_sta(src,offset,addr)

#define bgq_add(dst,lhs,rhs) \
	dst = vec_add(lhs, rhs)

#define bgq_sub(dst,lhs,rhs) \
	dst = vec_sub(lhs, rhs)

#define bgq_xxnpmadd(dst,a,b,c) \
	dst = vec_xxnpmadd(a,b,c)

#define bgq_iadd(dst,lhs,rhs) \
	bgq_xxnpmadd(dst,(vector4double)(1),rhs,lhs)

#define bgq_isub(dst,lhs,rhs) \
	bgq_xxcpnmadd(dst,rhs,(vector4double)(1),lhs)

#define bgq_merge2(dst, a23_to01, b01_to23) \
	dst = vec_perm(a23_to01, b01_to23, vec_gpci(02345))

#define bgq_xmul(dst,lhs,rhs) \
	dst = vec_xmul(lhs,rhs)

#define bgq_xmadd(dst,a,b,c) \
	dst = vec_xmadd(a,b,c)

#define bgq_mov(dst,src) \
	dst = src

#define bgq_const(dst,re,im) \
	dst = (vector4double){re,im,re,im}

#define bgq_xxcpnmadd(dst,a,b,c) \
	dst = vec_xxcpnmadd(a,b,c)

#define bgq_cconst(dst,re,im) \
	dst = (vector4double){re,im,re,im}

#endif

// vec_xmul(a, b)
// re =    a.re * b.re
// im =    a.re * b.im

// vec_xmul(b, a)
// re =    a.re * b.re
// im =    a.im * b.re

// vec_xmadd(a, b, c)
// re =    a.re * b.re + c.re
// im =    a.re * b.im + c.im

// vec_xmadd(b, a, c)
// re =    a.re * b.re + c.re
// im =    a.im * b.re + c.im

// vec_xxnpmadd(a, b, c)
// re =   -a.im * b.im + c.re
// im =    a.re * b.im + c.im

// vec_xxnpmadd(b, a, c)
// re =   -a.im * b.im + c.re
// im =    a.im * b.re + c.im

// vec_xxcpnmadd(a, b, c)
// re =    a.im * b.im + c.re
// im =   -a.re * b.im + c.im

// vec_xxcpnmadd(b, a, c)
// re =    a.im * b.im + c.re
// im =   -a.im * b.re + c.im

// vec_xxmadd(a, b, c)
// re =    a.im * b.im + c.re
// im =    a.re * b.im + c.im

// vec_xxmadd(b, a, c)
// re =    a.im * b.im + c.re
// im =    a.im * b.re + c.im

#define cvec_mul(a,b) vec_xxnpmadd(b,a,vec_xmul(a,b))
// vec_xxnpmadd(b,a,vec_xmul(a,b))
// vec_xxnpmadd(a,b,vec_xmul(b,a))
// re = a.re * b.re - a.im * b.im
// im = a.re * b.im + a.im * b.re
//    = a * b

#define bgq_cmul(dst,lhs,rhs)                                                                       \
	{                                                                                               \
		bgq_vector4double_decl(MAKENAME4(cmul,dst,lhs,rhs));                                         \
		bgq_xmul              (MAKENAME4(cmul,dst,lhs,rhs), lhs, rhs);                              \
		bgq_xxnpmadd          (dst                        , rhs, lhs, MAKENAME4(cmul,dst,lhs,rhs)); \
	}

#define cvec_madd(a,b,c) vec_xxnpmadd(b,a,vec_xmadd(a,b,c))
// vec_xxnpmadd(b,a,vec_xmadd(a,b,c))
// vec_xxnpmadd(a,b,vec_xmadd(b,a,c))
// re = a.re * b.re - a.im * b.im + c.re
// im = a.re * b.im + a.im * b.re + c.im
//    = a * b + c

#define bgq_cmadd(dst,a,b,c)                                                                                  \
	{                                                                                                         \
		bgq_vector4double_decl(MAKENAME5(cmadd,dst,a,b,c));                                           \
		bgq_xmadd             (MAKENAME5(cmadd,dst,a,b,c), a, b, c);                                  \
		bgq_xxnpmadd          (dst                       , b, a, MAKENAME5(cmadd,dst,a,b,c)); \
	}

// conjugated second argument
#define cvec_pcmul(a,b) vec_xxcpnmadd(a,b,vec_xmul(b,a))
// vec_xxcpnmadd(a,b,vec_xmul(b,a))
// re = a.re * b.re + a.im * b.im
// im = a.im * b.re - a.re * b.im
//    = a * conjugate(b)

#define bgq_ccmul(dst,lhs,rhs) \
	{                                  \
		bgq_vector4double_decl(tmq);   \
		bgq_xmul(tmq,rhs,lhs);         \
		bgq_xxnpmadd(dst,lhs,rhs,tmq); \
	}

// conjugated second argument
#define cvec_pcpmadd(a,b,c) vec_xxcpnmadd(a,b,vec_xmadd(b,a,c))
// vec_xxcpnmadd(a,b,vec_xmadd(b,a,c))
// vec_xmadd(b,a,vec_xxcpnmadd(a,b,c))
// re = a.re * b.re + a.im * b.im + c.re
// im = a.im * b.re - a.re * b.im + c.im
//    = a * conjugate(b) + c

#define bgq_ccmadd(dst,a,b,c)        \
	{                                \
		bgq_vector4double_decl(tmq); \
		bgq_xmadd(tmq,b,a,c);        \
		bgq_xxcpnmadd(dst,a,b,tmq);  \
	}

#define cvec_piadd(a,b) vec_xxnpmadd((vector4double)(1),b,a)
// vec_xxnpmadd((vector4double)(1),b,a)
// re = a.re - b.im
// im = a.im + b.re
//    = a + bi

#define cvec_pisub(a,b) vec_xxcpnmadd(b,(vector4double)(1),a)
// vec_xxcpnmadd(b,(vector4double)(1),a)
// re = a.re + b.im
// im = a.im - b.re
//    = a - bi

#define cvec_merge(a23_to01, b01_to23) vec_perm(a23_to01, b01_to23, vec_gpci(02345))
// [0] = a[2] = a.re
// [1] = a[3] = a.im
// [2] = b[0] = b.re
// [3] = b[1] = b.im

#define bgq_su3_vdecl(name)            \
	bgq_vector4double_decl(CONCAT(name,_c0)); \
	bgq_vector4double_decl(CONCAT(name,_c1)); \
	bgq_vector4double_decl(CONCAT(name,_c2))

#define bgq_su3_vdecl_leftonly(name)            \
	bgq_vector4double_decl_leftonly(CONCAT(name,_c0)); \
	bgq_vector4double_decl_leftonly(CONCAT(name,_c1)); \
	bgq_vector4double_decl_leftonly(CONCAT(name,_c2))

#define bgq_su3_vdecl_rightonly(name)            \
	bgq_vector4double_decl_rightonly(CONCAT(name,_c0)); \
	bgq_vector4double_decl_rightonly(CONCAT(name,_c1)); \
	bgq_vector4double_decl_rightonly(CONCAT(name,_c2))

#define bgq_su3_mdecl(name)             \
	bgq_vector4double_decl(name##_c00); \
	bgq_vector4double_decl(name##_c01); \
	bgq_vector4double_decl(name##_c02); \
	bgq_vector4double_decl(name##_c10); \
	bgq_vector4double_decl(name##_c11); \
	bgq_vector4double_decl(name##_c12); \
	bgq_vector4double_decl(name##_c20); \
	bgq_vector4double_decl(name##_c21); \
	bgq_vector4double_decl(name##_c22)

#define bgq_su3_spinor_decl(name) \
	bgq_su3_vdecl(name##_v0);     \
	bgq_su3_vdecl(name##_v1);     \
	bgq_su3_vdecl(name##_v2);     \
	bgq_su3_vdecl(name##_v3)

#define bgq_su3_spinor_decl_leftonly(name) \
	bgq_su3_vdecl_leftonly(name##_v0);     \
	bgq_su3_vdecl_leftonly(name##_v1);     \
	bgq_su3_vdecl_leftonly(name##_v2);     \
	bgq_su3_vdecl_leftonly(name##_v3)

#define bgq_su3_spinor_decl_rightonly(name) \
	bgq_su3_vdecl_rightonly(name##_v0);     \
	bgq_su3_vdecl_rightonly(name##_v1);     \
	bgq_su3_vdecl_rightonly(name##_v2);     \
	bgq_su3_vdecl_rightonly(name##_v3)

#define bgq_su3_weyl_decl(name) \
	bgq_su3_vdecl(name##_v0);	\
	bgq_su3_vdecl(name##_v1)

#define bgq_su3_spinor_double_load(dst, addr) \
	bgq_lda(NAME3(dst,v0,c0),   0, (double _Complex*)(addr));          \
	bgq_lda(NAME3(dst,v0,c1),  32, (double _Complex*)(addr));          \
	bgq_lda(NAME3(dst,v0,c2),  64, (double _Complex*)(addr));          \
	bgq_lda(NAME3(dst,v1,c0),  96, (double _Complex*)(addr));          \
	bgq_lda(NAME3(dst,v1,c1), 128, (double _Complex*)(addr));          \
	bgq_lda(NAME3(dst,v1,c2), 160, (double _Complex*)(addr));          \
	bgq_lda(NAME3(dst,v2,c0), 192, (double _Complex*)(addr));          \
	bgq_lda(NAME3(dst,v2,c1), 224, (double _Complex*)(addr));          \
	bgq_lda(NAME3(dst,v2,c2), 256, (double _Complex*)(addr));          \
	bgq_lda(NAME3(dst,v3,c0), 288, (double _Complex*)(addr));          \
	bgq_lda(NAME3(dst,v3,c1), 320, (double _Complex*)(addr));          \
	bgq_lda(NAME3(dst,v3,c2), 352, (double _Complex*)(addr))
// NOTE: qvlfdux is possibly more effective, but no compiler built-in exists

#define bgq_su3_weyl_double_load(dest, addr) \
	bgq_lda(dest##_v0_c0,   0, (double _Complex*)(addr));        \
	bgq_lda(dest##_v0_c1,  32, (double _Complex*)(addr));        \
	bgq_lda(dest##_v0_c2,  64, (double _Complex*)(addr));        \
	bgq_lda(dest##_v1_c0,  96, (double _Complex*)(addr));        \
	bgq_lda(dest##_v1_c1, 128, (double _Complex*)(addr));        \
	bgq_lda(dest##_v1_c2, 160, (double _Complex*)(addr))

#define bgq_su3_matrix_double_load(dest, addr)   \
	bgq_lda(dest##_c00,   0, (double _Complex*)(addr)); \
	bgq_lda(dest##_c01,  32, (double _Complex*)(addr)); \
	bgq_lda(dest##_c02,  64, (double _Complex*)(addr)); \
	bgq_lda(dest##_c10,  96, (double _Complex*)(addr)); \
	bgq_lda(dest##_c11, 128, (double _Complex*)(addr)); \
	bgq_lda(dest##_c12, 160, (double _Complex*)(addr)); \
	bgq_lda(dest##_c20, 192, (double _Complex*)(addr)); \
	bgq_lda(dest##_c21, 224, (double _Complex*)(addr)); \
	bgq_lda(dest##_c22, 256, (double _Complex*)(addr))

#define bgq_su3_matrix_double_load_left(dest, addr)   \
	dest##_c00 = vec_ld2a(  0, (double _Complex*)(addr)); \
	dest##_c01 = vec_ld2a( 32, (double _Complex*)(addr)); \
	dest##_c02 = vec_ld2a( 64, (double _Complex*)(addr)); \
	dest##_c10 = vec_ld2a( 96, (double _Complex*)(addr)); \
	dest##_c11 = vec_ld2a(128, (double _Complex*)(addr)); \
	dest##_c12 = vec_ld2a(160, (double _Complex*)(addr)); \
	dest##_c20 = vec_ld2a(192, (double _Complex*)(addr)); \
	dest##_c21 = vec_ld2a(224, (double _Complex*)(addr)); \
	dest##_c22 = vec_ld2a(256, (double _Complex*)(addr))

#define bgq_su3_matrix_double_load_right(dest, addr)   \
	dest##_c00 = vec_ld2a(16+  0, (double _Complex*)(addr)); \
	dest##_c01 = vec_ld2a(16+ 32, (double _Complex*)(addr)); \
	dest##_c02 = vec_ld2a(16+ 64, (double _Complex*)(addr)); \
	dest##_c10 = vec_ld2a(16+ 96, (double _Complex*)(addr)); \
	dest##_c11 = vec_ld2a(16+128, (double _Complex*)(addr)); \
	dest##_c12 = vec_ld2a(16+160, (double _Complex*)(addr)); \
	dest##_c20 = vec_ld2a(16+192, (double _Complex*)(addr)); \
	dest##_c21 = vec_ld2a(16+224, (double _Complex*)(addr)); \
	dest##_c22 = vec_ld2a(16+256, (double _Complex*)(addr))

#define bgq_su3_spinor_double_load_left(dst,addr) \
	bgq_ld2a(dst##_v0_c0,   0, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v0_c1,  32, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v0_c2,  64, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v1_c0,  96, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v1_c1, 128, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v1_c2, 160, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v2_c0, 192, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v2_c1, 224, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v2_c2, 256, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v3_c0, 288, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v3_c1, 320, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v3_c2, 352, (double _Complex*)(addr))

#define bgq_su3_spinor_double_load_left_toleftonly(dst,addr) \
	bgq_ld2a_leftonly(dst##_v0_c0,   0, (double*)(addr)); \
	bgq_ld2a_leftonly(dst##_v0_c1,  32, (double*)(addr)); \
	bgq_ld2a_leftonly(dst##_v0_c2,  64, (double*)(addr)); \
	bgq_ld2a_leftonly(dst##_v1_c0,  96, (double*)(addr)); \
	bgq_ld2a_leftonly(dst##_v1_c1, 128, (double*)(addr)); \
	bgq_ld2a_leftonly(dst##_v1_c2, 160, (double*)(addr)); \
	bgq_ld2a_leftonly(dst##_v2_c0, 192, (double*)(addr)); \
	bgq_ld2a_leftonly(dst##_v2_c1, 224, (double*)(addr)); \
	bgq_ld2a_leftonly(dst##_v2_c2, 256, (double*)(addr)); \
	bgq_ld2a_leftonly(dst##_v3_c0, 288, (double*)(addr)); \
	bgq_ld2a_leftonly(dst##_v3_c1, 320, (double*)(addr)); \
	bgq_ld2a_leftonly(dst##_v3_c2, 352, (double*)(addr))

#define bgq_su3_spinor_double_load_right(dst,addr) \
	bgq_ld2a(dst##_v0_c0,   0+16, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v0_c1,  32+16, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v0_c2,  64+16, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v1_c0,  96+16, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v1_c1, 128+16, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v1_c2, 160+16, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v2_c0, 192+16, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v2_c1, 224+16, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v2_c2, 256+16, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v3_c0, 288+16, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v3_c1, 320+16, (double _Complex*)(addr)); \
	bgq_ld2a(dst##_v3_c2, 352+16, (double _Complex*)(addr))

#define bgq_su3_spinor_double_load_right_torightonly(dst,addr) \
	bgq_ld2a_rightonly(dst##_v0_c0,   0+16, (double*)(addr)); \
	bgq_ld2a_rightonly(dst##_v0_c1,  32+16, (double*)(addr)); \
	bgq_ld2a_rightonly(dst##_v0_c2,  64+16, (double*)(addr)); \
	bgq_ld2a_rightonly(dst##_v1_c0,  96+16, (double*)(addr)); \
	bgq_ld2a_rightonly(dst##_v1_c1, 128+16, (double*)(addr)); \
	bgq_ld2a_rightonly(dst##_v1_c2, 160+16, (double*)(addr)); \
	bgq_ld2a_rightonly(dst##_v2_c0, 192+16, (double*)(addr)); \
	bgq_ld2a_rightonly(dst##_v2_c1, 224+16, (double*)(addr)); \
	bgq_ld2a_rightonly(dst##_v2_c2, 256+16, (double*)(addr)); \
	bgq_ld2a_rightonly(dst##_v3_c0, 288+16, (double*)(addr)); \
	bgq_ld2a_rightonly(dst##_v3_c1, 320+16, (double*)(addr)); \
	bgq_ld2a_rightonly(dst##_v3_c2, 352+16, (double*)(addr))

#define bgq_su3_spinor_double_store(addr,src) \
	bgq_sta(src##_v0_c0,   0, (double _Complex*)(addr));          \
	bgq_sta(src##_v0_c1,  32, (double _Complex*)(addr));          \
	bgq_sta(src##_v0_c2,  64, (double _Complex*)(addr));          \
	bgq_sta(src##_v1_c0,  96, (double _Complex*)(addr));          \
	bgq_sta(src##_v1_c1, 128, (double _Complex*)(addr));          \
	bgq_sta(src##_v1_c2, 160, (double _Complex*)(addr));          \
	bgq_sta(src##_v2_c0, 192, (double _Complex*)(addr));          \
	bgq_sta(src##_v2_c1, 224, (double _Complex*)(addr));          \
	bgq_sta(src##_v2_c2, 256, (double _Complex*)(addr));          \
	bgq_sta(src##_v3_c0, 288, (double _Complex*)(addr));          \
	bgq_sta(src##_v3_c1, 320, (double _Complex*)(addr));          \
	bgq_sta(src##_v3_c2, 352, (double _Complex*)(addr))

#define bgq_su3_weyl_double_store(addr,src) \
	bgq_sta(src##_v0_c0,   0, (double _Complex*)(addr));          \
	bgq_sta(src##_v0_c1,  32, (double _Complex*)(addr));          \
	bgq_sta(src##_v0_c2,  64, (double _Complex*)(addr));          \
	bgq_sta(src##_v1_c0,  96, (double _Complex*)(addr));          \
	bgq_sta(src##_v1_c1, 128, (double _Complex*)(addr));          \
	bgq_sta(src##_v1_c2, 160, (double _Complex*)(addr))

#define bgq_su3_spinor_merge(dst,a,b)         \
	bgq_su3_vmerge(dst##_v0, a##_v0, b##_v0); \
	bgq_su3_vmerge(dst##_v1, a##_v1, b##_v1); \
	bgq_su3_vmerge(dst##_v2, a##_v2, b##_v2); \
	bgq_su3_vmerge(dst##_v3, a##_v3, b##_v3)

#define bgq_su3_vmov(dst,src)    \
	bgq_mov(NAME2(dst,c0), NAME2(src,c0)); \
	bgq_mov(NAME2(dst,c1), NAME2(src,c1)); \
	bgq_mov(NAME2(dst,c2), NAME2(src,c2))

#define bgq_su3_vmerge(dst,a,b)           \
	bgq_merge2(dst##_c0, a##_c0, b##_c0); \
	bgq_merge2(dst##_c1, a##_c1, b##_c1); \
	bgq_merge2(dst##_c2, a##_c2, b##_c2)

#define bgq_su3_mmerge(dst,a,b)               \
	dst##_c00 = cvec_merge(a##_c00, b##_c00); \
	dst##_c01 = cvec_merge(a##_c01, b##_c01); \
	dst##_c02 = cvec_merge(a##_c02, b##_c02); \
	dst##_c10 = cvec_merge(a##_c10, b##_c10); \
	dst##_c11 = cvec_merge(a##_c11, b##_c11); \
	dst##_c12 = cvec_merge(a##_c12, b##_c12); \
	dst##_c20 = cvec_merge(a##_c20, b##_c20); \
	dst##_c21 = cvec_merge(a##_c21, b##_c21); \
	dst##_c22 = cvec_merge(a##_c22, b##_c22); \

#define bgq_su3_vadd(dst,v1,v2)          \
	bgq_add(dst##_c0, v1##_c0, v2##_c0); \
	bgq_add(dst##_c1, v1##_c1, v2##_c1); \
	bgq_add(dst##_c2, v1##_c2, v2##_c2)

#define bgq_su3_vsub(dst,v1,v2)         \
	bgq_sub(dst##_c0, v1##_c0, v2##_c0); \
	bgq_sub(dst##_c1, v1##_c1, v2##_c1); \
	bgq_sub(dst##_c2, v1##_c2, v2##_c2)

#define bgq_su3_vpiadd(dst,v1,v2)         \
	bgq_iadd(dst##_c0, v1##_c0, v2##_c0); \
	bgq_iadd(dst##_c1, v1##_c1, v2##_c1); \
	bgq_iadd(dst##_c2, v1##_c2, v2##_c2)

#define bgq_su3_vpisub(dst,v1,v2)         \
	bgq_isub(dst##_c0, v1##_c0, v2##_c0); \
	bgq_isub(dst##_c1, v1##_c1, v2##_c1); \
	bgq_isub(dst##_c2, v1##_c2, v2##_c2)

#define bgq_su3_cvmul(dst,c,v)     \
	bgq_cmul(CONCAT(dst,_c0), c, CONCAT(v,_c0)); \
	bgq_cmul(CONCAT(dst,_c1), c, CONCAT(v,_c1)); \
	bgq_cmul(CONCAT(dst,_c2), c, CONCAT(v,_c2))

#define bgq_su3_mvmul(dst,m,v)                                                      \
	{                                                                               \
		bgq_su3_vdecl(MAKENAME4(mvmul,dst,m,v));                           \
		bgq_cmul (MAKENAME5(mvmul,dst,m,v,c0), v##_c0, m##_c00);         \
		bgq_cmadd(MAKENAME5(mvmul,dst,m,v,c0), v##_c1, m##_c01, MAKENAME5(mvmul,dst,m,v,c0)); \
		bgq_cmadd(MAKENAME5(mvmul,dst,m,v,c0), v##_c2, m##_c02, MAKENAME5(mvmul,dst,m,v,c0)); \
		bgq_cmul (MAKENAME5(mvmul,dst,m,v,c1), v##_c0, m##_c10);         \
		bgq_cmadd(MAKENAME5(mvmul,dst,m,v,c1), v##_c1, m##_c11, MAKENAME5(mvmul,dst,m,v,c1)); \
		bgq_cmadd(MAKENAME5(mvmul,dst,m,v,c1), v##_c2, m##_c12, MAKENAME5(mvmul,dst,m,v,c1)); \
		bgq_cmul (MAKENAME5(mvmul,dst,m,v,c2), v##_c0, m##_c20);         \
		bgq_cmadd(MAKENAME5(mvmul,dst,m,v,c2), v##_c1, m##_c21, MAKENAME5(mvmul,dst,m,v,c2)); \
		bgq_cmadd(MAKENAME5(mvmul,dst,m,v,c2), v##_c2, m##_c22, MAKENAME5(mvmul,dst,m,v,c2)); \
		bgq_su3_vmov(dst, MAKENAME4(mvmul,dst,m,v));                       \
	}

#define bgq_su3_mvinvmul(dst,m,v)                    \
	{                                                \
		bgq_su3_vdecl(tmp);                          \
		bgq_ccmul (tmp_c0, v##_c0, m##_c00);         \
		bgq_ccmadd(tmp_c0, v##_c1, m##_c10, tmp_c0); \
		bgq_ccmadd(tmp_c0, v##_c2, m##_c20, tmp_c0); \
		bgq_ccmul (tmp_c1, v##_c0, m##_c01);         \
		bgq_ccmadd(tmp_c1, v##_c1, m##_c11, tmp_c1); \
		bgq_ccmadd(tmp_c1, v##_c2, m##_c21, tmp_c1); \
		bgq_ccmul (tmp_c2, v##_c0, m##_c02);         \
		bgq_ccmadd(tmp_c2, v##_c1, m##_c12, tmp_c2); \
		bgq_ccmadd(tmp_c2, v##_c2, m##_c22, tmp_c2); \
		bgq_su3_vmov(dst, tmp);                      \
	}

#define bgq_su3_spinor_mov(dst,src)  \
	bgq_su3_vmov(dst##_v0,src##_v0); \
	bgq_su3_vmov(dst##_v1,src##_v1); \
	bgq_su3_vmov(dst##_v2,src##_v2); \
	bgq_su3_vmov(dst##_v3,src##_v3)

#define bgq_su3_weyl_mov(dst,src)  \
	bgq_su3_vmov(dst##_v0,src##_v0); \
	bgq_su3_vmov(dst##_v1,src##_v1)


#define master_print(...)           \
	if (g_proc_id == 0)              \
		fprintf(stderr, __VA_ARGS__)

#define master_error(errcode, ...) \
	do {                            \
		master_print(__VA_ARGS__);  \
		exit(errcode);              \
	} while (0)

#endif /* BGQ_H_ */

