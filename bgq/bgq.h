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
#include <string.h>

#ifndef BGQ_QPX
#define BGQ_QPX 0
#endif


typedef struct {
	double q[4];
} v4d;

typedef struct {
	float q[4];
} v4f;

typedef struct {
	_Complex double q[2];
} v2cd;

typedef struct {
	_Complex float q[2];
} v2cf;

#ifdef XLC
	#define BGQ_VECTOR4DOUBLE_SUBSCRIPT(addr,idx) ( (*((vector4double*)(addr)))[(idx)] ) /* allows subscript like an array */
#else
	typedef v4d vector4double;
	#define BGQ_VECTOR4DOUBLE_SUBSCRIPT(addr,idx) ( ((v4d*)(addr))->q[(idx)] ) /* emulation using a struct */
#endif
#define BGQ_VECTOR4FLOAT_SUBSCRIPT(addr,idx) ( ((v4f*)(addr))->q[(idx)] ) /* emulation using a struct */




#if !BGQ_QPX
#define bgq_elem0(arg) \
	NAME2(arg,q0)

#define bgq_elem1(arg) \
	NAME2(arg,q1)

#define bgq_elem2(arg) \
	NAME2(arg,q2)

#define bgq_elem3(arg) \
	NAME2(arg,q3)

#define bgq_vars(name) \
	NAME2(name,q0),NAME2(name,q1),NAME2(name,q2),NAME2(name,q3)

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

#define bgq_lda_double(dst,offset,addr)                                 \
	assert( (((size_t)addr) + offset) % 32 == 0);                \
	NAME2(dst,q0) = BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 0); \
	NAME2(dst,q1) = BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 1); \
	NAME2(dst,q2) = BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 2); \
	NAME2(dst,q3) = BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 3)
//TODO: setting the 5 least significant bits of addr+offset to zero

#define bgq_lda_float(dst,offset,addr)                                 \
	assert( (((size_t)addr) + offset) % 16 == 0);                \
	NAME2(dst,q0) = BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 0); \
	NAME2(dst,q1) = BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 1); \
	NAME2(dst,q2) = BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 2); \
	NAME2(dst,q3) = BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 3)

#define bgq_ld2a_double(dst,offset,addr)                                      \
	assert( (((size_t)(addr)) + (offset)) % 16 == 0);                          \
	NAME2(dst,q0) = BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 0); \
	NAME2(dst,q1) = BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 1); \
	NAME2(dst,q2) = NAME2(dst,q0);                                             \
	NAME2(dst,q3) = NAME2(dst,q1)

#define bgq_ld2a_float(dst,offset,addr) \
	assert( (((size_t)(addr)) + (offset)) % 8 == 0);                         \
	NAME2(dst,q0) = BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 0); \
	NAME2(dst,q1) = BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 1); \
	NAME2(dst,q2) = NAME2(dst,q0);                                           \
	NAME2(dst,q3) = NAME2(dst,q1)

#define bgq_sta_double(src,offset,addr) \
	assert( (((size_t)(addr)) + (offset)) % 32 == 0);                         \
	BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 0) = NAME2(src,q0); \
	BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 1) = NAME2(src,q1); \
	BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 2) = NAME2(src,q2); \
	BGQ_VECTOR4DOUBLE_SUBSCRIPT((char*)(addr) + (offset), 3) = NAME2(src,q3)

#define bgq_sta_float(src,offset,addr) \
	assert( (((size_t)(addr)) + (offset)) % 16 == 0);                        \
	BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 0) = NAME2(src,q0); \
	BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 1) = NAME2(src,q1); \
	BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 2) = NAME2(src,q2); \
	BGQ_VECTOR4FLOAT_SUBSCRIPT((char*)(addr) + (offset), 3) = NAME2(src,q3)

#define bgq_neg(dst,arg)            \
	NAME2(dst,q0) = - NAME2(arg,q0); \
	NAME2(dst,q1) = - NAME2(arg,q1); \
	NAME2(dst,q2) = - NAME2(arg,q2); \
	NAME2(dst,q3) = - NAME2(arg,q3)

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

#define bgq_mul(dst,lhs,rhs)                                          \
	{                                                                  \
		bgq_vector4double_decl(MAKENAME4(mul,dst,lhs,rhs));            \
		MAKENAME5(mul,dst,lhs,rhs,q0) = NAME2(lhs,q0) * NAME2(rhs,q0); \
		MAKENAME5(mul,dst,lhs,rhs,q1) = NAME2(lhs,q1) * NAME2(rhs,q1); \
		MAKENAME5(mul,dst,lhs,rhs,q2) = NAME2(lhs,q2) * NAME2(rhs,q2); \
		MAKENAME5(mul,dst,lhs,rhs,q3) = NAME2(lhs,q3) * NAME2(rhs,q3); \
		bgq_mov(dst,MAKENAME4(mul,dst,lhs,rhs));                       \
	}

#define bgq_merge2(dst,a23_to01,b01_to23) \
	{                                       \
		bgq_vector4double_decl(MAKENAME4(merge2,dst,a23_to01,b01_to23));      \
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
		bgq_vector4double_decl(MAKENAME5(xmadd,dst,a,b,c));                            \
		MAKENAME6(xmadd,dst,a,b,c,q0) = CONCAT(a,_q0) * CONCAT(b,_q0) + CONCAT(c,_q0); \
		MAKENAME6(xmadd,dst,a,b,c,q1) = CONCAT(a,_q0) * CONCAT(b,_q1) + CONCAT(c,_q1); \
		MAKENAME6(xmadd,dst,a,b,c,q2) = CONCAT(a,_q2) * CONCAT(b,_q2) + CONCAT(c,_q2); \
		MAKENAME6(xmadd,dst,a,b,c,q3) = CONCAT(a,_q2) * CONCAT(b,_q3) + CONCAT(c,_q3); \
		bgq_mov(dst,MAKENAME5(xmadd,dst,a,b,c));                                       \
	}

#define bgq_madd(dst,a,b,c)                                                    \
	{                                                                           \
		bgq_vector4double_decl(MAKENAME5(madd,dst,a,b,c));                      \
		MAKENAME6(madd,dst,a,b,c,q0) = NAME2(a,q0) * NAME2(b,q0) + NAME2(c,q0); \
		MAKENAME6(madd,dst,a,b,c,q1) = NAME2(a,q1) * NAME2(b,q1) + NAME2(c,q1); \
		MAKENAME6(madd,dst,a,b,c,q2) = NAME2(a,q2) * NAME2(b,q2) + NAME2(c,q2); \
		MAKENAME6(madd,dst,a,b,c,q3) = NAME2(a,q3) * NAME2(b,q3) + NAME2(c,q3); \
		bgq_mov(dst,MAKENAME5(madd,dst,a,b,c));                                 \
	}

#define bgq_mov(dst,src) \
	NAME2(dst,q0) = NAME2(src,q0); \
	NAME2(dst,q1) = NAME2(src,q1); \
	NAME2(dst,q2) = NAME2(src,q2); \
	NAME2(dst,q3) = NAME2(src,q3)

#define bgq_xxcpnmadd(dst,a,b,c)                      \
	{                                                 \
		bgq_vector4double_decl(MAKENAME5(xmul,dst,a,b,c));   \
		MAKENAME6(xmul,dst,a,b,c,q0) =    CONCAT(a,_q1) * CONCAT(b,_q1) + CONCAT(c,_q0) ; \
		MAKENAME6(xmul,dst,a,b,c,q1) = - (CONCAT(a,_q0) * CONCAT(b,_q1) - CONCAT(c,_q1)); \
		MAKENAME6(xmul,dst,a,b,c,q2) =    CONCAT(a,_q3) * CONCAT(b,_q3) + CONCAT(c,_q2) ; \
		MAKENAME6(xmul,dst,a,b,c,q3) = - (CONCAT(a,_q2) * CONCAT(b,_q3) - CONCAT(c,_q3)); \
		bgq_mov(dst,MAKENAME5(xmul,dst,a,b,c));                                           \
	}

#define bgq_cconst(dst,re,im) \
	dst##_q0 = re;            \
	dst##_q1 = im;            \
	dst##_q2 = re;            \
	dst##_q3 = im

#define bgq_zero(dst) \
	NAME2(dst,q0) = 0;     \
	NAME2(dst,q1) = 0;     \
	NAME2(dst,q2) = 0;     \
	NAME2(dst,q3) = 0

#else
#define bgq_elem0(arg) \
	(arg)[0]

#define bgq_elem1(arg) \
	(arg)[1]

#define bgq_elem2(arg) \
	(arg)[2]

#define bgq_elem3(arg) \
	(arg)[3]

#define bgq_vars(name) \
	name

#define bgq_vector4double_decl(name) \
	vector4double name

#define bgq_lda_double(dst,offset,addr) \
	(dst) = vec_lda(offset,(double*)(addr))

#define bgq_lda_float(dst,offset,addr) \
	(dst) = vec_lda(offset,(float*)(addr))

#define bgq_ld2a_double(dst,offset,addr) \
	(dst) = vec_ld2a(offset, (double*)(addr))

#define bgq_ld2a_float(dst,offset,addr) \
	(dst) = vec_ld2a(offset, (float*)(addr))

#define bgq_sta_double(src,offset,addr) \
	vec_sta(src, offset, (double*)(addr))

#define bgq_sta_float(src,offset,addr) \
	vec_sta(src, offset, (float*)(addr))

#define bgq_neg(dst,arg) \
	(dst) = vec_neg(arg);

#define bgq_add(dst,lhs,rhs) \
	(dst) = vec_add(lhs, rhs)

#define bgq_sub(dst,lhs,rhs) \
	(dst) = vec_sub(lhs, rhs)

#define bgq_mul(dst,lhs,rhs) \
	(dst) = vec_mul(lhs,rhs)

#define bgq_xxnpmadd(dst,a,b,c) \
	(dst) = vec_xxnpmadd(a,b,c)

#define bgq_iadd(dst,lhs,rhs) \
	bgq_xxnpmadd(dst,rhs,(vector4double)(1),lhs)

#define bgq_isub(dst,lhs,rhs) \
	bgq_xxcpnmadd(dst,rhs,(vector4double)(1),lhs)

#define bgq_merge2(dst, a23_to01, b01_to23) \
	(dst) = vec_perm(a23_to01, b01_to23, vec_gpci(02345))

#define bgq_xmul(dst,lhs,rhs) \
	(dst) = vec_xmul(lhs,rhs)

#define bgq_xmadd(dst,a,b,c) \
	(dst) = vec_xmadd(a,b,c)

#define bgq_madd(dst,a,b,c) \
	(dst) = vec_madd(a,b,c)

#define bgq_mov(dst,src) \
	(dst) = (src)

#define bgq_xxcpnmadd(dst,a,b,c) \
	(dst) = vec_xxcpnmadd(a,b,c)

#define bgq_cconst(dst,re,im) \
	(dst) = (vector4double){re,im,re,im}

#define bgq_zero(dst) \
	(dst) = (vector4double)(0)
//	dst = (vector4double){0,0,0,0}
//	dst = vec_logical(dst, dst, 0);

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
#define bgq_rxxcpnmadd(dst, a, b, c) bgq_xxcpnmadd(dst, b, a, c)
// re =    a.im * b.im + c.re
// im =   -a.im * b.re + c.im

// vec_xxmadd(a, b, c)
// re =    a.im * b.im + c.re
// im =    a.re * b.im + c.im

// vec_xxmadd(b, a, c)
// re =    a.im * b.im + c.re
// im =    a.im * b.re + c.im



#define bgq_cmplxval1(name) \
	(bgq_elem0(name) + (bgq_elem1(name) * _Complex_I))

#define bgq_cmplxval2(name) \
	(bgq_elem2(name) + (bgq_elem3(name) * _Complex_I))


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
		bgq_vector4double_decl(MAKENAME4(ccmul,dst,lhs,rhs));   \
		bgq_xmul(MAKENAME4(ccmul,dst,lhs,rhs),rhs,lhs);         \
		bgq_xxcpnmadd(dst,lhs,rhs,MAKENAME4(ccmul,dst,lhs,rhs)); \
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
		bgq_vector4double_decl(MAKENAME5(ccmadd,dst,a,b,c)); \
		bgq_xmadd(MAKENAME5(ccmadd,dst,a,b,c),b,a,c);        \
		bgq_xxcpnmadd(dst,a,b,MAKENAME5(ccmadd,dst,a,b,c));  \
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
	bgq_su3_vdecl(NAME2(name,v0));     \
	bgq_su3_vdecl(NAME2(name,v1));     \
	bgq_su3_vdecl(NAME2(name,v2));     \
	bgq_su3_vdecl(NAME2(name,v3))

#define bgq_su3_weyl_decl(name) \
	bgq_su3_vdecl(name##_v0);	\
	bgq_su3_vdecl(name##_v1)

#define bgq_su3_vzero(dst) \
	bgq_zero(NAME2(dst,c0)); \
	bgq_zero(NAME2(dst,c1)); \
	bgq_zero(NAME2(dst,c2))

#define bgq_su3_spinor_zero(dst) \
	bgq_su3_vzero(NAME2(dst,v0));     \
	bgq_su3_vzero(NAME2(dst,v1));     \
	bgq_su3_vzero(NAME2(dst,v2));     \
	bgq_su3_vzero(NAME2(dst,v3))

#define bgq_su3_spinor_load NAME2(bgq_su3_spinor_load,PRECISION)
#define bgq_su3_spinor_load_double(dst, addr) \
	bgq_lda_double(NAME3(dst,v0,c0),   0, addr);          \
	bgq_lda_double(NAME3(dst,v0,c1),  32, addr);          \
	bgq_lda_double(NAME3(dst,v0,c2),  64, addr);          \
	bgq_lda_double(NAME3(dst,v1,c0),  96, addr);          \
	bgq_lda_double(NAME3(dst,v1,c1), 128, addr);          \
	bgq_lda_double(NAME3(dst,v1,c2), 160, addr);          \
	bgq_lda_double(NAME3(dst,v2,c0), 192, addr);          \
	bgq_lda_double(NAME3(dst,v2,c1), 224, addr);          \
	bgq_lda_double(NAME3(dst,v2,c2), 256, addr);          \
	bgq_lda_double(NAME3(dst,v3,c0), 288, addr);          \
	bgq_lda_double(NAME3(dst,v3,c1), 320, addr);          \
	bgq_lda_double(NAME3(dst,v3,c2), 352, addr)

#define bgq_su3_spinor_load_float(dst, addr) \
	bgq_lda_float(NAME3(dst,v0,c0),   0, addr);          \
	bgq_lda_float(NAME3(dst,v0,c1),  16, addr);          \
	bgq_lda_float(NAME3(dst,v0,c2),  32, addr);          \
	bgq_lda_float(NAME3(dst,v1,c0),  48, addr);          \
	bgq_lda_float(NAME3(dst,v1,c1),  64, addr);          \
	bgq_lda_float(NAME3(dst,v1,c2),  80, addr);          \
	bgq_lda_float(NAME3(dst,v2,c0),  96, addr);          \
	bgq_lda_float(NAME3(dst,v2,c1), 112, addr);          \
	bgq_lda_float(NAME3(dst,v2,c2), 128, addr);          \
	bgq_lda_float(NAME3(dst,v3,c0), 144, addr);          \
	bgq_lda_float(NAME3(dst,v3,c1), 160, addr);          \
	bgq_lda_float(NAME3(dst,v3,c2), 176, addr)

#define bgq_su3_weyl_load NAME2(bgq_su3_weyl_load,PRECISION)
#define bgq_su3_weyl_load_double(dest, addr) \
	bgq_lda_double(NAME3(dest,v0,c0),   0, addr);        \
	bgq_lda_double(NAME3(dest,v0,c1),  32, addr);        \
	bgq_lda_double(NAME3(dest,v0,c2),  64, addr);        \
	bgq_lda_double(NAME3(dest,v1,c0),  96, addr);        \
	bgq_lda_double(NAME3(dest,v1,c1), 128, addr);        \
	bgq_lda_double(NAME3(dest,v1,c2), 160, addr)

#define bgq_su3_weyl_load_float(dest, addr) \
	bgq_lda_float(NAME3(dest,v0,c0),   0, addr);        \
	bgq_lda_float(NAME3(dest,v0,c1),  16, addr);        \
	bgq_lda_float(NAME3(dest,v0,c2),  32, addr);        \
	bgq_lda_float(NAME3(dest,v1,c0),  48, addr);        \
	bgq_lda_float(NAME3(dest,v1,c1),  64, addr);        \
	bgq_lda_float(NAME3(dest,v1,c2),  80, addr)

#define bgq_su3_weyl_load_left NAME2(bgq_su3_weyl_load_left,PRECISION)
#define bgq_su3_weyl_load_left_double(dest, addr) \
	bgq_ld2a_double(NAME3(dest,v0,c0),   0, addr);        \
	bgq_ld2a_double(NAME3(dest,v0,c1),  32, addr);        \
	bgq_ld2a_double(NAME3(dest,v0,c2),  64, addr);        \
	bgq_ld2a_double(NAME3(dest,v1,c0),  96, addr);        \
	bgq_ld2a_double(NAME3(dest,v1,c1), 128, addr);        \
	bgq_ld2a_double(NAME3(dest,v1,c2), 160, addr)

#define bgq_su3_weyl_load_left_float(dest, addr) \
	bgq_ld2a_float(NAME3(dest,v0,c0),   0, addr);        \
	bgq_ld2a_float(NAME3(dest,v0,c1),  16, addr);        \
	bgq_ld2a_float(NAME3(dest,v0,c2),  32, addr);        \
	bgq_ld2a_float(NAME3(dest,v1,c0),  48, addr);        \
	bgq_ld2a_float(NAME3(dest,v1,c1),  64, addr);        \
	bgq_ld2a_float(NAME3(dest,v1,c2),  80, addr)

#define bgq_su3_weyl_load_right NAME2(bgq_su3_weyl_load_right,PRECISION)
#define bgq_su3_weyl_load_right_double(dest, addr) \
	bgq_ld2a_double(NAME3(dest,v0,c0),   0+16, addr);        \
	bgq_ld2a_double(NAME3(dest,v0,c1),  32+16, addr);        \
	bgq_ld2a_double(NAME3(dest,v0,c2),  64+16, addr);        \
	bgq_ld2a_double(NAME3(dest,v1,c0),  96+16, addr);        \
	bgq_ld2a_double(NAME3(dest,v1,c1), 128+16, addr);        \
	bgq_ld2a_double(NAME3(dest,v1,c2), 160+16, addr)

#define bgq_su3_weyl_load_right_float(dest, addr) \
	bgq_ld2a_float(NAME3(dest,v0,c0),   0+8, addr);        \
	bgq_ld2a_float(NAME3(dest,v0,c1),  16+8, addr);        \
	bgq_ld2a_float(NAME3(dest,v0,c2),  32+8, addr);        \
	bgq_ld2a_float(NAME3(dest,v1,c0),  48+8, addr);        \
	bgq_ld2a_float(NAME3(dest,v1,c1),  64+8, addr);        \
	bgq_ld2a_float(NAME3(dest,v1,c2),  80+8, addr)

#define bgq_su3_matrix_load NAME2(bgq_su3_matrix_load,PRECISION)
#define bgq_su3_matrix_load_double(dest, addr) \
	bgq_lda_double(NAME2(dest,c00),   0, addr); \
	bgq_lda_double(NAME2(dest,c01),  32, addr); \
	bgq_lda_double(NAME2(dest,c02),  64, addr); \
	bgq_lda_double(NAME2(dest,c10),  96, addr); \
	bgq_lda_double(NAME2(dest,c11), 128, addr); \
	bgq_lda_double(NAME2(dest,c12), 160, addr); \
	bgq_lda_double(NAME2(dest,c20), 192, addr); \
	bgq_lda_double(NAME2(dest,c21), 224, addr); \
	bgq_lda_double(NAME2(dest,c22), 256, addr)

#define bgq_su3_matrix_load_float(dest, addr) \
	bgq_lda_float(NAME2(dest,c00),   0, addr); \
	bgq_lda_float(NAME2(dest,c01),  16, addr); \
	bgq_lda_float(NAME2(dest,c02),  32, addr); \
	bgq_lda_float(NAME2(dest,c10),  48, addr); \
	bgq_lda_float(NAME2(dest,c11),  64, addr); \
	bgq_lda_float(NAME2(dest,c12),  80, addr); \
	bgq_lda_float(NAME2(dest,c20),  96, addr); \
	bgq_lda_float(NAME2(dest,c21), 112, addr); \
	bgq_lda_float(NAME2(dest,c22), 128, addr)

#define bgq_su3_spinor_load_left NAME2(bgq_su3_spinor_load_left,PRECISION)
#define bgq_su3_spinor_load_left_double(dst,addr) \
	bgq_ld2a_double(NAME3(dst,v0,c0),   0, addr); \
	bgq_ld2a_double(NAME3(dst,v0,c1),  32, addr); \
	bgq_ld2a_double(NAME3(dst,v0,c2),  64, addr); \
	bgq_ld2a_double(NAME3(dst,v1,c0),  96, addr); \
	bgq_ld2a_double(NAME3(dst,v1,c1), 128, addr); \
	bgq_ld2a_double(NAME3(dst,v1,c2), 160, addr); \
	bgq_ld2a_double(NAME3(dst,v2,c0), 192, addr); \
	bgq_ld2a_double(NAME3(dst,v2,c1), 224, addr); \
	bgq_ld2a_double(NAME3(dst,v2,c2), 256, addr); \
	bgq_ld2a_double(NAME3(dst,v3,c0), 288, addr); \
	bgq_ld2a_double(NAME3(dst,v3,c1), 320, addr); \
	bgq_ld2a_double(NAME3(dst,v3,c2), 352, addr)

#define bgq_su3_spinor_load_left_float(dst,addr) \
	bgq_ld2a_float(NAME3(dst,v0,c0),   0, addr); \
	bgq_ld2a_float(NAME3(dst,v0,c1),  16, addr); \
	bgq_ld2a_float(NAME3(dst,v0,c2),  32, addr); \
	bgq_ld2a_float(NAME3(dst,v1,c0),  48, addr); \
	bgq_ld2a_float(NAME3(dst,v1,c1),  64, addr); \
	bgq_ld2a_float(NAME3(dst,v1,c2),  80, addr); \
	bgq_ld2a_float(NAME3(dst,v2,c0),  96, addr); \
	bgq_ld2a_float(NAME3(dst,v2,c1), 112, addr); \
	bgq_ld2a_float(NAME3(dst,v2,c2), 128, addr); \
	bgq_ld2a_float(NAME3(dst,v3,c0), 144, addr); \
	bgq_ld2a_float(NAME3(dst,v3,c1), 160, addr); \
	bgq_ld2a_float(NAME3(dst,v3,c2), 176, addr)

#define bgq_su3_spinor_load_right NAME2(bgq_su3_spinor_load_right,PRECISION)
#define bgq_su3_spinor_load_right_double(dst,addr) \
	bgq_ld2a_double(NAME3(dst,v0,c0),   0+16, addr); \
	bgq_ld2a_double(NAME3(dst,v0,c1),  32+16, addr); \
	bgq_ld2a_double(NAME3(dst,v0,c2),  64+16, addr); \
	bgq_ld2a_double(NAME3(dst,v1,c0),  96+16, addr); \
	bgq_ld2a_double(NAME3(dst,v1,c1), 128+16, addr); \
	bgq_ld2a_double(NAME3(dst,v1,c2), 160+16, addr); \
	bgq_ld2a_double(NAME3(dst,v2,c0), 192+16, addr); \
	bgq_ld2a_double(NAME3(dst,v2,c1), 224+16, addr); \
	bgq_ld2a_double(NAME3(dst,v2,c2), 256+16, addr); \
	bgq_ld2a_double(NAME3(dst,v3,c0), 288+16, addr); \
	bgq_ld2a_double(NAME3(dst,v3,c1), 320+16, addr); \
	bgq_ld2a_double(NAME3(dst,v3,c2), 352+16, addr)

#define bgq_su3_spinor_load_right_float(dst,addr) \
	bgq_ld2a_float(NAME3(dst,v0,c0),   0+8, addr); \
	bgq_ld2a_float(NAME3(dst,v0,c1),  16+8, addr); \
	bgq_ld2a_float(NAME3(dst,v0,c2),  32+8, addr); \
	bgq_ld2a_float(NAME3(dst,v1,c0),  48+8, addr); \
	bgq_ld2a_float(NAME3(dst,v1,c1),  64+8, addr); \
	bgq_ld2a_float(NAME3(dst,v1,c2),  80+8, addr); \
	bgq_ld2a_float(NAME3(dst,v2,c0),  96+8, addr); \
	bgq_ld2a_float(NAME3(dst,v2,c1), 112+8, addr); \
	bgq_ld2a_float(NAME3(dst,v2,c2), 128+8, addr); \
	bgq_ld2a_float(NAME3(dst,v3,c0), 144+8, addr); \
	bgq_ld2a_float(NAME3(dst,v3,c1), 160+8, addr); \
	bgq_ld2a_float(NAME3(dst,v3,c2), 176+8, addr)

#define bgq_su3_spinor_store NAME2(bgq_su3_spinor_store,PRECISION)
#define bgq_su3_spinor_store_double(addr,src) \
	bgq_sta_double(NAME3(src,v0,c0),   0, addr);          \
	bgq_sta_double(NAME3(src,v0,c1),  32, addr);          \
	bgq_sta_double(NAME3(src,v0,c2),  64, addr);          \
	bgq_sta_double(NAME3(src,v1,c0),  96, addr);          \
	bgq_sta_double(NAME3(src,v1,c1), 128, addr);          \
	bgq_sta_double(NAME3(src,v1,c2), 160, addr);          \
	bgq_sta_double(NAME3(src,v2,c0), 192, addr);          \
	bgq_sta_double(NAME3(src,v2,c1), 224, addr);          \
	bgq_sta_double(NAME3(src,v2,c2), 256, addr);          \
	bgq_sta_double(NAME3(src,v3,c0), 288, addr);          \
	bgq_sta_double(NAME3(src,v3,c1), 320, addr);          \
	bgq_sta_double(NAME3(src,v3,c2), 352, addr)

#define bgq_su3_spinor_store_float(addr,src) \
	bgq_sta_float(NAME3(src,v0,c0),   0, addr);          \
	bgq_sta_float(NAME3(src,v0,c1),  16, addr);          \
	bgq_sta_float(NAME3(src,v0,c2),  32, addr);          \
	bgq_sta_float(NAME3(src,v1,c0),  48, addr);          \
	bgq_sta_float(NAME3(src,v1,c1),  64, addr);          \
	bgq_sta_float(NAME3(src,v1,c2),  80, addr);          \
	bgq_sta_float(NAME3(src,v2,c0),  96, addr);          \
	bgq_sta_float(NAME3(src,v2,c1), 112, addr);          \
	bgq_sta_float(NAME3(src,v2,c2), 128, addr);          \
	bgq_sta_float(NAME3(src,v3,c0), 144, addr);          \
	bgq_sta_float(NAME3(src,v3,c1), 160, addr);          \
	bgq_sta_float(NAME3(src,v3,c2), 176, addr)

#define bgq_su3_weyl_store NAME2(bgq_su3_weyl_store,PRECISION)
#define bgq_su3_weyl_store_double(addr,src) \
	bgq_sta_double(NAME3(src,v0,c0),   0, addr);          \
	bgq_sta_double(NAME3(src,v0,c1),  32, addr);          \
	bgq_sta_double(NAME3(src,v0,c2),  64, addr);          \
	bgq_sta_double(NAME3(src,v1,c0),  96, addr);          \
	bgq_sta_double(NAME3(src,v1,c1), 128, addr);          \
	bgq_sta_double(NAME3(src,v1,c2), 160, addr)

#define bgq_su3_weyl_store_float(addr,src) \
	bgq_sta_float(NAME3(src,v0,c0),   0, addr);          \
	bgq_sta_float(NAME3(src,v0,c1),  16, addr);          \
	bgq_sta_float(NAME3(src,v0,c2),  32, addr);          \
	bgq_sta_float(NAME3(src,v1,c0),  48, addr);          \
	bgq_sta_float(NAME3(src,v1,c1),  64, addr);          \
	bgq_sta_float(NAME3(src,v1,c2),  80, addr)


#define bgq_su3_spinor_merge(dst,a,b)         \
	bgq_su3_vmerge(dst##_v0, a##_v0, b##_v0); \
	bgq_su3_vmerge(dst##_v1, a##_v1, b##_v1); \
	bgq_su3_vmerge(dst##_v2, a##_v2, b##_v2); \
	bgq_su3_vmerge(dst##_v3, a##_v3, b##_v3)

#define bgq_su3_weyl_merge(dst,a,b)         \
	bgq_su3_vmerge(NAME2(dst,v0), NAME2(a,v0), NAME2(b,v0)); \
	bgq_su3_vmerge(NAME2(dst,v1), NAME2(a,v1), NAME2(b,v1))

#define bgq_su3_vmov(dst,src)    \
	bgq_mov(NAME2(dst,c0), NAME2(src,c0)); \
	bgq_mov(NAME2(dst,c1), NAME2(src,c1)); \
	bgq_mov(NAME2(dst,c2), NAME2(src,c2))

#define bgq_su3_vmerge(dst,a,b)           \
	bgq_merge2(NAME2(dst,c0), NAME2(a,c0), NAME2(b,c0)); \
	bgq_merge2(NAME2(dst,c1), NAME2(a,c1), NAME2(b,c1)); \
	bgq_merge2(NAME2(dst,c2), NAME2(a,c2), NAME2(b,c2))

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

#define bgq_su3_cvmul(dst,c,v)              \
	bgq_cmul(NAME2(dst,c0), c, NAME2(v,c0)); \
	bgq_cmul(NAME2(dst,c1), c, NAME2(v,c1)); \
	bgq_cmul(NAME2(dst,c2), c, NAME2(v,c2))

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
		bgq_su3_vdecl(MAKENAME4(mvmul,dst,m,v));                          \
		bgq_ccmul (MAKENAME5(mvmul,dst,m,v,c0), NAME2(v,c0), NAME2(m,c00));         \
		bgq_ccmadd(MAKENAME5(mvmul,dst,m,v,c0), NAME2(v,c1), NAME2(m,c10), MAKENAME5(mvmul,dst,m,v,c0)); \
		bgq_ccmadd(MAKENAME5(mvmul,dst,m,v,c0), NAME2(v,c2), NAME2(m,c20), MAKENAME5(mvmul,dst,m,v,c0)); \
		bgq_ccmul (MAKENAME5(mvmul,dst,m,v,c1), NAME2(v,c0), NAME2(m,c01));         \
		bgq_ccmadd(MAKENAME5(mvmul,dst,m,v,c1), NAME2(v,c1), NAME2(m,c11), MAKENAME5(mvmul,dst,m,v,c1)); \
		bgq_ccmadd(MAKENAME5(mvmul,dst,m,v,c1), NAME2(v,c2), NAME2(m,c21), MAKENAME5(mvmul,dst,m,v,c1)); \
		bgq_ccmul (MAKENAME5(mvmul,dst,m,v,c2), NAME2(v,c0), NAME2(m,c02));         \
		bgq_ccmadd(MAKENAME5(mvmul,dst,m,v,c2), NAME2(v,c1), NAME2(m,c12), MAKENAME5(mvmul,dst,m,v,c2)); \
		bgq_ccmadd(MAKENAME5(mvmul,dst,m,v,c2), NAME2(v,c2), NAME2(m,c22), MAKENAME5(mvmul,dst,m,v,c2)); \
		bgq_su3_vmov(dst, MAKENAME4(mvmul,dst,m,v));                      \
	}

#define bgq_su3_spinor_mov(dst,src)  \
	bgq_su3_vmov(dst##_v0,src##_v0); \
	bgq_su3_vmov(dst##_v1,src##_v1); \
	bgq_su3_vmov(dst##_v2,src##_v2); \
	bgq_su3_vmov(dst##_v3,src##_v3)

#define bgq_su3_weyl_mov(dst,src)  \
	bgq_su3_vmov(dst##_v0,src##_v0); \
	bgq_su3_vmov(dst##_v1,src##_v1)


#if !BGQ_QPX

// No semantic effects
#define bgq_prefetch(addr)
#define bgq_prefetchforwrite(addr)
#define bgq_prefetch_forward(addr)
#define bgq_prefetch_backward(addr)
#define bgq_flush(addr)

#define bgq_l1_zero(addr) \
	/*memset((addr),0,128)*/

#else

#if defined(XLC)
	#define bgq_prefetch(addr) \
		__dcbt(addr)
// __prefetch_by_load(addr) generates an lbz instruction, which is 'blocking'
	#define bgq_prefetchforwrite(addr) \
		__dcbtst(addr)
	#define bgq_prefetch_forward(addr) \
		__prefetch_by_stream(1/*forward*/,(addr))
	#define bgq_prefetch_backward(addr) \
		__prefetch_by_stream(3/*backward*/,(addr))
	#define bgq_l1_zero(addr) \
		__dcbz(addr) /* sets 128 bytes (L2 chache line size) to zero */
	#define bgq_flush(addr) \
		__dcbf(addr)
#elif defined(__GNUC__)
	#define bgq_prefetch(addr) \
		__builtin_prefetch((addr),0/*read*/)
	#define bgq_prefetchforwrite(addr) \
		__builtin_prefetch((addr),1/*write*/)
	#define bgq_prefetch_forward(addr) \
		bgq_prefetch(addr)
	#define bgq_prefetch_backward(addr) \
		bgq_prefetch(addr)
	#define bgq_l1_zero(addr)
	#define bgq_flush(addr)
#else
	#define bgq_prefetch(addr)
	#define bgq_prefetchforwrite(addr)
	#define bgq_prefetch_forward(addr)
	#define bgq_prefetch_backward(addr)
	#define bgq_l1_zero(addr)
	#define bgq_flush(addr)
#endif

#endif


#define bgq_su3_spinor_prefetch NAME2(bgq_su3_spinor_prefetch,PRECISION)
#define bgq_su3_spinor_prefetch_double(addr)     \
	bgq_prefetch((char*)(addr) +   0);    \
	bgq_prefetch((char*)(addr) +  64);    \
	bgq_prefetch((char*)(addr) + 128);    \
	bgq_prefetch((char*)(addr) + 192);    \
	bgq_prefetch((char*)(addr) + 256);    \
	bgq_prefetch((char*)(addr) + 320)

#define bgq_su3_spinor_prefetch_float(addr)     \
	bgq_prefetch((char*)(addr) +   0);    \
	bgq_prefetch((char*)(addr) +  64);    \
	bgq_prefetch((char*)(addr) + 128)


#define bgq_su3_weyl_prefetch NAME2(bgq_su3_weyl_prefetch,PRECISION)
#define bgq_su3_weyl_prefetch_double(addr)      \
	bgq_prefetch((char*)(addr) +   0);    \
	bgq_prefetch((char*)(addr) +  64);    \
	bgq_prefetch((char*)(addr) + 128)

#define bgq_su3_weyl_prefetch_float(addr)      \
	bgq_prefetch((char*)(addr) +   0);    \
	bgq_prefetch((char*)(addr) +  64) /* half a cacheline */


#define bgq_su3_matrix_prefetch NAME2(bgq_su3_matrix_prefetch,PRECISION)
#define bgq_su3_matrix_prefetch_double(addr)      \
	bgq_prefetch((char*)(addr) +   0);    \
	bgq_prefetch((char*)(addr) +  64);    \
	bgq_prefetch((char*)(addr) + 128);    \
	bgq_prefetch((char*)(addr) + 192);    \
	bgq_prefetch((char*)(addr) + 256)

#define bgq_su3_matrix_prefetch_float(addr)      \
	bgq_prefetch((char*)(addr) +   0);    \
	bgq_prefetch((char*)(addr) +  64);    \
	bgq_prefetch((char*)(addr) + 128) /* quarter of a cacheline */


#define bgq_su3_spinor_zeroload NAME2(bgq_su3_spinor_zeroload,PRECISION)
#define bgq_su3_spinor_zeroload_double(addr) \
	bgq_l1_zero((char*)(addr) +   0);    \
	bgq_l1_zero((char*)(addr) + 128);    \
	bgq_l1_zero((char*)(addr) + 256)
	// 384

#define bgq_su3_spinor_zeroload_float(addr) \
	bgq_l1_zero((char*)(addr) + 0);          \
	bgq_prefetchforwrite((char*)(addr) +  128)
	// 192


#define bgq_su3_weyl_zeroload NAME2(bgq_su3_weyl_zeroload,PRECISION)
#define bgq_su3_weyl_zeroload_double(addr) \
	bgq_l1_zero((char*)(addr) + 0);          \
	bgq_prefetchforwrite((char*)(addr) +  128)
	// 192

#define bgq_su3_weyl_zeroload_float(addr)    \
	bgq_prefetchforwrite((char*)(addr) +  0); \
	bgq_prefetchforwrite((char*)(addr) +  64)
	// 96



#define bgq_su3_spinor_flush NAME2(bgq_su3_spinor_flush,PRECISION)
#define bgq_su3_spinor_flush_double(addr) \
	bgq_flush((char*)(addr) +   0);    \
	bgq_flush((char*)(addr) +  64);    \
	bgq_flush((char*)(addr) + 128);    \
	bgq_flush((char*)(addr) + 192);    \
	bgq_flush((char*)(addr) + 256);    \
	bgq_flush((char*)(addr) + 320)
	// 384 bytes
#define bgq_su3_spinor_flush_float(addr) \
	bgq_flush((char*)(addr) +   0);    \
	bgq_flush((char*)(addr) +  64);    \
	bgq_flush((char*)(addr) + 128)
	// 192 bytes

#define bgq_su3_weyl_flush NAME2(bgq_su3_weyl_flush,PRECISION)
#define bgq_su3_weyl_flush_double(addr) \
	bgq_flush((char*)(addr) +   0);    \
	bgq_flush((char*)(addr) +  64);    \
	bgq_flush((char*)(addr) + 128)
	// 192 bytes

#define bgq_su3_weyl_flush_float(addr)    \
	bgq_flush((char*)(addr) +   0);    \
	bgq_flush((char*)(addr) +  64)
	// 96 bytes

#endif /* BGQ_H_ */

