/*
 * bgq.h
 *
 *  Created on: Jul 25, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_H_
#define BGQ_H_

typedef double vector4double[4];



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

#define cvec_madd(a,b,c) vec_xxnpmadd(b,a,vec_xmadd(a,b,c))
// vec_xxnpmadd(b,a,vec_xmadd(a,b,c))
// vec_xxnpmadd(a,b,vec_xmadd(b,a,c))
// re = a.re * b.re - a.im * b.im + c.re
// im = a.re * b.im + a.im * b.re + c.im
//    = a * b + c

// conjugated second argument
#define cvec_pcmul(a,b) vec_xxcpnmadd(a,b,vec_xmul(b,a))
// vec_xxcpnmadd(a,b,vec_xmul(b,a))
// re = a.re * b.re + a.im * b.im
// im = a.im * b.re - a.re * b.im
//    = a * conjugate(b)

// conjugated second argument
#define cvec_pcpmadd(a,b,c) vec_xxcpnmadd(a,b,vec_xmadd(b,a,c))
// vec_xxcpnmadd(a,b,vec_xmadd(b,a,c))
// vec_xmadd(b,a,vec_xxcpnmadd(a,b,c))
// re = a.re * b.re + a.im * b.im + c.re
// im = a.im * b.re - a.re * b.im + c.im
//    = a * conjugate(b) + c

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

#define bgq_su3_vdecl(name)  \
	vector4double name##_c0; \
	vector4double name##_c1; \
	vector4double name##_c2

#define bgq_su3_mdecl(name)   \
	vector4double name##_c00; \
	vector4double name##_c01; \
	vector4double name##_c02; \
	vector4double name##_c10; \
	vector4double name##_c11; \
	vector4double name##_c12; \
	vector4double name##_c20; \
	vector4double name##_c21; \
	vector4double name##_c22

#define bgq_su3_spinor_decl(name) \
	bgq_su3_vdecl(name##_v0);     \
	bgq_su3_vdecl(name##_v1);     \
	bgq_su3_vdecl(name##_v2);     \
	bgq_su3_vdecl(name##_v3)

#define bgq_su3_weyl_decl(name) \
	bgq_su3_vdecl(name##_v0);	\
	bgq_su3_vdecl(name##_v1)


#define bgq_su3_spinor_double_load(dest, addr)   \
	dest##_v0_c0 = vec_lda(  0, addr); \
	dest##_v0_c1 = vec_lda( 32, addr); \
	dest##_v0_c2 = vec_lda( 64, addr); \
	dest##_v1_c0 = vec_lda( 96, addr); \
	dest##_v1_c1 = vec_lda(128, addr); \
	dest##_v1_c2 = vec_lda(160, addr); \
	dest##_v2_c0 = vec_lda(192, addr); \
	dest##_v2_c1 = vec_lda(224, addr); \
	dest##_v2_c2 = vec_lda(256, addr); \
	dest##_v3_c0 = vec_lda(288, addr); \
	dest##_v3_c1 = vec_lda(320, addr); \
	dest##_v3_c2 = vec_lda(352, addr)
// NOTE: qvlfdux is more effective, but no compiler built-in exists

#define bgq_su3_weyl_double_load(dest, addr)   \
	dest##_v0_c0 = vec_lda(  0, addr); \
	dest##_v0_c1 = vec_lda( 32, addr); \
	dest##_v0_c2 = vec_lda( 64, addr); \
	dest##_v1_c0 = vec_lda( 96, addr); \
	dest##_v1_c1 = vec_lda(128, addr); \
	dest##_v1_c2 = vec_lda(160, addr)

#define bgq_su3_matrix_double_load(dest, addr)   \
	dest##_c00 = vec_lda(  0, addr); \
	dest##_c01 = vec_lda( 32, addr); \
	dest##_c02 = vec_lda( 64, addr); \
	dest##_c10 = vec_lda( 96, addr); \
	dest##_c11 = vec_lda(128, addr); \
	dest##_c12 = vec_lda(160, addr); \
	dest##_c20 = vec_lda(192, addr); \
	dest##_c21 = vec_lda(224, addr); \
	dest##_c22 = vec_lda(256, addr)

#define bgq_su3_matrix_double_load_left(dest, addr)   \
	dest##_c00 = vec_ld2a(  0, addr); \
	dest##_c01 = vec_ld2a( 32, addr); \
	dest##_c02 = vec_ld2a( 64, addr); \
	dest##_c10 = vec_ld2a( 96, addr); \
	dest##_c11 = vec_ld2a(128, addr); \
	dest##_c12 = vec_ld2a(160, addr); \
	dest##_c20 = vec_ld2a(192, addr); \
	dest##_c21 = vec_ld2a(224, addr); \
	dest##_c22 = vec_ld2a(256, addr)

#define bgq_su3_matrix_double_load_right(dest, addr)   \
	dest##_c00 = vec_ld2a(16+  0, addr); \
	dest##_c01 = vec_ld2a(16+ 32, addr); \
	dest##_c02 = vec_ld2a(16+ 64, addr); \
	dest##_c10 = vec_ld2a(16+ 96, addr); \
	dest##_c11 = vec_ld2a(16+128, addr); \
	dest##_c12 = vec_ld2a(16+160, addr); \
	dest##_c20 = vec_ld2a(16+192, addr); \
	dest##_c21 = vec_ld2a(16+224, addr); \
	dest##_c22 = vec_ld2a(16+256, addr)

#define bgq_su3_spinor_double_load_left(dst,addr) \
	dst##_v0_c0 = vec_ld2a(  0, addr); \
	dst##_v0_c1 = vec_ld2a( 32, addr); \
	dst##_v0_c2 = vec_ld2a( 64, addr); \
	dst##_v1_c0 = vec_ld2a( 96, addr); \
	dst##_v1_c1 = vec_ld2a(128, addr); \
	dst##_v1_c2 = vec_ld2a(160, addr); \
	dst##_v2_c0 = vec_ld2a(192, addr); \
	dst##_v2_c1 = vec_ld2a(224, addr); \
	dst##_v2_c2 = vec_ld2a(256, addr); \
	dst##_v3_c0 = vec_ld2a(288, addr); \
	dst##_v3_c1 = vec_ld2a(320, addr); \
	dst##_v3_c2 = vec_ld2a(352, addr)

#define bgq_su3_spinor_double_load_right(dst,addr) \
	dst##_v0_c0 = vec_ld2a(16+  0, addr); \
	dst##_v0_c1 = vec_ld2a(16+ 32, addr); \
	dst##_v0_c2 = vec_ld2a(16+ 64, addr); \
	dst##_v1_c0 = vec_ld2a(16+ 96, addr); \
	dst##_v1_c1 = vec_ld2a(16+128, addr); \
	dst##_v1_c2 = vec_ld2a(16+160, addr); \
	dst##_v2_c0 = vec_ld2a(16+192, addr); \
	dst##_v2_c1 = vec_ld2a(16+224, addr); \
	dst##_v2_c2 = vec_ld2a(16+256, addr); \
	dst##_v3_c0 = vec_ld2a(16+288, addr); \
	dst##_v3_c1 = vec_ld2a(16+320, addr); \
	dst##_v3_c2 = vec_ld2a(16+352, addr)

#define bgq_su3_spinor_double_store(addr,src) \
	vec_sta(src##_v0_c0,   0, addr);          \
	vec_sta(src##_v0_c1,  32, addr);          \
	vec_sta(src##_v0_c2,  64, addr);          \
	vec_sta(src##_v1_c0,  96, addr);          \
	vec_sta(src##_v1_c1, 128, addr);          \
	vec_sta(src##_v1_c2, 160, addr);          \
	vec_sta(src##_v2_c0, 192, addr);          \
	vec_sta(src##_v2_c1, 224, addr);          \
	vec_sta(src##_v2_c2, 256, addr);          \
	vec_sta(src##_v3_c0, 288, addr);          \
	vec_sta(src##_v3_c1, 320, addr);          \
	vec_sta(src##_v3_c2, 352, addr)

#define bgq_su3_spinor_merge(dst,a,b)        \
	bgq_su3_vmerge(dst##_v0, a##_v0, b##_v0) \
	bgq_su3_vmerge(dst##_v1, a##_v1, b##_v1) \
	bgq_su3_vmerge(dst##_v2, a##_v2, b##_v2)

#define bgq_su3_vmov(dst,src) \
	dst##_c0 = src##_c0;     \
	dst##_c1 = src##_c1;     \
	dst##_c2 = src##_c2

#define bgq_su3_vmerge(dst,a,b)            \
	dst##_c0 = cvec_merge(a##_c0, b##_c0); \
	dst##_c1 = cvec_merge(a##_c1, b##_c1); \
	dst##_c2 = cvec_merge(a##_c2, b##_c2)

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

#define bgq_su3_vadd(dest,v1,v2)           \
	dest##_c0 = vec_add(v1##_c0, v2##_c0); \
	dest##_c1 = vec_add(v1##_c1, v2##_c1); \
	dest##_c2 = vec_add(v1##_c2, v2##_c2)

#define bgq_su3_vsub(dest,v1,v2)           \
	dest##_c0 = vec_sub(v1##_c0, v2##_c0); \
	dest##_c1 = vec_sub(v1##_c1, v2##_c1); \
	dest##_c2 = vec_sub(v1##_c2, v2##_c2)

#define bgq_su3_vpiadd(dst,v1,v2)            \
	dst##_c0 = cvec_piadd(v1##_c0, v2##_c0); \
	dst##_c1 = cvec_piadd(v1##_c1, v2##_c1); \
	dst##_c2 = cvec_piadd(v1##_c2, v2##_c2)

#define bgq_su3_vpisub(dst,v1,v2)            \
	dst##_c0 = cvec_pisub(v1##_c0, v2##_c0); \
	dst##_c1 = cvec_pisub(v1##_c1, v2##_c1); \
	dst##_c2 = cvec_pisub(v1##_c2, v2##_c2)


#define bgq_su3_cvmul(dst,c,v)       \
	dst##_c0 = cvec_mul(c, v##_c0); \
	dst##_c1 = cvec_mul(c, v##_c1); \
	dst##_c2 = cvec_mul(c, v##_c2); \

#define bgq_su3_mvmul(dest,m,v)                        \
	dest##_c0 = cvec_mul (v##_c0, m##_c00);            \
	dest##_c0 = cvec_madd(v##_c1, m##_c01, dest##_c0); \
	dest##_c0 = cvec_madd(v##_c2, m##_c01, dest##_c0); \
	dest##_c1 = cvec_mul (v##_c0, m##_c10);            \
	dest##_c1 = cvec_madd(v##_c1, m##_c11, dest##_c1); \
	dest##_c1 = cvec_madd(v##_c2, m##_c11, dest##_c1); \
	dest##_c2 = cvec_mul (v##_c0, m##_c20);            \
	dest##_c2 = cvec_madd(v##_c1, m##_c21, dest##_c2); \
	dest##_c2 = cvec_madd(v##_c2, m##_c21, dest##_c2);

#define bgq_su3_mvinvmul(dest,m,v)                        \
	dest##_c0 = cvec_pcmul  (v##_c0, m##_c00);            \
	dest##_c0 = cvec_pcpmadd(v##_c1, m##_c10, dest##_c0); \
	dest##_c0 = cvec_pcpmadd(v##_c2, m##_c10, dest##_c0); \
	dest##_c1 = cvec_pcmul  (v##_c0, m##_c01);            \
	dest##_c1 = cvec_pcpmadd(v##_c1, m##_c11, dest##_c1); \
	dest##_c1 = cvec_pcpmadd(v##_c2, m##_c11, dest##_c1); \
	dest##_c2 = cvec_pcmul  (v##_c0, m##_c02);            \
	dest##_c2 = cvec_pcpmadd(v##_c1, m##_c12, dest##_c2); \
	dest##_c2 = cvec_pcpmadd(v##_c2, m##_c12, dest##_c2);

#endif /* BGQ_H_ */
