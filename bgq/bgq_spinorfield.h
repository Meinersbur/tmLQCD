/*
 * bgq_spinorfield.h
 *
 *  Created on: Oct 26, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_SPINORFIELD_H_
#define BGQ_SPINORFIELD_H_

#include "bgq_field.h"
#include "bgq_utils.h"
#include "bgq_qpx.h"

#include <stdbool.h>

#ifndef BGQ_SPINORFIELD_C_
#define EXTERN_INLINE EXTERN_INLINE_DECLARATION
#define EXTERN_FIELD extern
#define EXTERN_INIT(val)
#else
#define EXTERN_INLINE EXTERN_INLINE_DEFINITION
#define EXTERN_FIELD
#define EXTERN_INIT(val) = (val)
#endif


typedef struct {
	COMPLEX_PRECISION c[3];
} bgq_su3vector;

typedef struct {
	bgq_su3vector v[4];
} bgq_spinor;
typedef bgq_spinor bgq_spinor_nonvec;

typedef struct {
	COMPLEX_PRECISION s[2][3]; // 96 byte
} bgq_weyl_nonvec;

#define bgq_weyl_fromqpx(arg) bgq_weyl_fromqpx_raw(bgq_su3_weyl_vars(arg))
EXTERN_INLINE bgq_weyl_vec bgq_weyl_fromqpx_raw(bgq_su3_weyl_params(weyl)) {
	bgq_weyl_vec result;
	result.s[0][0][0] = bgq_cmplxval1(weyl_v0_c0);
	result.s[0][0][1] = bgq_cmplxval2(weyl_v0_c0);
	result.s[0][1][0] = bgq_cmplxval1(weyl_v0_c1);
	result.s[0][1][1] = bgq_cmplxval2(weyl_v0_c1);
	result.s[0][2][0] = bgq_cmplxval1(weyl_v0_c2);
	result.s[0][2][1] = bgq_cmplxval2(weyl_v0_c2);
	result.s[1][0][0] = bgq_cmplxval1(weyl_v1_c0);
	result.s[1][0][1] = bgq_cmplxval2(weyl_v1_c0);
	result.s[1][1][0] = bgq_cmplxval1(weyl_v1_c1);
	result.s[1][1][1] = bgq_cmplxval2(weyl_v1_c1);
	result.s[1][2][0] = bgq_cmplxval1(weyl_v1_c2);
	result.s[1][2][1] = bgq_cmplxval2(weyl_v1_c2);
	return result;
}


EXTERN_INLINE bgq_weyl_nonvec bgq_weyl_extractvec(bgq_weyl_vec weylvec, size_t k) {
	assert(0 <= k && k < PHYSICAL_LK);
	bgq_weyl_nonvec result;
	for (size_t v = 0; v < 2; v+=1) {
		for (size_t c = 0; c < 3; c+=1) {
			result.s[v][c] = weylvec.s[v][c][k];
		}
	}
	return result;
}


#define bgq_spinor_fromqpx(spinor,k) bgq_spinor_fromqpx_raw(bgq_su3_spinor_vars(spinor),k)
EXTERN_INLINE bgq_spinor bgq_spinor_fromqpx_raw(bgq_su3_spinor_params(spinor), ucoord k) {
	assert(0 <= k && k < PHYSICAL_LK);

	bgq_spinor result;
	result.v[0].c[0] = bgq_cmplxval(spinor_v0_c0,k);
	result.v[0].c[1] = bgq_cmplxval(spinor_v0_c1,k);
	result.v[0].c[2] = bgq_cmplxval(spinor_v0_c2,k);
	result.v[1].c[0] = bgq_cmplxval(spinor_v1_c0,k);
	result.v[1].c[1] = bgq_cmplxval(spinor_v1_c1,k);
	result.v[1].c[2] = bgq_cmplxval(spinor_v1_c2,k);
	result.v[2].c[0] = bgq_cmplxval(spinor_v2_c0,k);
	result.v[2].c[1] = bgq_cmplxval(spinor_v2_c1,k);
	result.v[2].c[2] = bgq_cmplxval(spinor_v2_c2,k);
	result.v[3].c[0] = bgq_cmplxval(spinor_v3_c0,k);
	result.v[3].c[1] = bgq_cmplxval(spinor_v3_c1,k);
	result.v[3].c[2] = bgq_cmplxval(spinor_v3_c2,k);

	return result;
}


EXTERN_INLINE bgq_spinor bgq_spinor_fromvec(bgq_spinor_vec spinorvec, ucoord k) {
	assert(0 <= k && k < PHYSICAL_LK);

	bgq_spinor result;
	for (ucoord v = 0; v < 4; v+=1) {
		for (ucoord c = 0; c < 3; c+=1) {
			result.v[v].c[c] = spinorvec.s[v][c][k];
		}
	}
	return result;
}



#define bgq_weyl_extractfromqpxvec(weyl, k) bgq_weyl_extractfromqpxvec_raw(bgq_su3_weyl_vars(weyl),k)
EXTERN_INLINE bgq_weyl_nonvec bgq_weyl_extractfromqpxvec_raw(bgq_su3_weyl_params(weyl), size_t k) {
	bgq_weyl_vec vec = bgq_weyl_fromqpx(weyl);
	return bgq_weyl_extractvec(vec,k);
}


void bgq_spinorfield_setup(bgq_weylfield_controlblock *field, bool isOdd, bool readFullspinor, bool writeFullspinor, bool readWeyl, bool writeWeyl);
void bgq_spinorfield_transfer(bool isOdd, bgq_weylfield_controlblock *targetfield, spinor* sourcefield);
double bgq_spinorfield_compare(bool isOdd, bgq_weylfield_controlblock *bgqfield, spinor *reffield, bool silent);





EXTERN_INLINE void bgq_spinor_expect(bgq_spinor spinor, scoord t,scoord x,scoord y,scoord z) {
#ifdef BGQ_COORDCHECK
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	ucoord t_global = bgq_local2global_t(t);
	ucoord x_global = bgq_local2global_x(x);
	ucoord y_global = bgq_local2global_y(y);
	ucoord z_global = bgq_local2global_z(z);

	assert(spinor.v[0].c[0] == t_global);
	assert(spinor.v[0].c[1] == x_global);
	assert(spinor.v[0].c[2] == y_global);
	assert(spinor.v[1].c[0] == z_global);
	assert(spinor.v[1].c[1] == 0);
	assert(spinor.v[1].c[2] == 0);
	assert(spinor.v[2].c[0] == 0);
	assert(spinor.v[2].c[1] == 0);
	assert(spinor.v[2].c[2] == 0);
	assert(spinor.v[3].c[0] == 0);
	assert(spinor.v[3].c[1] == 0);
	assert(spinor.v[3].c[2] == 0);
#endif
}


EXTERN_INLINE void bgq_spinorveck_expect(bgq_spinor_vec spinor, ucoord k, scoord t,scoord x,scoord y,scoord z) {
	assert(0 <= k && k < PHYSICAL_LK);
	bgq_spinor_expect(bgq_spinor_fromvec(spinor,k),t,x,y,z);
}

EXTERN_INLINE void bgq_spinorvec_expect(bgq_spinor_vec spinor, scoord t1, scoord t2, scoord x,scoord y,scoord z) {
	bgq_spinor_expect(bgq_spinor_fromvec(spinor,0),t1,x,y,z);
	bgq_spinor_expect(bgq_spinor_fromvec(spinor,1),t2,x,y,z);
}



#define bgq_spinorqpx_expect(spinor,t_left,t_right,x,y,z) bgq_spinorqpx_expect_raw(bgq_su3_spinor_vars(spinor),t_left,t_right,x,y,z)
EXTERN_INLINE void bgq_spinorqpx_expect_raw(bgq_su3_spinor_params(spinor),scoord t_left,scoord t_right,scoord x,scoord y,scoord z) {
	bgq_spinor_expect(bgq_spinor_fromqpx(spinor,0),t_left,x,y,z);
	bgq_spinor_expect(bgq_spinor_fromqpx(spinor,1),t_right,x,y,z);
}


#ifdef BGQ_COORDCHECK
#define BGQ_WEYLLEFT_EXPECT(weyl,t,x,y,z) \
		bgq_spinorfield_weyl_expect(bgq_weyl_extractfromqpxvec(weyl,0),t,x,y,z)
#define BGQ_WEYLRIGHT_EXPECT(weyl,t,x,y,z) \
		bgq_spinorfield_weyl_expect(bgq_weyl_extractfromqpxvec(weyl,1),t,x,y,z)
#define BGQ_WEYL_EXPECT(weyl,t1,t2,x,y,z) \
		BGQ_WEYLLEFT_EXPECT(weyl,t1,x,y,z); \
		BGQ_WEYLRIGHT_EXPECT(weyl,t2,x,y,z)
#else
#define BGQ_WEYLLEFT_EXPECT(weyl,t,x,y,z)
#define BGQ_WEYLRIGHT_EXPECT(weyl,t,x,y,z)
#define BGQ_WEYL_EXPECT(weyl,t,x,y,z)
#endif


EXTERN_INLINE bgq_weyl_nonvec bgq_weyl_fromvec(bgq_weyl_vec weylvec, ucoord k) {
	assert(0 <= k && k < PHYSICAL_LK);
	bgq_weyl_nonvec result;
	for (size_t v = 0; v < 2; v += 1) {
		for (size_t c = 0; c < 3; c += 1) {
			result.s[v][c] = weylvec.s[v][c][k];
		}
	}
	return result;
}


void bgq_weyl_expect(bgq_weyl_nonvec weyl, ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc);

EXTERN_INLINE void bgq_weylvec_expect(bgq_weyl_vec weyl, ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc) {
#ifdef BGQ_COORDCHECK
	bgq_weyl_expect(bgq_weyl_fromvec(weyl, 0), t1, x, y, z, d, isSrc);
	bgq_weyl_expect(bgq_weyl_fromvec(weyl, 1), t2, x, y, z, d, isSrc);
#endif
}




bgq_spinor bgq_spinorfield_getspinor(bgq_weylfield_controlblock *field, size_t t, size_t x, size_t y, size_t z) ;

void bgq_weylveck_written(bgq_weyl_vec *targetweyl, ucoord k, ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc);

EXTERN_INLINE void bgq_weylvec_written(bgq_weyl_vec *targetweyl, ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc) {
#ifdef BGQ_COORDCHECK
	bgq_weylveck_written(targetweyl, 0, t1, x,y,z,d,isSrc);
	bgq_weylveck_written(targetweyl, 1, t2, x,y,z,d,isSrc);
#endif
}



#undef EXTERN_INLINE
#undef EXTERN_FIELD
#undef EXTERN_INIT

#endif /* BGQ_SPINORFIELD_H_ */
