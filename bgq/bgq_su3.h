#ifndef BGQ_SU3_H_
#define BGQ_SU3_H_

#include "bgq_utils.h"

typedef struct {
	complex_double c[3];
} bgq_su3vector_double;
typedef bgq_su3vector_double bgq_su3vector;

typedef struct {
	bgq_su3vector v[4];
} bgq_spinor;
typedef bgq_spinor bgq_spinor_nonvec;

typedef struct {
	complex_double s[2][3]; // 96 byte
} bgq_weyl_nonvec_double;
typedef struct {
	complex_float s[2][3]; // 48 byte
} bgq_weyl_nonvec_float;
#define bgq_weyl_nonvec NAME2(bgq_weyl_nonvec,PRECISION)


#endif /* BGQ_SU3_H_ */
