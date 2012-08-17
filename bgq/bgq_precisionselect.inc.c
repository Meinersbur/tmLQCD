/*
 * bgq_precisionselect.inc.c
 *
 *  Created on: Aug 15, 2012
 *      Author: meinersbur
 */


#undef PRECISION
#undef COMPLEX_PRECISION
#undef PRECISION_DOUBLE
#undef PRECISION_FLOAT
#undef PRECISION_BITS

#undef bgq_su3_spinor_load
#undef bgq_su3_spinor_load_left
#undef bgq_su3_spinor_load_right
#undef bgq_su3_matrix_load
#undef bgq_su3_weyl_load
#undef bgq_su3_weyl_load_left
#undef bgq_su3_spinor_store
#undef bgq_su3_weyl_store
#undef bgq_su3_spinor_prefetch
#undef bgq_su3_matrix_prefetch
#undef bgq_su3_weyl_prefetch


#ifndef BGQ_PRECISION

// Just undef everything

#else

#if (BGQ_PRECISION==64)
// double precision
#define PRECISION double
#define COMPLEX_PRECISION complexdouble
#define PRECISION_DOUBLE 1
#define PRECISION_FLOAT 0
#define PRECISION_BITS 64

#define bgq_su3_spinor_load bgq_su3_spinor_load_double
#define bgq_su3_spinor_load_left bgq_su3_spinor_load_left_double
#define bgq_su3_spinor_load_right bgq_su3_spinor_load_right_double
#define bgq_su3_matrix_load bgq_su3_matrix_load_double
#define bgq_su3_weyl_load bgq_su3_weyl_load_double
#define bgq_su3_weyl_load_left bgq_su3_weyl_load_left_double

#define bgq_su3_spinor_store bgq_su3_spinor_store_double
#define bgq_su3_weyl_store bgq_su3_weyl_store_double

#define bgq_su3_spinor_prefetch bgq_su3_spinor_prefetch_double
#define bgq_su3_matrix_prefetch bgq_su3_matrix_prefetch_double
#define bgq_su3_weyl_prefetch bgq_su3_weyl_prefetch_double

#elif (BGQ_PRECISION==32)
// single precision
#define PRECISION float
#define COMPLEX_PRECISION complexfloat
#define PRECISION_DOUBLE 0
#define PRECISION_FLOAT 1
#define PRECISION_BITS 32

#define bgq_su3_spinor_load bgq_su3_spinor_load_float
#define bgq_su3_spinor_load_left bgq_su3_spinor_load_left_float
#define bgq_su3_spinor_load_right bgq_su3_spinor_load_right_float
#define bgq_su3_matrix_load bgq_su3_matrix_load_float
#define bgq_su3_weyl_load bgq_su3_weyl_load_float
#define bgq_su3_weyl_load_left bgq_su3_weyl_load_left_float

#define bgq_su3_spinor_store bgq_su3_spinor_store_float
#define bgq_su3_weyl_store bgq_su3_weyl_store_float

#define bgq_su3_spinor_prefetch bgq_su3_spinor_prefetch_float
#define bgq_su3_matrix_prefetch bgq_su3_matrix_prefetch_float
#define bgq_su3_weyl_prefetch bgq_su3_weyl_prefetch_float

#else

#warning Unsupported precision

#endif
#endif
