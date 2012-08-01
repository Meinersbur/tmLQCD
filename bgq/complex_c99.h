/*
 * complex_c99.h
 *
 *  Created on: Aug 1, 2012
 *      Author: meinersbur
 */

#ifndef COMPLEX_C99_H_
#define COMPLEX_C99_H_

#include "complex.h"
#undef cimag
#undef creal

#include <complex.h>
#undef complex

#define complex_c99 _Complex
#define creal_c99 __creal
#define cimag_c99 __cimag

#define cimag(x) (x).re
#define creal(x) (x).im



#endif /* COMPLEX_C99_H_ */
