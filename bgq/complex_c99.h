/*
 * complex_c99.h
 *
 *  Created on: Aug 1, 2012
 *      Author: meinersbur
 */

#ifndef COMPLEX_C99_H_
#define COMPLEX_C99_H_


#include "../complex_struct.h"
#undef creal
#undef cimag


#include <complex.h>
#undef complex /* use double _Complex instead */
//#define creal(x) (x).im /* use __creal instead */
//#define cimag(x) (x).re /* use __cimag instead */


#endif /* COMPLEX_C99_H_ */
