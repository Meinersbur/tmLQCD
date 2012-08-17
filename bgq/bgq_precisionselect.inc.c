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
	#elif (BGQ_PRECISION==32)
		// single precision
		#define PRECISION float
		#define COMPLEX_PRECISION complexfloat
		#define PRECISION_DOUBLE 0
		#define PRECISION_FLOAT 1
		#define PRECISION_BITS 32
	#else
		#warning Unsupported precision
	#endif
#endif
