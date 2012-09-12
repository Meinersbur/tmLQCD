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
#undef PRECISION_BYTES
#undef MPI_PRECISION

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
		#define PRECISION_BYTES 8
		#define MPI_PRECISION MPI_DOUBLE
	#elif (BGQ_PRECISION==32)
		// single precision
		#define PRECISION float
		#define COMPLEX_PRECISION complexfloat
		#define PRECISION_DOUBLE 0
		#define PRECISION_FLOAT 1
		#define PRECISION_BITS 32
		#define PRECISION_BYTES 4
		#define MPI_PRECISION MPI_FLOAT
	#else
		#warning Unsupported precision
	#endif
#endif

#define PRECISION_VECTOR_ALIGNMENT (PRECISION_BITS*4/8)
