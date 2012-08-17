/*
 * bgq_HoppingMatrix.c
 *
 *  Created on: Aug 16, 2012
 *      Author: meinersbur
 */

#include "bgq_HoppingMatrix.h"



//#ifdef HAVE_CONFIG_H
//#include <config.h>
//#endif

//#include "bgq_field.h"
//#include "bgq.h"
//#include "global.h"
//#include "boundary.h"
//#include <mpi.h>
//#include <omp.h>
//#include <stddef.h>


#define BGQ_PRECISION 64
#include "bgq_precisionselect.inc.c"


// Hopping_Matrix.c compatibility layer
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k) {
	const bool isOdd = (ieo!=0);
	bgq_spinorfield target = bgq_translate_spinorfield(l);
	bgq_spinorfield source = bgq_translate_spinorfield(k);

	bgq_HoppingMatrix_double(isOdd, target, source, g_gaugefield, hm_prefetchexplicit | hm_nokamul);
}


void Hopping_Matrix_nocom(const int ieo, spinor * const l, spinor * const k) {
	const bool isOdd = (ieo!=0);
	bgq_spinorfield target = bgq_translate_spinorfield(l);
	bgq_spinorfield source = bgq_translate_spinorfield(k);

	bgq_HoppingMatrix_double(isOdd, target, source, g_gaugefield, hm_nocom | hm_prefetchexplicit | hm_nokamul);
}


