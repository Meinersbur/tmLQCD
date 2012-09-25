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


#if BGQ_REPLACE
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
#endif

void bgq_HoppingMatrix_double(bool isOdd, bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bgq_hmflags opts) {
	const bool nokamul = opts & hm_nokamul;

	if (nokamul) {
		bgq_HoppingMatrix_nokamul_double(isOdd, targetfield, spinorfield, gaugefield, opts);
	} else {
		bgq_HoppingMatrix_kamul_double(isOdd, targetfield, spinorfield, gaugefield, opts);
	}

	bgq_spinorfield_setOdd_double(targetfield, isOdd, true);
	bgq_spinorfield_setOdd_double(spinorfield, isOdd, false);
}

void bgq_HoppingMatrix_float(bool isOdd, bgq_spinorfield_float targetfield, bgq_spinorfield_float spinorfield, bgq_gaugefield_float gaugefield, bgq_hmflags opts) {
	const bool nokamul = opts & hm_nokamul;

	if (nokamul) {
		bgq_HoppingMatrix_nokamul_float(isOdd, targetfield, spinorfield, gaugefield, opts);
	} else {
		bgq_HoppingMatrix_kamul_float(isOdd, targetfield, spinorfield, gaugefield, opts);
	}

	bgq_spinorfield_setOdd_float(targetfield, isOdd, true);
	bgq_spinorfield_setOdd_float(spinorfield, isOdd, false);
}

