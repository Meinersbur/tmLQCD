/*
 * bgq_HoppingMatrix.h
 *
 *  Created on: Aug 1, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_HOPPINGMATRIX_H_
#define BGQ_HOPPINGMATRIX_H_

#include "bgq_field_double.h"
#include "bgq_field_float.h"

typedef enum {
	hm_nocom = 1 << 0,
	hm_nooverlap = 1 << 1,
	hm_nokamul = 1 << 2,
	hm_prefetchlist = 1 << 3,
	hm_prefetchstream = 1 << 4,
	hm_prefetchexplicit = 1 << 5,
} bgq_hmflags;

void bgq_HoppingMatrix_double(bool isOdd, bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bgq_hmflags opts);
void bgq_HoppingMatrix_float(bool isOdd, bgq_spinorfield_float targetfield, bgq_spinorfield_float spinorfield, bgq_gaugefield_float gaugefield, bgq_hmflags opts);

#endif /* BGQ_HOPPINGMATRIX_H_ */
