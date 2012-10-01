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

	hm_noprefetchexplicit = 1 << 3,
	hm_noprefetchlist = 1 << 4,
	hm_noprefetchstream = 1 << 5,

	hm_noweylsend = 1 << 6,
	hm_nobody = 1 << 7,
	hm_nosurface = 1 << 8,

	hm_prefetchimplicitdisable = 1 << 9,
	hm_prefetchimplicitoptimistic = 2 << 9,
	hm_prefetchimplicitconfirmed = 3 << 9
} bgq_hmflags;

extern uint64_t fetch_depth[64];
extern uint64_t generate_depth[64];
extern L1P_Status_t l1plist_status[64];


#define bgq_HoppingMatrix NAME2(bgq_HoppingMatrix,PRECISION)
void bgq_HoppingMatrix_double(bool isOdd, bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bgq_hmflags opts);
void bgq_HoppingMatrix_float(bool isOdd, bgq_spinorfield_float targetfield, bgq_spinorfield_float spinorfield, bgq_gaugefield_float gaugefield, bgq_hmflags opts);

void bgq_HoppingMatrix_kamul_double(bool isOdd, bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bgq_hmflags opts);
void bgq_HoppingMatrix_kamul_float(bool isOdd, bgq_spinorfield_float targetfield, bgq_spinorfield_float spinorfield, bgq_gaugefield_float gaugefield, bgq_hmflags opts);
void bgq_HoppingMatrix_nokamul_double(bool isOdd, bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bgq_hmflags opts);
void bgq_HoppingMatrix_nokamul_float(bool isOdd, bgq_spinorfield_float targetfield, bgq_spinorfield_float spinorfield, bgq_gaugefield_float gaugefield, bgq_hmflags opts);

void bgq_HoppingMatrix_kamul_double_nodcbt(bool isOdd, bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bgq_hmflags opts);
void bgq_HoppingMatrix_kamul_float_nodcbt(bool isOdd, bgq_spinorfield_float targetfield, bgq_spinorfield_float spinorfield, bgq_gaugefield_float gaugefield, bgq_hmflags opts);
void bgq_HoppingMatrix_nokamul_double_nodcbt(bool isOdd, bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bgq_hmflags opts);
void bgq_HoppingMatrix_nokamul_float_nodcbt(bool isOdd, bgq_spinorfield_float targetfield, bgq_spinorfield_float spinorfield, bgq_gaugefield_float gaugefield, bgq_hmflags opts);


#endif /* BGQ_HOPPINGMATRIX_H_ */
