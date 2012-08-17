/*
 * bgq_HoppingMatrix_float.c
 *
 *  Created on: Aug 16, 2012
 *      Author: meinersbur
 */

#include "bgq_HoppingMatrix.h"
#include "bgq.h"
#include "../boundary.h"


#define BGQ_PRECISION 32
#include "bgq_precisionselect.inc.c"

#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wunused-variable"

void bgq_HoppingMatrix_float(bool isOdd, bgq_spinorfield targetfield, bgq_spinorfield spinorfield, bgq_gaugefield gaugefield, bgq_hmflags opts) {
	const bool nocom = opts & hm_nocom;
	const bool nooverlap = opts & hm_nooverlap;
	const bool nokamul = opts & hm_nokamul;
	const bool prefetchlist = opts & hm_prefetchlist;
	const bool prefetchstream = opts & hm_prefetchstream;
	const bool prefetchexplicit = opts & hm_prefetchexplicit;

	#define BGQ_HM_NOFUNC 1
	#define BGQ_HM_ZLINE_NOFUNC 1
	#define BGQ_HM_SITE_NOFUNC 1
	#define BGQ_HM_DIR_NOFUNC 1
	#include "bgq_HoppingMatrix.inc.c"
}
