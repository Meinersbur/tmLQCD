/*
 * bgq_HoppingMatrix_base.inc.c
 *
 *  Created on: Oct 5, 2012
 *      Author: meinersbur
 */

#undef bgq_HoppingMatrix
#define bgq_HoppingMatrix_impl NAME5(bgq_HoppingMatrix,KAMUL_NAME,PRECISION,DCBT_NAME,ODDNESS_NAME)

void bgq_HoppingMatrix_impl(bool isOdd, bgq_spinorfield targetfield, bgq_spinorfield spinorfield, bgq_gaugefield gaugefield, bgq_hmflags opts) {
#if !defined(ODDNESS)
#elif ODDNESS
	assert(!isOdd);
	isOdd = true;
#define isOdd true
#else
	assert(isOdd);
	isOdd = false;
#define isOdd false
#endif
	const bool nocom = opts & hm_nocom;
	const bool nooverlap = opts & hm_nooverlap;
#if !defined(BGQ_HM_NOKAMUL)
	const bool nokamul = opts & hm_nokamul;
#elif BGQ_HM_NOKAMUL
	assert(opts & hm_nokamul);
	const bool nokamul = true;
#else
	assert(!(opts & hm_nokamul));
	const bool nokamul = false;
#endif
	const bool noprefetchlist = opts & hm_noprefetchlist;
	const bool noprefetchstream = opts & hm_noprefetchstream;
#if !defined(BGQ_DCBT_DISABLE)
	const bool noprefetchexplicit = opts & hm_noprefetchexplicit;
#elif BGQ_DCBT_DISABLE
	assert(opts & hm_noprefetchexplicit);
	const bool noprefetchexplicit = true;
#else
	assert(!(opts & hm_noprefetchexplicit));
	const bool noprefetchexplicit = false;
#endif
	const bool noweylsend = opts & hm_noweylsend;
	const bool nobody = opts & hm_nobody;
	const bool nosurface = opts & hm_nosurface;

	#define BGQ_HM_NOFUNC 1
	#define BGQ_HM_ZLINE_NOFUNC 1
	#define BGQ_HM_SITE_NOFUNC 1
	#define BGQ_HM_DIR_NOFUNC 1
	#include "bgq_HoppingMatrix.inc.c"
	#undef BGQ_HM_NOFUNC
	#undef BGQ_HM_ZLINE_NOFUNC
	#undef BGQ_HM_SITE_NOFUNC
	#undef BGQ_HM_DIR_NOFUNC
}

#ifdef ODDNESS
#undef isOdd
#endif

