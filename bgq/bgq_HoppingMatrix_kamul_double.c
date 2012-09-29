
#if HAVE_CONFIG_H
#include "config.h"
#endif
#include "bgq_HoppingMatrix.h"
#include "bgq.h"
#include "../boundary.h"


#define BGQ_PRECISION 64
#include "bgq_precisionselect.inc.c"

#ifndef XLC
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#if !defined(BGQ_HM_NOKAMUL) || (BGQ_HM_NOKAMUL==0)
#define BGQ_HM_NOKAMUL 0

void bgq_HoppingMatrix_kamul_double(bool isOdd, bgq_spinorfield targetfield, bgq_spinorfield spinorfield, bgq_gaugefield gaugefield, bgq_hmflags opts) {
	const bool nocom = opts & hm_nocom;
	const bool nooverlap = opts & hm_nooverlap;
	const bool nokamul = false;
	const bool prefetchlist = opts & hm_prefetchlist;
	const bool prefetchstream = opts & hm_prefetchstream;
	const bool prefetchexplicit = opts & hm_prefetchexplicit;
	const bool noweylsend = opts & hm_noweylsend;
	const bool nobody = opts & hm_nobody;
	const bool nosurface = opts & hm_nosurface;
	const bool nol1plist = opts & hm_nol1plist;

	#define BGQ_HM_NOFUNC 1
	#define BGQ_HM_ZLINE_NOFUNC 1
	#define BGQ_HM_SITE_NOFUNC 1
	#define BGQ_HM_DIR_NOFUNC 1
	#include "bgq_HoppingMatrix.inc.c"
}

#endif

#if 0
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bgq_field.h"
#include "bgq.h"
#include "global.h"
#include "boundary.h"
#include <mpi.h>
#include <omp.h>
#include <stddef.h>

#if BGQ_PREFETCH_LIST
#include <l1p/pprefetch.h>
#endif
#if BGQ_PREFETCH_LIST
#include <l1p/sprefetch.h>
#endif



#define BGQ_HOPPINGMATRIX_C_
#include "bgq_HoppingMatrix.h"







#define BGQ_HM_BORDER_NOFUNC 1
#define BGQ_HM_BORDERDIST_NOFUNC 1

//#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
//#pragma GCC diagnostic ignored "-Wunused-variable"


#define HoppingMatrix bgq_HoppingMatrix_double
#include "bgq_HoppingMatrix.inc.c"
#undef HoppingMatrix
#endif








