/*
 * bgq_HoppingMatrix_double.c
 *
 *  Created on: Aug 16, 2012
 *      Author: meinersbur
 */

#if HAVE_CONFIG_H
#include "config.h"
#endif
#include "bgq_HoppingMatrix.h"
#include "bgq.h"
#include "../boundary.h"
#include "bgq_utils.h"


#define BGQ_PRECISION 64
#include "bgq_precisionselect.inc.c"

#ifndef XLC
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#if !defined(BGQ_HM_NOKAMUL) || (BGQ_HM_NOKAMUL==0)
#define BGQ_HM_NOKAMUL 0

#define KAMUL_NAME kamul
#include "bgq_HoppingMatrix_dcbt.inc.c"
#undef KAMUL_NAME

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








