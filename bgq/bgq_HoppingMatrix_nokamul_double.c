
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

#if !defined(BGQ_HM_NOKAMUL) || (BGQ_HM_NOKAMUL==1)
#define BGQ_HM_NOKAMUL 1

#define KAMUL_NAME nokamul
#include "bgq_HoppingMatrix_dcbt.inc.c"
#undef KAMUL_NAME

#endif




