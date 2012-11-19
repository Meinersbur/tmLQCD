

#include "bgq_HoppingMatrix.h"
#include <bgq_utils.h>

#define PRECISION double

void bgq_HoppingMatrix_worker(void * restrict arg, size_t tid, size_t threads) {
	bool const kamul = false;
	bool const readFulllayout = false;

#define BGQ_HOPPINGMATRIXWORKER_INC_ 1
#include "bgq_HoppingMatrixWorker.inc.c"
}


