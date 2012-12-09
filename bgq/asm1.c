
#include "bgq_HoppingMatrix.h"

#include "bgq_field.h"
#include "bgq_spinorfield.h"
#include "bgq_qpx.h"
#include "bgq_dispatch.h"
#include "bgq_comm.h"
#include "bgq_workers.h"
#include "bgq_gaugefield.h"

#include "../boundary.h"
#include "../update_backward_gauge.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

void bgq_HoppingMatrix_nokamul_worker_readWeyllayout_double(void *arg, size_t tid, size_t threads) {
	//bgq_HoppingMatrix_worker(arg,tid,threads,false,false);

	bool const kamul = false;
	bool const readFulllayout = false;

#define PRECISION double
#define BGQ_HOPPINGMATRIXWORKER_INC_ 1
#include "bgq_HoppingMatrixWorker.inc.c"
#undef PRECISION
}
