/*
 * bgq_workers.h
 *
 *  Created on: Nov 19, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_WORKERS_H_
#define BGQ_WORKERS_H_

#include "bgq_spinorfield.h"

#include <string.h>
#include <stdbool.h>

typedef struct {
	bool isOdd;
	bgq_weylfield_controlblock *field;
} bgq_unvectorize_workload;

void bgq_HoppingMatrix_unvectorize_double(void *arg_untyped, size_t tid, size_t threads);
void bgq_HoppingMatrix_unvectorize_float(void *arg_untyped, size_t tid, size_t threads);
#define bgq_HoppingMatrix_unvectorize NAME2(bgq_HoppingMatrix_unvectorize,PRECISION)

void bgq_HoppingMatrix_worker_datamove_double(void *arg_untyped, size_t tid, size_t threads);
void bgq_HoppingMatrix_worker_datamove_float(void *arg_untyped, size_t tid, size_t threads);
#define bgq_HoppingMatrix_worker_datamove NAME2(bgq_HoppingMatrix_worker_datamove,PRECISION)

void bgq_HoppingMatrix_datamovet_worker_double(void *arg_untyped, size_t tid, size_t threads);
void bgq_HoppingMatrix_datamovet_worker_float(void *arg_untyped, size_t tid, size_t threads);
#define bgq_HoppingMatrix_datamovet_worker NAME2(bgq_HoppingMatrix_datamovet_worker,PRECISION)

#endif /* BGQ_WORKERS_H_ */
