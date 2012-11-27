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


#define BGQ_SPINORFIELD_GENWORKER(name) \
	static void NAME2(name,readFull)(void *arg_untyped, size_t tid, size_t threads) { \
		name(arg_untyped, tid, threads, false, false, false, false); \
	} \
	static void NAME2(name,readWeyl)(void *arg_untyped, size_t tid, size_t threads) { \
		name(arg_untyped, tid, threads, true, false, false, false); \
	} \
	static void NAME2(name,readFullSloppy)(void *arg_untyped, size_t tid, size_t threads) { \
		name(arg_untyped, tid, threads, false, true, false, false); \
	} \
	static void NAME2(name,readWeylSloppy)(void *arg_untyped, size_t tid, size_t threads) { \
		name(arg_untyped, tid, threads, true, true, false, false); \
	} \
	\
	bgq_worker_func NAME3(g,name,list)[BGQ_SPINORFIELD_LAYOUT_COUNT] = { \
		&NAME2(name,readFull), &NAME2(name,readWeyl), &NAME2(name,readFullSloppy), &NAME2(name,readWeylSloppy) \
	};


#define PRECISION_ISSLOPPY_double false
#define PRECISION_ISSLOPPY_float true
#define PRECISION_ISSLOPPY NAME2(PRECISION_ISSLOPPY,PRECISION)

#define PRECISION_SIZEOF_double 8
#define PRECISION_SIZEOF_float 8
#define PRECISION_SIZEOF NAME2(PRECISION_SIZEOF,PRECISION)


typedef struct {
	bool isOdd;
	bgq_weylfield_controlblock *field;
} bgq_unvectorize_workload;

typedef struct {
	bool isOdd;
	bgq_weylfield_controlblock *field;
	spinor *target;
} bgq_copyToLegacy_workload;


typedef struct {
	bool isOdd;
	spinor *source;
	bgq_weylfield_controlblock *target;
} bgq_copyFromLegacy_workload;

void bgq_HoppingMatrix_unvectorize_double(void *arg_untyped, size_t tid, size_t threads);
void bgq_HoppingMatrix_unvectorize_float(void *arg_untyped, size_t tid, size_t threads);
#define bgq_HoppingMatrix_unvectorize NAME2(bgq_HoppingMatrix_unvectorize,PRECISION)

void bgq_HoppingMatrix_worker_datamove_double(void *arg_untyped, size_t tid, size_t threads);
void bgq_HoppingMatrix_worker_datamove_float(void *arg_untyped, size_t tid, size_t threads);
#define bgq_HoppingMatrix_worker_datamove NAME2(bgq_HoppingMatrix_worker_datamove,PRECISION)

void bgq_HoppingMatrix_datamovet_worker_double(void *arg_untyped, size_t tid, size_t threads);
void bgq_HoppingMatrix_datamovet_worker_float(void *arg_untyped, size_t tid, size_t threads);
#define bgq_HoppingMatrix_datamovet_worker NAME2(bgq_HoppingMatrix_datamovet_worker,PRECISION)


#if 0
void bgq_copyToLegacy_worker_fulllayout_double(void *arg_untyped, size_t tid, size_t threads);
void bgq_copyToLegacy_worker_fulllayout_float(void *arg_untyped, size_t tid, size_t threads);
#define bgq_copyToLegacy_worker_fulllayout NAME2(bgq_copyToLegacy_worker_fulllayout,PRECISION)

void bgq_copyToLegacy_worker_weyllayout_double(void *arg_untyped, size_t tid, size_t threads);
void bgq_copyToLegacy_worker_weyllayout_float(void *arg_untyped, size_t tid, size_t threads);
#define bgq_copyToLegacy_worker_weyllayout NAME2(bgq_copyToLegacy_worker_weyllayout,PRECISION)
#endif


void bgq_copyFromLegacy_worker_double(void *arg_untyped, size_t tid, size_t threads);
void bgq_copyFromLegacy_worker_float(void *arg_untyped, size_t tid, size_t threads);
#define bgq_copyFromLegacy_worker NAME2(bgq_copyFromLegacy_worker,PRECISION)

#endif /* BGQ_WORKERS_H_ */
