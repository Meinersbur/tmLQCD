/*
 * bgq_workers.h
 *
 *  Created on: Nov 19, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_WORKERS_H_
#define BGQ_WORKERS_H_

#include "bgq_spinorfield.h"
#include "bgq_dispatch.h"

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
	static void NAME2(name,readLegacy)(void *arg_untyped, size_t tid, size_t threads) { \
		name(arg_untyped, tid, threads, false, false, false, true); \
	} \
	\
	bgq_worker_func NAME3(g,name,list)[BGQ_SPINORFIELD_LAYOUT_COUNT] = { \
		&NAME2(name,readFull), &NAME2(name,readWeyl), &NAME2(name,readFullSloppy), &NAME2(name,readWeylSloppy), NULL,NULL,NULL,NULL, NAME2(name,readLegacy) \
	};


#define PRECISION_WEYLLAYOUT_double ly_weyl_double
#define PRECISION_WEYLLAYOUT_float ly_weyl_float
#define PRECISION_WEYLLAYOUT NAME2(PRECISION_WEYLLAYOUT,PRECISION)

typedef struct {
	bool isOdd;
	bgq_weylfield_controlblock *field;
	bgq_hmflags opts;
} bgq_unvectorize_workload;
#if 0
typedef struct {
	bool isOdd;
	bgq_weylfield_controlblock *field;
	spinor *target;
} bgq_copyToLegacy_workload;
#endif
typedef struct {
	bool isOdd;
	spinor *source;
	bgq_weylfield_controlblock *target;
} bgq_copyFromLegacy_workload;

typedef struct {
	bgq_weylfield_controlblock *field;
	bgq_spinor_vec_double *target_double;
	bgq_spinor_vec_float *target_float;
	spinor *target_legacy;
	bool isOdd;
} bgq_spinorfield_rewrite_work;
#define BGQ_TARGETFIELD NAME2(target,PRECISION)

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

extern bgq_worker_func g_bgq_spinorfield_rewrite_worker_double_list[BGQ_SPINORFIELD_LAYOUT_COUNT];
extern bgq_worker_func g_bgq_spinorfield_rewrite_worker_float_list[BGQ_SPINORFIELD_LAYOUT_COUNT];
#define bgq_spinorfield_rewrite_worker NAME2(bgq_spinorfield_rewrite_worker,PRECISION)

void bgq_spinorfield_hmfull_writeToSendbuf_double(void *arg_untyped, size_t tid, size_t threads);
void bgq_spinorfield_hmfull_writeToSendbuf_float(void *arg_untyped, size_t tid, size_t threads);
#define bgq_spinorfield_hmfull_writeToSendbuf NAME2(bgq_spinorfield_hmfull_writeToSendbuf,PRECISION)

#endif /* BGQ_WORKERS_H_ */
