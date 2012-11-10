/*
 * bgq_reduction.inc.c
 *
 *  Created on: Nov 10, 2012
 *      Author: meinersbur
 */

#include "bgq_spinorfield.h"
#include "bgq_dispatch.h"
#include "bgq_HoppingMatrix.h"
#include "bgq_qpx.h"
#include "bgq_utils.h"

#include <stdbool.h>


#define REDUCTION_NAME bgq_reduce
#define REDUCTION_ARGFIELDS 2
#define REDUCTION_EXTRAPARMS double factor
#define REDUCTION_EXTRAARGS factor

#define REDUCTION_RETURNTYPE(varname) double (varname) /* remember this can be a struct */

#define REDUCTION_SITEREDUCEFUNC bgq_sitereduce
#define REDUCTION_COMBINEFUNC bgq_reduction_combine  /* must be associative (commutative?) */



static inline double bgq_sitereduce(bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), double factor, ucoord tv, ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z) {
return 0;
}

static inline double bgq_reduction_combine(double sum1, double sum2) {
return 0;
}


#define REDUCTION_EXTRAPARMLIST IF_EMPTY_INPAREN((),(,),OPERATOR_EXTRAPARMS) REDUCTION_EXTRAPARMS
#define REDUCTION_EXTRAARGLIST IF_EMPTY_INPAREN((),(,),OPERATOR_EXTRAARGS) REDUCTION_EXTRAARGS




typedef struct {
	PREPROCESSOR_FOREACH( ,;,;, ,IDENTITY,REDUCTION_EXTRAPARMS)
} bgq_reduction_extraparm_t;

typedef struct {
	PREPROCESSOR_FOREACH( ,;,;, ,IDENTITY,REDUCTION_EXTRAPARMS)
} bgq_reduction_reduceparm_t;

typedef struct {
	bool isOdd;
	bgq_weylfield_controlblock *argfield1;
	bgq_weylfield_controlblock *argfield2;
	bgq_reduction_extraparm_t extra;
	REDUCTION_RETURNTYPE((threadresult[64]));
	bool hasResult[64];
} bgq_reduction_args_t;

#define REDUCTION_EXTRA_ASSIGN(varname) (varname) = ((arg->extra).varname);
#define REDUCTION_REF(varname) (&varname)
#define REDUCTION_SITE_DECL(varname)

static inline void bgq_reduction_worker(void *arg_untyped, size_t tid, size_t threads, bool readFulllayout1, bool readFulllayout2) {
	bgq_reduction_args_t *arg = (bgq_reduction_args_t*)arg_untyped;
	bool isOdd = arg->isOdd;
	bgq_weylfield_controlblock *argfield1 = arg->argfield1;
	bgq_weylfield_controlblock *argfield2 = arg->argfield2;
	PREPROCESSOR_FOREACH( ,;,;, ,IDENTITY,REDUCTION_EXTRAPARMS)
	PREPROCESSOR_FOREACH(,,,,REDUCTION_EXTRA_ASSIGN,REDUCTION_EXTRAARGS)

	size_t workload = PHYSICAL_VOLUME;
	size_t threadload = (workload+threads-1)/threads;
	ucoord beginj = tid*threadload;
	ucoord endj = min_sizet((tid+1)*threadload,workload);

	REDUCTION_RETURNTYPE(threadresult);
	bool hasResult = false;
	ucoord ic = beginj;
	if (ic < endj) {
		ucoord ih = bgq_collapsed2halfvolume(isOdd, ic);
		ucoord tv = bgq_halfvolume2tv(ih);
		ucoord t1 = bgq_halfvolume2t1(isOdd, ih);
		ucoord t2 = bgq_halfvolume2t2(isOdd, ih);
		ucoord x = bgq_halfvolume2x(ih);
		ucoord y = bgq_halfvolume2y(ih);
		ucoord z = bgq_halfvolume2z(ih);

		bgq_su3_spinor_decl(spinor1);
		if (readFulllayout1) {
			bgq_spinorsite *fulladdr1 = &argfield1->sec_fullspinor[ic];
			bgq_HoppingMatrix_loadFulllayout(spinor1, fulladdr1,t1,t2,x,y,z);
		} else {
			bgq_weylsite *weyladdr1 = &argfield1->sec_collapsed[ic];
			bgq_HoppingMatrix_loadWeyllayout(spinor1, weyladdr1,t1,t2,x,y,z);
		}

		bgq_su3_spinor_decl(spinor2);
		if (readFulllayout2) {
			bgq_spinorsite *fulladdr2 = &argfield2->sec_fullspinor[ic];
			bgq_HoppingMatrix_loadFulllayout(spinor2, fulladdr2,t1,t2,x,y,z);
		} else {
			bgq_weylsite *weyladdr2 = &argfield2->sec_collapsed[ic];
			bgq_HoppingMatrix_loadWeyllayout(spinor2, weyladdr2,t1,t2,x,y,z);
		}

		REDUCTION_RETURNTYPE(siteresult);
		siteresult = REDUCTION_SITEREDUCEFUNC(bgq_su3_spinor_vars(spinor1), bgq_su3_spinor_vars(spinor2) REDUCTION_EXTRAARGLIST, tv, t1, t2, x, y, z);

		threadresult = siteresult;
		hasResult = true;
	}

	ic += 1;
	for (; ic < endj; ic+=1) {
		ucoord ih = bgq_collapsed2halfvolume(isOdd, ic);
		ucoord tv = bgq_halfvolume2tv(ih);
		ucoord t1 = bgq_halfvolume2t1(isOdd, ih);
		ucoord t2 = bgq_halfvolume2t2(isOdd, ih);
		ucoord x = bgq_halfvolume2x(ih);
		ucoord y = bgq_halfvolume2y(ih);
		ucoord z = bgq_halfvolume2z(ih);

		bgq_su3_spinor_decl(spinor1);
		if (readFulllayout1) {
			bgq_spinorsite *fulladdr1 = &argfield1->sec_fullspinor[ic];
			bgq_HoppingMatrix_loadFulllayout(spinor1, fulladdr1,t1,t2,x,y,z);
		} else {
			bgq_weylsite *weyladdr1 = &argfield1->sec_collapsed[ic];
			bgq_HoppingMatrix_loadWeyllayout(spinor1, weyladdr1,t1,t2,x,y,z);
		}

		bgq_su3_spinor_decl(spinor2);
		if (readFulllayout2) {
			bgq_spinorsite *fulladdr2 = &argfield2->sec_fullspinor[ic];
			bgq_HoppingMatrix_loadFulllayout(spinor2, fulladdr2,t1,t2,x,y,z);
		} else {
			bgq_weylsite *weyladdr2 = &argfield2->sec_collapsed[ic];
			bgq_HoppingMatrix_loadWeyllayout(spinor2, weyladdr2,t1,t2,x,y,z);
		}

		REDUCTION_RETURNTYPE(siteresult);
		siteresult = REDUCTION_SITEREDUCEFUNC(bgq_su3_spinor_vars(spinor1), bgq_su3_spinor_vars(spinor2) REDUCTION_EXTRAARGLIST, tv, t1, t2, x, y, z);


		threadresult = REDUCTION_COMBINEFUNC(threadresult, siteresult);
	}

	arg->hasResult[tid] = hasResult;
	arg->threadresult[tid] = threadresult;
}


static void bgq_reduction_worker_readWeyllayout1_readWeyllayout2(void *arg_untyped, size_t tid, size_t threads) {
	bgq_reduction_worker(arg_untyped, tid, threads, false, false);
}
static void bgq_reduction_worker_readWeyllayout1_readFulllayout2(void *arg_untyped, size_t tid, size_t threads) {
	bgq_reduction_worker(arg_untyped, tid, threads, false, true);
}
static void bgq_reduction_worker_readFulllayout1_readWeyllayout2(void *arg_untyped, size_t tid, size_t threads) {
	bgq_reduction_worker(arg_untyped, tid, threads, true, false);
}
static void bgq_reduction_worker_readFulllayout1_readFulllayout2(void *arg_untyped, size_t tid, size_t threads) {
	bgq_reduction_worker(arg_untyped, tid, threads, true, true);
}

static bgq_worker_func bgq_worker_funcs[2][2] = {
		{&bgq_reduction_worker_readWeyllayout1_readWeyllayout2, &bgq_reduction_worker_readWeyllayout1_readFulllayout2},
		{&bgq_reduction_worker_readFulllayout1_readWeyllayout2, &bgq_reduction_worker_readFulllayout1_readFulllayout2}
};


REDUCTION_RETURNTYPE(bgq_reduction)(bgq_weylfield_controlblock *argfield1, bgq_weylfield_controlblock *argfield2 REDUCTION_EXTRAPARMLIST) {
	bool isOdd = argfield1->isOdd;

	bool useFulllayout1 = argfield1->hasFullspinorData;
	bool useWeyllayout1 = !useFulllayout1;
	bgq_spinorfield_setup(argfield1, isOdd, useFulllayout1, false, useWeyllayout1, false);

	bool useFulllayout2 = argfield2->hasFullspinorData;
	bool useWeyllayout2 = !useFulllayout2;
	bgq_spinorfield_setup(argfield2, isOdd, useFulllayout2, false, useWeyllayout2, false);

	bgq_master_sync();
	static bgq_reduction_args_t operator_args;
	bgq_reduction_args_t call_args = {
			.isOdd = isOdd,
			.argfield1 = argfield1,
			.argfield2 = argfield2,
			.extra = { REDUCTION_EXTRAARGS },
			.hasResult = {0}
	};
	operator_args = call_args;
	bgq_master_call(bgq_worker_funcs[useFulllayout1][useFulllayout2], &operator_args);

	int threads = omp_get_thread_num();
	size_t tid = 0;
	REDUCTION_RETURNTYPE(result);
	bgq_master_sync();
	for (; tid < threads;tid+=1) {
		if (operator_args.hasResult[tid]) {
			result = operator_args.threadresult[tid];
		}
	}
	assert(tid < threads);
	tid+=1;
	for (; tid < threads; tid+=1) {
		if (operator_args.hasResult[tid]) {
			REDUCTION_RETURNTYPE(threadresult);
			threadresult = operator_args.threadresult[tid];
			result = REDUCTION_COMBINEFUNC(result, threadresult);
		}
	}

	return result;
}



