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

#ifndef REDUCTION_NAME
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
#endif





#ifndef REDUCTION_EXTRAPARMS
#define REDUCTION_EXTRAPARMS
#endif
#ifndef REDUCTION_EXTRAARGS
#define REDUCTION_EXTRAARGS
#endif

#define REDUCTION_EXTRAPARMLIST IF_EMPTY_INPAREN((),(,),REDUCTION_EXTRAPARMS) REDUCTION_EXTRAPARMS
#define REDUCTION_EXTRAARGLIST IF_EMPTY_INPAREN((),(,),REDUCTION_EXTRAARGS) REDUCTION_EXTRAARGS

#if REDUCTION_ARGFIELDS>=1
#define IF1ARG(...) __VA_ARGS__
#else
#define IF1ARG(...)
#endif

#if REDUCTION_ARGFIELDS>=2
#define IF2ARG(...) __VA_ARGS__
#else
#define IF2ARG(...)
#endif


#if !ISEMPTY(REDUCTION_EXTRAPARMS)
typedef struct {
	PREPROCESSOR_FOREACH( ,;,;, ,IDENTITY,REDUCTION_EXTRAPARMS)
} bgq_reduction_extraparm_t;
#endif

typedef struct {
	bool isOdd;
#if REDUCTION_ARGFIELDS>=1
	bgq_weylfield_controlblock *argfield1;
#endif
#if REDUCTION_ARGFIELDS>=2
	bgq_weylfield_controlblock *argfield2;
#endif
#if !ISEMPTY(REDUCTION_EXTRAPARMS)
	bgq_reduction_extraparm_t extra;
#endif
	REDUCTION_RETURNTYPE((threadresult[64]));
	bool hasResult[64];
} NAME2(REDUCTION_NAME,args_t);

#define REDUCTION_EXTRA_ASSIGN(varname) (varname) = ((arg->extra).varname);
#define REDUCTION_REF(varname) (&varname)
#define REDUCTION_SITE_DECL(varname)

static inline void NAME2(REDUCTION_NAME,worker)(void *arg_untyped, size_t tid, size_t threads IF1ARG(, bool readFulllayout1) IF2ARG(, bool readFulllayout2)) {
	NAME2(REDUCTION_NAME,args_t) *arg = (NAME2(REDUCTION_NAME,args_t)*)arg_untyped;
	bool isOdd = arg->isOdd;
#if REDUCTION_ARGFIELDS>=1
	bgq_weylfield_controlblock *argfield1 = arg->argfield1;
#endif
#if REDUCTION_ARGFIELDS>=2
	bgq_weylfield_controlblock *argfield2 = arg->argfield2;
#endif
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

#if REDUCTION_ARGFIELDS>=1
		bgq_su3_spinor_decl(spinor1);
		if (readFulllayout1) {
			bgq_spinorsite *fulladdr1 = &argfield1->sec_fullspinor[ic];
			bgq_HoppingMatrix_loadFulllayout(spinor1, fulladdr1,t1,t2,x,y,z);
		} else {
			bgq_weylsite *weyladdr1 = &argfield1->sec_collapsed[ic];
			bgq_HoppingMatrix_loadWeyllayout(spinor1, weyladdr1,t1,t2,x,y,z);
		}
#endif

#if REDUCTION_ARGFIELDS>=2
		bgq_su3_spinor_decl(spinor2);
		if (readFulllayout2) {
			bgq_spinorsite *fulladdr2 = &argfield2->sec_fullspinor[ic];
			bgq_HoppingMatrix_loadFulllayout(spinor2, fulladdr2,t1,t2,x,y,z);
		} else {
			bgq_weylsite *weyladdr2 = &argfield2->sec_collapsed[ic];
			bgq_HoppingMatrix_loadWeyllayout(spinor2, weyladdr2,t1,t2,x,y,z);
		}
#endif

		REDUCTION_RETURNTYPE(siteresult);
		siteresult = REDUCTION_SITEREDUCEFUNC(
#if REDUCTION_ARGFIELDS>=1
				bgq_su3_spinor_vars(spinor1)
#endif
#if REDUCTION_ARGFIELDS>=2
				, bgq_su3_spinor_vars(spinor2)
#endif
				REDUCTION_EXTRAARGLIST, tv, t1, t2, x, y, z);

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

#if REDUCTION_ARGFIELDS>=1
		bgq_su3_spinor_decl(spinor1);
		if (readFulllayout1) {
			bgq_spinorsite *fulladdr1 = &argfield1->sec_fullspinor[ic];
			bgq_HoppingMatrix_loadFulllayout(spinor1, fulladdr1,t1,t2,x,y,z);
		} else {
			bgq_weylsite *weyladdr1 = &argfield1->sec_collapsed[ic];
			bgq_HoppingMatrix_loadWeyllayout(spinor1, weyladdr1,t1,t2,x,y,z);
		}
#endif

#if REDUCTION_ARGFIELDS>=2
		bgq_su3_spinor_decl(spinor2);
		if (readFulllayout2) {
			bgq_spinorsite *fulladdr2 = &argfield2->sec_fullspinor[ic];
			bgq_HoppingMatrix_loadFulllayout(spinor2, fulladdr2,t1,t2,x,y,z);
		} else {
			bgq_weylsite *weyladdr2 = &argfield2->sec_collapsed[ic];
			bgq_HoppingMatrix_loadWeyllayout(spinor2, weyladdr2,t1,t2,x,y,z);
		}
#endif

		REDUCTION_RETURNTYPE(siteresult);
		siteresult = REDUCTION_SITEREDUCEFUNC(
#if REDUCTION_ARGFIELDS>=1
				bgq_su3_spinor_vars(spinor1)
#endif
#if REDUCTION_ARGFIELDS>=2
				, bgq_su3_spinor_vars(spinor2)
#endif
				REDUCTION_EXTRAARGLIST, tv, t1, t2, x, y, z);


		threadresult = REDUCTION_COMBINEFUNC(threadresult, siteresult);
	}

	arg->hasResult[tid] = hasResult;
	arg->threadresult[tid] = threadresult;
}


#if REDUCTION_ARGFIELDS==0
static bgq_worker_func NAME2(REDUCTION_NAME,worker_funcs) = NAME2(REDUCTION_NAME,worker);
#endif


#if REDUCTION_ARGFIELDS==1
static void NAME3(REDUCTION_NAME,worker,readWeyllayout1)(void *arg_untyped, size_t tid, size_t threads) {
	NAME2(REDUCTION_NAME,worker)(arg_untyped, tid, threads, false);
}
static void NAME3(REDUCTION_NAME,worker,readFulllayout1)(void *arg_untyped, size_t tid, size_t threads) {
	NAME2(REDUCTION_NAME,worker)(arg_untyped, tid, threads, true);
}

static bgq_worker_func NAME2(REDUCTION_NAME,worker_funcs)[2] = {
		&NAME3(REDUCTION_NAME,worker,readWeyllayout1), &NAME3(REDUCTION_NAME,worker,readFulllayout1)
};
#endif


#if REDUCTION_ARGFIELDS==2
static void NAME4(REDUCTION_NAME,worker,readWeyllayout1,readWeyllayout2)(void *arg_untyped, size_t tid, size_t threads) {
	NAME2(REDUCTION_NAME,worker)(arg_untyped, tid, threads, false, false);
}
static void NAME4(REDUCTION_NAME,worker,readWeyllayout1,readFulllayout2)(void *arg_untyped, size_t tid, size_t threads) {
	NAME2(REDUCTION_NAME,worker)(arg_untyped, tid, threads, false, true);
}
static void NAME4(REDUCTION_NAME,worker,readFulllayout1,readWeyllayout2)(void *arg_untyped, size_t tid, size_t threads) {
	NAME2(REDUCTION_NAME,worker)(arg_untyped, tid, threads, true, false);
}
static void NAME4(REDUCTION_NAME,worker,readFulllayout1,readFulllayout2)(void *arg_untyped, size_t tid, size_t threads) {
	NAME2(REDUCTION_NAME,worker)(arg_untyped, tid, threads, true, true);
}

static bgq_worker_func NAME2(REDUCTION_NAME,worker_funcs)[2][2] = {
		{&NAME4(REDUCTION_NAME,worker,readWeyllayout1,readWeyllayout2), &NAME4(REDUCTION_NAME,worker,readWeyllayout1,readFulllayout2)},
		{&NAME4(REDUCTION_NAME,worker,readFulllayout1,readWeyllayout2), &NAME4(REDUCTION_NAME,worker,readFulllayout1,readFulllayout2)}
};
#endif


static inline REDUCTION_RETURNTYPE(REDUCTION_NAME)(
#if REDUCTION_ARGFIELDS>=1
		bgq_weylfield_controlblock *argfield1
#endif
#if REDUCTION_ARGFIELDS>=2
		, bgq_weylfield_controlblock *argfield2
#endif
		REDUCTION_EXTRAPARMLIST) {

	bool isOdd = argfield1->isOdd;


#if REDUCTION_ARGFIELDS>=1
	bool useFulllayout1 = argfield1->hasFullspinorData;
	bool useWeyllayout1 = !useFulllayout1;
	bgq_spinorfield_setup(argfield1, isOdd, useFulllayout1, false, useWeyllayout1, false);
#endif

#if REDUCTION_ARGFIELDS>=2
	bool useFulllayout2 = argfield2->hasFullspinorData;
	bool useWeyllayout2 = !useFulllayout2;
	bgq_spinorfield_setup(argfield2, isOdd, useFulllayout2, false, useWeyllayout2, false);
#endif

	bgq_master_sync();
	static NAME2(REDUCTION_NAME,args_t) operator_args;
	NAME2(REDUCTION_NAME,args_t) call_args = {
			.isOdd = isOdd,
#if REDUCTION_ARGFIELDS>=1
			.argfield1 = argfield1,
#endif
#if REDUCTION_ARGFIELDS>=2
			.argfield2 = argfield2,
#endif
#if !ISEMPTY(REDUCTION_EXTRAPARMS)
			.extra = { REDUCTION_EXTRAARGS },
#endif
			.hasResult = {0}
	};
	operator_args = call_args;
	bgq_master_call(NAME2(REDUCTION_NAME,worker_funcs)
#if REDUCTION_ARGFIELDS>=1
			[useFulllayout1]
#endif
#if REDUCTION_ARGFIELDS>=2
			 [useFulllayout2]
#endif
			  , &operator_args);

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


#undef IF1ARG
#undef IF2ARG

#undef REDUCTION_NAME
#undef REDUCTION_ARGFIELDS
#undef REDUCTION_RETURNTYPE
#undef REDUCTION_SITEREDUCEFUNC
#undef REDUCTION_COMBINEFUNC


