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
// Make this .inc.c file compilable as stand-alone
#define REDUCTION_NAME bgq_reduce
#define REDUCTION_ARGFIELDS 2
#define REDUCTION_EXTRAPARMS double factor
#define REDUCTION_EXTRAARGS factor
#define REDUCTION_RETURNTYPE(varname) double (varname) /* remember this can be a struct */

#define REDUCTION_VARINIT bgq_initzero
#define REDUCTION_SITEREDUCEFUNC bgq_sitereduce
#define REDUCTION_COMBINEFUNC bgq_reduction_combine  /* must be associative (commutative?) */

static inline double bgq_initzero() {
	return 0;
}

static inline void bgq_sitereduce(double *sum, bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), double factor, ucoord ic) {
	*sum += 1;
}

static inline void bgq_reduction_combine(double *sum, double intermediateSum) {
	*sum+=intermediateSum;
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
#define REDUCTION_EXTRA_ASSIGN(varname) (varname) = ((arg->extra).varname);
#define REDUCTION_REF(varname) (&varname)
#define REDUCTION_SITE_DECL(varname)

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

#if ISEMPTY(REDUCTION_EXTRAPARMS)
#define IFEXTRA(...)
#else
#define IFEXTRA(...) __VA_ARGS__
#endif

#define REDUCTION_LAYOUTNAME_ly_full_double Fulllayout
#define REDUCTION_LAYOUTNAME_ly_weyl_double Weyllayout
#define REDUCTION_LAYOUTNAME_ly_full_float FulllayoutSloppy
#define REDUCTION_LAYOUTNAME_ly_weyl_float WeyllayoutSloppy
#define REDUCTION_LAYOUTNAME_0 REDUCTION_LAYOUTNAME_ly_full_double
#define REDUCTION_LAYOUTNAME_1 REDUCTION_LAYOUTNAME_ly_weyl_double
#define REDUCTION_LAYOUTNAME_2 REDUCTION_LAYOUTNAME_ly_full_float
#define REDUCTION_LAYOUTNAME_3 REDUCTION_LAYOUTNAME_ly_weyl_float


#if !ISEMPTY(REDUCTION_EXTRAPARMS)
typedef struct {
	PREPROCESSOR_FOREACH( ,;,;, ,IDENTITY,REDUCTION_EXTRAPARMS)
} bgq_reduction_extraparm_t;
#endif

typedef struct {
	bool isOdd;

	IF1ARG(bgq_weylfield_controlblock *argfield1;)
	IF2ARG(bgq_weylfield_controlblock *argfield2;)
	IFEXTRA(bgq_reduction_extraparm_t extra;)
	REDUCTION_RETURNTYPE((threadresult[64]));
} NAME2(REDUCTION_NAME,args_t);



static inline void NAME2(REDUCTION_NAME,worker)(void *arg_untyped, size_t tid, size_t threads IF1ARG(, bool readWeyllayout1, bool sloppy1, bool mul1) IF2ARG(, bool readWeyllayout2, bool sloppy2, bool mul2)) {
	NAME2(REDUCTION_NAME,args_t) *arg = arg_untyped;
	bool isOdd = arg->isOdd;

	IF1ARG(bgq_weylfield_controlblock *argfield1 = arg->argfield1;)
	IF2ARG(bgq_weylfield_controlblock *argfield2 = arg->argfield2;)
	PREPROCESSOR_FOREACH( ,;,;, ,IDENTITY, REDUCTION_EXTRAPARMS)
	PREPROCESSOR_FOREACH(,,,,REDUCTION_EXTRA_ASSIGN, REDUCTION_EXTRAARGS)

	size_t workload = PHYSICAL_VOLUME;
	size_t threadload = (workload+threads-1)/threads;
	ucoord beginj = tid*threadload;
	ucoord endj = min_sizet((tid+1)*threadload,workload);
	REDUCTION_RETURNTYPE(threadresult) = REDUCTION_VARINIT();
	for (ucoord ic = beginj; ic < endj; ic+=1) {
#ifndef NDEBUG
		ucoord ih = bgq_collapsed2halfvolume(isOdd, ic);
		ucoord tv = bgq_halfvolume2tv(ih);
		ucoord t1 = bgq_halfvolume2t1(isOdd, ih);
		ucoord t2 = bgq_halfvolume2t2(isOdd, ih);
		ucoord x = bgq_halfvolume2x(ih);
		ucoord y = bgq_halfvolume2y(ih);
		ucoord z = bgq_halfvolume2z(ih);
#endif

#if REDUCTION_ARGFIELDS>=1
		bgq_su3_spinor_decl(spinor1);
		bgq_spinorfield_readSpinor(&spinor1, argfield1, ic, readWeyllayout1, sloppy1, mul1);
#endif

#if REDUCTION_ARGFIELDS>=2
		bgq_su3_spinor_decl(spinor2);
		bgq_spinorfield_readSpinor(&spinor2, argfield1, ic, readWeyllayout2, sloppy2, mul2);
#endif

		REDUCTION_SITEREDUCEFUNC(&threadresult IF1ARG(, bgq_su3_spinor_vars(spinor1)) IF2ARG(, bgq_su3_spinor_vars(spinor2)) REDUCTION_EXTRAARGLIST, ic);
	}

	arg->threadresult[tid] = threadresult;
}


#if REDUCTION_ARGFIELDS==0
static bgq_worker_func NAME2(REDUCTION_NAME,worker_funcs) = NAME2(REDUCTION_NAME,worker);


#elif REDUCTION_ARGFIELDS==1
#define REDUCTION_WORKERNAME(layout) NAME3(REDUCTION_NAME,worker,NAME2(REDUCTION_LAYOUTNAME,layout))
#define REDUCTION_WORKERINST(layout1) \
		static void REDUCTION_WORKERNAME(layout1)(void *arg_untyped, size_t tid, size_t threads) { \
			NAME2(REDUCTION_NAME,worker)(arg_untyped, tid, threads, (layout1)&ly_weyl, (layout1)&ly_sloppy, (layout1)&ly_mul); \
		}

REDUCTION_WORKERINST(ly_full_double)
REDUCTION_WORKERINST(ly_weyl_double)
REDUCTION_WORKERINST(ly_full_float)
REDUCTION_WORKERINST(ly_weyl_float)

static bgq_worker_func NAME2(REDUCTION_NAME,worker_funcs)[BGQ_SPINORFIELD_LAYOUT_COUNT] =
	   /* ly_full_double */                   /* ly_weyl_double */                   /* ly_full_float */                  /* ly_weyl_float */
 	 { &REDUCTION_WORKERNAME(ly_full_double), &REDUCTION_WORKERNAME(ly_weyl_double), &REDUCTION_WORKERNAME(ly_full_float), &REDUCTION_WORKERNAME(ly_weyl_float) };


#elif REDUCTION_ARGFIELDS==2
#define REDUCTION_WORKERNAME(layout1, layout2) NAME4(REDUCTION_NAME,worker,NAME2(REDUCTION_LAYOUTNAME,layout1), NAME2(REDUCTION_LAYOUTNAME,layout2))
#define REDUCTION_WORKERINST(layout1,layout2) \
		static void REDUCTION_WORKERNAME(layout1,layout2)(void *arg_untyped, size_t tid, size_t threads) { \
			NAME2(REDUCTION_NAME,worker)(arg_untyped, tid, threads, (layout1)&ly_weyl, (layout1)&ly_sloppy, (layout1)&ly_mul, (layout2)&ly_weyl, (layout2)&ly_sloppy, (layout2)&ly_mul); \
		}

REDUCTION_WORKERINST(ly_full_double,ly_full_double)
REDUCTION_WORKERINST(ly_weyl_double,ly_full_double)
REDUCTION_WORKERINST(ly_full_double,ly_weyl_double)
REDUCTION_WORKERINST(ly_weyl_double,ly_weyl_double)

REDUCTION_WORKERINST(ly_full_float,ly_full_float)
REDUCTION_WORKERINST(ly_weyl_float,ly_full_float)
REDUCTION_WORKERINST(ly_full_float,ly_weyl_float)
REDUCTION_WORKERINST(ly_weyl_float,ly_weyl_float)

static bgq_worker_func NAME2(REDUCTION_NAME,worker_funcs)[BGQ_SPINORFIELD_LAYOUT_COUNT][BGQ_SPINORFIELD_LAYOUT_COUNT] = {
		                /* ly_full_double */                                  /* ly_weyl_double */                                  /* ly_full_float */                   /* ly_weyl_float */
/* ly_full_double */  { &REDUCTION_WORKERNAME(ly_full_double,ly_full_double), &REDUCTION_WORKERNAME(ly_full_double,ly_weyl_double), NULL,                                               NULL },
/* ly_weyl_double */  { &REDUCTION_WORKERNAME(ly_weyl_double,ly_full_double), &REDUCTION_WORKERNAME(ly_weyl_double,ly_weyl_double), NULL,                                               NULL },
/* ly_full_float */   { NULL,                                                NULL,                                                  &REDUCTION_WORKERNAME(ly_full_float,ly_full_float), &REDUCTION_WORKERNAME(ly_full_float,ly_weyl_float) },
/* ly_weyl_float */   { NULL,                                                NULL,                                                  &REDUCTION_WORKERNAME(ly_weyl_float,ly_full_float), &REDUCTION_WORKERNAME(ly_weyl_float,ly_weyl_float) }
};


#endif


static inline REDUCTION_RETURNTYPE(REDUCTION_NAME)(
	IF1ARG(bgq_weylfield_controlblock *argfield1)
	IF2ARG(, bgq_weylfield_controlblock *argfield2)
	REDUCTION_EXTRAPARMLIST) {

	bool isOdd = argfield1->isOdd;

#if REDUCTION_ARGFIELDS==0
	bgq_worker_func workerfunc = NAME2(REDUCTION_NAME,worker_funcs);
#elif REDUCTION_ARGFIELDS==1
	bgq_spinorfield_layout layout1 = bgq_spinorfield_prepareRead(argfield1, isOdd, true, true, true, true);
	if (!NAME2(REDUCTION_NAME,worker_funcs)[layout1])
		layout1 = bgq_spinorfield_prepareRead(argfield1, isOdd, false, true, false, false); // Force a rewrite
	bgq_worker_func workerfunc = NAME2(REDUCTION_NAME,worker_funcs)[layout1];
#elif REDUCTION_ARGFIELDS==2
	bgq_spinorfield_layout layout1 = bgq_spinorfield_prepareRead(argfield1, isOdd, true, true, true, true);
	bgq_spinorfield_layout layout2 = bgq_spinorfield_prepareRead(argfield2, isOdd, true, true, true, true);
	if (!NAME2(REDUCTION_NAME,worker_funcs)[layout1][layout2]) {
		if (NAME2(REDUCTION_NAME,worker_funcs)[ly_full_double][layout2])
			layout1 = bgq_spinorfield_prepareRead(argfield1, isOdd, false, true, false, false);
		else if (NAME2(REDUCTION_NAME,worker_funcs)[layout1][ly_full_double])
			layout2 = bgq_spinorfield_prepareRead(argfield2, isOdd, false, true, false, false);
		else {
			// Force rewrite of both
			layout1 = bgq_spinorfield_prepareRead(argfield1, isOdd, false, true, false, false);
			layout2 = bgq_spinorfield_prepareRead(argfield2, isOdd, false, true, false, false);
		}
	}
	bgq_worker_func workerfunc = NAME2(REDUCTION_NAME,worker_funcs)[layout1][layout2];
#endif
	assert(workerfunc);

	bgq_master_sync();
	static NAME2(REDUCTION_NAME,args_t) reduction_args;
	NAME2(REDUCTION_NAME,args_t) call_args = {
			.isOdd = isOdd
			IF1ARG(, .argfield1 = argfield1)
			IF2ARG(, .argfield2 = argfield2)
			IFEXTRA(, .extra = { REDUCTION_EXTRAARGS })
	};
	reduction_args = call_args;
	bgq_master_call(workerfunc, &reduction_args);

	size_t threads = omp_get_num_threads();
	REDUCTION_RETURNTYPE(result) = REDUCTION_VARINIT();
	for (size_t tid = 0; tid < threads; tid+=1) {
		REDUCTION_RETURNTYPE(threadresult) = reduction_args.threadresult[tid];
		REDUCTION_COMBINEFUNC(&result, threadresult);
	}

	return result;
}


#undef IF1ARG
#undef IF2ARG

#undef REDUCTION_WORKERNAME
#undef REDUCTION_WORKERINST

#undef REDUCTION_NAME
#undef REDUCTION_ARGFIELDS
#undef REDUCTION_RETURNTYPE
#undef REDUCTION_SITEREDUCEFUNC
#undef REDUCTION_COMBINEFUNC
#undef REDUCTION_VARINIT

