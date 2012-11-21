/*
 * bgq_operator.inc.c
 *
 *  Created on: Nov 9, 2012
 *      Author: meinersbur
 */

#include "bgq_spinorfield.h"
#include "bgq_dispatch.h"
#include "bgq_HoppingMatrix.h"
#include "bgq_qpx.h"
#include "bgq_utils.h"

#include <stdbool.h>




#ifndef OPERATOR_NAME
#define OPERATOR_NAME bgq_operator
#define OPERATOR_ARGFIELDS 1
#define OPERATOR_EXTRAPARMS double factor
#define OPERATOR_EXTRAARGS factor
#define OPERATOR_VECSITEFUNC vecsitefunc_raw

static inline void vecsitefunc_raw(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor1), double factor, ucoord tv, ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z) {
}
#endif


#ifndef OPERATOR_EXTRAPARMS
#define OPERATOR_EXTRAPARMS
#define OPERATOR_EXTRAARGS
#endif




#if 0
#define OPERATOR_OUTPLACENAME bgq_operator_outplace
#define OPERATOR_INPLACENAME bgq_operator_inplace


static inline void vecsitefunc_raw(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), double factor, ucoord tv, ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z) {

}


#ifdef OPERATOR_NAME
#define OPERATOR_INCLUDED 0
#else
#define OPERATOR_INCLUDED 1

#define OPERATOR_INPLACENAME bgq_operator_inplace
#define OPERATOR_OUTPLACENAME bgq_operator_outplace




#endif


#define OPERATOR_EXTRAPARMLIST IF_EMPTY_INPAREN((),(,),OPERATOR_EXTRAPARMS) OPERATOR_EXTRAPARMS
#define OPERATOR_EXTRAARGLIST IF_EMPTY_INPAREN((),(,),OPERATOR_EXTRAARGS) OPERATOR_EXTRAARGS




typedef struct {
	PREPROCESSOR_FOREACH( ,;,;, ,IDENTITY,OPERATOR_EXTRAPARMS)
} bgq_operator_extraparm_t;

typedef struct {
	bool isOdd;
	bgq_weylfield_controlblock *targetfield;
	bgq_weylfield_controlblock *argfield1;
	bgq_weylfield_controlblock *argfield2;
	bgq_operator_extraparm_t extra;
} bgq_operator_args_t;





static void bgq_operator_worker_readWeyllayout1_readWeyllayout2(void *arg_untyped, size_t tid, size_t threads) {
	bgq_operator_worker(arg_untyped, tid, threads, false, false);
}
static void bgq_operator_worker_readWeyllayout1_readFulllayout2(void *arg_untyped, size_t tid, size_t threads) {
	bgq_operator_worker(arg_untyped, tid, threads, false, true);
}
static void bgq_operator_worker_readFulllayout1_readWeyllayout2(void *arg_untyped, size_t tid, size_t threads) {
	bgq_operator_worker(arg_untyped, tid, threads, true, false);
}
static void bgq_operator_worker_readFulllayout1_readFulllayout2(void *arg_untyped, size_t tid, size_t threads) {
	bgq_operator_worker(arg_untyped, tid, threads, true, true);
}

static bgq_worker_func bgq_worker_funcs[2][2] = {
		{&bgq_operator_worker_readWeyllayout1_readWeyllayout2, &bgq_operator_worker_readWeyllayout1_readFulllayout2},
		{&bgq_operator_worker_readFulllayout1_readWeyllayout2, &bgq_operator_worker_readFulllayout1_readFulllayout2}
};

#endif


#define OPERATOR_EXTRAPARMLIST IF_EMPTY_INPAREN((),(,),OPERATOR_EXTRAPARMS) OPERATOR_EXTRAPARMS
#define OPERATOR_EXTRAARGLIST IF_EMPTY_INPAREN((),(,),OPERATOR_EXTRAARGS) OPERATOR_EXTRAARGS

#if OPERATOR_ARGFIELDS>=1
#define IF1ARG(...) __VA_ARGS__
#else
#define IF1ARG(...)
#endif

#if OPERATOR_ARGFIELDS>=2
#define IF2ARG(...) __VA_ARGS__
#else
#define IF2ARG(...)
#endif

#if ISEMPTY(OPERATOR_EXTRAPARMS)
#define IFEXTRA(...)
#else
#define IFEXTRA(...) __VA_ARGS__
#endif


#if !ISEMPTY(OPERATOR_EXTRAPARMS)
typedef struct {
	PREPROCESSOR_FOREACH( ,;,;, ,IDENTITY,OPERATOR_EXTRAPARMS)
} bgq_operator_extraparm_t;
#endif

typedef struct {
	bool isOdd;
	bgq_weylfield_controlblock *targetfield;
	IF1ARG(bgq_weylfield_controlblock *argfield1;)
	IF2ARG(bgq_weylfield_controlblock *argfield2;)
	IFEXTRA(bgq_operator_extraparm_t extra;)
} NAME2(REDUCTION_NAME,args_t);


#define OPERATOR_EXTRA_ASSIGN(varname) (varname) = ((arg->extra).varname);


static inline void bgq_operator_worker(void *arg_untyped, size_t tid, size_t threads, bool writeSloppy IF1ARG(, bool readWeyllayout1, bool sloppy1, bool mul1) IF2ARG(, bool readWeyllayout2, bool sloppy2, bool mul2)) {
	NAME2(REDUCTION_NAME,args_t) *arg = (NAME2(REDUCTION_NAME,args_t)*)arg_untyped;
	bool isOdd = arg->isOdd;
	bgq_weylfield_controlblock *targetfield = arg->targetfield;
	IF1ARG(bgq_weylfield_controlblock *argfield1 = arg->argfield1;)
	IF2ARG(bgq_weylfield_controlblock *argfield2 = arg->argfield2;)
	PREPROCESSOR_FOREACH( ,;,;, ,IDENTITY,OPERATOR_EXTRAPARMS)
	PREPROCESSOR_FOREACH(,,,,OPERATOR_EXTRA_ASSIGN,OPERATOR_EXTRAARGS)

	size_t workload = PHYSICAL_VOLUME;
	size_t threadload = (workload+threads-1)/threads;
	ucoord beginj = tid*threadload;
	ucoord endj = min_sizet((tid+1)*threadload,workload);
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

#if OPERATOR_ARGFIELDS>=1
		bgq_su3_spinor_decl(spinor1);
		bgq_spinorfield_readSpinor(&spinor1, argfield1, ic, readWeyllayout1, sloppy1, mul1);
#endif

#if OPERATOR_ARGFIELDS>=2
		bgq_su3_spinor_decl(spinor2);
		bgq_spinorfield_readSpinor(&spinor2, argfield2, ic, readWeyllayout2, sloppy2, mul2);
#endif

		bgq_su3_spinor_decl(targetspinor);
		OPERATOR_VECSITEFUNC(bgq_su3_spinor_vars(&targetspinor) IF1ARG(, bgq_su3_spinor_vars(spinor1)) IF2ARG(, bgq_su3_spinor_vars(spinor2)) OPERATOR_EXTRAARGLIST, tv, t1, t2, x, y, z);

		// Warning: targetsite==argfield1 is possible and the inline assembler does not have memory barriers! If there is not data dependency, the compiler might arrange the stores before the loads!
		//asm volatile ("" : : : );
		if (writeSloppy) {
			bgq_spinorsite_float *targetsite = &targetfield->sec_fullspinor_float[ic];
			bgq_su3_spinor_store_float(targetsite, targetspinor);
		} else {
			bgq_spinorsite_double *targetsite = &targetfield->sec_fullspinor_double[ic];
			bgq_su3_spinor_store_double(targetsite, targetspinor);
		}
	}
}




#if 0
#define OPERATOR_MAKELIST0(MACRO) \
		MACRO()

#define OPERATOR_MAKELIST1(MACRO) \
		MACRO(0),MACRO(1)

#define OPERATOR_MAKELIST2(MACRO) \
		MACRO(0,0),MACRO(0,1),MACRO(1,0),MACRO(1,1)

#define OPERATOR_MAKELIST3(MACRO) \
		MACRO(0,0,0),MACRO(0,0,1),MACRO(0,1,0),MACRO(0,1,1),MACRO(1,0,0),MACRO(1,0,1),MACRO(1,1,0),MACRO(1,1,1)


#define OPERATOR_WORKERINST(layout1,layout2,PRECISION) \
	static void OPERATOR_WORKERNAME(layout1,layout2,PRECISION)(void *arg_untyped, size_t tid, size_t threads) { \
		NAME2(bgq_operator_worker,PRECISION)(arg_untyped, tid, threads, (layout1)&ly_weyl, (layout1)&ly_sloppy, (layout1)&ly_mul, (layout2)&ly_weyl, (layout2)&ly_sloppy, (layout2)&ly_mul); \
	}


#define OPERATOR_INSTWORKER(arg1,arg2) \
	OPERATOR_WORKERINST(NAME2(OPERATOR_),NAME2(),NAME2())

OPERATOR_MAKELIST2()




#define OPERATOR_WORKERLIST(prefix,postfix) \
		{ NAME2(prefix,CONCAT(readFulllayout,postfix)), NAME2(prefix,CONCAT(readWeyllayout,postfix)), NAME2(prefix,CONCAT(readFulllayoutSloppy,postfix)), NAME2(prefix,CONCAT(readWeyllayoutSloppy,postfix)), NULL }
#endif


#if OPERATOR_ARGFIELDS==1
static void bgq_operator_worker_full_double(void *arg_untyped, size_t tid, size_t threads) {
	bgq_operator_worker(arg_untyped, tid, threads, false, false, false, false );
}

static void bgq_operator_worker_weyl_double(void *arg_untyped, size_t tid, size_t threads) {
	bgq_operator_worker(arg_untyped, tid, threads, false, true, false, false );
}


static void bgq_operator_worker_full_float(void *arg_untyped, size_t tid, size_t threads) {
	bgq_operator_worker(arg_untyped, tid, threads, true, false, true, false );
}

static void bgq_operator_worker_weyl_float(void *arg_untyped, size_t tid, size_t threads) {
	bgq_operator_worker(arg_untyped, tid, threads, true, true, true, false );
}


static bgq_worker_func bgq_worker_funcs[BGQ_SPINORFIELD_LAYOUT_COUNT] =
	                       /* ly_full_double */              /* ly_full_float */              /* ly_weyl_double */              /* ly_weyl_float */
	/* ly_full_double */ { &bgq_operator_worker_full_double, &bgq_operator_worker_full_float, &bgq_operator_worker_weyl_double, &bgq_operator_worker_weyl_float };

#elif OPERATOR_ARGFIELDS==2

static void bgq_operator_worker_full_full_double(void *arg_untyped, size_t tid, size_t threads) {
	bgq_operator_worker(arg_untyped, tid, threads, false, false, false, false, false, false, false );
}

static void bgq_operator_worker_full_weyl_double(void *arg_untyped, size_t tid, size_t threads) {
	bgq_operator_worker(arg_untyped, tid, threads, false, true, false, false, false, false, false );
}

static void bgq_operator_worker_weyl_full_double(void *arg_untyped, size_t tid, size_t threads) {
	bgq_operator_worker(arg_untyped, tid, threads, false,  false, false, false, true, false, false );
}

static void bgq_operator_worker_weyl_weyl_double(void *arg_untyped, size_t tid, size_t threads) {
	bgq_operator_worker(arg_untyped, tid, threads, false,  true, false, false, true, false, false );
}


static void bgq_operator_worker_full_full_float(void *arg_untyped, size_t tid, size_t threads) {
	bgq_operator_worker(arg_untyped, tid, threads, true, false, true, false, false, true, false );
}

static void bgq_operator_worker_full_weyl_float(void *arg_untyped, size_t tid, size_t threads) {
	bgq_operator_worker(arg_untyped, tid, threads, true, true, true, false, false, true, false );
}

static void bgq_operator_worker_weyl_full_float(void *arg_untyped, size_t tid, size_t threads) {
	bgq_operator_worker(arg_untyped, tid, threads, true, false, true, false, true, true, false );
}

static void bgq_operator_worker_weyl_weyl_float(void *arg_untyped, size_t tid, size_t threads) {
	bgq_operator_worker(arg_untyped, tid, threads, true, true, true, false, true, true, false );
}


static bgq_worker_func bgq_worker_funcs[BGQ_SPINORFIELD_LAYOUT_COUNT][BGQ_SPINORFIELD_LAYOUT_COUNT] = {
	                       /* ly_full_double */                   /* ly_full_float */                    /* ly_weyl_double */                   /* ly_weyl_float */
	/* ly_full_double */ { &bgq_operator_worker_full_full_double, NULL,                                  &bgq_operator_worker_full_weyl_double, NULL                                 },
	/* ly_full_float  */ { NULL,                                  &bgq_operator_worker_full_full_float,  NULL,                                  &bgq_operator_worker_full_weyl_float },
	/* ly_weyl_double */ { &bgq_operator_worker_weyl_full_double, NULL,                                  &bgq_operator_worker_weyl_weyl_double, NULL                                 },
	/* ly_weyl_float  */ { NULL,                                  &bgq_operator_worker_weyl_full_float,  NULL,                                  &bgq_operator_worker_weyl_weyl_float }
};
#endif


static void NAME2(OPERATOR_NAME,selector)(bool sloppy, bgq_weylfield_controlblock *targetfield IF1ARG(, bgq_weylfield_controlblock *argfield1) IF2ARG(, bgq_weylfield_controlblock *argfield2) OPERATOR_EXTRAPARMLIST) {
	bool isOdd = argfield1->isOdd;

	IF1ARG(bgq_spinorfield_layout layout1 = bgq_spinorfield_prepareRead(argfield1, isOdd, true, !sloppy, sloppy, false);)
	IF2ARG(bgq_spinorfield_layout layout2 = bgq_spinorfield_prepareRead(argfield2, isOdd, true, !sloppy, sloppy, false);)
	bgq_spinorfield_prepareWrite(targetfield, isOdd, ly_full_double);

	bgq_master_sync();
	static NAME2(REDUCTION_NAME,args_t) operator_args;
	NAME2(REDUCTION_NAME,args_t) call_args = {
			.isOdd = isOdd,
			.targetfield = targetfield
			IF1ARG(, .argfield1 = argfield1)
			IF2ARG(, .argfield2 = argfield2)
			IFEXTRA(, .extra = { OPERATOR_EXTRAARGS })
	};
	operator_args = call_args;
	bgq_worker_func workerfunc = bgq_worker_funcs IF1ARG([layout1]) IF2ARG([layout2]);
	bgq_master_call(workerfunc, &operator_args);
}


void NAME2(OPERATOR_NAME,double)(bgq_weylfield_controlblock *targetfield IF1ARG(, bgq_weylfield_controlblock *argfield1) IF2ARG(, bgq_weylfield_controlblock *argfield2) OPERATOR_EXTRAPARMLIST) {
	NAME2(OPERATOR_NAME,selector)(false, targetfield IF1ARG(, argfield1) IF2ARG(, argfield2) OPERATOR_EXTRAARGLIST);
}


void NAME2(OPERATOR_NAME,float)(bgq_weylfield_controlblock *targetfield IF1ARG(, bgq_weylfield_controlblock *argfield1) IF2ARG(, bgq_weylfield_controlblock *argfield2) OPERATOR_EXTRAPARMLIST) {
	NAME2(OPERATOR_NAME,selector)(true, targetfield IF1ARG(, argfield1) IF2ARG(, argfield2) OPERATOR_EXTRAARGLIST);
}


#undef OPERATOR_NAME
#undef OPERATOR_OUTPLACENAME
#undef OPERATOR_INPLACENAME
#undef OPERATOR_ARGFIELDS
#undef OPERATOR_EXTRAPARMS
#undef OPERATOR_EXTRAARGS
#undef OPERATOR_VECSITEFUNC

#undef OPERATOR_EXTRAARGLIST
#undef OPERATOR_EXTRAPARMLIST
