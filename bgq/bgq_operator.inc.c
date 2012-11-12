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

#include <stdbool.h>

#if 0
#define OPERATOR_OUTPLACENAME bgq_operator_outplace
#define OPERATOR_INPLACENAME bgq_operator_inplace
#define OPERATOR_ARGFIELDS 1
#define OPERATOR_VECSITEFUNC vecsitefunc_raw


static void vecsitefunc_raw(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), double factor, ucoord tv, ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z) {

}


#ifdef OPERATOR_NAME
#define OPERATOR_INCLUDED 0
#else
#define OPERATOR_INCLUDED 1

#define OPERATOR_INPLACENAME bgq_operator_inplace
#define OPERATOR_OUTPLACENAME bgq_operator_outplace








#define OPERATOR_EXTRAPARMS , double factor
#define OPERATOR_EXTRAARGS factor
#endif

#if OPERATOR_ARGFIELDS==0
#define OPERATOR_ARGFIELDS_PARMS
#define OPERATOR_ARGFIELDS_ARGS
#elif OPERATOR_ARGFIELDS==1
#define OPERATOR_ARGFIELDS_PARMS , bgq_weylfield_controlblock *argfield1
#define OPERATOR_ARGFIELDS_ARGS argfield1
#elif OPERATOR_ARGFIELDS==2
#define OPERATOR_ARGFIELDS_PARMS , bgq_weylfield_controlblock *argfield1, bgq_weylfield_controlblock *argfield2
#define OPERATOR_ARGFIELDS_ARGS argfield1, argfield2
#else
#error Unsupported number of field arguments
#endif






typedef struct {
	bool isOdd;
	bgq_weylfield_controlblock *targetfield;
	bgq_weylfield_controlblock *argfield1;
	bgq_weylfield_controlblock *argfield2;
	double factor;
} bgq_operator_args_t;


static inline void bgq_operator_worker(void *arg_untyped, size_t tid, size_t threads, bool readFulllayout1, bool readFulllayout2) {
	bgq_operator_args_t *arg = (bgq_operator_args_t*)arg_untyped;
	bool isOdd = arg->isOdd;
	bgq_weylfield_controlblock *targetfield = arg->targetfield;
	bgq_weylfield_controlblock *argfield1 = arg->argfield1;
	bgq_weylfield_controlblock *argfield2 = arg->argfield2;
	double factor = arg->factor;

	size_t workload = PHYSICAL_VOLUME;
	size_t threadload = (workload+threads-1)/threads;
	ucoord beginj = tid*threadload;
	ucoord endj = min_sizet((tid+1)*threadload,workload);
	for (ucoord ic = beginj; ic < endj; ic+=1) {
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

		bgq_su3_spinor_decl(targetspinor);
		OPERATOR_VECSITEFUNC(bgq_su3_spinor_vars(&targetspinor), bgq_su3_spinor_vars(spinor1), bgq_su3_spinor_vars(spinor2), factor, tv, t1, t2, x, y, z);

		bgq_spinorsite *targetsite = &targetfield->sec_fullspinor[ic];
		bgq_su3_spinor_store_double(targetsite, targetspinor);
	}
}

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

#ifdef OPERATOR_OUTPLACENAME
void OPERATOR_OUTPLACENAME(bool isOdd, bgq_weylfield_controlblock *targetfield, bgq_weylfield_controlblock *argfield1, bgq_weylfield_controlblock *argfield2, double factor) {
	bgq_weylfield_controlblock *argfield[] = { argfield1, argfield2 };
	bool useFulllayout[2];
	bool useWeyllayout[2];

	bgq_spinorfield_setup(targetfield, isOdd , false, true, false, false);
	size_t nFieldargs = OPERATOR_ARGFIELDS;
	for (size_t i = 0; i < nFieldargs; i+=1) {
		useFulllayout[i] = argfield[i]->hasFullspinorData;
		useWeyllayout[i] = !useFulllayout[i];
		bgq_spinorfield_setup(argfield[i], isOdd, useFulllayout[i], false, useWeyllayout[i], false);
	}


	bgq_master_sync();
	static bgq_operator_args_t operator_args;
	bgq_operator_args_t call_args = {
			.targetfield = targetfield,
			.argfield1 = argfield1,
			.argfield2 = argfield2,
			.factor = factor
	};
	operator_args = call_args;
	bgq_master_call(&bgq_worker_funcs[useFulllayout[0]][useFulllayout[1]], &operator_args);
}
#endif


#ifdef OPERATOR_INPLACENAME
void bgq_operator_inplace(bgq_weylfield_controlblock *targetfield) {
	bool isOdd = targetfield->isOdd;

	bool useFulllayout = targetfield->hasFullspinorData;
	bool useWeyllayout = !useFulllayout;

	bgq_spinorfield_setup(targetfield, isOdd, useFulllayout, true, false, useWeyllayout);

	bgq_master_sync();
	bgq_master_call(&bgq_operator_worker[useFulllayout[0]][useFulllayout[1]], &operator_args);
}
#endif



#undef OPERATOR_INCLUDED
#undef OPERATOR_NAME
#endif
