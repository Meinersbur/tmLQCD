
#include "bgq_qpx.h"
#include "bgq_HoppingMatrix.h"
#include "bgq_comm.h"

#include <complex.h>


typedef struct {
	bool isOdd_src;
	bool isOdd_dst;
	bgq_weylfield_controlblock *targetfield;
	bgq_weylfield_controlblock *spinorfield;
	ucoord ic_begin;
	ucoord ic_end;
} bgq_HoppingMatrix_workload;


void *somewhere;
static inline void bgq_HoppingMatrix_worker_precomm_readFulllayout2(void * restrict arg, size_t tid, size_t threads, bool kamul, bool readFulllayout) {
	bgq_HoppingMatrix_workload *work = arg;
	bool isOdd = work->isOdd_src;
	bgq_weylfield_controlblock *spinorfield = work->spinorfield;
	bgq_weylfield_controlblock *targetfield = work->targetfield;
	ucoord ic_begin = work->ic_begin;
	ucoord ic_end = work->ic_end;

	bgq_vector4double_decl(qka0);
	bgq_complxval_splat(qka0,ka0);
	bgq_vector4double_decl(qka1);
	bgq_complxval_splat(qka1,ka1);
	bgq_vector4double_decl(qka2);
	bgq_complxval_splat(qka2,ka2);
	bgq_vector4double_decl(qka3);
	bgq_complxval_splat(qka3,ka3);

	const size_t workload = ic_end - ic_begin;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = ic_begin + tid*threadload;
	const size_t end = ic_end + min_sizet(workload, begin+threadload);
	for (ucoord ic = begin; ic<end; ic+=1) {
		//TODO: Check optaway
		ucoord ih = bgq_collapsed2halfvolume(isOdd,ic);
		ucoord t1 = bgq_halfvolume2t1(isOdd,ih);
		ucoord t2 = bgq_halfvolume2t2(isOdd,ih);
		ucoord x = bgq_halfvolume2x(ih);
		ucoord y = bgq_halfvolume2y(ih);
		ucoord z = bgq_halfvolume2z(ih);

		bgq_gaugesite *gaugesite = &g_bgq_gaugefield_fromCollapsed[isOdd][ic];
		bgq_weyl_ptr_t *destptrs = &targetfield->sendptr[ic];

		//TODO: prefetching
		//TODO: Check inlining
		bgq_su3_spinor_decl(spinor);
		if (readFulllayout) {
			bgq_spinorsite *spinorsite = &spinorfield->sec_fullspinor[ic];
			bgq_HoppingMatrix_loadFulllayout(spinor, spinorsite, t1, t2, x, y, z);
		} else {
			bgq_weylsite *weylsite = &spinorfield->sec_collapsed[ic];
			bgq_HoppingMatrix_loadWeyllayout(spinor, weylsite, t1, t2, x, y, z);
		}
		bgq_HoppingMatrix_compute_storeWeyllayout_alldir(destptrs, gaugesite, spinor, t1, t2, x, y, z, qka0,qka1,qka2,qka3,kamul);
	}
}

void bgq_HoppingMatrix_nokamul_worker_precomm_readFulllayout2(void *arg, size_t tid, size_t threads) {
	bgq_HoppingMatrix_worker_precomm_readFulllayout2(arg,tid,threads,false,true);
}
