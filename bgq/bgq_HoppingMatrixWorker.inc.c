
#ifndef BGQ_HOPPINGMATRIXWORKER_INC_
#include "bgq_utils.h"
#include "bgq_HoppingMatrix.h"

#include <stdbool.h>

void bgq_HoppingMatrix_worker(void *arg, size_t tid, size_t threads, bool kamul, bool readFulllayout)
#endif
{
	bgq_HoppingMatrix_workload *work = arg;
	bool isOdd = work->isOdd_src;
	bgq_weylfield_controlblock * restrict spinorfield = work->spinorfield;
	bgq_weylfield_controlblock * restrict targetfield = work->targetfield;
	ucoord ic_begin = work->ic_begin;
	ucoord ic_end = work->ic_end;
	bool noprefetchstream = work->noprefetchstream;

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
	const size_t end = min_sizet(ic_end, begin+threadload);

	if (!noprefetchstream) {
		bgq_prefetch_forward(&g_bgq_gaugefield_fromCollapsed[isOdd][begin]);
		if (readFulllayout) {
			bgq_prefetch_forward(&spinorfield->sec_fullspinor[begin]);
		} else {
			bgq_prefetch_forward(&spinorfield->sec_collapsed[begin]);
		}
		bgq_prefetch_forward(&targetfield->sendptr[begin]);
	}

	for (ucoord ic = begin; ic<end; ic+=1) {
		//TODO: Check optaway
		ucoord ih = bgq_collapsed2halfvolume(isOdd,ic);
		ucoord t1 = bgq_halfvolume2t1(isOdd,ih);
		ucoord t2 = bgq_halfvolume2t2(isOdd,ih);
		ucoord x = bgq_halfvolume2x(ih);
		ucoord y = bgq_halfvolume2y(ih);
		ucoord z = bgq_halfvolume2z(ih);


		bgq_gaugesite *gaugesite = &g_bgq_gaugefield_fromCollapsed[isOdd][ic];
#if 0
		bgq_weyl_ptr_t destptrsx = {
				&targetfield->sec_collapsed[0].d[TUP],
				&targetfield->sec_collapsed[0].d[TUP],
				&targetfield->sec_collapsed[0].d[TUP],
				&targetfield->sec_collapsed[0].d[TUP],
				&targetfield->sec_collapsed[0].d[TUP],
				&targetfield->sec_collapsed[0].d[TUP],
				&targetfield->sec_collapsed[0].d[TUP],
				&targetfield->sec_collapsed[0].d[TUP]
		};
		bgq_weyl_ptr_t * restrict destptrs = &destptrsx;
#elif 0
		bgq_weyl_ptr_t destptrsx = {
				&targetfield->sec_collapsed[ic].d[TUP],
				&targetfield->sec_collapsed[ic].d[TDOWN],
				&targetfield->sec_collapsed[ic].d[XUP],
				&targetfield->sec_collapsed[ic].d[XDOWN],
				&targetfield->sec_collapsed[ic].d[YUP],
				&targetfield->sec_collapsed[ic].d[YDOWN],
				&targetfield->sec_collapsed[ic].d[ZUP],
				&targetfield->sec_collapsed[ic].d[ZDOWN]
		};
		bgq_weyl_ptr_t * restrict destptrs = &destptrsx;
#elif 0
		bgq_weyl_ptr_t destptrsx = {
				&targetfield->sec_collapsed->d[TUP],
				&targetfield->sec_collapsed->d[TDOWN],
				&targetfield->sec_collapsed->d[XUP],
				&targetfield->sec_collapsed->d[XDOWN],
				&targetfield->sec_collapsed->d[YUP],
				&targetfield->sec_collapsed->d[YDOWN],
				&targetfield->sec_collapsed->d[ZUP],
				&targetfield->sec_collapsed->d[ZDOWN]
		};
		bgq_weyl_ptr_t * restrict destptrs = &destptrsx;
#else
		bgq_weyl_ptr_t * restrict targetptrs = &targetfield->sendptr[ic];
#endif

		//TODO: prefetching
		//TODO: Check inlining
		bgq_su3_spinor_decl(spinor);
		if (readFulllayout) {
			bgq_spinorsite *spinorsite = &spinorfield->sec_fullspinor[ic];
			assert(spinorsite->s[1][0][0]!=0);
			//bgq_su3_spinor_prefetch_double(&spinorfield->sec_fullspinor[ic+1]); // TODO: This prefetch is too early
			bgq_HoppingMatrix_loadFulllayout(spinor, spinorsite, t1, t2, x, y, z);
		} else {
			bgq_weylsite *weylsite = &spinorfield->sec_collapsed[ic];
			assert(weylsite->d[TUP].s[1][0][0]!=0);
			bgq_HoppingMatrix_loadWeyllayout(spinor, weylsite, t1, t2, x, y, z);
		}
		//bgq_HoppingMatrix_compute_storeWeyllayout_alldir(destptrs, gaugesite, spinor, t1, t2, x, y, z, qka0,qka1,qka2,qka3,kamul);

#define BGQ_COMPUTEWEYL_INC_
#include "bgq_ComputeWeyl.inc.c"
	}

}


#undef BGQ_HOPPINGMATRIXWORKER_INC_
