
#define HALOCONDITION_TUP (HALO_T && t2_dst==LOCAL_LT-1) /* ??? */
#define HALOCONDITION_TDOWN (HALO_T && t1_dst==0) /* ??? */
#define HALOCONDITION_XUP (HALO_X && x_dst==LOCAL_LX-1)
#define HALOCONDITION_XDOWN (HALO_X && x_dst==0)
#define HALOCONDITION_YUP (HALO_Y && y_dst==LOCAL_LY-1)
#define HALOCONDITION_YDOWN (HALO_Y && y_dst==0)
#define HALOCONDITION_ZUP (HALO_Z && z_dst==LOCAL_LZ-1)
#define HALOCONDITION_ZDOWN (HALO_Z && z_dst==0)

#ifndef BGQ_HOPPINGMATRIXWORKER_INC_
#include "bgq_utils.h"
#include "bgq_HoppingMatrix.h"
#include "bgq_gaugefield.h"
#include "bgq_comm.h"

#include "../boundary.h"

#include <stdbool.h>

#define PRECISION double

void bgq_HoppingMatrix_worker(void *arg, size_t tid, size_t threads, bool kamul, bool readFulllayout, bool writeFulllayout)
#endif

{
bgq_HoppingMatrix_workload *work = arg;
bool isOdd_src = work->isOdd_src;
bool isOdd_dst = work->isOdd_dst;
bgq_weylfield_controlblock * restrict spinorfield = work->spinorfield;
bgq_weylfield_controlblock * restrict targetfield = work->targetfield;
ucoord ic_begin = work->ic_begin;
ucoord ic_end = work->ic_end;
bool noprefetchstream = work->noprefetchstream;

assert(spinorfield->isOdd == isOdd_src);
assert(targetfield->isOdd == isOdd_dst);

bgq_vector4double_decl(qka0);
bgq_complxval_splat(qka0,ka0);
bgq_vector4double_decl(qka1);
bgq_complxval_splat(qka1,ka1);
bgq_vector4double_decl(qka2);
bgq_complxval_splat(qka2,ka2);
bgq_vector4double_decl(qka3);
bgq_complxval_splat(qka3,ka3);

if (writeFulllayout) {
	assert(readFulllayout);
	const ucoord workload = ic_end - ic_begin;
	const ucoord threadload = (workload+threads-1)/threads;
	const ucoord begin = ic_begin + tid*threadload;
	const ucoord end = min_sizet(ic_end, begin+threadload);

	bgq_gaugesite *gaugesite = &g_bgq_gaugefield_fromCollapsed_dst[isOdd_dst][begin];
#ifndef NDEBUG
	if (begin < end) {
		ucoord t1_dst = bgq_collapsed2t1(isOdd_dst, begin);
		ucoord t2_dst = bgq_collapsed2t2(isOdd_dst, begin);
		ucoord tv_dst = bgq_collapsed2tv(isOdd_dst, begin);
		ucoord x_dst = bgq_collapsed2x(isOdd_dst, begin);
		ucoord y_dst = bgq_collapsed2y(isOdd_dst, begin);
		ucoord z_dst = bgq_collapsed2z(isOdd_dst, begin);
		bgq_su3_mdecl(gauge);
		bgq_su3_matrix_load_double(gauge, gaugesite);
		bgq_gaugeqpx_expect(gauge, t1_dst, t2_dst, x_dst, y_dst, z_dst, TUP, false);
	}
#endif
	gaugesite = (bgq_gaugesite*)(((uint8_t*)gaugesite)-32);

	for (ucoord ic_dst = begin; ic_dst<end; ic_dst+=1) {
#ifndef NDEBUG
		ucoord t1_dst = bgq_collapsed2t1(isOdd_dst, ic_dst);
		ucoord t2_dst = bgq_collapsed2t2(isOdd_dst, ic_dst);
		ucoord tv_dst = bgq_collapsed2tv(isOdd_dst, ic_dst);
		ucoord x_dst = bgq_collapsed2x(isOdd_dst, ic_dst);
		ucoord y_dst = bgq_collapsed2y(isOdd_dst, ic_dst);
		ucoord z_dst = bgq_collapsed2z(isOdd_dst, ic_dst);
#endif

		bgq_spinor_vec *target = &targetfield->BGQ_SEC_FULLLAYOUT[ic_dst];
		const bool gaugemul = true;
		#define BGQ_HOPPINGMATRIXSTENCIL_INC_ 1
		#include "bgq_HoppingMatrixStencil.inc.c"
	}


#if 0
	for (ucoord ic_next = begin; ic_next<end; ) {
		if (bgq_collapsed2isSurface(ic_next)) {
			const ucoord begin_surface = ic_next;
			const ucoord end_surface = min_sizet(end, PHYSICAL_SURFACE);

			for (ucoord ic_dst = begin_surface; ic_dst < end_surface; ic_dst+=1) {
				assert(bgq_collapsed2isSurface(ic_dst));
			}
			ic_next = end_surface;
		} else if (bgq_collapsed2isOuter(ic_next)) {
			const ucoord begin_outer = ic_next;
			const ucoord end_outer = min_sizet(end, PHYSICAL_OUTER);

			for (ucoord ic_dst = begin_outer; ic_dst < end_outer; ic_dst+=1) {
				assert(bgq_collapsed2isOuterBody(ic_dst));
			}
			ic_next = end_outer;
		} else {
			const ucoord begin_inner = ic_next;
			const ucoord end_inner = end;

			for (ucoord ic_inner = begin_inner; ic_inner < end_inner; ) {
				const ucoord begin_line = ic_inner;
				ucoord z_dst = bgq_collapsed2z(isOdd_dst, begin_line);
				const ucoord end_line = min_sizet(end_inner, begin_line + (PHYSICAL_LZ - z_dst - 1));

				for (ucoord ic_dst = begin_line; ic_dst < end_line; ic_dst+=1) {
					assert(bgq_collapsed2isInnerBody(ic_dst));

				}
				ic_inner = end_line;
			}
			ic_next = end_inner;
		}
	}
#endif


} else {
	const size_t workload = ic_end - ic_begin;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = ic_begin + tid*threadload;
	const size_t end = min_sizet(ic_end, begin+threadload);

	if (!noprefetchstream) {
		bgq_prefetch_forward(&g_bgq_gaugefield_fromCollapsed_src[isOdd_src][begin]);
		if (readFulllayout) {
			bgq_prefetch_forward(&spinorfield->BGQ_SEC_FULLLAYOUT[begin]);
		} else {
			bgq_prefetch_forward(&spinorfield->BGQ_SEC_WEYLLAYOUT[begin]);
		}
		bgq_prefetch_forward(&targetfield->BGQ_SENDPTR[isOdd_dst][begin]);
	}

	const bool gaugemul = true;
	bgq_gaugesite *gaugesite = &g_bgq_gaugefield_fromCollapsed_src[isOdd_src][begin];
	gaugesite = (bgq_gaugesite*)(((uint8_t*)gaugesite)-32);
	bgq_weylsite *weylsite = &spinorfield->BGQ_SEC_WEYLLAYOUT[begin];
	weylsite = (bgq_weylsite*)(((uint8_t*)weylsite) - BGQ_VECTORSIZE);

	for (ucoord ic = begin; ic<end; ic+=1) {
#ifndef NDEBUG
		ucoord ih = bgq_collapsed2halfvolume(isOdd_src, ic);
		ucoord t1 = bgq_halfvolume2t1(isOdd_src, ih);
		ucoord t2 = bgq_halfvolume2t2(isOdd_src, ih);
		ucoord x = bgq_halfvolume2x(ih);
		ucoord y = bgq_halfvolume2y(ih);
		ucoord z = bgq_halfvolume2z(ih);
#endif
		bgq_weyl_ptr_t * restrict targetptrs = &targetfield->BGQ_SENDPTR[isOdd_dst][ic];

		//TODO: Check inlining
		bgq_su3_spinor_decl(spinor);

		if (readFulllayout) {
			bgq_su3_matrix_prefetch(gaugesite);
			bgq_spinorsite *spinorsite = &spinorfield->BGQ_SEC_FULLLAYOUT[ic];
			//assert(spinorsite->s[1][0][0]!=0);
			//bgq_su3_spinor_prefetch_double(&spinorfield->sec_fullspinor[ic+1]); // TODO: This prefetch is too early
			//bgq_HoppingMatrix_loadFulllayout(spinor, spinorsite, t1, t2, x, y, z);
			//bgq_su3_spinor_valgen(spinor);
			bgq_su3_spinor_load(spinor, spinorsite);
					bgq_spinorqpx_expect(spinor, t1, t2, x, y, z);

			#define BGQ_COMPUTEWEYL_INC_ 1
			#define BGQ_COMPUTEWEYL_INSERTPREFETCH bgq_su3_spinor_prefetch(weylsite);
			#include "bgq_ComputeWeyl.inc.c"
		} else {
			#define BGQ_READWEYLLAYOUT_INC_ 1
			#define BGQ_READWEYLLAYOUT_INSERTPREFETCH bgq_su3_matrix_prefetch_double(gaugesite);
			#include "bgq_ReadWeyllayout.inc.c"

			#define BGQ_COMPUTEWEYL_INC_ 1
			#define BGQ_COMPUTEWEYL_INSERTPREFETCH bgq_su3_weyl_prefetch(weylsite);
			#include "bgq_ComputeWeyl.inc.c"
		}
	}

}
}
#undef BGQ_HOPPINGMATRIXWORKER_INC_
