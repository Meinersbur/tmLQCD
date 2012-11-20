
#ifndef BGQ_HOPPINGMATRIXWORKER_INC_
#include "bgq_utils.h"
#include "bgq_HoppingMatrix.h"

#include <stdbool.h>

#define PRECISION double


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

	bgq_gaugesite *gaugesite = &g_bgq_gaugefield_fromCollapsed[isOdd][begin];
	gaugesite = (bgq_gaugesite*)(((uint8_t*)gaugesite)-32);
	bgq_weylsite *weylsite = &spinorfield->BGQ_SEC_WEYLLAYOUT[begin];
	weylsite = (bgq_weylsite*)(((uint8_t*)weylsite) - BGQ_VECTORSIZE);

	if (true || readFulllayout) {
	for (ucoord ic = begin; ic<end; ic+=1) {
#ifndef NDEBUG
		ucoord ih = bgq_collapsed2halfvolume(isOdd, ic);
		ucoord t1 = bgq_halfvolume2t1(isOdd, ih);
		ucoord t2 = bgq_halfvolume2t2(isOdd, ih);
		ucoord x = bgq_halfvolume2x(ih);
		ucoord y = bgq_halfvolume2y(ih);
		ucoord z = bgq_halfvolume2z(ih);
#endif
#if 0
		bgq_weyl_ptr_t destptrsx = {
				&targetfield->sec_collapsed[begin].d[TUP],
				&targetfield->sec_collapsed[begin].d[TUP],
				&targetfield->sec_collapsed[begin].d[TUP],
				&targetfield->sec_collapsed[begin].d[TUP],
				&targetfield->sec_collapsed[begin].d[TUP],
				&targetfield->sec_collapsed[begin].d[TUP],
				&targetfield->sec_collapsed[begin].d[TUP],
				&targetfield->sec_collapsed[begin].d[TUP]
		};
		bgq_weyl_ptr_t * restrict targetptrs = &destptrsx;
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
		bgq_weyl_ptr_t * restrict targetptrs = &targetfield->BGQ_SENDPTR[ic];
#endif

		//TODO: prefetching
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


		//bgq_su3_weylnext_prefetch_double(weylsite); //TODO: too late
	}
	} else {
#if 0
		ucoord ic_load = begin;
		bgq_weylsite *weylsite = (bgq_weylsite*)(((uintptr_t)&spinorfield->sec_collapsed[ic_load])-32);
		bgq_gaugesite *gaugesite = (bgq_gaugesite *)(((uintptr_t)&g_bgq_gaugefield_fromCollapsed[isOdd][ic_load])-32);
		bgq_su3_spinor_decl(spinornext);

		// Prologue
		{
			bgq_su3_weylnext_prefetch_double(weylsite);

			// TUP
			{
				bgq_su3_weyl_decl(weylnext_tup);
				bgq_qvlfduxa(weylnext_tup_v0_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_tup_v0_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_tup_v0_c2, weylsite, 32);
				bgq_qvlfduxa(weylnext_tup_v1_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_tup_v1_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_tup_v1_c2, weylsite, 32);
				bgq_su3_expand_weyl_tup(spinornext, weylnext_tup);
			}

			bgq_su3_weylnext_prefetch_double(weylsite);

			// TDOWN
			{
				bgq_su3_weyl_decl(weylnext_tdown);
				bgq_qvlfduxa(weylnext_tdown_v0_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_tdown_v0_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_tdown_v0_c2, weylsite, 32);
				bgq_qvlfduxa(weylnext_tdown_v1_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_tdown_v1_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_tdown_v1_c2, weylsite, 32);
				bgq_su3_accum_weyl_tdown(spinornext, weylnext_tdown);
			}

			bgq_su3_weylnext_prefetch_double(weylsite);

			// X+ /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weylnext_xup);
				bgq_qvlfduxa(weylnext_xup_v0_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_xup_v0_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_xup_v0_c2, weylsite, 32);
				bgq_qvlfduxa(weylnext_xup_v1_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_xup_v1_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_xup_v1_c2, weylsite, 32);
				bgq_su3_accum_weyl_xup(spinornext, weylnext_xup);
			}

			bgq_su3_weylnext_prefetch_double(weylsite);

			// X- /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weylnext_xdown);
				bgq_qvlfduxa(weylnext_xdown_v0_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_xdown_v0_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_xdown_v0_c2, weylsite, 32);
				bgq_qvlfduxa(weylnext_xdown_v1_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_xdown_v1_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_xdown_v1_c2, weylsite, 32);
				bgq_su3_accum_weyl_xdown(spinornext, weylnext_xdown);
			}

			bgq_su3_weylnext_prefetch_double(weylsite);

			// Y+ /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weylnext_yup);
				bgq_qvlfduxa(weylnext_yup_v0_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_yup_v0_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_yup_v0_c2, weylsite, 32);
				bgq_qvlfduxa(weylnext_yup_v1_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_yup_v1_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_yup_v1_c2, weylsite, 32);
				bgq_su3_accum_weyl_yup(spinornext, weylnext_yup);
			}

			bgq_su3_weylnext_prefetch_double(weylsite);

			// Y- /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weylnext_ydown);
				bgq_qvlfduxa(weylnext_ydown_v0_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_ydown_v0_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_ydown_v0_c2, weylsite, 32);
				bgq_qvlfduxa(weylnext_ydown_v1_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_ydown_v1_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_ydown_v1_c2, weylsite, 32);
				bgq_su3_accum_weyl_ydown(spinornext, weylnext_ydown);
			}

			bgq_su3_weylnext_prefetch_double(weylsite);

			// Z+ /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weylnext_zup);
				bgq_qvlfduxa(weylnext_zup_v0_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_zup_v0_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_zup_v0_c2, weylsite, 32);
				bgq_qvlfduxa(weylnext_zup_v1_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_zup_v1_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_zup_v1_c2, weylsite, 32);
				bgq_su3_accum_weyl_zup(spinornext, weylnext_zup);
			}

			bgq_su3_weylnext_prefetch_double(weylsite);
			bgq_su3_matrixnext_prefetch_double(gaugesite);

			// Z- /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weylnext_zdown);
				bgq_qvlfduxa(weylnext_zdown_v0_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_zdown_v0_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_zdown_v0_c2, weylsite, 32);
				bgq_qvlfduxa(weylnext_zdown_v1_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_zdown_v1_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_zdown_v1_c2, weylsite, 32);
				bgq_su3_accum_weyl_zdown(spinornext, weylnext_zdown);
			}
		}
		ucoord ic_comp = ic_load;
		ic_load+=1;
		for (; ic_load<end; ic_load+=1) {
			bgq_su3_spinor_decl(spinor);
			bgq_su3_spinor_mov(spinor, spinornext);
			//weylsite = (bgq_weylsite*)(((uintptr_t)&spinorfield->sec_collapsed[ic_load])-32);
			//gaugesite = (bgq_gaugesite *)(((uintptr_t)&g_bgq_gaugefield_fromCollapsed[isOdd][ic_comp])-32);
			bgq_weyl_ptr_t *targetptrs = &targetfield->sendptr[ic_comp];// 8*sizeof(bgq_weyl_vec*)=64 bytes

			bgq_su3_weylnext_prefetch_double(weylsite);
			bgq_su3_matrixnext_prefetch_double(gaugesite);

			// TUP
			{
				bgq_su3_weyl_decl(weylnext_tup);
				bgq_qvlfduxa(weylnext_tup_v0_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_tup_v0_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_tup_v0_c2, weylsite, 32);
				bgq_qvlfduxa(weylnext_tup_v1_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_tup_v1_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_tup_v1_c2, weylsite, 32);
				bgq_su3_expand_weyl_tup(spinornext, weylnext_tup);

				bgq_su3_weyl_decl(weyl_tup);
				bgq_su3_reduce_weyl_tdown(weyl_tup, spinor);

				bgq_su3_mdecl(gauge_tup);
				bgq_qvlfduxa(gauge_tup_c00, gaugesite, 32);
				bgq_qvlfduxa(gauge_tup_c01, gaugesite, 32);
				bgq_qvlfduxa(gauge_tup_c02, gaugesite, 32);
				bgq_qvlfduxa(gauge_tup_c10, gaugesite, 32);
				bgq_qvlfduxa(gauge_tup_c11, gaugesite, 32);
				bgq_qvlfduxa(gauge_tup_c12, gaugesite, 32);
				bgq_qvlfduxa(gauge_tup_c20, gaugesite, 32);
				bgq_qvlfduxa(gauge_tup_c21, gaugesite, 32);
				bgq_qvlfduxa(gauge_tup_c22, gaugesite, 32);

				bgq_su3_weyl_mvinvmul(weyl_tup, gauge_tup, weyl_tup);
				if (kamul) {
					bgq_su3_weyl_cmul(weyl_tup, qka0, weyl_tup);
				}

				bgq_su3_weyl_store_double(targetptrs->d[TUP], weyl_tup);
				bgq_weylvec_written(targetptrs->d[TUP], t1, t2, x, y, z, TUP, true);
			}

			bgq_su3_weylnext_prefetch_double(weylsite);
			bgq_su3_matrixnext_prefetch_double(gaugesite);

			// TDOWN
			{
				bgq_su3_weyl_decl(weylnext_tdown);
				bgq_qvlfduxa(weylnext_tdown_v0_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_tdown_v0_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_tdown_v0_c2, weylsite, 32);
				bgq_qvlfduxa(weylnext_tdown_v1_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_tdown_v1_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_tdown_v1_c2, weylsite, 32);
				bgq_su3_accum_weyl_tdown(spinornext, weylnext_tdown);

				bgq_su3_weyl_decl(weyl_tdown);
				bgq_su3_reduce_weyl_tup(weyl_tdown, spinor);

				bgq_su3_mdecl(gauge_tdown);
				bgq_qvlfduxa(gauge_tdown_c00, gaugesite, 32);
				bgq_qvlfduxa(gauge_tdown_c01, gaugesite, 32);
				bgq_qvlfduxa(gauge_tdown_c02, gaugesite, 32);
				bgq_qvlfduxa(gauge_tdown_c10, gaugesite, 32);
				bgq_qvlfduxa(gauge_tdown_c11, gaugesite, 32);
				bgq_qvlfduxa(gauge_tdown_c12, gaugesite, 32);
				bgq_qvlfduxa(gauge_tdown_c20, gaugesite, 32);
				bgq_qvlfduxa(gauge_tdown_c21, gaugesite, 32);
				bgq_qvlfduxa(gauge_tdown_c22, gaugesite, 32);

				bgq_su3_weyl_mvmul(weyl_tdown, gauge_tdown, weyl_tdown);
				if (kamul) {
					bgq_su3_weyl_cmul(weyl_tdown, qka0, weyl_tdown);
				}

				bgq_su3_weyl_store_double(targetptrs->d[TDOWN], weyl_tdown);
				bgq_weylvec_written(targetptrs->d[TDOWN], t1, t2, x, y, z, TDOWN, true);
			}

			// X+ /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weylnext_xup);
				bgq_qvlfduxa(weylnext_xup_v0_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_xup_v0_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_xup_v0_c2, weylsite, 32);
				bgq_qvlfduxa(weylnext_xup_v1_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_xup_v1_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_xup_v1_c2, weylsite, 32);
				bgq_su3_accum_weyl_xup(spinornext, weylnext_xup);

				bgq_su3_weyl_decl(weyl_xup);
				bgq_su3_reduce_weyl_xdown(weyl_xup, spinor);

				bgq_su3_mdecl(gauge_xup);
				bgq_qvlfduxa(gauge_xup_c00, gaugesite, 32);
				bgq_qvlfduxa(gauge_xup_c01, gaugesite, 32);
				bgq_qvlfduxa(gauge_xup_c02, gaugesite, 32);
				bgq_qvlfduxa(gauge_xup_c10, gaugesite, 32);
				bgq_qvlfduxa(gauge_xup_c11, gaugesite, 32);
				bgq_qvlfduxa(gauge_xup_c12, gaugesite, 32);
				bgq_qvlfduxa(gauge_xup_c20, gaugesite, 32);
				bgq_qvlfduxa(gauge_xup_c21, gaugesite, 32);
				bgq_qvlfduxa(gauge_xup_c22, gaugesite, 32);

				bgq_su3_weyl_mvinvmul(weyl_xup, gauge_xup, weyl_xup);
				if (kamul) {
					bgq_su3_weyl_cmul(weyl_xup, qka1, weyl_xup);
				}

				bgq_su3_weyl_store_double(targetptrs->d[XUP], weyl_xup);
				bgq_weylvec_written(targetptrs->d[XUP], t1, t2, x,y,z,XUP, true);
			}

			bgq_su3_weylnext_prefetch_double(weylsite);
			bgq_su3_matrixnext_prefetch_double(gaugesite);

			// X- /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weylnext_xdown);
				bgq_qvlfduxa(weylnext_xdown_v0_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_xdown_v0_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_xdown_v0_c2, weylsite, 32);
				bgq_qvlfduxa(weylnext_xdown_v1_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_xdown_v1_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_xdown_v1_c2, weylsite, 32);
				bgq_su3_accum_weyl_xdown(spinornext, weylnext_xdown);

				bgq_su3_weyl_decl(weyl_xdown);
				bgq_su3_reduce_weyl_xup(weyl_xdown, spinor);

				bgq_su3_mdecl(gauge_xdown);
				bgq_qvlfduxa(gauge_xdown_c00, gaugesite, 32);
				bgq_qvlfduxa(gauge_xdown_c01, gaugesite, 32);
				bgq_qvlfduxa(gauge_xdown_c02, gaugesite, 32);
				bgq_qvlfduxa(gauge_xdown_c10, gaugesite, 32);
				bgq_qvlfduxa(gauge_xdown_c11, gaugesite, 32);
				bgq_qvlfduxa(gauge_xdown_c12, gaugesite, 32);
				bgq_qvlfduxa(gauge_xdown_c20, gaugesite, 32);
				bgq_qvlfduxa(gauge_xdown_c21, gaugesite, 32);
				bgq_qvlfduxa(gauge_xdown_c22, gaugesite, 32);

				bgq_su3_weyl_mvmul(weyl_xdown, gauge_xdown, weyl_xdown);
				if (kamul) {
					bgq_su3_weyl_cmul(weyl_xdown, qka1, weyl_xdown);
				}

				bgq_su3_weyl_store_double(targetptrs->d[XDOWN], weyl_xdown);
				bgq_weylvec_written(targetptrs->d[XDOWN], t1, t2, x,y,z,XDOWN, true);
			}

			bgq_su3_weylnext_prefetch_double(weylsite);
			bgq_su3_matrixnext_prefetch_double(gaugesite);

			// Y+ /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weylnext_yup);
				bgq_qvlfduxa(weylnext_yup_v0_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_yup_v0_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_yup_v0_c2, weylsite, 32);
				bgq_qvlfduxa(weylnext_yup_v1_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_yup_v1_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_yup_v1_c2, weylsite, 32);
				bgq_su3_accum_weyl_yup(spinornext, weylnext_yup);

				bgq_su3_weyl_decl(weyl_yup);
				bgq_su3_reduce_weyl_ydown(weyl_yup, spinor);

				bgq_su3_mdecl(gauge_yup);
				bgq_qvlfduxa(gauge_yup_c00, gaugesite, 32);
				bgq_qvlfduxa(gauge_yup_c01, gaugesite, 32);
				bgq_qvlfduxa(gauge_yup_c02, gaugesite, 32);
				bgq_qvlfduxa(gauge_yup_c10, gaugesite, 32);
				bgq_qvlfduxa(gauge_yup_c11, gaugesite, 32);
				bgq_qvlfduxa(gauge_yup_c12, gaugesite, 32);
				bgq_qvlfduxa(gauge_yup_c20, gaugesite, 32);
				bgq_qvlfduxa(gauge_yup_c21, gaugesite, 32);
				bgq_qvlfduxa(gauge_yup_c22, gaugesite, 32);

				bgq_su3_weyl_mvinvmul(weyl_yup, gauge_yup, weyl_yup);
				if (kamul) {
					bgq_su3_weyl_cmul(weyl_yup, qka2, weyl_yup);
				}

				bgq_su3_weyl_store_double(targetptrs->d[YUP], weyl_yup);
				bgq_weylvec_written(targetptrs->d[YUP], t1, t2, x,y,z,YUP, true);
			}

			bgq_su3_weylnext_prefetch_double(weylsite);
			bgq_su3_matrixnext_prefetch_double(gaugesite);

			// Y- /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weylnext_ydown);
				bgq_qvlfduxa(weylnext_ydown_v0_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_ydown_v0_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_ydown_v0_c2, weylsite, 32);
				bgq_qvlfduxa(weylnext_ydown_v1_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_ydown_v1_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_ydown_v1_c2, weylsite, 32);
				bgq_su3_accum_weyl_ydown(spinornext, weylnext_ydown);

				bgq_su3_weyl_decl(weyl_ydown);
				bgq_su3_reduce_weyl_yup(weyl_ydown, spinor);

				bgq_su3_mdecl(gauge_ydown);
				bgq_qvlfduxa(gauge_ydown_c00, gaugesite, 32);
				bgq_qvlfduxa(gauge_ydown_c01, gaugesite, 32);
				bgq_qvlfduxa(gauge_ydown_c02, gaugesite, 32);
				bgq_qvlfduxa(gauge_ydown_c10, gaugesite, 32);
				bgq_qvlfduxa(gauge_ydown_c11, gaugesite, 32);
				bgq_qvlfduxa(gauge_ydown_c12, gaugesite, 32);
				bgq_qvlfduxa(gauge_ydown_c20, gaugesite, 32);
				bgq_qvlfduxa(gauge_ydown_c21, gaugesite, 32);
				bgq_qvlfduxa(gauge_ydown_c22, gaugesite, 32);

				bgq_su3_weyl_mvmul(weyl_ydown, gauge_ydown, weyl_ydown);
				if (kamul) {
					bgq_su3_weyl_cmul(weyl_ydown, qka2, weyl_ydown);
				}

				bgq_su3_weyl_store_double(targetptrs->d[YDOWN], weyl_ydown);
				bgq_weylvec_written(targetptrs->d[YDOWN], t1, t2, x,y,z,YDOWN, true);
			}

			bgq_su3_weylnext_prefetch_double(weylsite);
			bgq_su3_matrixnext_prefetch_double(gaugesite);

			// Z+ /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weylnext_zup);
				bgq_qvlfduxa(weylnext_zup_v0_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_zup_v0_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_zup_v0_c2, weylsite, 32);
				bgq_qvlfduxa(weylnext_zup_v1_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_zup_v1_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_zup_v1_c2, weylsite, 32);
				bgq_su3_accum_weyl_zup(spinornext, weylnext_zup);

				bgq_su3_weyl_decl(weyl_zup);
				bgq_su3_reduce_weyl_zdown(weyl_zup, spinor);

				bgq_su3_mdecl(gauge_zup);
				bgq_qvlfduxa(gauge_zup_c00, gaugesite, 32);
				bgq_qvlfduxa(gauge_zup_c01, gaugesite, 32);
				bgq_qvlfduxa(gauge_zup_c02, gaugesite, 32);
				bgq_qvlfduxa(gauge_zup_c10, gaugesite, 32);
				bgq_qvlfduxa(gauge_zup_c11, gaugesite, 32);
				bgq_qvlfduxa(gauge_zup_c12, gaugesite, 32);
				bgq_qvlfduxa(gauge_zup_c20, gaugesite, 32);
				bgq_qvlfduxa(gauge_zup_c21, gaugesite, 32);
				bgq_qvlfduxa(gauge_zup_c22, gaugesite, 32);

				bgq_su3_weyl_mvinvmul(weyl_zup, gauge_zup, weyl_zup);
				if (kamul) {
					bgq_su3_weyl_cmul(weyl_zup, qka3, weyl_zup);
				}

				bgq_su3_weyl_store_double(targetptrs->d[ZUP], weyl_zup);
				bgq_weylvec_written(targetptrs->d[ZUP], t1, t2, x,y,z,ZUP, true);
			}

			bgq_su3_weylnext_prefetch_double(weylsite);
			bgq_su3_matrixnext_prefetch_double(gaugesite);

			// Z- /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weylnext_zdown);
				bgq_qvlfduxa(weylnext_zdown_v0_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_zdown_v0_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_zdown_v0_c2, weylsite, 32);
				bgq_qvlfduxa(weylnext_zdown_v1_c0, weylsite, 32);
				bgq_qvlfduxa(weylnext_zdown_v1_c1, weylsite, 32);
				bgq_qvlfduxa(weylnext_zdown_v1_c2, weylsite, 32);
				bgq_su3_accum_weyl_zdown(spinornext, weylnext_zdown);

				bgq_su3_weyl_decl(weyl_zdown);
				bgq_su3_reduce_weyl_zup(weyl_zdown, spinor);

				bgq_su3_mdecl(gauge_zdown);
				bgq_qvlfduxa(gauge_zdown_c00, gaugesite, 32);
				bgq_qvlfduxa(gauge_zdown_c01, gaugesite, 32);
				bgq_qvlfduxa(gauge_zdown_c02, gaugesite, 32);
				bgq_qvlfduxa(gauge_zdown_c10, gaugesite, 32);
				bgq_qvlfduxa(gauge_zdown_c11, gaugesite, 32);
				bgq_qvlfduxa(gauge_zdown_c12, gaugesite, 32);
				bgq_qvlfduxa(gauge_zdown_c20, gaugesite, 32);
				bgq_qvlfduxa(gauge_zdown_c21, gaugesite, 32);
				bgq_qvlfduxa(gauge_zdown_c22, gaugesite, 32);

				bgq_su3_weyl_mvmul(weyl_zdown, gauge_zdown, weyl_zdown);
				if (kamul) {
					bgq_su3_weyl_cmul(weyl_zdown, qka3, weyl_zdown);
				}

				bgq_su3_weyl_store_double(targetptrs->d[ZDOWN], weyl_zdown);
				bgq_weylvec_written(targetptrs->d[ZDOWN], t1, t2, x,y,z,ZDOWN, true);
			}

			bgq_prefetch(targetptrs+1);
			ic_comp = ic_load;
		}

		// Epilogue
		{
			bgq_su3_spinor_decl(spinor);
			bgq_su3_spinor_mov(spinor, spinornext);
			bgq_weyl_ptr_t *targetptrs = &targetfield->sendptr[ic_comp];

			// tup
			{
				bgq_su3_weyl_decl(weyl_tup);
				bgq_su3_reduce_weyl_tdown(weyl_tup, spinor);

				bgq_su3_mdecl(gauge_tup);
				bgq_qvlfduxa(gauge_tup_c00, gaugesite, 32);
				bgq_qvlfduxa(gauge_tup_c01, gaugesite, 32);
				bgq_qvlfduxa(gauge_tup_c02, gaugesite, 32);
				bgq_qvlfduxa(gauge_tup_c10, gaugesite, 32);
				bgq_qvlfduxa(gauge_tup_c11, gaugesite, 32);
				bgq_qvlfduxa(gauge_tup_c12, gaugesite, 32);
				bgq_qvlfduxa(gauge_tup_c20, gaugesite, 32);
				bgq_qvlfduxa(gauge_tup_c21, gaugesite, 32);
				bgq_qvlfduxa(gauge_tup_c22, gaugesite, 32);

				bgq_su3_weyl_mvinvmul(weyl_tup, gauge_tup, weyl_tup);
				if (kamul) {
					bgq_su3_weyl_cmul(weyl_tup, qka0, weyl_tup);
				}

				bgq_su3_weyl_store_double(targetptrs->d[TUP], weyl_tup);
				bgq_weylvec_written(targetptrs->d[TUP], t1, t2, x, y, z, TUP, true);
			}

			// TDOWN
				{
					bgq_su3_weyl_decl(weyl_tdown);
					bgq_su3_reduce_weyl_tup(weyl_tdown, spinor);

					bgq_su3_mdecl(gauge_tdown);
					bgq_qvlfduxa(gauge_tdown_c00, gaugesite, 32);
					bgq_qvlfduxa(gauge_tdown_c01, gaugesite, 32);
					bgq_qvlfduxa(gauge_tdown_c02, gaugesite, 32);
					bgq_qvlfduxa(gauge_tdown_c10, gaugesite, 32);
					bgq_qvlfduxa(gauge_tdown_c11, gaugesite, 32);
					bgq_qvlfduxa(gauge_tdown_c12, gaugesite, 32);
					bgq_qvlfduxa(gauge_tdown_c20, gaugesite, 32);
					bgq_qvlfduxa(gauge_tdown_c21, gaugesite, 32);
					bgq_qvlfduxa(gauge_tdown_c22, gaugesite, 32);

					bgq_su3_weyl_mvmul(weyl_tdown, gauge_tdown, weyl_tdown);
					if (kamul) {
						bgq_su3_weyl_cmul(weyl_tdown, qka0, weyl_tdown);
					}

					bgq_su3_weyl_store_double(targetptrs->d[TDOWN], weyl_tdown);
					bgq_weylvec_written(targetptrs->d[TDOWN], t1, t2, x, y, z, TDOWN, true);
				}

				// X+ /////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_weyl_decl(weyl_xup);
					bgq_su3_reduce_weyl_xdown(weyl_xup, spinor);

					bgq_su3_mdecl(gauge_xup);
					bgq_qvlfduxa(gauge_xup_c00, gaugesite, 32);
					bgq_qvlfduxa(gauge_xup_c01, gaugesite, 32);
					bgq_qvlfduxa(gauge_xup_c02, gaugesite, 32);
					bgq_qvlfduxa(gauge_xup_c10, gaugesite, 32);
					bgq_qvlfduxa(gauge_xup_c11, gaugesite, 32);
					bgq_qvlfduxa(gauge_xup_c12, gaugesite, 32);
					bgq_qvlfduxa(gauge_xup_c20, gaugesite, 32);
					bgq_qvlfduxa(gauge_xup_c21, gaugesite, 32);
					bgq_qvlfduxa(gauge_xup_c22, gaugesite, 32);

					bgq_su3_weyl_mvinvmul(weyl_xup, gauge_xup, weyl_xup);
					if (kamul) {
						bgq_su3_weyl_cmul(weyl_xup, qka1, weyl_xup);
					}

					bgq_su3_weyl_store_double(targetptrs->d[XUP], weyl_xup);
					bgq_weylvec_written(targetptrs->d[XUP], t1, t2, x,y,z,XUP, true);
				}

				// X- /////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_weyl_decl(weyl_xdown);
					bgq_su3_reduce_weyl_xup(weyl_xdown, spinor);

					bgq_su3_mdecl(gauge_xdown);
					bgq_qvlfduxa(gauge_xdown_c00, gaugesite, 32);
					bgq_qvlfduxa(gauge_xdown_c01, gaugesite, 32);
					bgq_qvlfduxa(gauge_xdown_c02, gaugesite, 32);
					bgq_qvlfduxa(gauge_xdown_c10, gaugesite, 32);
					bgq_qvlfduxa(gauge_xdown_c11, gaugesite, 32);
					bgq_qvlfduxa(gauge_xdown_c12, gaugesite, 32);
					bgq_qvlfduxa(gauge_xdown_c20, gaugesite, 32);
					bgq_qvlfduxa(gauge_xdown_c21, gaugesite, 32);
					bgq_qvlfduxa(gauge_xdown_c22, gaugesite, 32);

					bgq_su3_weyl_mvmul(weyl_xdown, gauge_xdown, weyl_xdown);
					if (kamul) {
						bgq_su3_weyl_cmul(weyl_xdown, qka1, weyl_xdown);
					}

					bgq_su3_weyl_store_double(targetptrs->d[XDOWN], weyl_xdown);
					bgq_weylvec_written(targetptrs->d[XDOWN], t1, t2, x,y,z,XDOWN, true);
				}

				// Y+ /////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_weyl_decl(weyl_yup);
					bgq_su3_reduce_weyl_ydown(weyl_yup, spinor);

					bgq_su3_mdecl(gauge_yup);
					bgq_qvlfduxa(gauge_yup_c00, gaugesite, 32);
					bgq_qvlfduxa(gauge_yup_c01, gaugesite, 32);
					bgq_qvlfduxa(gauge_yup_c02, gaugesite, 32);
					bgq_qvlfduxa(gauge_yup_c10, gaugesite, 32);
					bgq_qvlfduxa(gauge_yup_c11, gaugesite, 32);
					bgq_qvlfduxa(gauge_yup_c12, gaugesite, 32);
					bgq_qvlfduxa(gauge_yup_c20, gaugesite, 32);
					bgq_qvlfduxa(gauge_yup_c21, gaugesite, 32);
					bgq_qvlfduxa(gauge_yup_c22, gaugesite, 32);

					bgq_su3_weyl_mvinvmul(weyl_yup, gauge_yup, weyl_yup);
					if (kamul) {
						bgq_su3_weyl_cmul(weyl_yup, qka2, weyl_yup);
					}

					bgq_su3_weyl_store_double(targetptrs->d[YUP], weyl_yup);
					bgq_weylvec_written(targetptrs->d[YUP], t1, t2, x,y,z,YUP, true);
				}

				// Y- /////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_weyl_decl(weyl_ydown);
					bgq_su3_reduce_weyl_yup(weyl_ydown, spinor);

					bgq_su3_mdecl(gauge_ydown);
					bgq_qvlfduxa(gauge_ydown_c00, gaugesite, 32);
					bgq_qvlfduxa(gauge_ydown_c01, gaugesite, 32);
					bgq_qvlfduxa(gauge_ydown_c02, gaugesite, 32);
					bgq_qvlfduxa(gauge_ydown_c10, gaugesite, 32);
					bgq_qvlfduxa(gauge_ydown_c11, gaugesite, 32);
					bgq_qvlfduxa(gauge_ydown_c12, gaugesite, 32);
					bgq_qvlfduxa(gauge_ydown_c20, gaugesite, 32);
					bgq_qvlfduxa(gauge_ydown_c21, gaugesite, 32);
					bgq_qvlfduxa(gauge_ydown_c22, gaugesite, 32);

					bgq_su3_weyl_mvmul(weyl_ydown, gauge_ydown, weyl_ydown);
					if (kamul) {
						bgq_su3_weyl_cmul(weyl_ydown, qka2, weyl_ydown);
					}

					bgq_su3_weyl_store_double(targetptrs->d[YDOWN], weyl_ydown);
					bgq_weylvec_written(targetptrs->d[YDOWN], t1, t2, x,y,z,YDOWN, true);
				}

				// Z+ /////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_weyl_decl(weyl_zup);
					bgq_su3_reduce_weyl_zdown(weyl_zup, spinor);

					bgq_su3_mdecl(gauge_zup);
					bgq_qvlfduxa(gauge_zup_c00, gaugesite, 32);
					bgq_qvlfduxa(gauge_zup_c01, gaugesite, 32);
					bgq_qvlfduxa(gauge_zup_c02, gaugesite, 32);
					bgq_qvlfduxa(gauge_zup_c10, gaugesite, 32);
					bgq_qvlfduxa(gauge_zup_c11, gaugesite, 32);
					bgq_qvlfduxa(gauge_zup_c12, gaugesite, 32);
					bgq_qvlfduxa(gauge_zup_c20, gaugesite, 32);
					bgq_qvlfduxa(gauge_zup_c21, gaugesite, 32);
					bgq_qvlfduxa(gauge_zup_c22, gaugesite, 32);

					bgq_su3_weyl_mvinvmul(weyl_zup, gauge_zup, weyl_zup);
					if (kamul) {
						bgq_su3_weyl_cmul(weyl_zup, qka3, weyl_zup);
					}

					bgq_su3_weyl_store_double(targetptrs->d[ZUP], weyl_zup);
					bgq_weylvec_written(targetptrs->d[ZUP], t1, t2, x,y,z,ZUP, true);
				}

				// Z- /////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_weyl_decl(weyl_zdown);
					bgq_su3_reduce_weyl_zup(weyl_zdown, spinor);

					bgq_su3_mdecl(gauge_zdown);
					bgq_qvlfduxa(gauge_zdown_c00, gaugesite, 32);
					bgq_qvlfduxa(gauge_zdown_c01, gaugesite, 32);
					bgq_qvlfduxa(gauge_zdown_c02, gaugesite, 32);
					bgq_qvlfduxa(gauge_zdown_c10, gaugesite, 32);
					bgq_qvlfduxa(gauge_zdown_c11, gaugesite, 32);
					bgq_qvlfduxa(gauge_zdown_c12, gaugesite, 32);
					bgq_qvlfduxa(gauge_zdown_c20, gaugesite, 32);
					bgq_qvlfduxa(gauge_zdown_c21, gaugesite, 32);
					bgq_qvlfduxa(gauge_zdown_c22, gaugesite, 32);

					bgq_su3_weyl_mvmul(weyl_zdown, gauge_zdown, weyl_zdown);
					if (kamul) {
						bgq_su3_weyl_cmul(weyl_zdown, qka3, weyl_zdown);
					}

					bgq_su3_weyl_store_double(targetptrs->d[ZDOWN], weyl_zdown);
					bgq_weylvec_written(targetptrs->d[ZDOWN], t1, t2, x,y,z,ZDOWN, true);
				}
		}
#endif
	}
}


#undef BGQ_HOPPINGMATRIXWORKER_INC_
