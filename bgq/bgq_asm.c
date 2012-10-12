/*
 * bgq_asm.c
 *
 *  Created on: Oct 8, 2012
 *      Author: meinersbur
 */

#include <stdbool.h>
#include <omp.h>

#include "bgq.h"
#include "bgq_utils.h"
#include "bgq_HoppingMatrix.h"

#include "bgq_field_double.h"

#define BGQ_PRECISION 64
#include "bgq_precisionselect.inc.c"

#define BGQ_HM_NOKAMUL 1
#define BGQ_HM_NOFUNC 1
#define BGQ_HM_ZLINE_NOFUNC 1
#define BGQ_HM_SITE_NOFUNC 1
#define BGQ_HM_DIR_NOFUNC 1

#define BGQ_DCBT_DISABLE 1

#if 1
#undef BGQ_SPINORSITE
#define BGQ_SPINORSITE(spinorfield, isOdd, tv, x, y, z, t1, t2, isRead, isWrite) (bgq_spinorsite*)(isRead ? /*read*/(tmp=addr,addr+=2*4*3*2,tmp) : (isWrite ? /*write*/(tmp=writeaddr,writeaddr+=2*4*3*2,tmp) : /*prefetch*/(addr) ) )
#undef BGQ_SPINORSITE_LEFT
#define BGQ_SPINORSITE_LEFT BGQ_SPINORSITE
#undef BGQ_SPINORSITE_RIGHT
#define BGQ_SPINORSITE_RIGHT BGQ_SPINORSITE

#undef BGQ_WEYLSITE_T
#define BGQ_WEYLSITE_T(weylfield, isOdd, t, xv, y, z, x1, x2, isRead, isWrite) (bgq_weylsite*)(isRead ? /*read*/(tmp=addr,addr+=2*2*3*2,tmp) : (isWrite ? /*write*/(tmp=writeaddr,writeaddr+=2*2*3*2,tmp) : /*prefetch*/(addr) ) )
#undef BGQ_WEYLSITE_X
#define BGQ_WEYLSITE_X BGQ_WEYLSITE_T
#undef BGQ_WEYLSITE_Y
#define BGQ_WEYLSITE_Y BGQ_WEYLSITE_T
#undef BGQ_WEYLSITE_Z
#define BGQ_WEYLSITE_Z BGQ_WEYLSITE_T
#undef BGQ_WEYLSITE_T_LEFT
#define BGQ_WEYLSITE_T_LEFT BGQ_WEYLSITE_T
#undef BGQ_WEYLSITE_T_RIGHT
#define BGQ_WEYLSITE_T_RIGHT BGQ_WEYLSITE_T

#undef BGQ_GAUGESITE
#define BGQ_GAUGESITE(gaugefield,isOdd,tv,x,y,z,dir,t1,t2,isRead,isWrite) (bgq_gaugesite*)(isRead ? /*read*/(tmp=addr,addr+=2*3*3*2,tmp) : (isWrite ? /*write*/(tmp=writeaddr,writeaddr+=2*3*3*2,tmp) : /*prefetch*/(addr) ) )
#undef BGQ_GAUGESITE_LEFT
#define BGQ_GAUGESITE_LEFT BGQ_GAUGESITE
#undef BGQ_GAUGESITE_RIGHT
#define BGQ_GAUGESITE_RIGHT BGQ_GAUGESITE
#endif

// To disable dcbt
//#define bgq_su3_spinor_prefetch(x)
//#define bgq_su3_matrix_prefetch(x)

void bgq_HoppingMatrix_asm(bool isOdd, bgq_spinorfield restrict targetfield, bgq_spinorfield restrict spinorfield, bgq_gaugefield restrict gaugefield, bgq_hmflags opts) {
	isOdd = false;
	const bool nocom = true;
	const bool nooverlap = true;
	const bool nokamul = true;
	const bool noprefetchlist = false;
	const bool noprefetchstream = true;
	const bool noprefetchexplicit = false;
	const bool noweylsend = true;
	const bool nobody = false;
	const bool nosurface = true;

//#pragma disjoint(targetfield, spinorfield, gaugefield) // Unnecessary due to restrict keyword
__alignx(BGQ_PRECISION,targetfield);
__alignx(BGQ_PRECISION,spinorfield);
__alignx(BGQ_PRECISION,gaugefield);

#pragma omp parallel
	{
		const int threads = omp_get_num_threads();
		const int tid = omp_get_thread_num();

		PRECISION *tmp;
		PRECISION *addr = ((PRECISION*) spinorfield) + (size_t) tid * READTOTLENGTH / threads;
		addr = (PRECISION*) ((size_t) addr & ~(BGQ_PRECISION - 1)); // alignment

		PRECISION *writeaddr = ((PRECISION*) targetfield) + (size_t) tid * WRITETOTLENGTH / threads;
		writeaddr = (PRECISION*) ((size_t) writeaddr & ~(BGQ_PRECISION - 1)); // alignment

#pragma omp for schedule(static)
		for (int txy = 0; txy < VOLUME/2; txy += 2) {
			const int tv = 0;
			const int x = 0;
			const int y = 0;
			const int z = 0;

			int t1 = tv * PHYSICAL_LP * PHYSICAL_LK + ((isOdd + x + y + 0) & 1);
			int t2 = t1 + PHYSICAL_LP;

			{
				bgq_su3_spinor_decl(result);
				bgq_su3_spinor_zero(result);

				bgq_su3_spinor_prefetch(addr + 2*4*3*2 + 2*3*3*2);
				bgq_su3_matrix_prefetch(addr + 2*4*3*2 + 2*3*3*2 + 2*4*3*2);

				// T+ //////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_spinor_decl(spinor_tup);
					bgq_su3_weyl_decl(weyl_tup);
					bgq_su3_mdecl(gauge_tup);

					bgq_spinorsite *spinorsite_tup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z, t1+1, t2+1, true,false);
					bgq_su3_spinor_load(spinor_tup, spinorsite_tup);
					//bgq_su3_spinor_flush(spinorsite_tup);
					bgq_su3_reduce_weyl_tup(weyl_tup, spinor_tup);

					bgq_gaugesite *gaugesite_tup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, TUP, t1, t2, true,false);
					bgq_su3_matrix_load(gauge_tup, gaugesite_tup);
					//bgq_su3_matrix_flush(gaugesite_tup);
					bgq_su3_weyl_mvmul(weyl_tup,gauge_tup,weyl_tup);

					bgq_su3_expand_weyl_tup(result,weyl_tup);
				}

				bgq_su3_spinor_prefetch(addr + 2*4*3*2 + 2*3*3*2);
				bgq_su3_matrix_prefetch(addr + 2*4*3*2 + 2*3*3*2 + 2*4*3*2);

				// T- //////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_spinor_decl(spinor_tdown);
					bgq_su3_weyl_decl(weyl_tdown);
					bgq_su3_mdecl(gauge_tdown);

					bgq_spinorsite *spinorsite_tdown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z, t1-1, t2-1, true,false);
					bgq_su3_spinor_load(spinor_tdown, spinorsite_tdown);
					//bgq_su3_spinor_flush(spinorsite_tdown);
					bgq_su3_reduce_weyl_tdown(weyl_tdown, spinor_tdown);

					bgq_gaugesite *gaugesite_tdown = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, TUP, t1-1, t2-1, true,false);
					bgq_su3_matrix_load(gauge_tdown, gaugesite_tdown);
					//bgq_su3_matrix_flush(gaugesite_tdown);
					bgq_su3_weyl_mvinvmul(weyl_tdown, gauge_tdown, weyl_tdown);

					bgq_su3_accum_weyl_tdown(result, weyl_tdown);
				}

				bgq_su3_spinor_prefetch(addr + 2*4*3*2 + 2*3*3*2);
				bgq_su3_matrix_prefetch(addr + 2*4*3*2 + 2*3*3*2 + 2*4*3*2);

				// X+ //////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_spinor_decl(spinor_xup);
					bgq_su3_weyl_decl(weyl_xup);
					bgq_su3_mdecl(gauge_xup);

					bgq_spinorsite *spinorsite_xup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x+1, y, z, t1, t2, true,false);
					bgq_su3_spinor_load(spinor_xup, spinorsite_xup);
					//bgq_su3_spinor_flush(spinorsite_xup);
					bgq_su3_reduce_weyl_xup(weyl_xup, spinor_xup);

					bgq_gaugesite *gaugesite_xup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, XUP, t1, t2, true,false);
					bgq_su3_matrix_load(gauge_xup, gaugesite_xup);
					//bgq_su3_matrix_flush(gaugesite_xup);
					bgq_su3_weyl_mvmul(weyl_xup, gauge_xup, weyl_xup);

					bgq_su3_accum_weyl_xup(result, weyl_xup);
				}

				bgq_su3_spinor_prefetch(addr + 2*4*3*2 + 2*3*3*2);
				bgq_su3_matrix_prefetch(addr + 2*4*3*2 + 2*3*3*2 + 2*4*3*2);

				// X- //////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_spinor_decl(spinor_xdown);
					bgq_su3_weyl_decl(weyl_xdown);
					bgq_su3_mdecl(gauge_xdown);

					bgq_spinorsite *spinorsite_xdown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x-1, y, z, t1, t2, true,false);
					bgq_su3_spinor_load(spinor_xdown, spinorsite_xdown);
					//bgq_su3_spinor_flush(spinorsite_xdown);
					bgq_su3_reduce_weyl_xdown(weyl_xdown, spinor_xdown);

					bgq_gaugesite *gaugesite_xdown = BGQ_GAUGESITE(gaugefield, isOdd, tv, x-1, y, z, XUP, t1, t2, true,false);
					bgq_su3_matrix_load(gauge_xdown, gaugesite_xdown);
					//bgq_su3_matrix_flush(gaugesite_xdown);
					bgq_su3_weyl_mvinvmul(weyl_xdown, gauge_xdown, weyl_xdown);

					bgq_su3_accum_weyl_xdown(result, weyl_xdown);
				}

				bgq_su3_spinor_prefetch(addr + 2*4*3*2 + 2*3*3*2);
				bgq_su3_matrix_prefetch(addr + 2*4*3*2 + 2*3*3*2 + 2*4*3*2);

				// Y+ //////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_spinor_decl(spinor_yup);
					bgq_su3_weyl_decl(weyl_yup);
					bgq_su3_mdecl(gauge_yup);

					bgq_spinorsite *spinorsite_yup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y+1, z, t1, t2, true,false);
					bgq_su3_spinor_load(spinor_yup, spinorsite_yup);
					//bgq_su3_spinor_flush(spinorsite_yup);
					bgq_su3_reduce_weyl_yup(weyl_yup, spinor_yup);

					bgq_gaugesite *gaugesite_yup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, YUP, t1, t2, true,false);
					bgq_su3_matrix_load(gauge_yup, gaugesite_yup);
					//bgq_su3_matrix_flush(gaugesite_yup);
					bgq_su3_weyl_mvmul(weyl_yup, gauge_yup, weyl_yup);

					bgq_su3_accum_weyl_yup(result, weyl_yup);
				}

				bgq_su3_spinor_prefetch(addr + 2*4*3*2 + 2*3*3*2);
				bgq_su3_matrix_prefetch(addr + 2*4*3*2 + 2*3*3*2 + 2*4*3*2);

				// Y- //////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_spinor_decl(spinor_ydown);
					bgq_su3_weyl_decl(weyl_ydown);
					bgq_su3_mdecl(gauge_ydown);

					bgq_spinorsite *spinorsite_ydown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y-1, z, t1, t2, true,false);
					bgq_su3_spinor_load(spinor_ydown, spinorsite_ydown);
					//bgq_su3_spinor_flush(spinorsite_ydown);
					bgq_su3_reduce_weyl_ydown(weyl_ydown, spinor_ydown);

					bgq_gaugesite *gaugesite_ydown = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y-1, z, YUP, t1, t2, true,false);
					bgq_su3_matrix_load(gauge_ydown, gaugesite_ydown);
					//bgq_su3_matrix_flush(gaugesite_ydown);
					bgq_su3_weyl_mvinvmul(weyl_ydown, gauge_ydown, weyl_ydown);

					bgq_su3_accum_weyl_ydown(result, weyl_ydown);
				}

				bgq_su3_spinor_prefetch(addr + 2*4*3*2 + 2*3*3*2);
				bgq_su3_matrix_prefetch(addr + 2*4*3*2 + 2*3*3*2 + 2*4*3*2);

				// Z+ //////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_spinor_decl(spinor_zup);
					bgq_su3_weyl_decl(weyl_zup);
					bgq_su3_mdecl(gauge_zup);

					bgq_spinorsite *spinorsite_zup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z+1, t1, t2, true,false);
					bgq_su3_spinor_load(spinor_zup, spinorsite_zup);
					//bgq_su3_spinor_flush(spinorsite_zup);
					bgq_su3_reduce_weyl_zup(weyl_zup, spinor_zup);

					bgq_gaugesite *gaugesite_zup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, ZUP, t1, t2, true,false);
					bgq_su3_matrix_load(gauge_zup, gaugesite_zup);
					//bgq_su3_matrix_flush(gaugesite_zup);
					bgq_su3_weyl_mvmul(weyl_zup, gauge_zup, weyl_zup);

					bgq_su3_accum_weyl_zup(result, weyl_zup);
				}

				bgq_su3_spinor_prefetch(addr + 2*4*3*2 + 2*3*3*2);
				bgq_su3_matrix_prefetch(addr + 2*4*3*2 + 2*3*3*2 + 2*4*3*2);

				// Z- //////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_spinor_decl(spinor_zdown);
					bgq_su3_weyl_decl(weyl_zdown);
					bgq_su3_mdecl(gauge_zdown);

					bgq_spinorsite *spinorsite_zdown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z-1, t1, t2, true, false);
					bgq_su3_spinor_load(spinor_zdown, spinorsite_zdown);
					//bgq_su3_spinor_flush(spinorsite_zdown);
					bgq_su3_reduce_weyl_zdown(weyl_zdown, spinor_zdown);

					bgq_gaugesite *gaugesite_zdown = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z-1, ZUP, t1, t2, true, false);
					bgq_su3_matrix_load(gauge_zdown, gaugesite_zdown);
					//bgq_su3_matrix_flush(gaugesite_zdown);
					bgq_su3_weyl_mvinvmul(weyl_zdown, gauge_zdown, weyl_zdown);

					bgq_su3_accum_weyl_zdown(result, weyl_zdown);
				}

				// Store the result
				bgq_spinorsite *targetsite = BGQ_SPINORSITE(targetfield, isOdd, tv, x, y, z, t1,t2, false, true);
				bgq_su3_spinor_zeroload(targetsite);
				bgq_su3_spinor_store(targetsite, result);
				bgq_su3_spinor_flush(targetsite);
			}

		}

	}
}






void bgq_HoppingMatrix_asm_parallel(bool isOdd, bgq_spinorfield restrict targetfield, bgq_spinorfield restrict spinorfield, bgq_gaugefield restrict gaugefield, bgq_hmflags opts, int tid, int tot_threads) {
	isOdd = false;
	const bool nocom = true;
	const bool nooverlap = true;
	const bool nokamul = true;
	const bool noprefetchlist = false;
	const bool noprefetchstream = true;
	const bool noprefetchexplicit = false;
	const bool noweylsend = true;
	const bool nobody = false;
	const bool nosurface = true;

//#pragma disjoint(targetfield, spinorfield, gaugefield) // Unnecessary due to restrict keyword
__alignx(BGQ_PRECISION,targetfield);
__alignx(BGQ_PRECISION,spinorfield);
__alignx(BGQ_PRECISION,gaugefield);

//#pragma omp parallel
	{
		const int threads = omp_get_num_threads();
		const int tid = omp_get_thread_num();

		PRECISION *tmp;
		PRECISION *addr = ((PRECISION*) spinorfield) + (size_t) tid * READTOTLENGTH / tot_threads;
		addr = (PRECISION*) ((size_t) addr & ~(BGQ_PRECISION - 1)); // alignment

		PRECISION *writeaddr = ((PRECISION*) targetfield) + (size_t) tid * WRITETOTLENGTH / tot_threads;
		writeaddr = (PRECISION*) ((size_t) writeaddr & ~(BGQ_PRECISION - 1)); // alignment

//#pragma omp for schedule(static)
		size_t count = VOLUME/2;
		for (size_t txy = 0; txy < count/tot_threads; txy += 2) {
			const int tv = 0;
			const int x = 0;
			const int y = 0;
			const int z = 0;

			int t1 = tv * PHYSICAL_LP * PHYSICAL_LK + ((isOdd + x + y + 0) & 1);
			int t2 = t1 + PHYSICAL_LP;

			{
				bgq_su3_spinor_decl(result);
				bgq_su3_spinor_zero(result);

				bgq_su3_spinor_prefetch(addr + 2*4*3*2 + 2*3*3*2);
				bgq_su3_matrix_prefetch(addr + 2*4*3*2 + 2*3*3*2 + 2*4*3*2);

				// T+ //////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_spinor_decl(spinor_tup);
					bgq_su3_weyl_decl(weyl_tup);
					bgq_su3_mdecl(gauge_tup);

					bgq_spinorsite *spinorsite_tup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z, t1+1, t2+1, true,false);
					bgq_su3_spinor_load(spinor_tup, spinorsite_tup);
					//bgq_su3_spinor_flush(spinorsite_tup);
					bgq_su3_reduce_weyl_tup(weyl_tup, spinor_tup);

					bgq_gaugesite *gaugesite_tup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, TUP, t1, t2, true,false);
					bgq_su3_matrix_load(gauge_tup, gaugesite_tup);
					//bgq_su3_matrix_flush(gaugesite_tup);
					bgq_su3_weyl_mvmul(weyl_tup,gauge_tup,weyl_tup);

					bgq_su3_expand_weyl_tup(result,weyl_tup);
				}

				bgq_su3_spinor_prefetch(addr + 2*4*3*2 + 2*3*3*2);
				bgq_su3_matrix_prefetch(addr + 2*4*3*2 + 2*3*3*2 + 2*4*3*2);

				// T- //////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_spinor_decl(spinor_tdown);
					bgq_su3_weyl_decl(weyl_tdown);
					bgq_su3_mdecl(gauge_tdown);

					bgq_spinorsite *spinorsite_tdown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z, t1-1, t2-1, true,false);
					bgq_su3_spinor_load(spinor_tdown, spinorsite_tdown);
					//bgq_su3_spinor_flush(spinorsite_tdown);
					bgq_su3_reduce_weyl_tdown(weyl_tdown, spinor_tdown);

					bgq_gaugesite *gaugesite_tdown = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, TUP, t1-1, t2-1, true,false);
					bgq_su3_matrix_load(gauge_tdown, gaugesite_tdown);
					//bgq_su3_matrix_flush(gaugesite_tdown);
					bgq_su3_weyl_mvinvmul(weyl_tdown, gauge_tdown, weyl_tdown);

					bgq_su3_accum_weyl_tdown(result, weyl_tdown);
				}

				bgq_su3_spinor_prefetch(addr + 2*4*3*2 + 2*3*3*2);
				bgq_su3_matrix_prefetch(addr + 2*4*3*2 + 2*3*3*2 + 2*4*3*2);

				// X+ //////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_spinor_decl(spinor_xup);
					bgq_su3_weyl_decl(weyl_xup);
					bgq_su3_mdecl(gauge_xup);

					bgq_spinorsite *spinorsite_xup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x+1, y, z, t1, t2, true,false);
					bgq_su3_spinor_load(spinor_xup, spinorsite_xup);
					//bgq_su3_spinor_flush(spinorsite_xup);
					bgq_su3_reduce_weyl_xup(weyl_xup, spinor_xup);

					bgq_gaugesite *gaugesite_xup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, XUP, t1, t2, true,false);
					bgq_su3_matrix_load(gauge_xup, gaugesite_xup);
					//bgq_su3_matrix_flush(gaugesite_xup);
					bgq_su3_weyl_mvmul(weyl_xup, gauge_xup, weyl_xup);

					bgq_su3_accum_weyl_xup(result, weyl_xup);
				}

				bgq_su3_spinor_prefetch(addr + 2*4*3*2 + 2*3*3*2);
				bgq_su3_matrix_prefetch(addr + 2*4*3*2 + 2*3*3*2 + 2*4*3*2);

				// X- //////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_spinor_decl(spinor_xdown);
					bgq_su3_weyl_decl(weyl_xdown);
					bgq_su3_mdecl(gauge_xdown);

					bgq_spinorsite *spinorsite_xdown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x-1, y, z, t1, t2, true,false);
					bgq_su3_spinor_load(spinor_xdown, spinorsite_xdown);
					//bgq_su3_spinor_flush(spinorsite_xdown);
					bgq_su3_reduce_weyl_xdown(weyl_xdown, spinor_xdown);

					bgq_gaugesite *gaugesite_xdown = BGQ_GAUGESITE(gaugefield, isOdd, tv, x-1, y, z, XUP, t1, t2, true,false);
					bgq_su3_matrix_load(gauge_xdown, gaugesite_xdown);
					//bgq_su3_matrix_flush(gaugesite_xdown);
					bgq_su3_weyl_mvinvmul(weyl_xdown, gauge_xdown, weyl_xdown);

					bgq_su3_accum_weyl_xdown(result, weyl_xdown);
				}

				bgq_su3_spinor_prefetch(addr + 2*4*3*2 + 2*3*3*2);
				bgq_su3_matrix_prefetch(addr + 2*4*3*2 + 2*3*3*2 + 2*4*3*2);

				// Y+ //////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_spinor_decl(spinor_yup);
					bgq_su3_weyl_decl(weyl_yup);
					bgq_su3_mdecl(gauge_yup);

					bgq_spinorsite *spinorsite_yup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y+1, z, t1, t2, true,false);
					bgq_su3_spinor_load(spinor_yup, spinorsite_yup);
					//bgq_su3_spinor_flush(spinorsite_yup);
					bgq_su3_reduce_weyl_yup(weyl_yup, spinor_yup);

					bgq_gaugesite *gaugesite_yup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, YUP, t1, t2, true,false);
					bgq_su3_matrix_load(gauge_yup, gaugesite_yup);
					//bgq_su3_matrix_flush(gaugesite_yup);
					bgq_su3_weyl_mvmul(weyl_yup, gauge_yup, weyl_yup);

					bgq_su3_accum_weyl_yup(result, weyl_yup);
				}

				bgq_su3_spinor_prefetch(addr + 2*4*3*2 + 2*3*3*2);
				bgq_su3_matrix_prefetch(addr + 2*4*3*2 + 2*3*3*2 + 2*4*3*2);

				// Y- //////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_spinor_decl(spinor_ydown);
					bgq_su3_weyl_decl(weyl_ydown);
					bgq_su3_mdecl(gauge_ydown);

					bgq_spinorsite *spinorsite_ydown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y-1, z, t1, t2, true,false);
					bgq_su3_spinor_load(spinor_ydown, spinorsite_ydown);
					//bgq_su3_spinor_flush(spinorsite_ydown);
					bgq_su3_reduce_weyl_ydown(weyl_ydown, spinor_ydown);

					bgq_gaugesite *gaugesite_ydown = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y-1, z, YUP, t1, t2, true,false);
					bgq_su3_matrix_load(gauge_ydown, gaugesite_ydown);
					//bgq_su3_matrix_flush(gaugesite_ydown);
					bgq_su3_weyl_mvinvmul(weyl_ydown, gauge_ydown, weyl_ydown);

					bgq_su3_accum_weyl_ydown(result, weyl_ydown);
				}

				bgq_su3_spinor_prefetch(addr + 2*4*3*2 + 2*3*3*2);
				bgq_su3_matrix_prefetch(addr + 2*4*3*2 + 2*3*3*2 + 2*4*3*2);

				// Z+ //////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_spinor_decl(spinor_zup);
					bgq_su3_weyl_decl(weyl_zup);
					bgq_su3_mdecl(gauge_zup);

					bgq_spinorsite *spinorsite_zup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z+1, t1, t2, true,false);
					bgq_su3_spinor_load(spinor_zup, spinorsite_zup);
					//bgq_su3_spinor_flush(spinorsite_zup);
					bgq_su3_reduce_weyl_zup(weyl_zup, spinor_zup);

					bgq_gaugesite *gaugesite_zup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, ZUP, t1, t2, true,false);
					bgq_su3_matrix_load(gauge_zup, gaugesite_zup);
					//bgq_su3_matrix_flush(gaugesite_zup);
					bgq_su3_weyl_mvmul(weyl_zup, gauge_zup, weyl_zup);

					bgq_su3_accum_weyl_zup(result, weyl_zup);
				}

				bgq_su3_spinor_prefetch(addr + 2*4*3*2 + 2*3*3*2);
				bgq_su3_matrix_prefetch(addr + 2*4*3*2 + 2*3*3*2 + 2*4*3*2);

				// Z- //////////////////////////////////////////////////////////////////////////
				{
					bgq_su3_spinor_decl(spinor_zdown);
					bgq_su3_weyl_decl(weyl_zdown);
					bgq_su3_mdecl(gauge_zdown);

					bgq_spinorsite *spinorsite_zdown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z-1, t1, t2, true, false);
					bgq_su3_spinor_load(spinor_zdown, spinorsite_zdown);
					//bgq_su3_spinor_flush(spinorsite_zdown);
					bgq_su3_reduce_weyl_zdown(weyl_zdown, spinor_zdown);

					bgq_gaugesite *gaugesite_zdown = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z-1, ZUP, t1, t2, true, false);
					bgq_su3_matrix_load(gauge_zdown, gaugesite_zdown);
					//bgq_su3_matrix_flush(gaugesite_zdown);
					bgq_su3_weyl_mvinvmul(weyl_zdown, gauge_zdown, weyl_zdown);

					bgq_su3_accum_weyl_zdown(result, weyl_zdown);
				}

				// Store the result
				//for (int i = 0; i < 4; i+=1) {
				bgq_spinorsite *targetsite = BGQ_SPINORSITE(targetfield, isOdd, tv, x, y, z, t1,t2, false, true);
				bgq_su3_spinor_zeroload(targetsite); /* no impact on performance */
				bgq_su3_spinor_store(targetsite, result);
				//bgq_su3_spinor_flush(targetsite); /* reduces performance */
				//}
			}

		}

	}
}
