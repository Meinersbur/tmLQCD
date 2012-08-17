
#ifndef BGQ_HM_TUP_FIRST
#define BGQ_HM_TUP_FIRST 0
#endif

#ifndef BGQ_HM_TUP_TLINEINDENT
#define BGQ_HM_TUP_TLINEINDENT -1
#endif

#ifndef BGQ_HM_TUP_PREFETCH
#define BGQ_HM_TUP_PREFETCH 0
#endif

#ifndef BGQ_HM_TUP_WEYLREAD
dfg
#define BGQ_HM_TUP_WEYLREAD 0
#endif

#ifndef BGQ_HM_TUP_COMPUTE
#define BGQ_HM_TUP_COMPUTE 0
#endif

#ifndef BGQ_HM_TUP_ACCUMULATE
#define BGQ_HM_TUP_ACCUMULATE 0
#endif


#define BGQ_LOADORPREFETCH_PREFETCH BGQ_HM_TUP_PREFETCH
#define BGQ_LOADORPREFETCH_LOAD !BGQ_HM_TUP_PREFETCH
#include "bgq_loadorprefetch.inc.c"


#ifndef BGQ_HM_DIR_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void HoppingMatrix_tup(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int tv, int x, int y, int z, int t1, int t2) {
	bgq_su3_spinor_decl(result);
#else
{
#endif


	bgq_su3_weyl_decl(weyl_tup);
	#if (BGQ_HM_TUP_WEYLREAD==-1)
	if ( tv != PHYSICAL_LTV-1 ) {
	#endif
	#if (BGQ_HM_TUP_WEYLREAD==-1) || (BGQ_HM_TUP_WEYLREAD==0)
		// Load the input spinor
		bgq_su3_spinor_decl(spinor_tup);
		#if (BGQ_HM_TUP_TLINEINDENT==-1)
		if (((tv+x+y+z)&1) == isOdd) {
		#endif
		#if (BGQ_HM_TUP_TLINEINDENT==-1) || (BGQ_HM_TUP_TLINEINDENT==0)
			bgq_spinorsite *spinorsite_tup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z, t1+1, t2+1, !BGQ_HM_TUP_PREFETCH,false);
			bgq_su3_spinor_loadorprefetch(spinor_tup, spinorsite_tup);
		#endif
		#if (BGQ_HM_TUP_TLINEINDENT==-1)
		} else {
		#endif
		#if (BGQ_HM_TUP_TLINEINDENT==-1) || (BGQ_HM_TUP_TLINEINDENT==1)
			bgq_spinorsite *spinorsite_tup_mid = BGQ_SPINORSITE_RIGHT(spinorfield, !isOdd, tv, x, y, z, t1-1, t1+1, !BGQ_HM_TUP_PREFETCH,false);
			bgq_su3_spinor_decl(spinor_tup_mid);
			bgq_su3_spinor_loadorprefetch_right(spinor_tup_mid, spinorsite_tup_mid);

			bgq_spinorsite *spinorsite_tup_right = BGQ_SPINORSITE_LEFT(spinorfield, !isOdd, tv+1, x, y, z, t2+1, t2+3, !BGQ_HM_TUP_PREFETCH,false);
			bgq_su3_spinor_decl(spinor_tup_right);
			bgq_su3_spinor_loadorprefetch_left(spinor_tup_right, spinorsite_tup_right);

			bgq_su3_spinor_merge(spinor_tup, spinor_tup_mid, spinor_tup_right);
		#endif
		#if (BGQ_HM_TUP_TLINEINDENT==-1)
		}
		#endif

		// Compute its halfspinor
		bgq_su3_vadd(weyl_tup_v0, spinor_tup_v0, spinor_tup_v2);
		bgq_su3_vadd(weyl_tup_v1, spinor_tup_v1, spinor_tup_v3);
	#endif
	#if (BGQ_HM_TUP_WEYLREAD==-1)
	} else {
	#endif
	#if (BGQ_HM_TUP_WEYLREAD==-1) || (BGQ_HM_TUP_WEYLREAD==1)
		#if (BGQ_HM_TUP_TLINEINDENT==-1)
		if ( ((tv+x+y+z)&1) == isOdd ) {
		#endif
		#if (BGQ_HM_TUP_TLINEINDENT==-1) || (BGQ_HM_TUP_TLINEINDENT==0)
			assert(t2+2 == LOCAL_LT);
			// Read as normal
			bgq_spinorsite *spinorsite_tup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z, t1+1, t2+1, !BGQ_HM_TUP_PREFETCH,false);
			bgq_su3_spinor_decl(spinor_tup);
			bgq_su3_spinor_loadorprefetch(spinor_tup, spinorsite_tup);

			// Compute its halfspinor
			bgq_su3_vadd(weyl_tup_v0, spinor_tup_v0, spinor_tup_v2);
			bgq_su3_vadd(weyl_tup_v1, spinor_tup_v1, spinor_tup_v3);
		#endif
		#if (BGQ_HM_TUP_TLINEINDENT==-1)
		} else {
		#endif
		#if (BGQ_HM_TUP_TLINEINDENT==-1) || (BGQ_HM_TUP_TLINEINDENT==1)
			assert(t2+1 == LOCAL_LT);
			// i.e. the rightmost value must read tup from the weyl received data

			// Read in the left part of the site like normal
			bgq_spinorsite *spinorsite_tup_mid = BGQ_SPINORSITE_RIGHT(spinorfield, !isOdd, tv, x, y, z, t1-1, t1+1, !BGQ_HM_TUP_PREFETCH,false);
			bgq_su3_spinor_decl(spinor_tup_mid);
			bgq_su3_spinor_loadorprefetch_left(spinor_tup_mid, spinorsite_tup_mid);

			// Also compute its halfspinor, but this we have only a right part, hence wasting a few flops
			bgq_su3_weyl_decl(weyl_tup_mid);
			bgq_su3_vadd(weyl_tup_mid_v0, spinor_tup_mid_v0, spinor_tup_mid_v2);
			bgq_su3_vadd(weyl_tup_mid_v1, spinor_tup_mid_v1, spinor_tup_mid_v3);

			// Read the other halfspinor from what the neighbor node sent us
			const int xeo = x / PHYSICAL_LP;
			const int xv = xeo / PHYSICAL_LK;
			const int kx = mod(xeo, PHYSICAL_LK);

			bgq_weylsite *weylsite_tup_right = BGQ_WEYLSITE_T(weylxchange_recv[TUP], !isOdd, t2+1, xv, y, z, kx==0 ? x : x-2, kx==1 ? x : x+2, !BGQ_HM_TUP_PREFETCH, false);
			weylsite_tup_right = (bgq_weylsite*) ((char*)weylsite_tup_right + kx*sizeof(COMPLEX_PRECISION)); // Some trick: if we are supposed to read k=1, shift the pointer to the right to match k=0, so we avoid some conditional
			bgq_su3_weyl_decl(weyl_tup_right);
			bgq_su3_weyl_loadorprefetch_left(weyl_tup_right, weylsite_tup_right);

			// Merge both weyls to the final result
			bgq_su3_weyl_merge(weyl_tup, weyl_tup_mid, weyl_tup_right);
		#endif
		#if (BGQ_HM_TUP_TLINEINDENT==-1)
		}
		#endif

	#endif
	#if (BGQ_HM_TUP_WEYLREAD==-1)
	}
	#endif


	#if BGQ_HM_TUP_COMPUTE
		bgq_gaugesite *gaugesite_tup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, TUP, t1, t2, !BGQ_HM_TUP_PREFETCH,false);
		bgq_su3_mdecl(gauge_tup);
		bgq_su3_matrix_loadorprefetch(gauge_tup, gaugesite_tup);

		bgq_su3_mvmul(weyl_tup_v0, gauge_tup, weyl_tup_v0);
		bgq_su3_mvmul(weyl_tup_v1, gauge_tup, weyl_tup_v1);

		#ifndef BGQ_HM_NOKAMUL
			bgq_su3_cvmul(weyl_tup_v0, qka0, weyl_tup_v0);
			bgq_su3_cvmul(weyl_tup_v1, qka0, weyl_tup_v1);
		#endif
	#endif


	#if BGQ_HM_TUP_ACCUMULATE
		#if BGQ_HM_TUP_FIRST
			bgq_su3_vmov(result_v0, weyl_tup_v0);
			bgq_su3_vmov(result_v1, weyl_tup_v1);
			bgq_su3_vmov(result_v2, weyl_tup_v0);
			bgq_su3_vmov(result_v3, weyl_tup_v1);
		#else
			bgq_su3_vadd(result_v0, result_v0, weyl_tup_v0);
			bgq_su3_vadd(result_v1, result_v1, weyl_tup_v1);
			bgq_su3_vadd(result_v2, result_v2, weyl_tup_v0);
			bgq_su3_vadd(result_v3, result_v3, weyl_tup_v1);
		#endif
	#endif


}


#include "bgq_loadorprefetch.inc.c"

#undef BGQ_HM_TUP_FIRST
#undef BGQ_HM_TUP_TLINEINDENT
#undef BGQ_HM_TUP_PREFETCH
#undef BGQ_HM_TUP_WEYLREAD
#undef BGQ_HM_TUP_COMPUTE
#undef BGQ_HM_TUP_ACCUMULATE

