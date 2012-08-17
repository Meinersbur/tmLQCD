
#ifndef BGQ_HM_TDOWN_TLINEINDENT
#define BGQ_HM_TDOWN_TLINEINDENT -1
#endif

#ifndef BGQ_HM_TDOWN_PREFETCH
#define BGQ_HM_TDOWN_PREFETCH 0
#endif

#ifndef BGQ_HM_TDOWN_WEYLREAD
#define BGQ_HM_TDOWN_WEYLREAD 0
#endif

#ifndef BGQ_HM_TDOWN_COMPUTE
#define BGQ_HM_TDOWN_COMPUTE 0
#endif

#ifndef BGQ_HM_TDOWN_ACCUMULATE
#define BGQ_HM_TDOWN_ACCUMULATE 0
#endif


#define BGQ_LOADORPREFETCH_PREFETCH BGQ_HM_TDOWN_PREFETCH
#define BGQ_LOADORPREFETCH_LOAD !BGQ_HM_TDOWN_PREFETCH
#include "bgq_loadorprefetch.inc.c"


#ifndef BGQ_HM_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void bgq_HoppingMatrix_tdown(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int tv, int x, int y, int z, int t1, int t2) {
	bgq_su3_spinor_decl(result);
	bgq_vector4double_decl(qka0);
#else
{
#endif


	bgq_su3_weyl_decl(weyl_tdown);
	#if (BGQ_HM_TDOWN_WEYLREAD==-1)
	if ( tv != 0 ) {
	#endif
	#if (BGQ_HM_TDOWN_WEYLREAD==-1) || (BGQ_HM_TDOWN_WEYLREAD==0)
		// Load the input spinor
		bgq_su3_spinor_decl(spinor_tdown);
		#if (BGQ_HM_TDOWN_TLINEINDENT==-1)
		if ( ((x+y+z)&1) == isOdd ) {
		#endif
		#if (BGQ_HM_TDOWN_TLINEINDENT==-1) || (BGQ_HM_TDOWN_TLINEINDENT==0)
			bgq_spinorsite *spinorsite_tdown_left = BGQ_SPINORSITE_RIGHT(spinorfield, !isOdd, tv-1, x, y, z, t1-3,t1-1, !BGQ_HM_TDOWN_PREFETCH,false);
			bgq_su3_spinor_decl(spinor_tdown_left);
			bgq_su3_spinor_loadorprefetch_right(spinor_tdown_left, spinorsite_tdown_left);

			bgq_spinorsite *spinorsite_tdown_mid = BGQ_SPINORSITE_LEFT(spinorfield, !isOdd, tv, x, y, z, t2-1,t2+1, !BGQ_HM_TDOWN_PREFETCH,false);
			bgq_su3_spinor_decl(spinor_tdown_mid);
			bgq_su3_spinor_loadorprefetch_left(spinor_tdown_mid, spinorsite_tdown_mid);

			#if !BGQ_HM_TDOWN_PREFETCH
				bgq_su3_spinor_merge(spinor_tdown, spinor_tdown_left, spinor_tdown_mid);
			#endif
		#endif
		#if (BGQ_HM_TDOWN_TLINEINDENT==-1)
		} else {
		#endif
		#if (BGQ_HM_TDOWN_TLINEINDENT==-1) || (BGQ_HM_TDOWN_TLINEINDENT==1)
			bgq_spinorsite *spinorsite_tdown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z, t1-1, t2-1, !BGQ_HM_TDOWN_PREFETCH,false);
			bgq_su3_spinor_loadorprefetch(spinor_tdown, spinorsite_tdown);
		#endif
		#if (BGQ_HM_TDOWN_TLINEINDENT==-1)
		}
		#endif

		#if !BGQ_HM_TDOWN_PREFETCH
			// Compute its halfspinor
			bgq_su3_vsub(weyl_tdown_v0, spinor_tdown_v0, spinor_tdown_v2);
			bgq_su3_vsub(weyl_tdown_v1, spinor_tdown_v1, spinor_tdown_v3);
		#endif
	#endif
	#if (BGQ_HM_TDOWN_WEYLREAD==-1)
	} else {
	#endif
	#if (BGQ_HM_TDOWN_WEYLREAD==-1) || (BGQ_HM_TDOWN_WEYLREAD==1)
		#if (BGQ_HM_TDOWN_TLINEINDENT==-1)
		if ( ((tv+x+y+z)&1) == isOdd ) {
		#endif
		#if (BGQ_HM_TDOWN_TLINEINDENT==-1) || (BGQ_HM_TDOWN_TLINEINDENT==0)
			assert( t1 == 0 );

			// Read the left component from the weyl receive region
			const int xeo = x / PHYSICAL_LP;
			const int xv = xeo / PHYSICAL_LK;
			const int kx = mod(xeo, PHYSICAL_LK);

			bgq_weylsite *weylsite_tdown_left = BGQ_WEYLSITE_T(weylxchange_recv[TDOWN], !isOdd, t1-1, xv, y, z, kx==0 ? x : x-2, kx==1 ? x : x+2, !BGQ_HM_TDOWN_PREFETCH, false);
			#if !BGQ_HM_TDOWN_PREFETCH
				weylsite_tdown_left = (bgq_weylsite_double*)((char*)weylsite_tdown_left + kx*sizeof(_Complex double)); // Some trick: if we are supposed to read kx=1, shift the pointer to the right to match k=0, so we avoid some conditional
			#endif
			bgq_su3_weyl_decl(weyl_tdown_left);
			bgq_su3_weyl_loadorprefetch_left(weyl_tdown_left, weylsite_tdown_left);

			// Read the second component as normal
			bgq_spinorsite *spinorsite_tdown_mid = BGQ_SPINORSITE_LEFT(spinorfield, !isOdd, tv, x, y, z, t2-1,t2+1, !BGQ_HM_TDOWN_PREFETCH, false);
			bgq_su3_spinor_decl(spinor_tdown_mid);
			bgq_su3_spinor_loadorprefetch_left(spinor_tdown_mid, spinorsite_tdown_mid);

			#if !BGQ_HM_TDOWN_PREFETCH
				// Compute the halfspinor og the left component, i.e. the result of the right component is to be ignored (it will be the same)
				bgq_su3_weyl_decl(weyl_tdown_mid);
				bgq_su3_vsub(weyl_tdown_mid_v0, spinor_tdown_mid_v0, spinor_tdown_mid_v2);
				bgq_su3_vsub(weyl_tdown_mid_v1, spinor_tdown_mid_v1, spinor_tdown_mid_v3);

				// Merge both weyls to get the weyl result
				bgq_su3_weyl_merge(weyl_tdown, weyl_tdown_left, spinor_tdown_mid);
			#endif
		#endif
		#if (BGQ_HM_TDOWN_TLINEINDENT==-1)
		} else {
		#endif
		#if (BGQ_HM_TDOWN_TLINEINDENT==-1) || (BGQ_HM_TDOWN_TLINEINDENT==1)
			assert( t1 == 1 );
			// Read spinor as normal
			bgq_spinorsite *spinorsite_tdown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z, t1-1, t2-1, !BGQ_HM_TDOWN_PREFETCH,false);
			bgq_su3_spinor_decl(spinor_tdown);
			bgq_su3_spinor_loadorprefetch(spinor_tdown, spinorsite_tdown);


			#if !BGQ_HM_TDOWN_PREFETCH
				// Compute its halfspinor
				bgq_su3_vsub(weyl_tdown_v0, spinor_tdown_v0, spinor_tdown_v2);
				bgq_su3_vsub(weyl_tdown_v1, spinor_tdown_v1, spinor_tdown_v3);
			#endif
		#endif
		#if (BGQ_HM_TDOWN_TLINEINDENT==-1)
		}
		#endif
	#endif
	#if (BGQ_HM_TDOWN_WEYLREAD==-1)
	}
	#endif


	#if BGQ_HM_TDOWN_COMPUTE
		bgq_gaugesite *gaugesite_tdown;
		#if (BGQ_HM_TDOWN_TLINEINDENT==-1)
		if ( ((tv+x+y+z)&1) == isOdd ) {
		#endif
		#if (BGQ_HM_TDOWN_TLINEINDENT==-1) || (BGQ_HM_TDOWN_TLINEINDENT==0)
			gaugesite_tdown = BGQ_GAUGESITE(gaugefield, !isOdd, tv, x, y, z, TUP_SHIFT, t1-1, t2-1, !BGQ_HM_TDOWN_PREFETCH,false);
		#endif
		#if (BGQ_HM_TDOWN_TLINEINDENT==-1)
		} else {
		#endif
		#if (BGQ_HM_TDOWN_TLINEINDENT==-1) || (BGQ_HM_TDOWN_TLINEINDENT==1)
			gaugesite_tdown = BGQ_GAUGESITE(gaugefield, !isOdd, tv, x, y, z, TUP, t1-1, t2-1, !BGQ_HM_TDOWN_PREFETCH,false);
		#endif
		#if (BGQ_HM_TDOWN_TLINEINDENT==-1)
		}
		#endif

		bgq_su3_mdecl(gauge_tdown);
		bgq_su3_matrix_loadorprefetch(gauge_tdown, gaugesite_tdown);

		bgq_su3_mvinvmul(weyl_tdown_v0, gauge_tdown, weyl_tdown_v0);
		bgq_su3_mvinvmul(weyl_tdown_v1, gauge_tdown, weyl_tdown_v1);

		#ifndef BGQ_HM_NOKAMUL
			bgq_su3_cvmul(weyl_tdown_v0, qka0, weyl_tdown_v0);
			bgq_su3_cvmul(weyl_tdown_v1, qka0, weyl_tdown_v1);
		#endif
	#endif


	#if BGQ_HM_TDOWN_ACCUMULATE
		bgq_su3_vadd(result_v0, result_v0, weyl_tdown_v0);
		bgq_su3_vadd(result_v1, result_v1, weyl_tdown_v1);
		bgq_su3_vsub(result_v2, result_v2, weyl_tdown_v0);
		bgq_su3_vsub(result_v3, result_v3, weyl_tdown_v1);
	#endif


}


#include "bgq_loadorprefetch.inc.c"

#undef BGQ_HM_TDOWN_TLINEINDENT
#undef BGQ_HM_TDOWN_PREFETCH
#undef BGQ_HM_TDOWN_WEYLREAD
#undef BGQ_HM_TDOWN_COMPUTE
#undef BGQ_HM_TDOWN_ACCUMULATE

