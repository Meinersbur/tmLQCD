/*
 * bgq_HoppingMatrix_carry.inc.c
 *
 *  Created on: Aug 20, 2012
 *      Author: meinersbur
 */

#define BGQ_LOADORPREFETCH_PREFETCH BGQ_HM_PREFETCH
#define BGQ_LOADORPREFETCH_LOAD !BGQ_HM_PREFETCH
#include "bgq_loadorprefetch.inc.c"


#ifndef BGQ_HM_DIR_NOFUNC
#include "bgq.h"
#include "bgq_field.h"
#include "bgq_field_double.h"

#define BGQ_PRECISION 64
#include "bgq_precisionselect.inc.c"

#define BGQ_LOADORPREFETCH_LOAD 1
#include "bgq_loadorprefetch.inc.c"

#define BGQ_HM_PREFETCH 0

void HoppingMatrix_carry(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int tv, int x, int y, int z, int t1, int t2) {
	bgq_su3_weyl_decl(weyl_tup);
	bgq_su3_weyl_decl(weyl_tdown);
	bgq_su3_weyl_decl(weyl_zup);
	bgq_su3_weyl_decl(weyl_zdown);

	bgq_su3_weyl_decl(weyl_znext_zdown);
	bgq_su3_weyl_decl(weyl_znext_tdown);
	bgq_su3_weyl_decl(weyl_znext_tup);

	bgq_su3_weyl_decl(weyl_znextnext_zdown);

	bgq_su3_spinor_decl(result);
#else
{
#endif


	//const int tv_left = tv - 1;
	//const int tv_mid = tv;
	//const int tv_right = tv + 1;

	const int z_up = mod(z+1, PHYSICAL_LZ); // wraparound for z==PHYSICAL_LZ


#ifndef BGQ_HM_CARRY_NOSHIFT
	// Shift registers
	//NOTE: There won't be enough registers for all these value, hence the compiler has to spill most of them
	// The hope is, that the total load traffic is still lower because weyls are half the size of spinors (note: this does not apply to the single-precision version)
	// The for-z loop gets 2-unrolled, hence the compiler should allocate space on stack such that no explicit move-instructions are necessary
	bgq_su3_weyl_mov(weyl_tup, weyl_znext_tup);
	bgq_su3_weyl_mov(weyl_tdown, weyl_znext_tdown);
	bgq_su3_weyl_mov(weyl_zdown, weyl_znext_zdown);

	bgq_su3_weyl_mov(weyl_znext_zdown, weyl_znextnext_zdown);
#endif


	bgq_spinorsite *spinorsite_tmid = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z_up, t1, t2, !BGQ_HM_PREFETCH,false);
	bgq_su3_spinor_decl(spinor_tmid);
	bgq_su3_spinor_loadorprefetch(spinor_tmid, spinorsite_tmid);


#ifndef BGQ_HM_CARRY_NO_ZUP
	// Prepare the weyls we need for the next iterations
	// Iteration z (this iteration): z_up
	//bgq_su3_weyl_decl(weyl_zup);
	bgq_su3_vpiadd(weyl_zup_v0, spinor_tmid_v0, spinor_tmid_v2);
	bgq_su3_vpisub(weyl_zup_v1, spinor_tmid_v1, spinor_tmid_v3);
#endif

#ifndef BGQ_HM_CARRY_NO_T
	#if (BGQ_HM_TLINEINDENT==-1)
	if ( ((x+y+z)&1) == isOdd ) {
	#endif
	#if (BGQ_HM_TLINEINDENT==-1) || (BGQ_HM_TLINEINDENT==0)
		// This line is flush, i.e. z_up is ragged

		// Prepare the weyls we need for the next iterations
		// Iteration z+1 (next iteration): t_up, t_down

		//bgq_su3_weyl_decl(weyl_znext_tup);
		#if (BGQ_HM_TUP_WEYLREAD==-1)
		if ( tv != PHYSICAL_LTV-1 ) {
		#endif
		#if (BGQ_HM_TUP_WEYLREAD==-1) || (BGQ_HM_TUP_WEYLREAD==0)
			bgq_spinorsite *spinorsite_tright = BGQ_SPINORSITE_LEFT(spinorfield, !isOdd, tv+1, x, y, z_up, t2+2, t2+4, !BGQ_HM_PREFETCH,false);
			bgq_su3_spinor_decl(spinor_tright);
			bgq_su3_spinor_loadorprefetch_left(spinor_tright, spinorsite_tright);

			bgq_su3_spinor_decl(spinor_tmerge);
			bgq_su3_spinor_merge(spinor_tmerge, spinor_tmid, spinor_tright);

			bgq_su3_vadd(weyl_znext_tup_v0, spinor_tmerge_v0, spinor_tmerge_v2);
			bgq_su3_vadd(weyl_znext_tup_v1, spinor_tmerge_v1, spinor_tmerge_v3);
		#endif
		#if (BGQ_HM_TUP_WEYLREAD==-1)
		} else {
		#endif
		#if (BGQ_HM_TUP_WEYLREAD==-1) || (BGQ_HM_TUP_WEYLREAD==1)
			// Also compute its halfspinor, but we have only a right part, hence wasting a few flops
			bgq_su3_weyl_decl(weyl_znext_tup_mid);
			bgq_su3_vadd(weyl_znext_tup_mid_v0, spinor_tmid_v0, spinor_tmid_v2);
			bgq_su3_vadd(weyl_znext_tup_mid_v1, spinor_tmid_v1, spinor_tmid_v3);

			// Read the other halfspinor from what the neighbor node sent us
			const int xeo = x / PHYSICAL_LP;
			const int xv = xeo / PHYSICAL_LK;
			const int kx = mod(xeo, PHYSICAL_LK);

			bgq_weylsite *weylsite_znext_tup_right = BGQ_WEYLSITE_T(weylxchange_recv[TUP], !isOdd, t2+2, xv, y, z_up, (kx==0) ? x : x-2, (kx==1) ? x : x+2, !BGQ_HM_PREFETCH, false);
			weylsite_znext_tup_right = (bgq_weylsite*) ((char*)weylsite_znext_tup_right + kx*sizeof(COMPLEX_PRECISION)); // Some trick: if we are supposed to read k=1, shift the pointer to the right to match k=0, so we avoid some conditional
			bgq_su3_weyl_decl(weyl_znext_tup_right);
			bgq_su3_weyl_loadorprefetch_left(weyl_znext_tup_right, weylsite_znext_tup_right);

			// Merge both weyls to the final result
			bgq_su3_weyl_merge(weyl_znext_tup, weyl_znext_tup_mid, weyl_znext_tup_right);
		#endif
		#if (BGQ_HM_TUP_WEYLREAD==-1)
		}
		#endif


		//bgq_su3_weyl_decl(weyl_znext_tdown);
		bgq_su3_vsub(weyl_znext_tdown_v0, spinor_tmid_v0, spinor_tmid_v2);
		bgq_su3_vsub(weyl_znext_tdown_v1, spinor_tmid_v1, spinor_tmid_v3);
	#endif
	#if (BGQ_HM_TLINEINDENT==-1)
	} else {
	#endif
	#if (BGQ_HM_TLINEINDENT==-1) || (BGQ_HM_TLINEINDENT==1)
		// This line is ragged, i.e. z_up is flush

		// Prepare the weyls we need for the next iterations
		// Iteration z+1 (next iteration): t_up, t_down
		//bgq_su3_weyl_decl(weyl_znext_tup);
		#if ( BGQ_HM_TUP_WEYLREAD == -1 )
		if ( tv != 0 ) {
		#endif
		#if (BGQ_HM_TUP_WEYLREAD==-1) || (BGQ_HM_TUP_WEYLREAD==0)
			bgq_spinorsite *spinorsite_tleft = BGQ_SPINORSITE_RIGHT(spinorfield, !isOdd, tv-1, x, y, z_up, t1-4, t1-2, !BGQ_HM_PREFETCH,false);
			bgq_su3_spinor_decl(spinor_tleft);
			bgq_su3_spinor_loadorprefetch_right(spinor_tleft, spinorsite_tleft);

			bgq_su3_spinor_decl(spinor_tleftmerge);
			bgq_su3_spinor_merge(spinor_tleftmerge, spinor_tleft, spinor_tmid);

			//bgq_su3_weyl_decl(weyl_znext_tdown);
			bgq_su3_vsub(weyl_znext_tdown_v0, spinor_tleftmerge_v0, spinor_tleftmerge_v2);
			bgq_su3_vsub(weyl_znext_tdown_v1, spinor_tleftmerge_v1, spinor_tleftmerge_v3);
		#endif
		#if (BGQ_HM_TUP_WEYLREAD==-1)
		} else {
		#endif
		#if (BGQ_HM_TUP_WEYLREAD==-1) || (BGQ_HM_TUP_WEYLREAD==1)
			// Compute the halfspinor of the left component, i.e. the result of the right component is to be ignored (it will be the same)
			bgq_su3_weyl_decl(weyl_znext_tdown_mid);
			bgq_su3_vsub(weyl_znext_tdown_mid_v0, spinor_tmid_v0, spinor_tmid_v2);
			bgq_su3_vsub(weyl_znext_tdown_mid_v1, spinor_tmid_v1, spinor_tmid_v3);

			// Read the left component from the weyl receive region
			const int xeo = x / PHYSICAL_LP;
			const int xv = xeo / PHYSICAL_LK;
			const int kx = mod(xeo, PHYSICAL_LK);

			bgq_weylsite *weylsite_znext_tdown_left = BGQ_WEYLSITE_T(weylxchange_recv[TDOWN], !isOdd, t1-2, xv, y, z_up, kx==0 ? x : x-2, kx==1 ? x : x+2, !BGQ_HM_PREFETCH, false);
			weylsite_znext_tdown_left = (bgq_weylsite*)((char*)weylsite_znext_tdown_left + kx*sizeof(COMPLEX_PRECISION)); // Some trick: if we are supposed to read kx=1, shift the pointer to the right to match k=0, so we avoid some conditional
			bgq_su3_weyl_decl(weyl_znext_tdown_left);
			bgq_su3_weyl_loadorprefetch_left(weyl_znext_tdown_left, weylsite_znext_tdown_left);

			// Merge both weyls to get the weyl result
			bgq_su3_weyl_merge(weyl_znext_tdown, weyl_znext_tdown_left, weyl_znext_tdown_mid);
		#endif
		#if (BGQ_HM_TUP_WEYLREAD==-1)
		}
		#endif


		bgq_su3_vadd(weyl_znext_tup_v0, spinor_tmid_v0, spinor_tmid_v2);
		bgq_su3_vadd(weyl_znext_tup_v1, spinor_tmid_v1, spinor_tmid_v3);
	#endif
	#if (BGQ_HM_TLINEINDENT==-1)
	}
	#endif
#endif


#ifndef BGQ_HM_CARRY_NO_ZDOWN
	// Prepare the weyls we need for the next iterations
	// Iteration z+2: z_down
	//bgq_su3_weyl_decl(weyl_znextnext_zdown);
	bgq_su3_vpisub(weyl_znextnext_zdown_v0, spinor_tmid_v0, spinor_tmid_v2);
	bgq_su3_vpiadd(weyl_znextnext_zdown_v1, spinor_tmid_v1, spinor_tmid_v3);
#endif
}


#include "bgq_loadorprefetch.inc.c"

