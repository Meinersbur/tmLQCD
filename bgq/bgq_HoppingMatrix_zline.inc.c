/*
 * bgq_HoppingMatrix_tline.inc.c
 *
 *  Created on: Jul 27, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_HM_ZLINE_STARTINDENT
#define BGQ_HM_ZLINE_STARTINDENT -1
#endif

#ifndef BGQ_HM_ZLINE_CARRY
#define BGQ_HM_ZLINE_CARRY 1
#endif

#ifndef BGQ_HM_ZLINE_ID
#error Need some unique string to identify goto labels
#define BGQ_HM_ZLINE_ID SOMESTRING
#endif
#define STARTFLUSH NAME2(STARTFLUSH,BGQ_HM_ZLINE_ID)
#define STARTRAGGED NAME2(STARTRAGGED,BGQ_HM_ZLINE_ID)

#ifndef BGQ_HM_ZLINE_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

#define BGQ_PRECISION 64
#include "bgq_precisionselect.inc.c"

#define BGQ_LOADORPREFETCH_LOAD 1
#include "bgq_loadorprefetch.inc.c"

#define BGQ_HM_SITE_NOFUNC 1
#define BGQ_HM_DIR_NOFUNC 1
void bgq_HoppingMatrix_zline(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int tv, int x, int y, int t1, int t2) {
#else
{
#endif

	int z = 0;

	#if BGQ_PREFETCH_STREAM
	if (!noprefetchstream) {
		// Lots of streams (14), still without the border weyl ones

		bgq_spinorsite *spinorsite_t_left = BGQ_SPINORSITE(spinorfield, !isOdd, tv-1, x, y, z, t1, t2, false,false);
		bgq_prefetch_forward(spinorsite_t_left);
		bgq_spinorsite *spinorsite_t_mid = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z, t1, t2, false,false);
		bgq_prefetch_forward(spinorsite_t_mid);
		bgq_spinorsite *spinorsite_t_right = BGQ_SPINORSITE(spinorfield, !isOdd, tv-1, x, y, z, t1, t2, false,false);
		bgq_prefetch_forward(spinorsite_t_right);

		bgq_gaugesite *gaugesite_tup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, TUP, t1, t2, false,false);
		bgq_prefetch_forward(gaugesite_tup);
		bgq_gaugesite *gaugesite_tdown_flush = BGQ_GAUGESITE(gaugefield, !isOdd, tv, x, y, z, TUP_SHIFT, t1-1, t2-1, false,false);
		bgq_prefetch_forward(gaugesite_tdown_flush);
		bgq_gaugesite *gaugesite_tdown_ragged = BGQ_GAUGESITE(gaugefield, !isOdd, tv, x, y, z, TUP, t1-1, t2-1, false,false);
		bgq_prefetch_forward(gaugesite_tdown_ragged);

		bgq_spinorsite *spinorsite_xup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x+1, y, z, t1, t2, false,false);
		bgq_prefetch_forward(spinorsite_xup);
		bgq_spinorsite *spinorsite_xdown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x-1, y, z, t1, t2, false,false);
		bgq_prefetch_forward(spinorsite_xdown);

		bgq_gaugesite *gaugesite_xup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, XUP, t1, t2, false,false);
		bgq_prefetch_forward(gaugesite_xup);
		bgq_gaugesite *gaugesite_xdown = BGQ_GAUGESITE(gaugefield, isOdd, tv, x-1, y, z, XUP, t1, t2, false,false);
		bgq_prefetch_forward(gaugesite_xdown);

		bgq_spinorsite *spinorsite_yup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y+1, z, t1, t2, false,false);
		bgq_prefetch_forward(spinorsite_yup);
		bgq_spinorsite *spinorsite_ydown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y-1, z, t1, t2, false,false);
		bgq_prefetch_forward(spinorsite_ydown);

		bgq_gaugesite *gaugesite_yup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, YUP, t1, t2, false,false);
		bgq_prefetch_forward(gaugesite_yup);
		bgq_gaugesite *gaugesite_ydown = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y-1, z, YUP, t1, t2, false,false);
		bgq_prefetch_forward(gaugesite_ydown);
	}
	#endif


#if BGQ_HM_CARRY
	bgq_su3_spinor_decl(weyl_zup);
	bgq_su3_spinor_decl(weyl_znext_tup);
	bgq_su3_spinor_decl(weyl_znext_tdown);
	bgq_su3_spinor_decl(weyl_znextnext_zdown);

	bgq_su3_spinor_decl(weyl_tup);
	bgq_su3_spinor_decl(weyl_tdown);
	bgq_su3_spinor_decl(weyl_znext_zdown);

	bgq_su3_spinor_decl(weyl_zdown);



	// Prologue for loading the carry values
	{
		// z==-1, carry value for z=0,z_down
		const int z = PHYSICAL_LZ-2; // such that z+1==-1==PHYSICAL_LZ-1 (z_up) is the site to load
		#define BGQ_HM_TLINEINDENT BGQ_HM_ZLINE_STARTINDENT
		#define BGQ_HM_CARRY_NOSHIFT 1
		#define BGQ_HM_CARRY_NO_ZUP 1
		#define BGQ_HM_CARRY_NO_T 1
		#define BGQ_HM_PREFETCH 0
		#include "bgq_HoppingMatrix_carry.inc.c"
		#undef BGQ_HM_TLINEINDENT
		#undef BGQ_HM_CARRY_NOSHIFT
		#undef BGQ_HM_CARRY_NO_ZUP
		#undef BGQ_HM_CARRY_NO_T
		#undef BGQ_HM_PREFETCH
	}


	const int t1_shadowed = t1;
	const int t2_shadowed = t2;
	{
		// shift
		bgq_su3_weyl_mov(weyl_znext_zdown, weyl_znextnext_zdown);

		const int t1 = t1_shadowed + 1-(t1_shadowed&1)*2;
		const int t2 = t2_shadowed + 1-(t2_shadowed&1)*2;
		const int z = -1; // such that z+1==0 (z_up) is the site to load

		// z==0, carry value for z=0,t_up,t_down z=1,z_down
		#if (BGQ_HM_ZLINE_STARTINDENT==-1)
		#define BGQ_HM_TLINEINDENT -1
		#else
		#define BGQ_HM_TLINEINDENT !BGQ_HM_ZLINE_STARTINDENT
		#endif
		#define BGQ_HM_CARRY_NOSHIFT 1
		#define BGQ_HM_CARRY_NO_ZUP 1
		#define BGQ_HM_PREFETCH 0
		#include "bgq_HoppingMatrix_carry.inc.c"
		#undef BGQ_HM_TLINEINDENT
		#undef BGQ_HM_CARRY_NOSHIFT
		#undef BGQ_HM_CARRY_NO_ZUP
		#undef BGQ_HM_PREFETCH
	}

#endif


	#if (BGQ_HM_ZLINE_STARTINDENT==-1)
	if ( ((x+y)&1) == isOdd ) {
	#endif
	#if (BGQ_HM_ZLINE_STARTINDENT==-1) || (BGQ_HM_ZLINE_STARTINDENT==0)
		goto STARTFLUSH;
	#endif
	#if (BGQ_HM_ZLINE_STARTINDENT==-1)
	} else {
	#endif
	#if (BGQ_HM_ZLINE_STARTINDENT==-1) || (BGQ_HM_ZLINE_STARTINDENT==1)
		goto STARTRAGGED;
	#endif
	#if (BGQ_HM_ZLINE_STARTINDENT==-1)
	}
	#endif

	// A 2-unrolled for-loop
	// How does the compiler cope with the jump into the middle? Does it still do strength reduction of variables dependent on z?
	// If there is a single jump target  (BGQ_HM_ZLINE_STARTINDENT!=-1), does it recognize the loop?
	// NOTE: This jump-into-the middle saves us dublication of the loop
	while (true) {
		// Begin with flush line
		STARTFLUSH:
		if (z >= PHYSICAL_LZ)
			break;
		{
			#define BGQ_HM_TLINEINDENT 0
			#include "bgq_HoppingMatrix_site.inc.c"
			#undef BGQ_HM_TLINEINDENT
		}
		t1 += 1;
		t2 += 1;
		z += 1;

		// There is always a ragged line following
		STARTRAGGED:
		if (z >= PHYSICAL_LZ)
			break;
		{
			#define BGQ_HM_TLINEINDENT 1
			#include "bgq_HoppingMatrix_site.inc.c"
			#undef BGQ_HM_TLINEINDENT
		}
		t1 -= 1;
		t2 -= 1;
		z += 1;
	}
}


#ifndef BGQ_HM_ZLINE_NOFUNC
#undef BGQ_HM_SITE_NOFUNC
#undef BGQ_HM_DIR_NOFUNC
#endif


#undef BGQ_HM_ZLINE_CARRY
#undef BGQ_HM_ZLINE_STARTINDENT
#undef BGQ_HM_ZLINE_ID
#undef STARTFLUSH
#undef STARTRAGGED

