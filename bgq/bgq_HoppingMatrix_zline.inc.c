/*
 * bgq_HoppingMatrix_tline.inc.c
 *
 *  Created on: Jul 27, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_HM_ZLINE_STARTINDENT
#define BGQ_HM_ZLINE_STARTINDENT -1
#endif

#ifndef BGQ_HM_ZLINE_TUP_WEYLREAD
#define BGQ_HM_ZLINE_TUP_WEYLREAD 0
#endif

#ifndef BGQ_HM_ZLINE_TDOWN_WEYLREAD
#define BGQ_HM_ZLINE_TDOWN_WEYLREAD 0
#endif

#ifndef BGQ_HM_ZLINE_XUP_WEYLREAD
#define BGQ_HM_ZLINE_XUP_WEYLREAD 0
#endif

#ifndef BGQ_HM_ZLINE_XDOWN_WEYLREAD
#define BGQ_HM_ZLINE_XDOWN_WEYLREAD 0
#endif

#ifndef BGQ_HM_ZLINE_YUP_WEYLREAD
#define BGQ_HM_ZLINE_YUP_WEYLREAD 0
#endif

#ifndef BGQ_HM_ZLINE_YDOWN_WEYLREAD
#define BGQ_HM_ZLINE_YDOWN_WEYLREAD 0
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

#define BGQ_HM_SITE_NOFUNC 1
#define BGQ_HM_DIR_NOFUNC 1
void bgq_HoppingMatrix_zline(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int t, int x, int y) {
	int t1 = 0;
	int t2 = 2;
#else
{
#endif

	int z = 0;

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
		if (z >= PHYSICAL_LX)
			break;
		{
			#define BGQ_HM_TUP_TLINEINDENT 0
			#define BGQ_HM_TDOWN_TLINEINDENT 0
			#define BGQ_HM_TUP_WEYLREAD BGQ_HM_ZLINE_TUP_WEYLREAD
			#define BGQ_HM_TDOWN_WEYLREAD BGQ_HM_ZLINE_TDOWN_WEYLREAD
			#define BGQ_HM_XUP_WEYLREAD BGQ_HM_ZLINE_XUP_WEYLREAD
			#define BGQ_HM_XDOWN_WEYLREAD BGQ_HM_ZLINE_XDOWN_WEYLREAD
			#define BGQ_HM_YUP_WEYLREAD BGQ_HM_ZLINE_YUP_WEYLREAD
			#define BGQ_HM_YDOWN_WEYLREAD BGQ_HM_ZLINE_YDOWN_WEYLREAD
			#include "bgq_HoppingMatrix_site.inc.c"
		}
		t1 += 1;
		t2 += 1;
		z += 1;

		// There is always a ragged line following
		STARTRAGGED:
		if (z >= PHYSICAL_LX)
			break;
		{
			#define BGQ_HM_TUP_TLINEINDENT 1
			#define BGQ_HM_TDOWN_TLINEINDENT 1
			#define BGQ_HM_TUP_WEYLREAD BGQ_HM_ZLINE_TUP_WEYLREAD
			#define BGQ_HM_TDOWN_WEYLREAD BGQ_HM_ZLINE_TDOWN_WEYLREAD
			#define BGQ_HM_XUP_WEYLREAD BGQ_HM_ZLINE_XUP_WEYLREAD
			#define BGQ_HM_XDOWN_WEYLREAD BGQ_HM_ZLINE_XDOWN_WEYLREAD
			#define BGQ_HM_YUP_WEYLREAD BGQ_HM_ZLINE_YUP_WEYLREAD
			#define BGQ_HM_YDOWN_WEYLREAD BGQ_HM_ZLINE_YDOWN_WEYLREAD
			#include "bgq_HoppingMatrix_site.inc.c"
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


#undef BGQ_HM_ZLINE_STARTINDENT
#undef BGQ_HM_ZLINE_TUP_WEYLREAD
#undef BGQ_HM_ZLINE_TDOWN_WEYLREAD
#undef BGQ_HM_ZLINE_XUP_WEYLREAD
#undef BGQ_HM_ZLINE_XDOWN_WEYLREAD
#undef BGQ_HM_ZLINE_YUP_WEYLREAD
#undef BGQ_HM_ZLINE_YDOWN_WEYLREAD

#undef BGQ_HM_ZLINE_ID
#undef STARTFLUSH
#undef STARTRAGGED

