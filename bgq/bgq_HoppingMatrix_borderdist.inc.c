/*
 * bgq_HoppingMatrix_borderdist.inc.c
 *
 *  Created on: Jul 30, 2012
 *      Author: meinersbur
 */


#ifndef BGQ_HM_BORDERDIST_NOFUNC
#include "bgq.h"
#include "bgq_field.h"
void bgq_HoppingMatrix_borderdist(bool isOdd)
{
#endif
	{
#if !BGQ_HM_BORDER_VERTEX
		SECTION_DECL;
#endif


#if BGQ_HM_BORDER_ZUP
#define BGQ_HM_ZUP_WEYLREAD 1
		const int z = PHYSICAL_LZ-1;
#elif BGQ_HM_BORDER_ZDOWN
#define BGQ_HM_ZDOWN_WEYLREAD 1
		const int z = 0;
#else
		const int z = SECTION_SLICE(PHYSICAL_LZ-2) + 1;
#endif

#if BGQ_HM_BORDER_YUP
#define BGQ_HM_YUP_WEYLREAD 1
		const int y = PHYSICAL_LY-1;
#elif BGQ_HM_BORDER_YDOWN
#define BGQ_HM_YDOWN_WEYLREAD 1
		const int y = 0;
#else
		const int y = SECTION_SLICE(PHYSICAL_LY-2) + 1;
#endif

#if BGQ_HM_BORDER_XUP
#define BGQ_HM_XUP_WEYLREAD 1
		const int x = PHYSICAL_LX-1;
#elif BGQ_HM_BORDER_XDOWN
#define BGQ_HM_XDOWN_WEYLREAD 1
		const int x = 0;
#else
		const int x = SECTION_SLICE(PHYSICAL_LX-2) + 1;
#endif

		assert(xyz == 0);

#include "bgq_HoppingMatrix_tline.inc.c"

	}
#ifndef BGQ_HM_BORDERDIST_NOFUNC
}
#endif


