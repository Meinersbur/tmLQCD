/*
 * bgq_HoppingMatrix_border.inc.c
 *
 *  Created on: Jul 30, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_HM_BORDER_XUP
#define BGQ_HM_BORDER_XUP 0
#endif

#ifndef BGQ_HM_BORDER_XDOWN
#define BGQ_HM_BORDER_XDOWN 0
#endif

#if BGQ_HM_BORDER_XUP || BGQ_HM_BORDER_XDOWN
#define BGQ_HM_BORDER_X 1
#else
#define BGQ_HM_BORDER_X 0
#endif


#ifndef BGQ_HM_BORDER_YUP
#define BGQ_HM_BORDER_YUP 0
#endif

#ifndef BGQ_HM_BORDER_YDOWN
#define BGQ_HM_BORDER_YDOWN 0
#endif

#if BGQ_HM_BORDER_YUP || BGQ_HM_BORDER_YDOWN
#define BGQ_HM_BORDER_Y 1
#else
#define BGQ_HM_BORDER_Y 0
#endif


#ifndef BGQ_HM_BORDER_ZUP
#define BGQ_HM_BORDER_ZUP 0
#endif

#ifndef BGQ_HM_BORDER_ZDOWN
#define BGQ_HM_BORDER_ZDOWN 0
#endif

#if BGQ_HM_BORDER_ZUP || BGQ_HM_BORDER_ZDOWN
#define BGQ_HM_BORDER_Z 1
#else
#define BGQ_HM_BORDER_Z 0
#endif

#if !BGQ_HM_BORDER_Z
#define BGQ_HM_BORDER_Z_INNERMOST 1
#elif !BGQ_HM_BORDER_Y
#define BGQ_HM_BORDER_Y_INNERMOST 1
#elif !BGQ_HM_BORDER_X
#define BGQ_HM_BORDER_X_INNERMOST 1
#else
#define BGQ_HM_BORDER_VERTEX 1
#endif




#ifndef BGQ_HM_BORDER_NOFUNC
#include "bgq.h"
#include "bgq_field.h"
void bgq_HoppingMatrix_border(bool isOdd)
{
#define BGQ_HM_BORDERDIST_NOFUNC 1
#define BGQ_HM_TLINE_NOFUNC 1
#define BGQ_HM_SITE_NOFUNC 1
#define BGQ_HM_DIR_NOFUNC 1
#endif
	{

	int volume = 1;
#if !BGQ_HM_BORDER_X
	volume *= (PHYSICAL_LX-2);
#endif
#if !BGQ_HM_BORDER_Y
	volume *= (PHYSICAL_LY-2);
#endif
#if !BGQ_HM_BORDER_Z
	volume *= (PHYSICAL_LZ-2);
#endif



#if BGQ_HM_BORDER_VERTEX


#if BGQ_HM_BORDER_XDOWN && BGQ_HM_BORDER_YDOWN && BGQ_HM_BORDER_XDOWN
// 3 downs
#define BGQ_HM_VERTEX_RAGGEDLINE 1
#elif !BGQ_HM_BORDER_XDOWN && !BGQ_HM_BORDER_YDOWN && !BGQ_HM_BORDER_XDOWN
// 0 downs
#define BGQ_HM_VERTEX_FLUSHLINE 1
#elif (BGQ_HM_BORDER_YDOWN && BGQ_HM_BORDER_ZDOWN) || (BGQ_HM_BORDER_XDOWN && BGQ_HM_BORDER_ZDOWN) || (BGQ_HM_BORDER_XDOWN && BGQ_HM_BORDER_YDOWN)
// 2 downs
#define BGQ_HM_VERTEX_FLUSHLINE 1
#else
// 1 down
#define BGQ_HM_VERTEX_RAGGEDLINE 1
#endif

#if BGQ_HM_VERTEX_FLUSHLINE
#define BGQ_HM_TLINE_FLUSHLINE 1
#endif
#if BGQ_HM_VERTEX_RAGGEDLINE
#define BGQ_HM_TLINE_RAGGEDLINE 1
#endif

	{
	int xyz = 0;
#include "bgq_HoppingMatrix_borderdist.inc.c"
	}


#else


#pragma omp parallel for schedule(static,1)
	for (int xyz = isOdd; xyz < volume; xyz += 2) {
#define BGQ_HM_TLINE_RAGGEDLINE 1
#include "bgq_HoppingMatrix_borderdist.inc.c"
	}

#pragma omp parallel for schedule(static,1)
	for (int xyz = !isOdd; xyz < volume; xyz += 2) {
#define BGQ_HM_TLINE_FLUSHLINE 1
#include "bgq_HoppingMatrix_borderdist.inc.c"
	}


#endif


	}
#ifndef BGQ_HM_BORDER_NOFUNC
}
#endif


#undef BGQ_HM_BORDER_XUP
#undef BGQ_HM_BORDER_XDOWN
#undef BGQ_HM_BORDER_X
#undef BGQ_HM_BORDER_X_INNERMOST

#undef BGQ_HM_BORDER_YUP
#undef BGQ_HM_BORDER_YDOWN
#undef BGQ_HM_BORDER_Y
#undef BGQ_HM_BORDER_Y_INNERMOST

#undef BGQ_HM_BORDER_ZUP
#undef BGQ_HM_BORDER_ZDOWN
#undef BGQ_HM_BORDER_Z
#undef BGQ_HM_BORDER_Z_INNERMOST

#undef BGQ_HM_BORDER_VERTEX
#undef BGQ_HM_VERTEX_RAGGEDLINE
#undef BGQ_HM_VERTEX_FLUSHLINE
