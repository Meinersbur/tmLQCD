/*
 * bgq_HoppingMatrix_tline.inc.c
 *
 *  Created on: Jul 27, 2012
 *      Author: meinersbur
 */

#if BGQ_HM_TLINE_FLUSHLINE
#define BGQ_HM_TLINE_RAGGEDLINE 0
#elif BGQ_HM_TLINE_RAGGEDLINE
#define BGQ_HM_TLINE_FLUSHLINE 0
#else
#error Need to define flush or ragged line
#endif


#ifndef BGQ_HM_TLINE_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void bgq_HoppingMatrix_tline(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int x, int y, int z) {
	bgq_su3_spinor_decl(result);
#define BGQ_HM_SITE_NOFUNC 1
#define BGQ_HM_DIR_NOFUNC 1
#endif

	{

#if BGQ_HM_TLINE_FLUSHLINE
		// Flush line (distance 0 from local left border)
#define BGQ_HM_TUP_FLUSHLINE 1
#define BGQ_HM_TDOWN_FLUSHLINE 1
#elif BGQ_HM_TLINE_RAGGEDLINE
		// Ragged line (distance 1 from local left border)
#define BGQ_HM_TUP_RAGGEDLINE 1
#define BGQ_HM_TDOWN_RAGGEDLINE 1
#else
#error need to define flush- or ragged line
#endif

		bgq_su3_spinor_decl(spinor_tcarry);

		// Prologue
		{
			int tv = 0;
#define BGQ_HM_TDOWN_LEFTWRAPAROUND 1
#define BGQ_HM_TUP_WRITECARRYSPINOR 1
#include "bgq_HoppingMatrix_site.inc.c"
		}

		for (int tv = 1; tv < PHYSICAL_LTV-1; tv+=1) {
#define BGQ_HM_TDOWN_READCARRYSPINOR 1
#define BGQ_HM_TUP_WRITECARRYSPINOR 1
			#include "bgq_HoppingMatrix_site.inc.c"
		}

		// Epilogue
		{
			int tv = PHYSICAL_LTV-1;
#define BGQ_HM_TDOWN_READCARRYSPINOR 1
#define BGQ_HM_TUP_RIGHTWRAPAROUND 1
#include "bgq_HoppingMatrix_site.inc.c"
		}

	}
#ifndef BGQ_HM_TLINE_NOFUNC
}
#undef BGQ_HM_SITE_NOFUNC
#endif

#undef BGQ_HM_TLINE_FLUSHLINE
#undef BGQ_HM_TLINE_RAGGEDLINE
