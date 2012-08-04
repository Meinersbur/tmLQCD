/*
 * bgq_HoppingMatrix_tline.inc.c
 *
 *  Created on: Jul 27, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_HM_ZLINE_ZLINEINDENT
#error Must define the line indention (0, 1 or -1 for runtime-conditional)
#define BGQ_HM_ZLINE_ZLINEINDENT -1
#endif


#ifndef BGQ_HM_ZLINE_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void bgq_HoppingMatrix_zline(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int t, int x, int y) {
	bgq_su3_spinor_decl(result);
#define BGQ_HM_SITE_NOFUNC 1
#define BGQ_HM_DIR_NOFUNC 1
#endif
	{


		int z1;
		int z2;
#if (BGQ_HM_ZLINE_ZLINEINDENT==-1)
		if (((t + x + y)&1) == isOdd) {
#endif
#if (BGQ_HM_ZLINE_ZLINEINDENT==-1) || (BGQ_HM_ZLINE_ZLINEINDENT==0)
			assert(((t+x+y)&1)==isOdd);
			// zline indention == 0
			// iterate from zbottom to ztop
			z1 = 0;
			z2 = 2;
#endif
#if (BGQ_HM_ZLINE_ZLINEINDENT==-1)
			} else {
#endif
#if (BGQ_HM_ZLINE_ZLINEINDENT==-1) || (BGQ_HM_ZLINE_ZLINEINDENT==1)
				assert(((t+x+y)&1)==!isOdd);
				// zline indention == 1
				z1 = 1;
				z2 = 3;
#endif
#if (BGQ_HM_ZLINE_ZLINEINDENT==-1)
			}
#endif


		for (int zv = 0; zv < PHYSICAL_LZV-1; zv+=1) {
			#define BGQ_HM_ZUP_ZLINEINDENT BGQ_HM_ZLINE_ZLINEINDENT
			#define BGQ_HM_ZDOWN_ZLINEINDENT BGQ_HM_ZLINE_ZLINEINDENT
			#include "bgq_HoppingMatrix_site.inc.c"

			z1 += PHYSICAL_LP * PHYSICAL_LK;
			z2 += PHYSICAL_LP * PHYSICAL_LK;
		}



	}
#ifndef BGQ_HM_ZLINE_NOFUNC
}
#undef BGQ_HM_SITE_NOFUNC
#endif


#undef BGQ_HM_ZLINE_ZLINEINDENT
