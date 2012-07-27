/*
 * bgq_HoppingMatrix_tline.inc.c
 *
 *  Created on: Jul 27, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_HM_TLINE_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void bgq_HoppingMatrix_tline(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefieldeo_double gaugefield, bool isOdd, int x, int y, int z) {
	bgq_su3_spinor_decl(result);
#endif

#define BGQ_HM_SITE_NOFUNC


	if ((x+y+z)%2==isOdd) {
		// Flush line (distance 0 from local left border)
		bgq_su3_spinor_decl(spinor_tdown_flush);
		bgq_spinorsite_double *spinorsite_flush = BGQ_SPINORSITE(spinorfield,isOdd,x,y,z,0);
		bgq_su3_spinor_double_load_left(spinor_tdown_flush, spinorsite_flush);

		bgq_su3_spinor_decl(spinorsite_wraparound);
		bgq_spinorsite_double *spinorsite_wraparound = BGQ_SPINORSITE(spinorfield,isOdd,x,y,z,PHYSICAL_LTV-1);
		bgq_su3_spinor_double_load_right(spinorsite_wraparound, spinorsite_wraparound);
		bgq_su3_spinor_merge(spinor_tdown, spinorsite_wraparound, spinor_tdown_flush);


		for (int tv = 0; tv < PHYSICAL_LTV; tv+=1) {
			#include "bgq_HoppingMatrix_site.inc.c"
		}

	} else {
		// Ragged line (distance 1 from left border)

		bgq_spinorsite_double *spinorsite = BGQ_SPINORSITE(spinorfield,isOdd,x,y,z,0);
				bgq_su3_spinor_double_load(spinor_tdown, spinorsite);
				tlinesize = PHYSICAL_LTV-1;

		for (int tv = 0; tv < PHYSICAL_LTV-1; tv+=1) {
			#include "bgq_HoppingMatrix_site.inc.c"
		}

		// Epilogue
		 {
				int tv = PHYSICAL_LTV-1;

				bgq_su3_spinor_decl(spinor_tup);

				bgq_su3_spinor_decl(spinor_tup_flush);
				bgq_spinorsite_double *spinorsite_tup_flush = BGQ_SPINORSITE(spinorfield,isOdd,x,y,z,tv);
				bgq_su3_spinor_double_load_right(spinor_tup_flush, spinorsite_tup_flush);

				bgq_su3_spinor_decl(spinor_tup_wraparound);
				bgq_spinorsite_double *spinorsite_tup_wraparound = BGQ_SPINORSITE(spinorfield,isOdd,x,y,z,0);
				bgq_su3_spinor_double_load_left(spinor_tup_wraparound, spinorsite_tup_wraparound);
				bgq_su3_spinor_merge(spinor_tdown, spinor_tdown_flush, spinor_tup_wraparound);

		#include "bgq_HoppingMatrix_site.inc.c"
			}
	}

#undef BGQ_HM_NOFUNC

#ifndef BGQ_HM_TLINE_NOFUNC
}
#endif


