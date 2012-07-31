/*
 * bgq_HoppingMatrix_tline.inc.c
 *
 *  Created on: Jul 27, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_HM_TLINE_TLINEINDENT
#error Must define the line indention (0, 1 or -1 for runtime-conditional)
#define BGQ_HM_TLINE_TLINEINDENT -1
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

		// Load constants
		bgq_vector4double_decl(qka0);
		bgq_cconst(qka0,creal(ka0),cimag(ka0));
		bgq_vector4double_decl(qka1);
		bgq_cconst(qka1,creal(ka1),cimag(ka1));
		bgq_vector4double_decl(qka2);
		bgq_cconst(qka2,creal(ka2),cimag(ka2));
		bgq_vector4double_decl(qka3);
		bgq_cconst(qka3,creal(ka3),cimag(ka3));


		bgq_su3_spinor_decl(spinor_tcarry);

		// Prologue
		{
			const int tv = 0;
#define BGQ_HM_TUP_TLINEINDENT BGQ_HM_TLINE_TLINEINDENT
#define BGQ_HM_TDOWN_TLINEINDENT BGQ_HM_TLINE_TLINEINDENT
#define BGQ_HM_TDOWN_LEFTWRAPAROUND 1
#define BGQ_HM_TUP_WRITECARRYSPINOR 1
#include "bgq_HoppingMatrix_site.inc.c"
		}

		for (int tv = 1; tv < PHYSICAL_LTV-1; tv+=1) {
#define BGQ_HM_TUP_TLINEINDENT BGQ_HM_TLINE_TLINEINDENT
#define BGQ_HM_TDOWN_TLINEINDENT BGQ_HM_TLINE_TLINEINDENT
#define BGQ_HM_TDOWN_READCARRYSPINOR 1
#define BGQ_HM_TUP_WRITECARRYSPINOR 1
			#include "bgq_HoppingMatrix_site.inc.c"
		}

		// Epilogue
		{
			const int tv = PHYSICAL_LTV-1;
#define BGQ_HM_TUP_TLINEINDENT BGQ_HM_TLINE_TLINEINDENT
#define BGQ_HM_TDOWN_TLINEINDENT BGQ_HM_TLINE_TLINEINDENT
#define BGQ_HM_TDOWN_READCARRYSPINOR 1
#define BGQ_HM_TUP_RIGHTWRAPAROUND 1
#include "bgq_HoppingMatrix_site.inc.c"
		}

	}




#ifndef BGQ_HM_TLINE_NOFUNC
}
#undef BGQ_HM_SITE_NOFUNC
#endif


#undef BGQ_HM_TLINE_TLINEINDENT
