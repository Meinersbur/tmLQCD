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

		// Load constants
		bgq_vector4double_decl(qka0);
		// t
		bgq_cconst(qka0, creal(ka0), cimag(ka0));
		bgq_vector4double_decl(qka1);
		// x
		bgq_cconst(qka1, creal(ka1), cimag(ka1));
		bgq_vector4double_decl(qka2);
		// y
		bgq_cconst(qka2, creal(ka2), cimag(ka2));
		bgq_vector4double_decl(qka3);
		// z
		bgq_cconst(qka3, creal(ka3), cimag(ka3));

		bgq_su3_weyl_decl(weyl_zcarry);

#if (BGQ_HM_ZLINE_ZLINEINDENT==-1)
		if (((t + x + y)&1) == isOdd) {
#endif
#if (BGQ_HM_ZLINE_ZLINEINDENT==-1) || (BGQ_HM_ZLINE_ZLINEINDENT==0)
			assert(((t+x+y)&1)==isOdd);
			// zline indention == 0
			// iterate from zbottom to ztop


			// Prologue
			{
				// We need the zcarry from the wraparound
				// We fake a site update from the zup-border such that it computes the carry for the first iteration
				const int zv = PHYSICAL_LZV - 1; // wraparound
				const int z1 = LOCAL_LZ-4;
				const int z2 = LOCAL_LZ-2;

				// Read the zcarry value for the first iteration
				#define BGQ_HM_ZUP_ZLINEINDENT 0
				#define BGQ_HM_ZUP_WRITECARRY 1
				#define BGQ_HM_ZUP_RIGHTWRAPAROUND 1
				#include "bgq_HoppingMatrix_zup.inc.c"
			}

			{
				int z1 = 0;
				int z2 = 2;

				for (int zv = 1; zv < PHYSICAL_LZV-1; zv+=1) {
					#define BGQ_HM_ZUP_ZLINEINDENT 0
					#define BGQ_HM_ZDOWN_ZLINEINDENT 0
					#define BGQ_HM_TUP_COMPUTE 1
					#define BGQ_HM_TDOWN_COMPUTE 1
					#define BGQ_HM_XUP_COMPUTE 1
					#define BGQ_HM_XDOWN_COMPUTE 1
					#define BGQ_HM_YUP_COMPUTE 1
					#define BGQ_HM_YDOWN_COMPUTE 1
					#define BGQ_HM_ZUP_COMPUTE 1
					#define BGQ_HM_ZDOWN_COMPUTE 1
					#define BGQ_HM_ZDOWN_READCARRYSPINOR 1
					#define BGQ_HM_ZUP_WRITECARRYSPINOR 1
					#include "bgq_HoppingMatrix_site.inc.c"

					z1 += PHYSICAL_LP * PHYSICAL_LK;
					z2 += PHYSICAL_LP * PHYSICAL_LK;
				}
			}

#endif
#if (BGQ_HM_ZLINE_ZLINEINDENT==-1)
			} else {
#endif
#if (BGQ_HM_ZLINE_ZLINEINDENT==-1) || (BGQ_HM_ZLINE_ZLINEINDENT==1)
				assert(((t+x+y)&1)==!isOdd);
				// zline indention == 1
				// iterate backwards from ztop to zbottom
				int z1 = LOCAL_LZ-1-2;
				int z2 = LOCAL_LZ-1;

				// Prologue
				{
					const int zv = 0; // wraparound

					// Read the zcarry value for the first iteration
					#define BGQ_HM_ZDOWN_ZLINEINDENT 1
					#define BGQ_HM_ZDOWN_WRITECARRY 1
					#define BGQ_HM_ZDOWN_RIGHTWRAPAROUND 1
					#include "bgq_HoppingMatrix_zdown.inc.c"
				}

				for (int zv = 1; zv < PHYSICAL_LZV; zv+=1) {
					#define BGQ_HM_ZUP_ZLINEINDENT 1
					#define BGQ_HM_ZDOWN_ZLINEINDENT 1
					#define BGQ_HM_TUP_COMPUTE 1
					#define BGQ_HM_TDOWN_COMPUTE 1
					#define BGQ_HM_XUP_COMPUTE 1
					#define BGQ_HM_XDOWN_COMPUTE 1
					#define BGQ_HM_YUP_COMPUTE 1
					#define BGQ_HM_YDOWN_COMPUTE 1
					#define BGQ_HM_ZUP_COMPUTE 1
					#define BGQ_HM_ZDOWN_COMPUTE 1
					#define BGQ_HM_ZDOWN_READCARRYSPINOR 1
					#define BGQ_HM_ZUP_WRITECARRYSPINOR 1
					#include "bgq_HoppingMatrix_site.inc.c"

					z1 -= PHYSICAL_LP * PHYSICAL_LK;
					z2 -= PHYSICAL_LP * PHYSICAL_LK;
				}

#endif
#if (BGQ_HM_ZLINE_ZLINEINDENT==-1)
			}
#endif

		}
#ifndef BGQ_HM_ZLINE_NOFUNC
	}
#undef BGQ_HM_SITE_NOFUNC
#endif

#undef BGQ_HM_ZLINE_ZLINEINDENT

