
#ifndef BGQ_HM_YUP_PREFETCH
#define BGQ_HM_YUP_PREFETCH 0
#endif

#ifndef BGQ_HM_YUP_COMPUTE
#define BGQ_HM_YUP_COMPUTE 0
#endif

#ifndef BGQ_HM_YUP_WEYL_SEND
#define BGQ_HM_YUP_WEYL_SEND 0
#endif

#ifndef BGQ_HM_YUP_ACCUMULATE
#define BGQ_HM_YUP_ACCUMULATE 0
#endif


#define BGQ_LOADORPREFETCH_PREFETCH BGQ_HM_YUP_PREFETCH
#define BGQ_LOADORPREFETCH_LOAD !BGQ_HM_YUP_PREFETCH
#include "bgq_loadorprefetch.inc.c"


#ifndef BGQ_HM_SITE_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void bgq_HoppingMatrix_yup(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int tv, int x, int y, int z, int t1, int t2) {
	bgq_su3_spinor_decl(result);
#else
{
#endif


	bgq_su3_weyl_decl(weyl_yup);
	#if BGQ_HM_YUP_WEYLREAD==-1
	if (y==PHYSICAL_LY-1) {
	#endif
	#if (BGQ_HM_YUP_WEYLREAD==-1) || (BGQ_HM_YUP_WEYLREAD==1)
		bgq_weylsite *weylsite_yup = BGQ_WEYLSITE_Y(weylxchange_recv[YUP], !isOdd, tv, x, y+1, z, t1, t2, !BGQ_HM_YUP_PREFETCH,false);
		bgq_su3_weyl_loadorprefetch(weyl_yup, weylsite_yup);
	#endif
	#if (BGQ_HM_YUP_WEYLREAD==-1)
	} else {
	#endif
	#if (BGQ_HM_YUP_WEYLREAD==-1) || (BGQ_HM_YUP_WEYLREAD==0)
		// Load the input spinor
		bgq_su3_spinor_decl(spinor_yup);
		bgq_spinorsite *spinorsite_yup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y+1, z, t1, t2, !BGQ_HM_YUP_PREFETCH,false);
		bgq_su3_spinor_loadorprefetch(spinor_yup, spinorsite_yup);

		#if !BGQ_HM_YUP_PREFETCH
			// Compute its halfspinor
			bgq_su3_vadd(weyl_yup_v0, spinor_yup_v0, spinor_yup_v3);
			bgq_su3_vsub(weyl_yup_v1, spinor_yup_v1, spinor_yup_v2);
		#endif
	#endif
	#if BGQ_HM_YUP_WEYLREAD==-1
	}
	#endif


	#if BGQ_HM_YUP_COMPUTE
		// Load the interaction matrix between the lattice sites
		bgq_su3_mdecl(gauge_yup);
		bgq_gaugesite *gaugesite_yup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, YUP, t1, t2, !BGQ_HM_YUP_PREFETCH,false);
		bgq_su3_matrix_loadorprefetch(gauge_yup, gaugesite_yup);

		// Multiply the halfspinor with the matrix
		bgq_su3_mvmul(weyl_yup_v0, gauge_yup, weyl_yup_v0);
		bgq_su3_mvmul(weyl_yup_v1, gauge_yup, weyl_yup_v1);

		#if !BGQ_HM_NOKAMUL
			// Multiply with custom constant
			bgq_su3_cvmul(weyl_yup_v0, qka2, weyl_yup_v0);
			bgq_su3_cvmul(weyl_yup_v1, qka2, weyl_yup_v1);
		#endif
	#endif


	#if BGQ_HM_YUP_WEYL_SEND
		// Store the halfspinor to be transfered to the neighbor node
		bgq_weylsite *weylsite_yup = BGQ_WEYLSITE_Y(weylxchange_send[YDOWN/*!!!*/], !isOdd, tv, x, y+1, z, t1,t2, false,true);
		bgq_su3_weyl_zeroload(weylsite_yup);
		bgq_su3_weyl_store(weylsite_yup, weyl_yup);
		bgq_su3_weyl_flush(weylsite_yup);
	#endif


	#if BGQ_HM_YUP_ACCUMULATE
		// Add up at the output lattice site
		bgq_su3_vadd(result_v0, result_v0, weyl_yup_v0);
		bgq_su3_vadd(result_v1, result_v1, weyl_yup_v1);
		bgq_su3_vsub(result_v2, result_v2, weyl_yup_v1);
		bgq_su3_vadd(result_v3, result_v3, weyl_yup_v0);

		bgq_setbgqvalue(t1, x, y, z, BGQREF_YUP, bgq_cmplxval1(weyl_yup_v1_c0), "weyl_yup");
		bgq_setbgqvalue(t2, x, y, z, BGQREF_YUP, bgq_cmplxval2(weyl_yup_v1_c0), "weyl_yup");
	#endif


}


#include "bgq_loadorprefetch.inc.c"

#undef BGQ_HM_YUP_PREFETCH
#undef BGQ_HM_YUP_COMPUTE
#undef BGQ_HM_YUP_WEYL_SEND
#undef BGQ_HM_YUP_ACCUMULATE
