
#ifndef BGQ_HM_YDOWN_PREFETCH
#define BGQ_HM_YDOWN_PREFETCH 0
#endif

#ifndef BGQ_HM_YDOWN_COMPUTE
#define BGQ_HM_YDOWN_COMPUTE 0
#endif

#ifndef BGQ_HM_YDOWN_WEYL_SEND
#define BGQ_HM_YDOWN_WEYL_SEND 0
#endif

#ifndef BGQ_HM_YDOWN_ACCUMULATE
#define BGQ_HM_YDOWN_ACCUMULATE 0
#endif


#define BGQ_LOADORPREFETCH_PREFETCH BGQ_HM_YDOWN_PREFETCH
#define BGQ_LOADORPREFETCH_LOAD !BGQ_HM_YDOWN_PREFETCH
#include "bgq_loadorprefetch.inc.c"


#ifndef BGQ_HM_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void bgq_HoppingMatrix_ydown(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int tv, int x, int y, int z, int t1, int t2) {
	bgq_su3_spinor_decl(result);
#else
{
#endif


	bgq_su3_weyl_decl(weyl_ydown);
	#if (BGQ_HM_YDOWN_WEYLREAD==-1)
	if (y==0) {
	#endif
	#if (BGQ_HM_YDOWN_WEYLREAD==-1) || (BGQ_HM_YDOWN_WEYLREAD==1)
		bgq_weylsite *weylsite_ydown = BGQ_WEYLSITE_Y(weylxchange_recv[YDOWN], !isOdd, tv, x, y-1, z, t1, t2, !BGQ_HM_YDOWN_PREFETCH,false);
		bgq_su3_weyl_loadorprefetch(weyl_ydown, weylsite_ydown);
	#endif
	#if (BGQ_HM_YDOWN_WEYLREAD==-1)
	} else {
	#endif
	#if (BGQ_HM_YDOWN_WEYLREAD==-1) || (BGQ_HM_YDOWN_WEYLREAD==0)
		// Load the input spinor
		bgq_su3_spinor_decl(spinor_ydown);
		bgq_spinorsite *spinorsite_ydown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y-1, z, t1,t2, !BGQ_HM_YDOWN_PREFETCH,false);
		bgq_su3_spinor_loadorprefetch(spinor_ydown, spinorsite_ydown);

		#if !BGQ_HM_YDOWN_PREFETCH
			// Compute its halfspinor
			bgq_su3_vsub(weyl_ydown_v0, spinor_ydown_v0, spinor_ydown_v3);
			bgq_su3_vadd(weyl_ydown_v1, spinor_ydown_v1, spinor_ydown_v2);
		#endif
	#endif
	#if (BGQ_HM_YDOWN_WEYLREAD==-1)
	}
	#endif



#if BGQ_HM_YDOWN_COMPUTE
		bgq_su3_mdecl(gauge_ydown);
		bgq_gaugesite *gaugesite_ydown = BGQ_GAUGESITE(gaugefield, !isOdd, tv, x, y-1, z, YUP, t1, t2, !BGQ_HM_YDOWN_PREFETCH,false);
		bgq_su3_matrix_loadorprefetch(gauge_ydown, gaugesite_ydown);

		bgq_su3_mvinvmul(weyl_ydown_v0, gauge_ydown, weyl_ydown_v0);
		bgq_su3_mvinvmul(weyl_ydown_v1, gauge_ydown, weyl_ydown_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(weyl_ydown_v0, qka2, weyl_ydown_v0);
		bgq_su3_cvmul(weyl_ydown_v1, qka2, weyl_ydown_v1);
#endif
#endif


#if BGQ_HM_YDOWN_WEYL_SEND
		// Store the halfspinor to be transfered to the neighbor node
		bgq_weylsite *weylsite_ydown = BGQ_WEYLSITE_Y(weylxchange_send[YUP/*!!!*/], !isOdd, tv, x, y-1, z, t1,t2, false,true);
		bgq_su3_weyl_zeroload(weylsite_ydown);
		bgq_su3_weyl_store(weylsite_ydown, weyl_ydown);
#endif


#if BGQ_HM_YDOWN_ACCUMULATE
		bgq_su3_vadd(result_v0, result_v0, weyl_ydown_v0);
		bgq_su3_vadd(result_v1, result_v1, weyl_ydown_v1);
		bgq_su3_vadd(result_v2, result_v2, weyl_ydown_v1);
		bgq_su3_vsub(result_v3, result_v3, weyl_ydown_v0);
#endif


}


#include "bgq_loadorprefetch.inc.c"

#undef BGQ_HM_YDOWN_PREFETCH
#undef BGQ_HM_YDOWN_COMPUTE
#undef BGQ_HM_YDOWN_WEYL_SEND
#undef BGQ_HM_YDOWN_ACCUMULATE
