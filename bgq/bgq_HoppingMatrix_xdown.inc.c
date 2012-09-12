
#ifndef BGQ_HM_XDOWN_PREFETCH
#define BGQ_HM_XDOWN_PREFETCH 0
#endif

#ifndef BGQ_HM_XDOWN_COMPUTE
#define BGQ_HM_XDOWN_COMPUTE 0
#endif

#ifndef BGQ_HM_XDOWN_WEYL_SEND
#define BGQ_HM_XDOWN_WEYL_SEND 0
#endif

#ifndef BGQ_HM_XDOWN_ACCUMULATE
#define BGQ_HM_XDOWN_ACCUMULATE 0
#endif


#define BGQ_LOADORPREFETCH_PREFETCH BGQ_HM_XDOWN_PREFETCH
#define BGQ_LOADORPREFETCH_LOAD !BGQ_HM_XDOWN_PREFETCH
#include "bgq_loadorprefetch.inc.c"


#ifndef BGQ_HM_DIR_NOFUNC
#include "bgq.h"
#include "bgq_field.h"
void bgq_HoppingMatrix_xdown(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int tv, int x, int y, int z, int k, int t1, int 22) {
	bgq_su3_spinor_decl(result);
#else
{
#endif


		bgq_su3_weyl_decl(weyl_xdown);
#if BGQ_HM_XDOWN_WEYLREAD==-1
		if (x==0) {
#endif
#if (BGQ_HM_XDOWN_WEYLREAD==-1) || (BGQ_HM_XDOWN_WEYLREAD==1)
		bgq_weylsite *weylsite_xdown = BGQ_WEYLSITE_X(weylxchange_recv[XDOWN], !isOdd, tv, x-1, y, z, t1, t2, !BGQ_HM_XDOWN_PREFETCH,false);
		bgq_su3_weyl_loadorprefetch(weyl_xdown, weylsite_xdown);

		#if !BGQ_HM_XDOWN_PREFETCH
			bgq_setbgqvalue(t1, x, y, z, BGQREF_XDOWN_WEYLREAD, bgq_cmplxval1(weyl_xdown_v0_c0), "weyl_recv_xdown");
			bgq_setbgqvalue(t2, x, y, z, BGQREF_XDOWN_WEYLREAD, bgq_cmplxval2(weyl_xdown_v0_c0), "weyl_recv_xdown");
		#endif
#endif
#if BGQ_HM_XDOWN_WEYLREAD==-1
		} else {
#endif
#if (BGQ_HM_XDOWN_WEYLREAD==-1) || (BGQ_HM_XDOWN_WEYLREAD==0)
		// Load the input spinor
		bgq_su3_spinor_decl(spinor_xdown);
		bgq_spinorsite *spinorsite_xdown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x-1, y, z, t1, t2, !BGQ_HM_XDOWN_PREFETCH,false);
		bgq_su3_spinor_loadorprefetch(spinor_xdown, spinorsite_xdown);

		#if !BGQ_HM_XDOWN_PREFETCH
			// Compute its halfspinor
			bgq_su3_vpisub(weyl_xdown_v0, spinor_xdown_v0, spinor_xdown_v3);
			bgq_su3_vpisub(weyl_xdown_v1, spinor_xdown_v1, spinor_xdown_v2);
		#endif
#endif
#if BGQ_HM_XDOWN_WEYLREAD==-1
		}
#endif


#if BGQ_HM_XDOWN_COMPUTE
		// Load the between sites-interaction matrix
		bgq_su3_mdecl(gauge_xdown);
		bgq_gaugesite *gaugesite_xdown = BGQ_GAUGESITE(gaugefield, !isOdd, tv, x-1, y, z, XUP, t1, t2, !BGQ_HM_XDOWN_PREFETCH,false);
		bgq_su3_matrix_loadorprefetch(gauge_xdown, gaugesite_xdown);

		// Multiply the halfspinor with the matrix
		bgq_su3_mvinvmul(weyl_xdown_v0, gauge_xdown, weyl_xdown_v0);
		bgq_su3_mvinvmul(weyl_xdown_v1, gauge_xdown, weyl_xdown_v1);

#if !BGQ_HM_NOKAMUL
		bgq_su3_cvmul(weyl_xdown_v0, qka1, weyl_xdown_v0);
		bgq_su3_cvmul(weyl_xdown_v1, qka1, weyl_xdown_v1);
#endif
#endif


#if BGQ_HM_XDOWN_WEYL_SEND
		// Store the halfspinor to be transfered to the neighbor node
		bgq_weylsite *weylsite_xdown = BGQ_WEYLSITE_X(weylxchange_send[XUP/*!!!*/], !isOdd, tv, x-1, y, z, t1,t2, false,true);
		bgq_su3_weyl_zeroload(weylsite_xdown);
		bgq_su3_weyl_store(weylsite_xdown, weyl_xdown);
#endif


#if BGQ_HM_XDOWN_ACCUMULATE
		bgq_su3_vadd(result_v0, result_v0, weyl_xdown_v0);
		bgq_su3_vadd(result_v1, result_v1, weyl_xdown_v1);
		bgq_su3_vpiadd(result_v2, result_v2, weyl_xdown_v1);
		bgq_su3_vpiadd(result_v3, result_v3, weyl_xdown_v0);

		bgq_setbgqvalue(t1, x, y, z, BGQREF_XDOWN, bgq_cmplxval1(weyl_xdown_v0_c0), "weyl_xdown");
		bgq_setbgqvalue(t2, x, y, z, BGQREF_XDOWN, bgq_cmplxval2(weyl_xdown_v0_c0), "weyl_xdown");
#endif


}


#include "bgq_loadorprefetch.inc.c"

#undef BGQ_HM_XDOWN_PREFETCH
#undef BGQ_HM_XDOWN_COMPUTE
#undef BGQ_HM_XDOWN_WEYL_SEND
#undef BGQ_HM_XDOWN_ACCUMULATE
