
#ifndef BGQ_HM_XUP_WEYLREAD
#define BGQ_HM_XUP_WEYLREAD 0
#endif

#ifndef BGQ_HM_XUP_COMPUTE
#define BGQ_HM_XUP_COMPUTE 0
#endif

#ifndef BGQ_HM_XUP_WEYL_SEND
#define BGQ_HM_XUP_WEYL_SEND 0
#endif

#ifndef BGQ_HM_XUP_ACCUMULATE
#define BGQ_HM_XUP_ACCUMULATE 0
#endif

#ifndef BGQ_HM_DIR_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void bgq_HoppingMatrix_site_xup(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int x, int y, int z, int tv, int k) {
	bgq_su3_spinor_decl(result);
#endif
	{


		bgq_su3_weyl_decl(weyl_xup);
#if BGQ_HM_XUP_WEYLREAD==-1
		if (x==PHYSICAL_LX-1) {
#endif
#if (BGQ_HM_XUP_WEYLREAD==-1) || (BGQ_HM_XUP_WEYLREAD==1)
		bgq_weylsite_double *weylsite_xup = BGQ_WEYLSITE_X(weylxchange_recv_double[X_UP], !isOdd, t, x+1, y, zv);
		bgq_su3_weyl_double_load(weyl_xup, weylsite_xup);
#endif
#if BGQ_HM_XUP_WEYLREAD==-1
		} else {
#endif
#if (BGQ_HM_XUP_WEYLREAD==-1) || (BGQ_HM_XUP_WEYLREAD==0)
		// Load the input spinor
		bgq_su3_spinor_decl(spinor_xup);
		bgq_spinorsite_double *spinorsite_xup = BGQ_SPINORSITE(spinorfield, !isOdd, t, x+1, y, zv);
		bgq_su3_spinor_double_load(spinor_xup, spinorsite_xup);

		// Compute its halfspinor
		bgq_su3_vadd(weyl_xup_v0, spinor_xup_v0, spinor_xup_v2);
		bgq_su3_vadd(weyl_xup_v1, spinor_xup_v1, spinor_xup_v3);
#endif
#if BGQ_HM_XUP_WEYLREAD==-1
		}
#endif


#if BGQ_HM_XUP_COMPUTE
		// Load the interaction matrix between the lattice sites
		bgq_su3_mdecl(gauge_xup);
		bgq_gaugesite_double *gaugesite_xup = BGQ_GAUGESITE(gaugefield, isOdd, t, x, y, zv, X_UP);
		bgq_su3_matrix_double_load(gauge_xup, gaugesite_xup);

		// Multiply the halfspinor with the matrix
		bgq_su3_mvmul(weyl_xup_v0, gauge_xup, weyl_xup_v0);
		bgq_su3_mvmul(weyl_xup_v1, gauge_xup, weyl_xup_v1);

#ifndef BGQ_HM_NOKAMUL
		// Multiply with custom constant
		bgq_su3_cvmul(weyl_xup_v0, qka0, weyl_xup_v0);
		bgq_su3_cvmul(weyl_xup_v1, qka0, weyl_xup_v1);
#endif
#endif


#if BGQ_HM_XUP_WEYL_SEND
		// Store the halfspinor to be transfered to the neighbor node
		bgq_su3_weyl_double_store(weylxchange_send_double[X_UP], weyl_xup);
#endif


#if BGQ_HM_XUP_ACCUMULATE
		// Add up at the output lattice site
		bgq_su3_vmov(result_v0, weyl_xup_v0);
		bgq_su3_vmov(result_v1, weyl_xup_v1);
		bgq_su3_vmov(result_v2, weyl_xup_v0);
		bgq_su3_vmov(result_v3, weyl_xup_v1);
#endif


	}
#ifndef BGQ_HM_DIR_NOFUNC
}
#endif

#undef BGQ_HM_XUP_WEYLREAD
#undef BGQ_HM_XUP_COMPUTE
#undef BGQ_HM_XUP_WEYL_SEND
#undef BGQ_HM_XUP_ACCUMULATE

