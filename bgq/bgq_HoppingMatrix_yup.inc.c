
#ifndef BGQ_HM_YUP_WEYLREAD
#define BGQ_HM_YUP_WEYLREAD 0
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

#ifndef BGQ_HM_SITE_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void bgq_HoppingMatrix_yup(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int tv, int x, int y, int z, int t1, int t2) {
	bgq_su3_spinor_decl(result);
#endif
	{


		bgq_su3_weyl_decl(weyl_yup);
#if BGQ_HM_YUP_WEYLREAD==-1
		if (y==PHYSICAL_LY-1) {
#endif
#if (BGQ_HM_YUP_WEYLREAD==-1) || (BGQ_HM_YUP_WEYLREAD==1)
		bgq_weylsite_double *weylsite_yup = BGQ_WEYLSITE_Y(weylxchange_recv_double[YUP], !isOdd, tv, x, y+1, z, t1, t2);
		bgq_su3_weyl_double_load(weyl_yup, weylsite_yup);
#endif
#if BGQ_HM_YUP_WEYLREAD==-1
	} else {
#endif
#if (BGQ_HM_YUP_WEYLREAD==-1) || (BGQ_HM_YUP_WEYLREAD==0)
		// Load the input spinor
		bgq_su3_spinor_decl(spinor_yup);
		bgq_spinorsite_double *spinorsite_yup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y+1, z, t1,t2, true,false);
		bgq_su3_spinor_double_load(spinor_yup, spinorsite_yup);

		// Compute its halfspinor
		bgq_su3_vadd(weyl_yup_v0, spinor_yup_v0, spinor_yup_v3);
		bgq_su3_vsub(weyl_yup_v1, spinor_yup_v1, spinor_yup_v2);
#endif
#if BGQ_HM_YUP_WEYLREAD==-1
	}
#endif


#if BGQ_HM_YUP_COMPUTE
		// Load the interaction matrix between the lattice sites
		bgq_su3_mdecl(gauge_yup);
		bgq_gaugesite_double *gaugesite_yup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, YUP, t1, t2, true,false);
		bgq_su3_matrix_double_load(gauge_yup, gaugesite_yup);

		// Multiply the halfspinor with the matrix
		bgq_su3_mvmul(weyl_yup_v0, gauge_yup, weyl_yup_v0);
		bgq_su3_mvmul(weyl_yup_v1, gauge_yup, weyl_yup_v1);

#ifndef BGQ_HM_NOKAMUL
		bgq_su3_cvmul(weyl_yup_v0, qka2, weyl_yup_v0);
		bgq_su3_cvmul(weyl_yup_v1, qka2, weyl_yup_v1);
#endif
#endif


#if BGQ_HM_YUP_WEYL_SEND
		// Store the halfspinor to be transfered to the neighbor node
		bgq_weylsite_double *weylsite_yup = BGQ_WEYLSITE_Y(weylxchange_send_double[YDOWN/*!!!*/], isOdd, tv, x, y, z, t1,t2);
		bgq_su3_weyl_double_store(weylsite_yup, weyl_yup);
#endif


#if BGQ_HM_YUP_ACCUMULATE
		// Add up at the output lattice site
		bgq_su3_vadd(result_v0, result_v0, weyl_yup_v0);
		bgq_su3_vadd(result_v1, result_v1, weyl_yup_v1);
		bgq_su3_vsub(result_v2, result_v2, weyl_yup_v1);
		bgq_su3_vadd(result_v3, result_v3, weyl_yup_v0);
#endif


	}
#ifndef BGQ_HM_SITE_NOFUNC
}
#endif


#undef BGQ_HM_YUP_WEYLREAD
#undef BGQ_HM_YUP_COMPUTE
#undef BGQ_HM_YUP_WEYL_SEND
#undef BGQ_HM_YUP_ACCUMULATE
