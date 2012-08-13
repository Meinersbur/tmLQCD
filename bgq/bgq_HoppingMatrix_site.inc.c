
#ifndef BGQ_HM_SITE_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void HoppingMatrix_site(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool isOdd, int t, int x, int y, int z1, int z2. int tv, int k) {
#define BGQ_HM_DIR_NOFUNC 1
#else
{
#endif
	bgq_su3_spinor_decl(result);
	bgq_su3_spinor_zero(result);


// direction T_UP /////////////////////////////////////////////////////////////
#define BGQ_HM_TUP_FIRST 1
#define BGQ_HM_TUP_COMPUTE 1
#define BGQ_HM_TUP_ACCUMULATE 1
#include "bgq_HoppingMatrix_tup.inc.c"

// direction T_DOWN ///////////////////////////////////////////////////////////
#define BGQ_HM_TDOWN_COMPUTE 1
#define BGQ_HM_TDOWN_ACCUMULATE 1
#include "bgq_HoppingMatrix_tdown.inc.c"

// direction X_UP /////////////////////////////////////////////////////////////
#define BGQ_HM_XUP_COMPUTE 1
#define BGQ_HM_XUP_ACCUMULATE 1
#include "bgq_HoppingMatrix_xup.inc.c"

// direction X_DOWN /////////////////////////////////////////////////////////////
#define BGQ_HM_XDOWN_COMPUTE 1
#define BGQ_HM_XDOWN_ACCUMULATE 1
#include "bgq_HoppingMatrix_xdown.inc.c"

// direction Y_UP /////////////////////////////////////////////////////////////
#define BGQ_HM_YUP_COMPUTE 1
#define BGQ_HM_YUP_ACCUMULATE 1
#include "bgq_HoppingMatrix_yup.inc.c"

// direction Y_DOWN /////////////////////////////////////////////////////////////
#define BGQ_HM_YDOWN_COMPUTE 1
#define BGQ_HM_YDOWN_ACCUMULATE 1
#include "bgq_HoppingMatrix_ydown.inc.c"

// direction Z_UP /////////////////////////////////////////////////////////////
#define BGQ_HM_ZUP_COMPUTE 1
#define BGQ_HM_ZUP_ACCUMULATE 1
#include "bgq_HoppingMatrix_zup.inc.c"

// direction Z_DOWN /////////////////////////////////////////////////////////////
#define BGQ_HM_ZDOWN_COMPUTE 1
#define BGQ_HM_ZDOWN_ACCUMULATE 1
#include "bgq_HoppingMatrix_zdown.inc.c"

///////////////////////////////////////////////////////////////////////////////
// Store the spinor

	bgq_spinorsite_double *targetsite = BGQ_SPINORSITE(targetfield, isOdd, tv, x, y, z, t1,t2, false,true);
	bgq_su3_spinor_double_store(targetsite, result);


}


#ifndef BGQ_HM_SITE_NOFUNC
#undef BGQ_HM_DIR_NOFUNC 1
#endif
