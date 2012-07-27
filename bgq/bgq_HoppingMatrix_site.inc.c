#ifndef BGQ_HM_NOFUNC
#include "bgq.h"
#include "bgq_field.h"

void HoppingMatrix_site(bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefieldeo_double gaugefield, bool isOdd, int x, int y, int z, int tv, int k) {
#endif

	bgq_su3_spinor_decl(result);

#define BGQ_HM_NOFUNC

// direction X_UP /////////////////////////////////////////////////////////////
#include "bgq_HoppingMatrix_xup.inc.c"

// direction X_DOWN /////////////////////////////////////////////////////////////
#include "bgq_HoppingMatrix_xdown.inc.c"

// direction Y_UP /////////////////////////////////////////////////////////////
#include "bgq_HoppingMatrix_yup.inc.c"

// direction Y_DOWN /////////////////////////////////////////////////////////////
#include "bgq_HoppingMatrix_ydown.inc.c"

// direction Z_UP /////////////////////////////////////////////////////////////
#include "bgq_HoppingMatrix_zup.inc.c"

	// direction Z_DOWN /////////////////////////////////////////////////////////////
#include "bgq_HoppingMatrix_zdown.inc.c"

	// direction T_UP /////////////////////////////////////////////////////////////
#include "bgq_HoppingMatrix_tup.inc.c"

	// direction T_DOWN /////////////////////////////////////////////////////////////
#include "bgq_HoppingMatrix_tdown.inc.c"

///////////////////////////////////////////////////////////////////////////////
// Store the spinor

bgq_spinorsite_double
		*
spinor_target
=
		bgq_spinorsite_double_physical_pointer(
				s
			pin
orfield,
isOdd,
x,
y,
z			, tv);
	bgq_su3_spinor_double_store(spinor_target, result);

#ifndef BGQ_HM_NOFUNC
}
#endif

