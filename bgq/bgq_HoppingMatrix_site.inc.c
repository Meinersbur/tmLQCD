
#ifndef BGQ_HM_ZLINE_TUP_WEYLREAD
#define BGQ_HM_ZLINE_TUP_WEYLREAD 0
#endif

#ifndef BGQ_HM_ZLINE_TDOWN_WEYLREAD
#define BGQ_HM_ZLINE_TDOWN_WEYLREAD 0
#endif

#ifndef BGQ_HM_ZLINE_XUP_WEYLREAD
#define BGQ_HM_ZLINE_XUP_WEYLREAD 0
#endif

#ifndef BGQ_HM_ZLINE_XDOWN_WEYLREAD
#define BGQ_HM_ZLINE_XDOWN_WEYLREAD 0
#endif

#ifndef BGQ_HM_ZLINE_YUP_WEYLREAD
#define BGQ_HM_ZLINE_YUP_WEYLREAD 0
#endif

#ifndef BGQ_HM_ZLINE_YDOWN_WEYLREAD
#define BGQ_HM_ZLINE_YDOWN_WEYLREAD 0
#endif


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


	#if BGQ_PREFETCH_EXPLICIT
		#define BGQ_HM_TDOWN_PREFETCH 1
		#include "bgq_HoppingMatrix_tdown.inc.c"
	#endif

// direction T_UP /////////////////////////////////////////////////////////////
#define BGQ_HM_TUP_FIRST 1
#define BGQ_HM_TUP_COMPUTE 1
#define BGQ_HM_TUP_ACCUMULATE 1
#include "bgq_HoppingMatrix_tup.inc.c"


	#if BGQ_PREFETCH_EXPLICIT
		#define BGQ_HM_XUP_PREFETCH 1
		#include "bgq_HoppingMatrix_xup.inc.c"
	#endif

// direction T_DOWN ///////////////////////////////////////////////////////////
#define BGQ_HM_TDOWN_COMPUTE 1
#define BGQ_HM_TDOWN_ACCUMULATE 1
#include "bgq_HoppingMatrix_tdown.inc.c"


	#if BGQ_PREFETCH_EXPLICIT
		#define BGQ_HM_XDOWN_PREFETCH 1
		#include "bgq_HoppingMatrix_xdown.inc.c"
	#endif

// direction X_UP /////////////////////////////////////////////////////////////
#define BGQ_HM_XUP_COMPUTE 1
#define BGQ_HM_XUP_ACCUMULATE 1
#include "bgq_HoppingMatrix_xup.inc.c"


	#if BGQ_PREFETCH_EXPLICIT
		#define BGQ_HM_YUP_PREFETCH 1
		#include "bgq_HoppingMatrix_yup.inc.c"
	#endif

// direction X_DOWN /////////////////////////////////////////////////////////////
#define BGQ_HM_XDOWN_COMPUTE 1
#define BGQ_HM_XDOWN_ACCUMULATE 1
#include "bgq_HoppingMatrix_xdown.inc.c"


	#if BGQ_PREFETCH_EXPLICIT
		#define BGQ_HM_YDOWN_PREFETCH 1
		#include "bgq_HoppingMatrix_ydown.inc.c"
	#endif

// direction Y_UP /////////////////////////////////////////////////////////////
#define BGQ_HM_YUP_COMPUTE 1
#define BGQ_HM_YUP_ACCUMULATE 1
#include "bgq_HoppingMatrix_yup.inc.c"


	#if BGQ_PREFETCH_EXPLICIT
		#define BGQ_HM_ZUP_PREFETCH 1
		#include "bgq_HoppingMatrix_zup.inc.c"
	#endif

// direction Y_DOWN /////////////////////////////////////////////////////////////
#define BGQ_HM_YDOWN_COMPUTE 1
#define BGQ_HM_YDOWN_ACCUMULATE 1
#include "bgq_HoppingMatrix_ydown.inc.c"


	#if BGQ_PREFETCH_EXPLICIT
		#define BGQ_HM_ZDOWN_PREFETCH 1
		#include "bgq_HoppingMatrix_zdown.inc.c"
	#endif

// direction Z_UP /////////////////////////////////////////////////////////////
#define BGQ_HM_ZUP_COMPUTE 1
#define BGQ_HM_ZUP_ACCUMULATE 1
#include "bgq_HoppingMatrix_zup.inc.c"


	#if BGQ_PREFETCH_EXPLICIT
		// prefetch for the next iteration
		const int z_shadowed = z;
		const int t1_shadowed = t1;
		const int t2_shadowed = t2;
		{
			const int z = z_shadowed + 1;
			#if (BGQ_HM_TLINEINDENT==-1)
			if ( (x+y+z)&1 == isOdd )
			#endif
			#if (BGQ_HM_TLINEINDENT==-1) || (BGQ_HM_TLINEINDENT==0)
				const int t1 = t1_shadowed + 1;
				const int t2 = t2_shadowed + 1;
			#endif
			#if (BGQ_HM_TLINEINDENT==-1)
			} else {
			#endif
			#if (BGQ_HM_TLINEINDENT==-1) || (BGQ_HM_TLINEINDENT==1)
				const int t1 = t1_shadowed - 1;
				const int t2 = t2_shadowed - 1;
			#endif
			#if (BGQ_HM_TLINEINDENT==-1)
			}
			#endif

			// Invert TLINEINDENT
			#if (BGQ_HM_TLINEINDENT==0)
			#undef BGQ_HM_TLINEINDENT
			#define BGQ_HM_TLINEINDENT 1
			#elif (BGQ_HM_TLINEINDENT==1)
			#undef BGQ_HM_TLINEINDENT
			#define BGQ_HM_TLINEINDENT 0
			#endif

			#define BGQ_HM_TDOWN_PREFETCH 1
			#include "bgq_HoppingMatrix_tdown.inc.c"

			// Revert TLINEINDENT
			#if (BGQ_HM_TLINEINDENT==0)
			#undef BGQ_HM_TLINEINDENT
			#define BGQ_HM_TLINEINDENT 1
			#elif (BGQ_HM_TLINEINDENT==1)
			#undef BGQ_HM_TLINEINDENT
			#define BGQ_HM_TLINEINDENT 0
			#endif
		}
	#endif

// direction Z_DOWN /////////////////////////////////////////////////////////////
#define BGQ_HM_ZDOWN_COMPUTE 1
#define BGQ_HM_ZDOWN_ACCUMULATE 1
#include "bgq_HoppingMatrix_zdown.inc.c"

///////////////////////////////////////////////////////////////////////////////
// Store the spinor


	bgq_spinorsite *targetsite = BGQ_SPINORSITE(targetfield, isOdd, tv, x, y, z, t1,t2, false,true);
	bgq_su3_spinor_zeroload(targetsite);
	bgq_su3_spinor_store(targetsite, result);

}


#ifndef BGQ_HM_SITE_NOFUNC
#undef BGQ_HM_DIR_NOFUNC 1
#endif


#undef BGQ_HM_ZLINE_TUP_WEYLREAD
#undef BGQ_HM_ZLINE_TDOWN_WEYLREAD
#undef BGQ_HM_ZLINE_XUP_WEYLREAD
#undef BGQ_HM_ZLINE_XDOWN_WEYLREAD
#undef BGQ_HM_ZLINE_YUP_WEYLREAD
#undef BGQ_HM_ZLINE_YDOWN_WEYLREAD

