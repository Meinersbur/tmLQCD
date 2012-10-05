/*
 * bgq_HoppingMatrix_fixedvol.inc.c
 *
 *  Created on: Oct 5, 2012
 *      Author: meinersbur
 */


#define ODDNESS_NAME evenodd
#undef ODDNESS
#include "bgq_HoppingMatrix_base.inc.c"
#undef ODDNESS_NAME


#define ODDNESS_NAME even
#define ODDNESS 0
#include "bgq_HoppingMatrix_base.inc.c"
#undef ODDNESS_NAME
#undef ODDNESS


#define ODDNESS_NAME odd
#define ODDNESS 1
#include "bgq_HoppingMatrix_base.inc.c"
#undef ODDNESS_NAME
#undef ODDNESS
