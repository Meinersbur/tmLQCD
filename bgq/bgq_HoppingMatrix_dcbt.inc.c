/*
 * bgq_HoppingMatrix_dcbt.inc.c
 *
 *  Created on: Oct 5, 2012
 *      Author: meinersbur
 */



#define DCBT_NAME dcbt
#define BGQ_DCBT_DISABLE 0
#include "bgq_HoppingMatrix_oddness.inc.c"
#undef DCBT_NAME
#undef BGQ_DCBT_DISABLE



#define DCBT_NAME nodcbt
#define BGQ_DCBT_DISABLE 1
#include "bgq_HoppingMatrix_oddness.inc.c"
#undef DCBT_NAME
#undef BGQ_DCBT_DISABLE
