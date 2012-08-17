/*
 * bgq_field_float.c
 *
 *  Created on: Aug 15, 2012
 *      Author: meinersbur
 */

#define BGQ_FIELD_FLOAT_C_

#include "../global.h"
#include "bgq_field_float.h"


#define BGQ_PRECISION 32
#include "bgq_precisionselect.inc.c"

#include "bgq_field.inc.c"
