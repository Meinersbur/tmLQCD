/*
 * bgq_field_double.h
 *
 *  Created on: Aug 15, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_FIELD_DOUBLE_H_
#define BGQ_FIELD_DOUBLE_H_

#include <stdint.h>


#define BGQ_PRECISION 64
#include "bgq_precisionselect.inc.c"

#include "bgq_field.inc.h"

#undef BGQ_PRECISION
#include "bgq_precisionselect.inc.c"

#endif /* BGQ_FIELD_DOUBLE_H_ */
