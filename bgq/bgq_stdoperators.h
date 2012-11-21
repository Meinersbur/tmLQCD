/*
 * bgq_stdoperators.h
 *
 *  Created on: Nov 21, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_STDOPERATORS_H_
#define BGQ_STDOPERATORS_H_

#include "bgq_spinorfield.h"

void bgq_spinorfield_diff_double(bgq_weylfield_controlblock *targetfield, bgq_weylfield_controlblock *argfield1, bgq_weylfield_controlblock *argfield2);
void bgq_spinorfield_diff_float(bgq_weylfield_controlblock *targetfield, bgq_weylfield_controlblock *argfield1, bgq_weylfield_controlblock *argfield2);

#endif /* BGQ_STDOPERATORS_H_ */
