/*
 * bgq_stdoperators.h
 *
 *  Created on: Nov 21, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_STDOPERATORS_H_
#define BGQ_STDOPERATORS_H_

#include "bgq_spinorfield.h"


void bgq_spinorfield_diff_double(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *argfield1, bgq_weylfield_controlblock *argfield2);
void bgq_spinorfield_diff_float(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *argfield1, bgq_weylfield_controlblock *argfield2);

void bgq_spinorfield_copy(bgq_weylfield_controlblock *target, bgq_spinorfield_layout targetLayout, bgq_weylfield_controlblock *source);

void bgq_spinorfield_cmul_double(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield, double c);
void bgq_spinorfield_cmul_float(bgq_weylfield_controlblock *targetfield, tristate isOdd, bgq_weylfield_controlblock *sourcefield, double c);

#endif /* BGQ_STDOPERATORS_H_ */
