/*
 * bgq_stdoperators.c
 *
 *  Created on: Nov 21, 2012
 *      Author: meinersbur
 */

#include "bgq_stdoperators.h"

#include "bgq_field.h"
#include "bgq_qpx.h"



static inline void bgq_site_diff(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor1), bgq_su3_spinor_params(spinor2), ucoord ic) {
	bgq_su3_spinor_decl(result);

	bgq_su3_vsub(result_v0, spinor1_v0, spinor2_v0);
	bgq_su3_vsub(result_v1, spinor1_v1, spinor2_v1);
	bgq_su3_vsub(result_v2, spinor1_v2, spinor2_v2);
	bgq_su3_vsub(result_v3, spinor1_v3, spinor2_v3);

	bgq_su3_spinor_mov(*target, result);
}


static inline void bgq_site_set(bgq_su3_spinor_params(*target), bgq_su3_spinor_params(spinor), ucoord ic) {
	bgq_su3_spinor_decl(result);

	bgq_su3_spinor_mov(result, spinor);

	bgq_su3_spinor_mov(*target, result);
}


#define OPERATOR_NAME bgq_spinorfield_diff
#define OPERATOR_ARGFIELDS 2
#define OPERATOR_VECSITEFUNC bgq_site_diff
#include "bgq_operator.inc.c"


#define OPERATOR_NAME bgq_spinorfield_set
#define OPERATOR_ARGFIELDS 0
#define OPERATOR_VECSITEFUNC bgq_site_set
#define OPERATOR_EXTRAPARMS bgq_su3_spinor_params(spinor)
#define OPERATOR_EXTRAARGS bgq_su3_spinor_vars(spinor)
#include "bgq_operator.inc.c"

