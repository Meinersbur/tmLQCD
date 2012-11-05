
#include "bgq_HoppingMatrix.h"
#include "bgq_comm.h"

void *somewhere;
bgq_spinor bgq_spinorfield_getspinor(bgq_weylfield_controlblock *field, ucoord t, ucoord x, ucoord y, ucoord z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);

	assert(bgq_local2isOdd(t,x,y,z)==field->isOdd);
	//ucoord ic = bgq_local2collapsed(t,x,y,z);
	ucoord k = bgq_local2k(t,x,y,z);
	for (ucoord ic = 0; ic < t; ic+=1) {
		bgq_spinorfield_setup(field, field->isOdd, false, false, true, false);
		bgq_su3_spinor_decl(spinor);
		bgq_HoppingMatrix_loadWeyllayout(spinor,&field->sec_collapsed[ic], bgq_t2t(t,0), bgq_t2t(t,1), x, y, z);
		//return bgq_spinor_fromqpx(spinor,k);
		bgq_su3_spinor_store_double(somewhere,spinor);
	}
}


