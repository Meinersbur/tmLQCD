/*
 * bgq_field.c
 *
 *  Created on: Aug 15, 2012
 *      Author: meinersbur
 */

#include "bgq_field.h"
#include "bgq_field_double.h"
#include "bgq_field_float.h"


void bgq_init_gaugefield_allprec() {
	bgq_init_gaugefield_double();
	bgq_init_gaugefield_float();
}


void bgq_free_gaugefield_allprec() {
	bgq_free_gaugefield_double();
	bgq_free_gaugefield_float();
}

void bgq_init_spinorfields_allprec(int count) {
	bgq_init_spinorfields_double(count);
	bgq_init_spinorfields_float(count);
}


void bgq_free_spinorfields_allprec() {
	bgq_free_spinofields_double();
	bgq_free_spinofields_float();
}


void bgq_hm_init_allprec() {
	bgq_hm_init_double();
	bgq_hm_init_float();
}


void bgq_hm_free_allprec() {
	bgq_hm_free_double();
	bgq_hm_free_float();
}

