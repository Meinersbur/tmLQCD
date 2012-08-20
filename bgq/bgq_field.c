/*
 * bgq_field.c
 *
 *  Created on: Aug 15, 2012
 *      Author: meinersbur
 */

#include "bgq_field.h"
#include "bgq_field_double.h"
#include "bgq_field_float.h"

#include "../update_backward_gauge.h"


typedef struct {
	unsigned char tv;
	unsigned char x;
	unsigned char y;
	unsigned char z;
} bgq_spinor_coord;

bgq_spinor_coord *g_spinor_body_zline_order;
bgq_spinor_coord *g_spinor_surface_zline_order;


void bgq_init_gaugefield_allprec() {
	bgq_init_gaugefield_double();
	bgq_init_gaugefield_float();
}


void bgq_free_gaugefield_allprec() {
	bgq_free_gaugefield_double();
	bgq_free_gaugefield_float();
}


void bgq_init_spinorfields_allprec(int count) {
	g_spinor_body_zline_order = malloc_aligned(sizeof(bgq_spinor_coord) * BODY_ZLINES, 128);
	g_spinor_surface_zline_order = malloc_aligned(sizeof(bgq_spinor_coord) * SURFACE_ZLINES, 128);

	bgq_init_spinorfields_double(count);
	bgq_init_spinorfields_float(count);
}


void bgq_free_spinorfields_allprec() {
	bgq_free_spinofields_double();
	bgq_free_spinofields_float();

	free(g_spinor_body_zline_order);
	g_spinor_body_zline_order = NULL;
	free(g_spinor_surface_zline_order);
	g_spinor_surface_zline_order = NULL;
}


void bgq_hm_init_allprec() {
	bgq_hm_init_double();
	bgq_hm_init_float();
}


void bgq_hm_free_allprec() {
	bgq_hm_free_double();
	bgq_hm_free_float();
}


void bgq_update_backward_gauge() {
	if (!g_update_gauge_copy)
		return;

	if (g_gaugefield_double)
		bgq_transfer_gaugefield_double(g_gaugefield_double, g_gauge_field);
	if (g_gaugefield_float)
		bgq_transfer_gaugefield_float(g_gaugefield_float, g_gauge_field);
	update_backward_gauge();

	g_update_gauge_copy = 0;
}

