/*
 * bgq_field.c
 *
 *  Created on: Aug 15, 2012
 *      Author: meinersbur
 */

#define BGQ_FIELD_C_H_

#include "bgq_field.h"
#include "bgq_field_double.h"
#include "bgq_field_float.h"

#include "../update_backward_gauge.h"

#include <sys/stat.h>

typedef struct {
	unsigned char tv;
	unsigned char x;
	unsigned char y;
	unsigned char z;
} bgq_spinor_coord;

bgq_spinor_coord *g_spinor_body_zline_order;
bgq_spinor_coord *g_spinor_surface_zline_order;


char *(g_idxdesc[BGQREF_count]);
complexdouble *g_bgqvalue = NULL;
complexdouble *g_refvalue = NULL;


void bgq_initbgqref() {
	int datasize = sizeof(complexdouble) * VOLUME * lengthof(g_idxdesc);
	if (g_refvalue == NULL) {
		g_bgqvalue = malloc_aligned(datasize, 128);
		g_refvalue = malloc_aligned(datasize, 128);
	}
	memset(g_bgqvalue, 0xFF, datasize);
	memset(g_refvalue, 0xFF, datasize);

	for (int idx = 0; idx <  lengthof(g_idxdesc); idx+=1) {
		g_idxdesc[idx] = NULL;
	}
}


void bgq_setrefvalue(int t, int x, int y, int z, int idx, complexdouble val, char *desc) {
	if (t < 0)
		t = 0;
	if (t >= LOCAL_LT)
		t = LOCAL_LT-1;
	if (x < 0)
		x = 0;
	if (x >= LOCAL_LX)
		x = LOCAL_LX-1;
	if (y < 0)
		y = 0;
	if (y >= LOCAL_LY)
		y = LOCAL_LY-1;
	if (z < 0)
		z = 0;
	if (z >= LOCAL_LZ)
		z = LOCAL_LZ-1;
	g_refvalue[(((idx*LOCAL_LT + t)*LOCAL_LX + x)*LOCAL_LY + y)*LOCAL_LZ + z] = val;
	g_idxdesc[idx] = desc;
}
void bgq_setbgqvalue(int t, int x, int y, int z, int idx, complexdouble val, char *desc) {
	if (t < 0)
		t = 0;
	if (t >= LOCAL_LT)
		t = LOCAL_LT-1;
	if (x < 0)
		x = 0;
	if (x >= LOCAL_LX)
		x = LOCAL_LX-1;
	if (y < 0)
		y = 0;
	if (y >= LOCAL_LY)
		y = LOCAL_LY-1;
	if (z < 0)
		z = 0;
	if (z >= LOCAL_LZ)
		z = LOCAL_LZ-1;
	g_bgqvalue[(((idx*LOCAL_LT + t)*LOCAL_LX + x)*LOCAL_LY + y)*LOCAL_LZ + z] = val;
	g_idxdesc[idx] = desc;
}

void bgq_savebgqref() {
	if (g_proc_id != 0)
		return;

	int i = 0;
	while (true) {
		char filename[100];
		snprintf(filename, sizeof(filename)-1, "cmp_%d.txt", i);

		struct stat buf;
		if (stat(filename, &buf) != -1) {
			master_print("MK file %s already exists\n", filename);
			i += 1;
			continue;
		}

		// Create marker file
		FILE *file =  fopen(filename, "w");
		fclose(file);

		break;
	}

	for (int idx = 0; idx <  lengthof(g_idxdesc); idx+=1) {
		if (!g_idxdesc[idx])
			continue;

		char reffilename[100];
		snprintf(reffilename, sizeof(reffilename)-1, "cmp_%d_idx%d_%s_ref.txt", i, idx, g_idxdesc[idx]);
		char bgqfilename[100];
		snprintf(bgqfilename, sizeof(bgqfilename)-1, "cmp_%d_idx%d_%s_bgq.txt", i, idx, g_idxdesc[idx]);
		master_print("Cmp Going to write to %s and %s\n", reffilename, bgqfilename);

		FILE *reffile = fopen(reffilename, "w");
		FILE *bgqfile = fopen(bgqfilename, "w");

		fprintf(reffile, "%s\n\n", g_idxdesc[idx]);
		fprintf(bgqfile, "%s\n\n", g_idxdesc[idx]);

		for (int t = 0; t < LOCAL_LT; t += 1) {
			for (int x = 0; x < LOCAL_LX; x += 1) {
				for (int y = 0; y < LOCAL_LY; y += 1) {
					fprintf(reffile, "t=%d x=%d y=%d: ", t,x,y);
					fprintf(bgqfile, "t=%d x=%d y=%d: ", t,x,y);
					for (int z = 0; z < LOCAL_LZ; z += 1) {
						complexdouble refval = g_refvalue[(((idx*LOCAL_LT + t)*LOCAL_LX + x)*LOCAL_LY + y)*LOCAL_LZ + z];
						complexdouble bgqval = g_bgqvalue[(((idx*LOCAL_LT + t)*LOCAL_LX + x)*LOCAL_LY + y)*LOCAL_LZ + z];

						fprintf(reffile, "%8f + %8fi	", creal(refval), cimag(refval));
						fprintf(bgqfile, "%8f + %8fi	", creal(bgqval), cimag(bgqval));
					}
					fprintf(reffile, "\n");
					fprintf(bgqfile, "\n");
				}
			}
		}

		fclose(reffile);
		fclose(bgqfile);

		master_print("Cmp data written to %s and %s\n", reffilename, bgqfilename);
	}
}


void bgq_transfer_spinorfield_allprec(const bool isOdd, int targetindex, spinor * const sourcefield) {
	bgq_transfer_spinorfield_double(isOdd, g_spinorfields_double[targetindex], sourcefield);
	bgq_transfer_spinorfield_float(isOdd, g_spinorfields_float[targetindex], sourcefield);
}


void bgq_init_gaugefield_allprec() {
	bgq_init_gaugefield_double();
	bgq_init_gaugefield_float();
}


void bgq_free_gaugefield_allprec() {
	bgq_free_gaugefield_double();
	bgq_free_gaugefield_float();
}


void bgq_init_spinorfields_allprec(int count, int chi_count) {
	g_spinor_body_zline_order = malloc_aligned(sizeof(bgq_spinor_coord) * BODY_ZLINES, 128);
	g_spinor_surface_zline_order = malloc_aligned(sizeof(bgq_spinor_coord) * SURFACE_ZLINES, 128);

	bgq_init_spinorfields_double(count, chi_count);
	bgq_init_spinorfields_float(count, chi_count);
}


void bgq_free_spinorfields_allprec() {
	bgq_free_spinofields_double();
	bgq_free_spinofields_float();

	free(g_spinor_body_zline_order);
	g_spinor_body_zline_order = NULL;
	free(g_spinor_surface_zline_order);
	g_spinor_surface_zline_order = NULL;
}


//double recvbuf;
//MPI_Request recvrequest;
//double sendbuf;
//MPI_Request sendrequest;

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

