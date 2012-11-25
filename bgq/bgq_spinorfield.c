/*
 * bgq_spinorfield.c
 *
 *  Created on: Oct 26, 2012
 *      Author: meinersbur
 */

#define BGQ_SPINORFIELD_C_
#include "bgq_spinorfield.h"

#include "bgq_HoppingMatrix.h"
#include "bgq_qpx.h"
#include "bgq_dispatch.h"
#include "bgq_comm.h"
#include "bgq_workers.h"

#include "../geometry_eo.h"

#include <mpi.h>
#include <sys/stat.h>
#include <stddef.h>



typedef struct {
	double _Complex c[3];
} su3_vector_array64;

typedef struct {
	su3_vector_array64 v[4];
} spinor_array64;

double bgq_spinorfield_compare(bool isOdd, bgq_weylfield_controlblock *bgqfield, spinor *reffield, bool silent) {
	assert(bgqfield);
	assert(reffield);

	bool readFulllayout = bgqfield->hasFullspinorData;
	bgq_spinorfield_setup(bgqfield, isOdd, readFulllayout, false, !readFulllayout, false, false);
	bgq_master_sync(); // Necessary after bgq_spinorfield_setup if field is accessed without bgq_master_call (which does this implicitely)

	double diff_max = 0;
	size_t count = 0;
	for (size_t z = 0; z < LOCAL_LZ ; z += 1) {
		for (size_t y = 0; y < LOCAL_LY ; y += 1) {
			for (size_t x = 0; x < LOCAL_LX ; x += 1) {
				for (size_t t = 0; t < LOCAL_LT ; t += 1) {
					if (bgq_local2isOdd(t, x, y, z) != isOdd)
						continue;

					const int ix = g_ipt[t][x][y][z]; /* lexic coordinate */
					assert(ix == Index(t,x,y,z));
					int iy = g_lexic2eo[ix]; /* even/odd coordinate (even and odd sites in two different fields of size VOLUME/2, first even field followed by odd) */
					assert(0 <= iy && iy < (VOLUME+RAND));
					int icx = g_lexic2eosub[ix]; /*  even/odd coordinate relative to field base */
					assert(0 <= icx && icx < VOLUME/2);
					assert(icx == iy - (isOdd ? (VOLUME+RAND)/2 : 0));
					spinor_array64 *sp = (spinor_array64*) &reffield[icx];

					size_t ic = bgq_local2collapsed(t, x, y, z);
					size_t k = bgq_local2k(t, x, y, z);
					bgq_spinor bgqspinor = bgq_spinorfield_getspinor(bgqfield, t,x,y,z);

					bool first = true;
					for (size_t v = 0; v < 4; v += 1) {
						for (size_t c = 0; c < 3; c += 1) {
							complexdouble bgqvalue = bgqspinor.v[v].c[c];
							complexdouble refvalue = sp->v[v].c[c];

							double diff = cabs(bgqvalue - refvalue);

							if (diff > 0.01) {
								if (!silent && first)
									master_print("Coordinate (%zu,%zu,%zu,%zu)(%zu,%zu): ref=(%f + %fi) != bgq=(%f + %fi) off by %f\n", t, x, y, z, v, c, creal(refvalue), cimag(refvalue), creal(bgqvalue), cimag(bgqvalue), diff);
								if (first)
									count += 1;
								first = false;
							}
							if (diff > diff_max) {
								diff_max = diff;
							}
						}
					}
				}
			}
		}
	}

	double global_diff_max;
	MPI_Allreduce(&diff_max, &global_diff_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	if (count > 0) {
		if (!silent)
			master_print("%zu sites of %d wrong\n", count, VOLUME/2);
	}

	return global_diff_max;
}


void bgq_weyl_expect(bgq_weyl_nonvec weyl, ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc) {
	size_t t_global = bgq_local2global_t(t);
	size_t x_global = bgq_local2global_x(x);
	size_t y_global = bgq_local2global_y(y);
	size_t z_global = bgq_local2global_z(z);
	if (isSrc) {
		bgq_direction_move_global(&t_global, &x_global, &y_global, &z_global, d);
		d = bgq_direction_revert(d);
	}

	assert(weyl.s[0][0] == t_global);
	assert(weyl.s[0][1] == x_global);
	assert(weyl.s[0][2] == y_global);
	assert(weyl.s[1][0] == z_global);
	assert(weyl.s[1][1] == d);
	assert(weyl.s[1][2] == 0.125);
}


static bgq_weyl_vec_double *bgq_offset2pointer_double(bgq_weylfield_controlblock *field, size_t offset) {
	assert(field);
	assert(offset % sizeof(bgq_weyl_vec_double) == 0);
	assert(bgq_weyl_section_offset(0) <= offset && offset < bgq_weyl_section_offset(sec_end));

	bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
	bgq_weyl_vec_double *baseptr = bgq_section_baseptr_double(field, sec);
	assert(baseptr);
	size_t baseoffset = bgq_weyl_section_offset(sec);
	assert(offset >= baseoffset);
	size_t reloffset = offset - baseoffset;
	return (bgq_weyl_vec_double*)((uint8_t*)baseptr + reloffset);
}


static bgq_weyl_vec_float *bgq_offset2pointer_float(bgq_weylfield_controlblock *field, size_t offset) {
	assert(field);
	assert(offset % sizeof(bgq_weyl_vec_double) == 0);
	assert(bgq_weyl_section_offset(0) <= offset && offset < bgq_weyl_section_offset(sec_end));

	bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
	bgq_weyl_vec_float *baseptr = bgq_section_baseptr_float(field, sec);
	assert(baseptr);
	size_t baseoffset = bgq_weyl_section_offset(sec);
	assert(offset >= baseoffset);
	size_t reloffset = offset - baseoffset;
	return (bgq_weyl_vec_float*)((uint8_t*)baseptr + reloffset/2);
}


static bgq_weyl_vec_double *bgq_index2pointer_double(bgq_weylfield_controlblock *field, ucoord index) {
	return bgq_offset2pointer_double(field, bgq_index2offset(index));
}

static bgq_weyl_vec_float *bgq_index2pointer_float(bgq_weylfield_controlblock *field, ucoord index) {
	return bgq_offset2pointer_float(field, bgq_index2offset(index));
}


static size_t bgq_weylfield_bufferoffset2consecutiveoffset(bool isOdd, size_t offset, size_t k) {
	assert(offset);
	assert(offset % sizeof(bgq_weyl_vec_double) == 0);
	bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
	size_t index = bgq_offset2index(offset);

	//bgq_direction d = bgq_section2direction(sec);
	//assert(d == bgq_offset2ddst(offset));
	bgq_direction d = bgq_offset2ddst(offset);
	ucoord ic = g_bgq_index2collapsed[isOdd][index];
	assert(ic != -1);
	return bgq_collapsed2consecutiveoffset(ic, d);
}


static bgq_spinor_nonvec bgq_spinor_coord_encode(scoord t, scoord x, scoord y, scoord z) {
	ucoord t_global = bgq_local2global_t(t);
	ucoord x_global = bgq_local2global_x(x);
	ucoord y_global = bgq_local2global_y(y);
	ucoord z_global = bgq_local2global_z(z);

	bgq_spinor_nonvec result = {{{{0}}}};
	result.v[0].c[0] = t_global;
	result.v[0].c[1] = x_global;
	result.v[0].c[2] = y_global;
	result.v[1].c[0] = z_global;

	return result;
}


static void bgq_spinorveck_write_double(bgq_spinorsite_double *target, ucoord k, bgq_spinor data) {
	for (ucoord i = 0; i < 4; i+=1) {
		for (ucoord l = 0; l < 3; l+=1) {
			target->s[i][l][k] = data.v[i].c[l];
		}
	}
}


static void bgq_spinorveck_write_float(bgq_spinorsite_float *target, ucoord k, bgq_spinor data) {
	for (ucoord i = 0; i < 4; i+=1) {
		for (ucoord l = 0; l < 3; l+=1) {
			target->s[i][l][k] = data.v[i].c[l];
		}
	}
}


void bgq_spinorveck_written_double(bgq_spinorsite_double *targetspinor, ucoord k, ucoord t, ucoord x, ucoord y, ucoord z) {
	assert(targetspinor->s[0][0][1] == targetspinor->s[0][0][1]); // For valgrind (or to catch NaN)
#ifdef BGQ_COORDCHECK
	bgq_spinor coord = bgq_spinor_coord_encode(t,x,y,z);
	bgq_spinorveck_write_double(targetspinor, k, coord);
#endif
}

void bgq_spinorveck_written_float(bgq_spinorsite_float *targetspinor, ucoord k, ucoord t, ucoord x, ucoord y, ucoord z) {
	if (targetspinor->s[0][0][0] == 0)
		assert(targetspinor->s[0][0][1] != 1); // For valgrind
#ifdef BGQ_COORDCHECK
	bgq_spinor coord = bgq_spinor_coord_encode(t,x,y,z);
	bgq_spinorveck_write_float(targetspinor, k, coord);
#endif
}

static bgq_weyl_nonvec bgq_weyl_coord_encode(ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc) {
	ucoord t_global = bgq_local2global_t(t);
	ucoord x_global = bgq_local2global_x(x);
	ucoord y_global = bgq_local2global_y(y);
	ucoord z_global = bgq_local2global_z(z);
	if (isSrc) {
		bgq_direction_move_global(&t_global, &x_global, &y_global, &z_global, d);
		d = bgq_direction_revert(d);
	}

	bgq_weyl_nonvec result = {{{0}}};
	result.s[0][0] = t_global;
	result.s[0][1] = x_global;
	result.s[0][2] = y_global;
	result.s[1][0] = z_global;
	result.s[1][1] = d;
	result.s[1][2] = 0.125;
	return result;
}

static void bgq_weylveck_write_double(bgq_weyl_vec_double *target, ucoord k, bgq_weyl_nonvec data) {
	for (ucoord i = 0; i < 2; i += 1) {
		for (ucoord l = 0; l < 3; l += 1) {
			target->s[i][l][k] = data.s[i][l];
		}
	}
}


static void bgq_weylveck_write_float(bgq_weyl_vec_float *target, ucoord k, bgq_weyl_nonvec data) {
	for (ucoord i = 0; i < 2; i += 1) {
		for (ucoord l = 0; l < 3; l += 1) {
			target->s[i][l][k] = data.s[i][l];
		}
	}
}

#define bgq_weylveck_written NAME2(bgq_weylveck_written,PRECISION)
void bgq_weylveck_written_double(bgq_weyl_vec_double *targetweyl, ucoord k, ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc) {
#ifdef BGQ_COORDCHECK
	bgq_weyl_nonvec coord = bgq_weyl_coord_encode(t,x,y,z,d,isSrc);
	bgq_weylveck_write_double(targetweyl, k, coord);
#endif
}


void bgq_weylveck_written_float(bgq_weyl_vec_float *targetweyl, ucoord k, ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc) {
#ifdef BGQ_COORDCHECK
	bgq_weyl_nonvec coord = bgq_weyl_coord_encode(t,x,y,z,d,isSrc);
	bgq_weylveck_write_float(targetweyl, k, coord);
#endif
}


static size_t bgq_physical_halo_sites(bgq_dimension dim) {
	switch (dim) {
	case DIM_T:
		return HALO_T * LOCAL_HALO_T/PHYSICAL_LP;
	case DIM_X:
		return PHYSICAL_HALO_X;
	case DIM_Y:
		return PHYSICAL_HALO_Y;
	case DIM_Z:
		return PHYSICAL_HALO_Z;
	}
	UNREACHABLE
	return -1;
}


static size_t bgq_physical_halo_vecs(bgq_dimension dim) {
	switch (dim) {
	case DIM_T:
		return HALO_T * LOCAL_HALO_T/(PHYSICAL_LP*PHYSICAL_LK);
	case DIM_X:
		return PHYSICAL_HALO_X;
	case DIM_Y:
		return PHYSICAL_HALO_Y;
	case DIM_Z:
		return PHYSICAL_HALO_Z;
	}
	UNREACHABLE
	return -1;
}


static bgq_weylfield_section bgq_direction2writesec(bgq_direction d) {
	switch (d) {
	case TUP:
		if (BGQ_UNVECTORIZE || !COMM_T) {
			return sec_temp_tup;
		} else {
			return sec_send_tup;
		}
	case TDOWN:
		if (BGQ_UNVECTORIZE || !COMM_T) {
			return sec_temp_tdown;
		} else {
			return sec_send_tdown;
		}
	case XUP:
		return sec_send_xup;
	case XDOWN:
		return sec_send_xdown;
	case YUP:
		return sec_send_yup;
	case YDOWN:
		return sec_send_ydown;
	case ZUP:
		return sec_send_zup;
	case ZDOWN:
		return sec_send_zdown;
	}
	UNREACHABLE
	return -1;
}


void bgq_spinorfield_setup(bgq_weylfield_controlblock *field, bool isOdd, bool readFullspinor, bool writeFullspinor, bool readWeyl, bool writeWeyl, bool writeFloat) {
	assert(field);
	assert(readFullspinor || writeFullspinor || readWeyl || writeWeyl);
	// Do something

	bool fullspinorAvailable;
	bool weylAvailable;
	if (field->isInitialized) {
		fullspinorAvailable = field->hasFullspinorData;
		weylAvailable = field->hasWeylfieldData || field->waitingForRecv;
	} else {
		fullspinorAvailable = false;
		weylAvailable = false;
		field->sec_weyl = NULL;
		field->sec_collapsed_double = NULL;
		field->sec_collapsed_float = NULL;
		field->sec_fullspinor_double = NULL;
		field->sec_fullspinor_float = NULL;
		field->hasFullspinorData = false;
		field->hasWeylfieldData = false;
		field->waitingForRecv = false;
		field->isOdd = isOdd;
		field->isFulllayoutSloppy = false;
		field->isWeyllayoutSloppy = false;
		field->pendingDatamove = false;

		field->isInitialized = true;
	}

	// Possible actions
	bool actionAllocFulllayout= false;
	bool actionInitFulllayout = false;
	bool actionAllocWeyllayout = false;
	bool actionInitWeylPtrs = false;
	bool actionWaitForRecv = false;
	bool actionDatamove = false;

	if (readFullspinor) {
		if (!fullspinorAvailable)
			master_error(1, "No fullspinor data written here\n");
		//TODO: We may implement a field data translation if the calling func does not support this layout
	}
	if (writeFullspinor) {
		if (field->sec_fullspinor_double==NULL) {
			actionAllocFulllayout = true;
			actionInitFulllayout = true;
		} else if (field->isOdd != isOdd) {
			actionInitFulllayout = true;
		}
	}
	if (readWeyl) {
		if (!weylAvailable)
			master_error(1, "No weyl data written here\n");
		//TODO: We may implement a field data translation if the calling func does not support this layout

		if (field->pendingDatamove) {
			actionDatamove = true;
		}
		if (field->waitingForRecv) { //TODO: make this a global property, there can be at most one communication going on!
			actionWaitForRecv = true;
		}
	}
	if (writeWeyl) {
		if (field->waitingForRecv) {
			// The means we are overwriting data that has never been read
			//master_print("PERFORMANCE WARNING: Overwriting data not yet received from remote node, i.e. has never been used yet\n");
			actionWaitForRecv = true;
		}

		if (field->sec_weyl==NULL) {
			actionAllocWeyllayout = true;
			actionInitWeylPtrs = true;
		} else if (field->isOdd != isOdd) {
			master_print("PERFORMANCE WARNING: Performance loss by reuse of spinorfield with different oddness\n");
			actionInitWeylPtrs = true;
			field->isOdd = isOdd;
		}
	}
	assert(!field->waitingForRecv || !(actionAllocWeyllayout || actionInitWeylPtrs)); // Do not change field while we are receiving


	if (actionWaitForRecv) {
		bool nospi = field->hmflags & hm_nospi;

		// 4. Wait for the communication to finish
		bgq_comm_wait(nospi, field->isWeyllayoutSloppy);
		field->waitingForRecv = false;
	}
	if (actionDatamove) {
		// 5. Move received to correct location
		bgq_master_sync();
		static bgq_work_datamove work;
		work.spinorfield = field;
		work.opts = field->hmflags;
		if (field->isWeyllayoutSloppy)
			bgq_master_call(&bgq_HoppingMatrix_worker_datamove_float, &work);
		else
			bgq_master_call(&bgq_HoppingMatrix_worker_datamove_double, &work);
		field->pendingDatamove = false;
	}


	if (actionAllocFulllayout) {
		size_t fullfieldsize = PHYSICAL_VOLUME * sizeof(bgq_spinorsite_double);
		void *spinorsite = malloc_aligned(fullfieldsize, BGQ_ALIGNMENT_L2);

		field->sec_fullspinor_double = spinorsite;
		field->sec_fullspinor_float = spinorsite;
		//field->sec_fullspinor_surface = spinorsite;
		spinorsite = (bgq_spinorsite_double*)spinorsite + PHYSICAL_SURFACE;

		//field->sec_fullspinor_body = spinorsite;
	}
	if (actionInitFulllayout) {
		//bgq_spinorfield_fulllayout_clear(field);
	}


	if (actionAllocWeyllayout) {
		size_t weylfieldsize = PHYSICAL_VOLUME * sizeof(bgq_weylsite_double);
		uint8_t *weylbase = ((uint8_t*)malloc_aligned(weylfieldsize, BGQ_ALIGNMENT_L2));

		field->sec_weyl = weylbase;
		field->sec_collapsed_double = (bgq_weylsite_double*) (weylbase + bgq_weyl_section_offset(sec_collapsed));
		field->sec_collapsed_float = (bgq_weylsite_float*)field->sec_collapsed_double;
		//field->sec_surface = (bgq_weylsite*) (weylbase + bgq_weyl_section_offset(sec_surface));
		//field->sec_body = (bgq_weylsite*) (weylbase + bgq_weyl_section_offset(sec_body));
		field->sec_end = weylbase + bgq_weyl_section_offset(sec_end);

		field->sendptr_double = malloc_aligned(PHYSICAL_VOLUME * sizeof(*field->sendptr_double), BGQ_ALIGNMENT_L2);
		field->sendptr_float = malloc_aligned(PHYSICAL_VOLUME * sizeof(*field->sendptr_float), BGQ_ALIGNMENT_L2);

		for (size_t d = 0; d < PHYSICAL_LD; d+=1) {
			bgq_dimension dim = bgq_direction2dimension(d);
			field->consptr_double[d] = malloc_aligned(bgq_physical_halo_sites(dim) * sizeof(*field->consptr_double[d]), BGQ_ALIGNMENT_L2);
			field->consptr_float[d] = malloc_aligned(bgq_physical_halo_sites(dim) * sizeof(*field->consptr_float[d]), BGQ_ALIGNMENT_L2);
		}
	}

	if (actionInitWeylPtrs) {
		uint8_t *weylbase = field->sec_weyl;
		assert(weylbase);

		// For 1st phase (distribute)
		for (ucoord ic_src = 0; ic_src < PHYSICAL_VOLUME; ic_src += 1) {
			for (size_t d_dst = 0; d_dst < PHYSICAL_LD; d_dst += 1) {
				ucoord index = g_bgq_collapsed2indexsend[isOdd/*_dst*/][ic_src].d[d_dst];

				{
					bgq_weyl_vec_double *ptr = bgq_index2pointer_double(field, index);
					field->isWeyllayoutSloppy = false;
					assert(bgq_pointer2offset(field, ptr) == bgq_index2offset(index));
					field->sendptr_double[ic_src].d[d_dst] = ptr;
				}

				{
					bgq_weyl_vec_float *ptr = bgq_index2pointer_float(field, index);
					field->isWeyllayoutSloppy = true;
					assert(bgq_pointer2offset(field, ptr) == bgq_index2offset(index));
					field->sendptr_float[ic_src].d[d_dst] = ptr;
				}
			}
		}


		// For 5th phase (datamove)
		for (ucoord d_dst = 0; d_dst < PHYSICAL_LD; d_dst+=1) {
			bgq_direction d_src = bgq_direction_revert(d_dst);
			bgq_dimension dim = bgq_direction2dimension(d_dst);
			ucoord sites = bgq_physical_halo_sites(dim);
			for (ucoord j = 0; j < sites; j+=1) {
				bgq_weylfield_section sec = bgq_direction2writesec(d_src);
				size_t baseoffset = bgq_weyl_section_offset(sec);
				ucoord baseindex = bgq_offset2index(baseoffset);
				ucoord index = baseindex + j;
				ucoord ic_dst = g_bgq_index2collapsed[isOdd][index]; // What did the previous step write here?
				size_t offset_cons = bgq_collapsed2consecutiveoffset(ic_dst, d_dst);
				{
					bgq_weyl_vec_double *ptr = bgq_offset2pointer_double(field, offset_cons);
					field->consptr_double[d_dst][j] = ptr;
				}

				{
					bgq_weyl_vec_float *ptr = bgq_offset2pointer_float(field, offset_cons);
					field->consptr_float[d_dst][j] = ptr;
				}
			}
		}
	}


	if (writeWeyl || writeFullspinor) {
		// Data is going to be written, so what actually is stored here changes
		field->hasWeylfieldData = writeWeyl;
		field->isWeyllayoutSloppy = writeFloat;
		field->hasFullspinorData = writeFullspinor;
		field->isFulllayoutSloppy = writeFloat;
	}
}


void bgq_spinorfield_setup_float(bgq_weylfield_controlblock *field, bool isOdd, bool readFullspinor, bool writeFullspinor, bool readWeyl, bool writeWeyl) {
	bgq_spinorfield_setup(field, isOdd, readFullspinor, writeFullspinor, readWeyl, writeWeyl, true);
}


void bgq_spinorfield_setup_double(bgq_weylfield_controlblock *field, bool isOdd, bool readFullspinor, bool writeFullspinor, bool readWeyl, bool writeWeyl) {
	bgq_spinorfield_setup(field, isOdd, readFullspinor, writeFullspinor, readWeyl, writeWeyl, true);
}


typedef struct {
	bgq_weylfield_controlblock *field;
} bgq_conversion_args;


void bgq_spinorfield_transfer(bool isOdd, bgq_weylfield_controlblock *targetfield, spinor *sourcefield) {
	bgq_spinorfield_setup(targetfield, isOdd, false, true, false, false, false);
	size_t ioff = isOdd ? (VOLUME+RAND)/2 : 0;

	for (size_t i_eosub = 0; i_eosub < VOLUME/2; i_eosub+=1) {
		size_t i_eo = i_eosub + ioff;
		size_t i_lexic = g_eo2lexic[i_eo];
		int t = g_coord[i_lexic][0] - g_proc_coords[0]*T;
		int x = g_coord[i_lexic][1] - g_proc_coords[1]*LX;
		int y = g_coord[i_lexic][2] - g_proc_coords[2]*LY;
		int z = g_coord[i_lexic][3] - g_proc_coords[3]*LZ;
		spinor_array64 *sourcespinor = (spinor_array64*)&sourcefield[i_eosub];
		assert(bgq_local2isOdd(t,x,y,z)==isOdd);

		ucoord ih = bgq_local2halfvolume(t,x,y,z);
		ucoord ic = bgq_local2collapsed(t,x,y,z);
		ucoord k = bgq_local2k(t,x,y,z);
		bgq_spinorsite_double *targetspinor = &targetfield->sec_fullspinor_double[ic];

		for (int v = 0; v < 4; v+=1) {
			for (int c = 0; c < 3; c+=1) {
				targetspinor->s[v][c][k] = sourcespinor->v[v].c[c];
			}
		}
		bgq_spinorveck_written_double(targetspinor, k, t, x, y, z);
	}
}


bgq_spinor bgq_legacy_getspinor(spinor *spinor, ucoord t, ucoord x, ucoord y, ucoord z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);

	int ix = g_ipt[t][x][y][z]; /* lexic coordinate */
	assert(ix == Index(t,x,y,z));
	int iy = g_lexic2eo[ix]; /* even/odd coordinate (even and odd sites in two different fields of size VOLUME/2, first even field followed by odd) */
	assert(0 <= iy && iy < (VOLUME+RAND));
	int icx = g_lexic2eosub[ix]; /*  even/odd coordinate relative to field base */
	assert(0 <= icx && icx < VOLUME/2);

	bgq_spinor *result = (bgq_spinor*)&spinor[icx];
	return *result;
}


bgq_spinor bgq_spinorfield_getspinor(bgq_weylfield_controlblock *field, ucoord t, ucoord x, ucoord y, ucoord z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);

	assert(bgq_local2isOdd(t, x, y, z) == field->isOdd);
	ucoord ic = bgq_local2collapsed(t, x, y, z);
	ucoord k = bgq_local2k(t, x, y, z);

	bgq_spinorfield_layout layout = bgq_spinorfield_prepareRead(field, field->isOdd, true, true, true, false);
	//bgq_spinorfield_setup(field, field->isOdd, !(layout & ly_weyl), false, (layout & ly_weyl), false, false);
	bgq_su3_spinor_decl(spinor);
	bgq_spinorfield_readSpinor(&spinor, field, ic, layout&ly_weyl, layout&ly_sloppy, layout&ly_mul);
	return bgq_spinor_fromqpx(spinor, k);

#if 0
	if (field->hasFullspinorData) {
		bgq_spinorfield_setup(field, field->isOdd, true, false, false, false, false);
		bgq_spinorsite spinor = field->sec_fullspinor[ic];
		return bgq_spinor_fromvec(spinor,k);
	} else if (field->hasWeylfieldData) {
		bgq_spinorfield_setup(field, field->isOdd, false, false, true, false, false);
		bgq_su3_spinor_decl(spinor);
		size_t offset = bgq_pointer2offset(field, &field->sec_collapsed[ic].d[XUP]);
		ucoord index = bgq_offset2index(offset);
		bgq_HoppingMatrix_loadWeyllayout(spinor,&field->sec_collapsed[ic], bgq_t2t(t,0), bgq_t2t(t,1), x, y, z);
		return bgq_spinor_fromqpx(spinor,k);
	} else {
		printf("Field contains no data\n");
		abort();
		bgq_spinor dummy;
		return dummy;
	}
#endif
}



char *(g_idxdesc[BGQREF_count]);
complexdouble *g_bgqvalue = NULL;
complexdouble *g_refvalue = NULL;


void bgq_initbgqref_impl() {
	size_t datasize = sizeof(complexdouble) * VOLUME * lengthof(g_idxdesc);
	if (g_refvalue == NULL) {
		g_bgqvalue = malloc_aligned(datasize, BGQ_ALIGNMENT_L2);
		g_refvalue = malloc_aligned(datasize, BGQ_ALIGNMENT_L2);
	}
	memset(g_bgqvalue, 0xFF, datasize);
	memset(g_refvalue, 0xFF, datasize);

	for (int idx = 0; idx <  lengthof(g_idxdesc); idx+=1) {
		g_idxdesc[idx] = NULL;
	}
}


void bgq_setdesc_impl(int idx, char *desc){
	g_idxdesc[idx] = desc;
}


void bgq_setrefvalue_impl(int t, int x, int y, int z, bgqref idx, complexdouble val) {
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
	if (!g_idxdesc[idx])
		g_idxdesc[idx] = "";
}


void bgq_setbgqvalue_impl(int t, int x, int y, int z, bgqref idx, complex_double val) {
	if (idx==BGQREF_TUP && t==0 && x==0 && y==0 && z==0) {
		int a = 0;
	}
	if (val == 0) {
		if (val == 1) {
			//assert(false);
		}
	}

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
	if (!g_idxdesc[idx])
		g_idxdesc[idx] = "";
}


void bgq_setbgqvalue_src_impl(ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bgqref idx, complexdouble val) {
	bgq_direction_move_local(&t, &x, &y, &z, d);
	bgq_setbgqvalue(t, x, y, z, idx, val);
}


void bgq_savebgqref_impl() {
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
		char refdir[100];
		snprintf(refdir, sizeof(refdir)-1, "cmpref");
		snprintf(reffilename, sizeof(reffilename)-1, "%s/cmp_idx%02d_%s.txt", refdir, idx, g_idxdesc[idx]);

		char bgqdir[100];
		char bgqfilename[100];
		snprintf(bgqdir, sizeof(bgqdir)-1, "cmpbgq");
		snprintf(bgqfilename, sizeof(bgqfilename)-1, "%s/cmp_idx%02d_%s.txt", bgqdir, idx, g_idxdesc[idx]);
		master_print("Cmp Going to write to %s and %s\n", reffilename, bgqfilename);

		mkdir(refdir, S_IRWXU | S_IRWXG | S_IRWXO);
		mkdir(bgqdir, S_IRWXU | S_IRWXG | S_IRWXO);

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

	bgq_initbgqref();
}


size_t bgq_fieldpointer2offset(void *ptr) {
	for (size_t i = 0; i < g_bgq_spinorfields_count; i+=1) {
		bgq_weylfield_controlblock *field = &g_bgq_spinorfields[i];
		if (!field->isInitialized)
			continue;
		if (!field->sec_weyl)
			continue;

		size_t result = bgq_pointer2offset_raw(field, ptr, false);
		if (result!=-1)
			return result;
	}

	assert(!"Pointer does not point to weyllayout");
	return -1;
}


bgq_weyl_vec_double *bgq_section_baseptr_double(bgq_weylfield_controlblock *field, bgq_weylfield_section section) {
	assert(field);
	assert(field->isInitialized);
	bgq_weyl_vec_double *result = NULL;

	switch (section) {
	case sec_surface:
		if (PHYSICAL_SURFACE==0)
			return NULL;
		result = (bgq_weyl_vec_double*)&field->sec_collapsed_double[bgq_surface2collapsed(0)];
		break;
	case sec_body:
		if (PHYSICAL_BODY==0)
			return NULL;
		result = (bgq_weyl_vec_double*)&field->sec_collapsed_double[bgq_body2collapsed(0)];
		break;
	case sec_send_tup:
		result = g_bgq_sec_send_double[TUP];
		break;
	case sec_send_tdown:
		result = g_bgq_sec_send_double[TDOWN];
		break;
	case sec_send_xup:
		result = g_bgq_sec_send_double[XUP];
		break;
	case sec_send_xdown:
		result = g_bgq_sec_send_double[XDOWN];
		break;
	case sec_send_yup:
		result = g_bgq_sec_send_double[YUP];
		break;
	case sec_send_ydown:
		result = g_bgq_sec_send_double[YDOWN];
		break;
	case sec_send_zup:
		result = g_bgq_sec_send_double[ZUP];
		break;
	case sec_send_zdown:
		result = g_bgq_sec_send_double[ZDOWN];
		break;
	case sec_recv_tup:
		result = g_bgq_sec_recv_double[TUP];
		break;
	case sec_recv_tdown:
		result = g_bgq_sec_recv_double[TDOWN];
		break;
	case sec_recv_xup:
		result = g_bgq_sec_recv_double[XUP];
		break;
	case sec_recv_xdown:
		result = g_bgq_sec_recv_double[XDOWN];
		break;
	case sec_recv_yup:
		result = g_bgq_sec_recv_double[YUP];
		break;
	case sec_recv_ydown:
		result = g_bgq_sec_recv_double[YDOWN];
		break;
	case sec_recv_zup:
		result = g_bgq_sec_recv_double[ZUP];
		break;
	case sec_recv_zdown:
		result = g_bgq_sec_recv_double[ZDOWN];
		break;
	case sec_temp_tup:
		result = g_bgq_sec_temp_tup_double;
		break;
	case sec_temp_tdown:
		result = g_bgq_sec_temp_tdown_double;
		break;
	case sec_vrecv_tdown:
	case sec_vrecv_tup:
		return NULL;
	default:
		assert(!"Section does not physically exist");
		UNREACHABLE
	}

	return result;
}


bgq_weyl_vec_float *bgq_section_baseptr_float(bgq_weylfield_controlblock *field, bgq_weylfield_section section) {
	assert(field);
	assert(field->isInitialized);
	bgq_weyl_vec_float *result = NULL;

	switch (section) {
	case sec_surface:
		if (PHYSICAL_SURFACE==0)
			return NULL;
		result = (bgq_weyl_vec_float*)&field->sec_collapsed_float[bgq_surface2collapsed(0)];
		break;
	case sec_body:
		if (PHYSICAL_BODY==0)
			return NULL;
		result = (bgq_weyl_vec_float*)&field->sec_collapsed_float[bgq_body2collapsed(0)];
		break;
	case sec_send_tup:
		result = g_bgq_sec_send_float[TUP];
		break;
	case sec_send_tdown:
		result = g_bgq_sec_send_float[TDOWN];
		break;
	case sec_send_xup:
		result = g_bgq_sec_send_float[XUP];
		break;
	case sec_send_xdown:
		result = g_bgq_sec_send_float[XDOWN];
		break;
	case sec_send_yup:
		result = g_bgq_sec_send_float[YUP];
		break;
	case sec_send_ydown:
		result = g_bgq_sec_send_float[YDOWN];
		break;
	case sec_send_zup:
		result = g_bgq_sec_send_float[ZUP];
		break;
	case sec_send_zdown:
		result = g_bgq_sec_send_float[ZDOWN];
		break;
	case sec_recv_tup:
		result = g_bgq_sec_recv_float[TUP];
		break;
	case sec_recv_tdown:
		result = g_bgq_sec_recv_float[TDOWN];
		break;
	case sec_recv_xup:
		result = g_bgq_sec_recv_float[XUP];
		break;
	case sec_recv_xdown:
		result = g_bgq_sec_recv_float[XDOWN];
		break;
	case sec_recv_yup:
		result = g_bgq_sec_recv_float[YUP];
		break;
	case sec_recv_ydown:
		result = g_bgq_sec_recv_float[YDOWN];
		break;
	case sec_recv_zup:
		result = g_bgq_sec_recv_float[ZUP];
		break;
	case sec_recv_zdown:
		result = g_bgq_sec_recv_float[ZDOWN];
		break;
	case sec_temp_tup:
		result = g_bgq_sec_temp_tup_float;
		break;
	case sec_temp_tdown:
		result = g_bgq_sec_temp_tdown_float;
		break;
	case sec_vrecv_tdown:
	case sec_vrecv_tup:
		return NULL;
	default:
		assert(!"Section does not physically exist");
		UNREACHABLE
	}

	return result;
}


void bgq_spinorfield_prepareWrite(bgq_weylfield_controlblock *field, bool isOdd, bgq_spinorfield_layout layout) {
	assert(field);

	bgq_spinorfield_setup(field, isOdd, false, !(layout & ly_weyl), false, (layout & ly_weyl), layout & ly_sloppy);
}



typedef struct {
	bgq_weylfield_controlblock *field;
	bool sloppy;
	bgq_spinorfield_layout layout;
} bgq_spinorfield_rewrite_work;

static void bgq_spinorfield_rewrite_worker(void *arg_untyped, size_t tid, size_t threads) {
	bgq_spinorfield_rewrite_work *arg = arg_untyped;
	bgq_weylfield_controlblock *field  = arg->field;
	bool sloppy = arg->sloppy;
	bgq_spinorfield_layout layout = arg->layout;

	const size_t workload = PHYSICAL_VOLUME;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);
	for (ucoord ic = begin; ic < end; ic+=1) {
#ifndef NDEBUG
		ucoord ih = bgq_collapsed2halfvolume(field->isOdd, ic);
		ucoord t1 = bgq_halfvolume2t1(field->isOdd, ih);
		ucoord t2 = bgq_halfvolume2t2(field->isOdd, ih);
		ucoord tv = bgq_halfvolume2tv(ih);
		ucoord x = bgq_halfvolume2x(ih);
		ucoord y = bgq_halfvolume2y(ih);
		ucoord z = bgq_halfvolume2z(ih);
#endif

		bgq_su3_spinor_decl(spinor);
		bgq_spinorfield_readSpinor(&spinor, field, ic, layout&ly_weyl, layout&ly_sloppy, layout&ly_weyl);

		if (sloppy) {
			bgq_spinor_vec_float *fulladdr = &field->sec_fullspinor_float[ic];
			bgq_su3_spinor_store_float(fulladdr, spinor);
					bgq_spinorvec_written_float(fulladdr, t1, t2, x, y, z);
		} else {
			bgq_spinor_vec_double *fulladdr = &field->sec_fullspinor_double[ic];
			bgq_su3_spinor_store_double(fulladdr, spinor);
					bgq_spinorvec_written_double(fulladdr, t1, t2, x, y, z);
		}
	}
}


bgq_spinorfield_layout bgq_spinorfield_prepareRead(bgq_weylfield_controlblock *field, bool isOdd, bool acceptWeyl, bool acceptDouble, bool acceptSloppy, bool acceptMul) {
	assert(field);
	assert(field->isInitialized);
	assert(field->hasWeylfieldData || field->hasFullspinorData);
	assert(field->isOdd == isOdd);

	// to wait for comm if necessary
	if (field->waitingForRecv)
		bgq_spinorfield_setup(field, isOdd, false, false, true, false, false);

	bool actionRewrite = false;
	if (field->hasFullspinorData) {
		if (field->isFulllayoutSloppy) {
			if (!acceptSloppy)
				actionRewrite = true;
		} else {
			if (!acceptDouble)
				actionRewrite = true;
		}
	} else if (field->hasWeylfieldData) {
		if (!acceptWeyl)
			actionRewrite = true;
		if (field->isWeyllayoutSloppy) {
			if (!acceptSloppy)
				actionRewrite = true;
		} else {
			if (!acceptDouble)
				actionRewrite = true;
		}
	} else {
		UNREACHABLE
	}


	bgq_spinorfield_layout result = -1;

	//TODO: This is not meant to be fast; If you need something fast, special-case it (i.e. accept more inputs)
	if (actionRewrite) {
		bgq_spinorfield_layout layout = bgq_spinorfield_bestLayout(field);
		bgq_spinorfield_setup(field, isOdd, false, true, false, false, !acceptDouble);

		if (!(layout&ly_weyl) && (void*)field->sec_fullspinor_float==(void*)field->sec_fullspinor_double) {
			// This is bad: we are going to overwrite the data we need to read
			// Solution: alloc a new memory area
			// Not that only half of the mem allocated for sec_fullspinor_float will stay in use
			//TODO: find something better, like doing it in reverse
			if (acceptDouble)
				field->sec_fullspinor_double = malloc_aligned(PHYSICAL_VOLUME * sizeof(*field->sec_fullspinor_double), BGQ_ALIGNMENT_L2);
			else
				field->sec_fullspinor_float = malloc_aligned(PHYSICAL_VOLUME * sizeof(*field->sec_fullspinor_float), BGQ_ALIGNMENT_L2);
		}

		bgq_master_sync();
		static bgq_spinorfield_rewrite_work work;
		work.field = field;
		work.sloppy = !acceptDouble;
		work.layout = layout;
		bgq_master_call(&bgq_spinorfield_rewrite_worker, &work);

		field->hasFullspinorData = true;
		field->isFulllayoutSloppy = !acceptDouble;

		result = acceptDouble ? ly_full_double : ly_full_float;
	} else {
		result = bgq_spinorfield_bestLayout(field);
	}

	assert(result!=-1);
	assert(acceptWeyl || !(result&ly_weyl));
	assert(acceptDouble || result&ly_sloppy);
	assert(acceptSloppy || !(result&ly_sloppy));
	assert(acceptMul || !(result&ly_mul));

	return result;
}
