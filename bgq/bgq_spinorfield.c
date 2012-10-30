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

#include "../geometry_eo.h"

#include <mpi.h>




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
	bgq_spinorfield_setup(bgqfield, isOdd, readFulllayout, false, !readFulllayout, false);

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

					for (size_t v = 0; v < 4; v += 1) {
						for (size_t c = 0; c < 3; c += 1) {
							complexdouble bgqvalue = bgqspinor.v[v].c[c];
							complexdouble refvalue = sp->v[v].c[c];

							double diff = cabs(bgqvalue - refvalue);

							if (diff > 0.01) {
								if (!silent)
									master_print("Coordinate (%zu,%zu,%zu,%zu)(%zu,%zu): ref=(%f + %fi) != bgb=(%f + %fi) off by %f\n", t, x, y, z, v, c, creal(refvalue), cimag(refvalue), creal(bgqvalue), cimag(bgqvalue), diff);
								count += 1;
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
	MPI_Allreduce(&diff_max, &global_diff_max, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	return diff_max;

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
	assert(weyl.s[1][2] == 0);
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
	return result;
}


static void bgq_weylveck_write(bgq_weyl_vec *target, ucoord k, bgq_weyl_nonvec data) {
	for (ucoord i = 0; i < 2; i += 1) {
		for (ucoord l = 0; l < 3; l += 1) {
			target->s[i][l][k] = data.s[i][l];
		}
	}
}


static void bgq_spinorveck_write(bgq_spinorsite *target, ucoord k, bgq_spinor data) {
	for (ucoord i = 0; i < 4; i+=1) {
		for (ucoord l = 0; l < 3; l+=1) {
			target->s[i][l][k] = data.v[i].c[l];
		}
	}
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


static bgq_weyl_vec bgq_weyl_mergevec(bgq_weyl_nonvec s_left, bgq_weyl_nonvec s_right) {
	bgq_weyl_vec result;
	for (size_t v = 0; v < 2; v+=1) {
		for (size_t c = 0; c < 3; c+=1) {
			result.s[v][c][0] = s_left.s[v][c];
			result.s[v][c][1] = s_right.s[v][c];
		}
	}
	return result;
}


static bgq_spinor_vec bgq_spinor_mergevec(bgq_spinor_nonvec s_left, bgq_spinor_nonvec s_right) {
	bgq_spinor_vec result;
	for (ucoord v = 0; v < 4; v+=1) {
		for (ucoord c = 0; c < 3; c+=1) {
			result.s[v][c][0] = s_left.v[v].c[c];
			result.s[v][c][1] = s_right.v[v].c[c];
		}
	}
	return result;
}




void bgq_HoppingMatrix_worker_datamove(void *argptr, size_t tid, size_t threads) {
	const bgq_weylfield_controlblock *spinorfield = argptr;

	const size_t workload_recvt = 2*(COMM_T ? 2*LOCAL_HALO_T/PHYSICAL_LP : LOCAL_HALO_T/PHYSICAL_LP);
	const size_t workload_recv = 2*PHYSICAL_HALO_X + 2*PHYSICAL_HALO_Y + 2*PHYSICAL_HALO_Z;
	const size_t workload = workload_recvt + workload_recv;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min(workload, begin+threadload);
	for (size_t i = begin; i<end; ) {
		WORKLOAD_DECL(i,workload);

		if (WORKLOAD_SPLIT(workload_recvt)) {
			// Do T-dimension
			(void)WORKLOAD_PARAM(2); // Count an T-iteration twice; better have few underloaded threads (so remaining SMT-threads have some more ressources) then few overloaded threads (so the master thread has to wait for them)

#if COMM_T
			size_t beginj = WORKLOAD_PARAM(2*LOCAL_HALO_T/PHYSICAL_LP);
			size_t endj = min(2*LOCAL_HALO_T/PHYSICAL_LP,beginj+threadload/2);
			for (size_t j = beginj; j < endj; j+=1) {
				//TODO: Check strength reduction
				//TODO: Prefetch
				//TODO: Inline assembler
				bgq_weyl_vec *weyladdr_left = &spinorfield->sec_recv[TDOWN][j]; // Note: Overlaps into sec_send_tdown
				bgq_weyl_vec *weyladdr_right = &spinorfield->sec_send[TUP][j]; // Note: Overlaps into sec_recv_tup
				bgq_weyl_vec *weyladdr_dst = spinorfield->destptrFromTRecv[j];

				bgq_su3_weyl_decl(weyl_left);
				bgq_su3_weyl_load_left_double(weyl_left, weyladdr_left);
				bgq_su3_weyl_decl(weyl_right);
				bgq_su3_weyl_load_right_double(weyl_right, weyladdr_right);

				bgq_su3_weyl_decl(weyl);
				bgq_su3_weyl_merge2(weyl, weyl_left, weyl_right);

				bgq_su3_weyl_store_double(weyladdr_dst, weyl);
			}
			i += 2*(endj - beginj);
#else
			assert(!"yet implemented");
#endif

		} else {
			// Do other dimensions

			size_t beginj = WORKLOAD_PARAM(workload_recv);
			size_t endj = min(workload_recv,beginj+threadload);
			for (size_t j = beginj; j < endj; j+=1) {
				//TODO: Check strength reduction
				//TODO: Prefetch
				//TODO: Inline assembler
				bgq_weyl_vec *weyladdr_src= &spinorfield->sec_recv[XUP][j]; // Note: overlaps into following sections
				bgq_weyl_vec *weyladdr_dst = spinorfield->destptrFromRecv[j];
				assert((bgq_weyl_vec *)spinorfield->sec_weyl <= weyladdr_dst && weyladdr_dst <  (bgq_weyl_vec *)spinorfield->sec_end);

				bgq_su3_weyl_decl(weyl);
				bgq_su3_weyl_load_double(weyl, weyladdr_src);
				bgq_su3_weyl_store_double(weyladdr_dst, weyl);
			}
			i += (endj - beginj);
		}
		WORKLOAD_CHECK
	}
}


static bgq_weyl_vec *bgq_offset2pointer(uint8_t *weylbase, size_t offset) {
	assert(weylbase);
	//assert(offset);
	assert(offset % sizeof(bgq_weyl_vec) == 0);
	assert(bgq_weyl_section_offset(0) <= offset && offset < bgq_weyl_section_offset(sec_end));

	if (offset >= bgq_weyl_section_offset(sec_comm)) {
		// Communication buffers
		size_t offset_comm = offset - bgq_weyl_section_offset(sec_comm);
		return (bgq_weyl_vec *)(g_bgq_sec_comm + offset_comm);
	} else {
		// Main field
		return (bgq_weyl_vec *)(weylbase + offset);
	}
}


static bgq_weyl_vec *bgq_index2pointer(void *weylbase, ucoord index) {
	return bgq_offset2pointer(weylbase, bgq_index2offset(index));
}


static bgq_weyl_vec *bgq_encodedoffset2pointer(uint8_t *weylbase, size_t code) {
	size_t offset = bgq_decode_offset(code);
	assert(offset % sizeof(bgq_weyl_vec) == 0);
	assert(bgq_weyl_section_offset(0) <= offset && offset < bgq_weyl_section_offset(sec_end));
	return bgq_offset2pointer(weylbase,offset);
}








static size_t bgq_weylfield_bufferoffset2consecutiveoffset(bool isOdd, size_t offset, size_t k) {
	assert(offset);
	assert(offset % sizeof(bgq_weyl_vec)==0);
	bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
	size_t index = offset/sizeof(bgq_weyl_vec);
	assert(index);
	bool isOdd_src = isOdd;
	bool isOdd_dst = !isOdd;

	switch (sec) {
	case sec_send_tup:
	case sec_send_tdown: {
		// These are interesting cases: since there are meant to be sent, there is no consecutive offset for them
		// However, in T-direction only one k is useful on the other side, the other is reused locally
		// If COMM_T==0, both are reused locally, but at different locations

		assert((sec==sec_send_tup && k!=1) || (sec==sec_send_tdown && k!=0)); // Values here are meant to be sent to another node
		bgq_direction d_src = bgq_section2direction(sec); // Data is sent into this direction
		bgq_direction d_dst = bgq_direction_revert(d_src); // Therefore, from the receiver's perspective, it comes from this direction
		assert(d_src == bgq_offset2dsrc(offset));

		// Who wrote them here?
		size_t ih_src = g_bgq_index2ih_src[isOdd][index];
		assert(d_src == g_bgq_index2d_src[isOdd][index]);
		size_t tv_src = bgq_halfvolume2tv(ih_src);

		size_t t_src = bgq_halfvolume2t(isOdd_src, ih_src, k);
		size_t x_src = bgq_halfvolume2x(ih_src);
		size_t y_src = bgq_halfvolume2y(ih_src);
		size_t z_src = bgq_halfvolume2z(ih_src);
		assert(bgq_halfvolume2isSurface(isOdd_src, ih_src));

		// To whom is it intended?
		size_t ih_dst = bgq_src2ih_dst(t_src, x_src, y_src, z_src, d_src);
		size_t k_dst = bgq_src2k_dst(t_src, x_src, y_src, z_src, d_src);
		assert((sec==sec_send_tup && k_dst==1) || (sec==sec_send_tdown && k_dst==0));

		size_t result = bgq_weyllayout_halfvolume2consecutiveoffset(isOdd_dst, ih_dst, d_dst);
		assert(d_src == bgq_offset2dsrc(result));
		return result;
	}

	case sec_recv_tup:
	case sec_recv_tdown:
	case sec_recv_xup:
	case sec_recv_xdown:
	case sec_recv_yup:
	case sec_recv_ydown:
	case sec_recv_zup:
	case sec_recv_zdown: {
		// Standard case

		// What is stored at the given offset?
		bgq_direction d_dst = bgq_section2direction(sec);
		bgq_direction d_src = bgq_direction_revert(d_dst);
		assert(d_src == bgq_offset2dsrc(offset));
		size_t ih_dst = g_bgq_index2ih_dst[isOdd][index];
		assert(d_dst == g_bgq_index2d_dst[isOdd][index]);

		if (sec==sec_recv_tup) {
			assert(k!=1); // The value stored here is garbage (useful to sending node only)
		}
		if (sec==sec_recv_tdown) {
			assert(k!=0); // The value stored here is garbage
		}

		// Convert to consecutive
		size_t result = bgq_weyllayout_halfvolume2consecutiveoffset(isOdd_dst, ih_dst, d_dst);
		assert(d_dst == bgq_offset2ddst(result));
		return result;
	}

	case sec_surface:
	case sec_body:
		return offset;// surface and body are already consecutive
	case sec_send_xup:
	case sec_send_xdown:
	case sec_send_yup:
	case sec_send_ydown:
	case sec_send_zup:
	case sec_send_zdown:
		// Send buffers have no alternative offset
	default:
		UNREACHABLE
	}
}


static void bgq_spinorveck_written(bgq_spinorsite *targetspinor, ucoord k, ucoord t, ucoord x, ucoord y, ucoord z) {
#ifdef BGQ_COORDCHECK
	bgq_spinor coord = bgq_spinor_coord_encode(t,x,y,z);
	bgq_spinorveck_write(targetspinor, k, coord);
#endif
}


void bgq_weylveck_written(bgq_weyl_vec *targetweyl, ucoord k, ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc) {
#ifdef BGQ_COORDCHECK
	bgq_weyl_nonvec coord = bgq_weyl_coord_encode(t,x,y,z,d,isSrc);
	bgq_weylveck_write(targetweyl, k, coord);
#endif
}


static void bgq_spinorfield_fulllayout_clear(bgq_weylfield_controlblock *field) {
#ifdef BGQ_COORDCHECK
	bool isOdd = field->isOdd;
	for (ucoord z = 0; z < LOCAL_LZ ; z += 1) {
		for (ucoord y = 0; y < LOCAL_LY ; y += 1) {
			for (ucoord x = 0; x < LOCAL_LX ; x += 1) {
				for (ucoord t = 0; t < LOCAL_LT ; t += 1) {
					if (bgq_local2isOdd(t, x, y, z) != isOdd)
						continue;
					ucoord ic = bgq_local2collapsed(t, x, y, z);
					ucoord k = bgq_local2k(t, x, y, z);
					bgq_spinorveck_written(&field->sec_fullspinor[ic], k, t,x,y,z);
				}
			}
		}
	}

#if 0
	bool isOdd = field->isOdd;
	for (ucoord ic = 0; ic < PHYSICAL_VOLUME; ic += 1) {
		ucoord t1 = bgq_collapsed2t(isOdd, ic, 0);
		ucoord t2 = bgq_collapsed2t(isOdd, ic, 1);
		ucoord x = bgq_collapsed2x(isOdd, ic);
		ucoord y = bgq_collapsed2y(isOdd, ic);
		ucoord z = bgq_collapsed2z(isOdd, ic);
		field->sec_fullspinor[ic] = bgq_spinor_mergevec(bgq_spinor_coord_encode(t1, x, y, z), bgq_spinor_coord_encode(t2, x, y, z));
	}
#endif
#endif
}


static void bgq_spinorfield_weyllayout_clear(bgq_weylfield_controlblock *field) {
#ifdef BGQ_COORDCHECK
	//ucoord indices = (bgq_weyl_section_offset(sec_end) - bgq_weyl_section_offset(0))/sizeof ()

	bool isOdd = field->isOdd;
	for (size_t z = 0; z < LOCAL_LZ ; z += 1) {
		for (size_t y = 0; y < LOCAL_LY ; y += 1) {
			for (size_t x = 0; x < LOCAL_LX ; x += 1) {
				for (size_t t = 0; t < LOCAL_LT ; t += 1) {
					bool isOdd_txyz = bgq_local2isOdd(t, x, y, z);
					ucoord ic = bgq_local2collapsed(t, x, y, z);
					ucoord k = bgq_local2k(t, x, y, z);



					for (ucoord d = 0; d < PHYSICAL_LD; d += 1) {
						if (isOdd == isOdd_txyz) {
							bgq_weyl_nonvec coord = bgq_weyl_coord_encode(t, x, y, z, d, false);
							bgq_weylveck_write(&field->sec_collapsed[ic].d[d], k, coord);

							ucoord index_recv = g_bgq_collapsed2indexrecv[isOdd][ic].d[d];
							bgq_weylveck_write(bgq_index2pointer(field->sec_weyl, index_recv), k, coord);
						} else {
							bgq_weyl_nonvec coord_send = bgq_weyl_coord_encode(t, x, y, z, d, true);
							bgq_direction d_dst = bgq_direction_revert(d);
							ucoord index_send = g_bgq_collapsed2indexsend[isOdd][ic].d[d];
							bgq_weylveck_write(bgq_index2pointer(field->sec_weyl, index_send), k, coord_send);
						}
					}
				}
			}
		}
	}
#endif
}


void bgq_spinorfield_setup(bgq_weylfield_controlblock *field, bool isOdd, bool readFullspinor, bool writeFullspinor, bool readWeyl, bool writeWeyl) {
	assert(field);
	assert(readFullspinor || writeFullspinor || readWeyl || writeWeyl);
	// Do something

	bool fullspinorAvailable;
	bool weylAvailable;
	if (field->isInitinialized) {
		fullspinorAvailable = field->hasFullspinorData;
		weylAvailable = field->hasWeylfieldData || field->waitingForRecv;
	} else {
		fullspinorAvailable = false;
		weylAvailable = false;
		field->sec_weyl = NULL;
		field->sec_fullspinor = NULL;
		field->hasFullspinorData = false;
		field->hasWeylfieldData = false;
		field->waitingForRecv = false;
		field->isOdd = isOdd;
		field->isSloppy = false;

		field->isInitinialized = true;
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
		if (field->sec_fullspinor==NULL) {
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

		if (field->waitingForRecv) {
			actionWaitForRecv = true;
			actionDatamove = true;
		}
	}
	if (writeWeyl) {
		if (field->waitingForRecv) {
			// The means we are overwriting data that has never been read
			master_print("PERFORMANCE WARNING: Overwriting data not yet received from remote node, i.e. has never been used yet\n");
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
		// 4. Wait for the communication to finish
		bgq_comm_wait();
		field->waitingForRecv = false;
	}


	if (actionDatamove) {
		// 5. Move received to correct location
		bgq_master_call(&bgq_HoppingMatrix_worker_datamove, field);
	}


	if (actionAllocFulllayout) {
		size_t fullfieldsize = PHYSICAL_VOLUME * sizeof(bgq_spinorsite);
		bgq_spinorsite *spinorsite = malloc_aligned(fullfieldsize, BGQ_ALIGNMENT_L2);

		field->sec_fullspinor = spinorsite;
		field->sec_fullspinor_surface = spinorsite;
		spinorsite += PHYSICAL_SURFACE;

		field->sec_fullspinor_body = spinorsite;
	}
	if (actionInitFulllayout) {
		bgq_spinorfield_fulllayout_clear(field);
	}


	if (actionAllocWeyllayout) {
		size_t weylfieldsize = PHYSICAL_VOLUME * sizeof(bgq_weylsite);
		uint8_t *weylbase = ((uint8_t*)malloc_aligned(weylfieldsize, BGQ_ALIGNMENT_L2));

		field->sec_weyl = weylbase;
		field->sec_index = (bgq_weyl_vec*)weylbase;
		for (size_t d = 0; d < PHYSICAL_LD; d += 1) {
			field->sec_send[d] = g_bgq_sec_send[d];
			field->sec_recv[d] = g_bgq_sec_recv[d];
		}
		field->sec_collapsed = (bgq_weylsite*) (weylbase + bgq_weyl_section_offset(sec_collapsed));
		field->sec_surface = (bgq_weylsite*) (weylbase + bgq_weyl_section_offset(sec_surface));
		field->sec_body = (bgq_weylsite*) (weylbase + bgq_weyl_section_offset(sec_body));
		field->sec_end = weylbase + bgq_weyl_section_offset(sec_end);

		field->sendptr = malloc(PHYSICAL_VOLUME * sizeof(*field->sendptr));
		field->destptrFromHalfvolume = malloc(PHYSICAL_VOLUME * sizeof(*field->destptrFromHalfvolume));
		field->destptrFromSurface = malloc(PHYSICAL_SURFACE * sizeof(*field->destptrFromSurface));
		field->destptrFromBody = malloc(PHYSICAL_SURFACE * sizeof(*field->destptrFromBody));

		field->destptrFromTRecv = malloc((LOCAL_HALO_T/PHYSICAL_LP + 0 + LOCAL_HALO_T/PHYSICAL_LP) * sizeof(bgq_weyl_vec*)); // COMM_T==1

		size_t offset_begin = bgq_weyl_section_offset(sec_recv_xup);
		size_t offset_end = bgq_weyl_section_offset(sec_recv_zdown + 1);
		size_t weylCount = (offset_end - offset_begin) / sizeof(bgq_weyl_vec);
		field->destptrFromRecv = malloc(weylCount * sizeof(bgq_weyl_vec*));
	}

	if (actionInitWeylPtrs) {
		uint8_t *weylbase = field->sec_weyl;
		assert(weylbase);

		// For 1st phase (distribute)
		for (ucoord ic_src=0; ic_src < PHYSICAL_VOLUME; ic_src+=1) {
			for (size_t d_src = 0; d_src < PHYSICAL_LD; d_src += 1) {
				field->sendptr[ic_src].d[d_src] = bgq_index2pointer(weylbase, g_bgq_collapsed2indexsend[isOdd][ic_src].d[d_src]);
			}
		}

		for (size_t ih_src = 0; ih_src < PHYSICAL_VOLUME; ih_src += 1) {
			for (size_t d_src = 0; d_src < PHYSICAL_LD; d_src += 1) {
				field->destptrFromHalfvolume[ih_src].d[d_src] = bgq_encodedoffset2pointer(weylbase, g_bgq_ihsrc2offsetwrite[isOdd][ih_src].d[d_src]);
			}
		}


		for (size_t is_src = 0; is_src < PHYSICAL_SURFACE; is_src += 1) {
			for (size_t d_src = 0; d_src < PHYSICAL_LD; d_src += 1) {
				field->destptrFromSurface[is_src].d[d_src] = bgq_encodedoffset2pointer(weylbase, g_bgq_issrc2offsetwrite[isOdd][is_src].d[d_src]);
				assert(field->destptrFromSurface[is_src].d[d_src] == bgq_encodedoffset2pointer(weylbase, g_bgq_ihsrc2offsetwrite[isOdd][bgq_surface2halfvolume(isOdd, is_src)].d[d_src]));
			}
		}

		for (size_t ib_src = 0; ib_src < PHYSICAL_BODY; ib_src += 1) {
			for (size_t d_src = 0; d_src < PHYSICAL_LD; d_src += 1) {
				field->destptrFromBody[ib_src].d[d_src] = bgq_encodedoffset2pointer(weylbase, g_bgq_ibsrc2offsetwrite[isOdd][ib_src].d[d_src]);
				assert(field->destptrFromBody[ib_src].d[d_src] == bgq_encodedoffset2pointer(weylbase, g_bgq_ihsrc2offsetwrite[isOdd][bgq_body2halfvolume(isOdd, ib_src)].d[d_src]));
			}
		}

		// For 5th phase (datamove)
		if (COMM_T) {
			size_t beginj = 0;
			size_t endj = 2*LOCAL_HALO_T/PHYSICAL_LP; // There should be no padding between sections
			assert(endj == (bgq_weyl_section_offset(sec_recv_tup+1) - bgq_weyl_section_offset(sec_send_tup))/sizeof(bgq_weyl_vec));
			for (int j = beginj; j < endj; j+=1) {
				size_t offset_left = bgq_weyl_section_offset(sec_send_tup) + j * sizeof(bgq_weyl_vec); // Overlap into sec_recv_tup
				bgq_weylfield_section sec_left = bgq_sectionOfOffset(offset_left);
				size_t offset_right = bgq_weyl_section_offset(sec_recv_tdown) + j * sizeof(bgq_weyl_vec); // Overlap into sec_send_tdown
				bgq_weylfield_section sec_right = bgq_sectionOfOffset(offset_right);
				assert(((sec_left==sec_send_tup)&&(sec_right==sec_recv_tdown)) || ((sec_left==sec_recv_tup)&&(sec_right==sec_send_tdown)));

				size_t consecutive_left = bgq_weylfield_bufferoffset2consecutiveoffset(isOdd, offset_left, 0);
				size_t consecutive_right = bgq_weylfield_bufferoffset2consecutiveoffset(isOdd, offset_right, 1);
				assert(consecutive_left == consecutive_right);

				assert(consecutive_left);
				field->destptrFromTRecv[j] = bgq_offset2pointer(weylbase, consecutive_left);
			}
		} else {
			assert(!"Not yet implemented");
			exit(0);
		}

		size_t offset_begin = bgq_weyl_section_offset(sec_recv_xup);
		size_t offset_end = bgq_weyl_section_offset(sec_recv_zdown + 1);
		size_t weylCount = (offset_end - offset_begin) / sizeof(bgq_weyl_vec);
		for (size_t j = 0; j < weylCount; j+=1) {
			size_t offset = offset_begin + j*sizeof(bgq_weyl_vec); // Overlaps into following recv sections
			bgq_weylfield_section sec = bgq_sectionOfOffset(offset);

			// To check that this is a recv buffer
			switch (sec) {
			case sec_recv_xup:
				case sec_recv_xdown:
				case sec_recv_yup:
				case sec_recv_ydown:
				case sec_recv_zup:
				case sec_recv_zdown:
				break;
			default:
				assert(!"The sections  must be consecutive");
				UNREACHABLE
			}

			size_t offset_consecutive = bgq_weylfield_bufferoffset2consecutiveoffset(isOdd, offset, -1/*doesn't matter*/);
			field->destptrFromRecv[j] = bgq_offset2pointer(weylbase, offset_consecutive);
		}

		bgq_spinorfield_weyllayout_clear(field);
	}


	if (writeWeyl || writeFullspinor) {
		// Data is going to be written, so what actually is stored here changes
		field->hasWeylfieldData = writeWeyl;
		field->hasFullspinorData = writeFullspinor;
	}
}

typedef struct {
	bgq_weylfield_controlblock *field;
} bgq_conversion_args;


static void bgq_spinorfield_fullspinor2weyl_worker(void *args_untyped, size_t tid, size_t threads) {
	bgq_conversion_args *args = args_untyped;
	bgq_weylfield_controlblock *field = args->field;
	bool isOdd = field->isOdd;

	const size_t workload = PHYSICAL_VOLUME;
	const size_t threadload = (workload + threads - 1) / threads;
	const size_t begin = tid * threadload;
	const size_t end = min(workload, begin + threadload);
	for (size_t ih = begin; ih < end; ih += 1) {
		bgq_spinorsite *fullspinoraddr = &field->sec_fullspinor[ih];
		bgq_su3_spinor_decl(fullspinor);
		bgq_su3_spinor_load_double(fullspinor, fullspinoraddr);

		bgq_spinorfield_weyl_store_fromHalfvolume(field, isOdd, ih, bgq_su3_spinor_vars(fullspinor));
	}
}


static void bgq_spinorfield_fullspinor2weyl(bgq_weylfield_controlblock *field) {
	assert(field->isInitinialized);
	assert(field->hasFullspinorData);

	bgq_spinorfield_setup(field, field->isOdd, true, false, false, true);
	bgq_conversion_args args = {
			.field = field
	};
	bgq_master_call(&bgq_spinorfield_fullspinor2weyl_worker, &args); field->hasWeylfieldData = true;
}





void bgq_spinorfield_transfer(bool isOdd, bgq_weylfield_controlblock *targetfield, spinor *sourcefield) {
	bgq_spinorfield_setup(targetfield, isOdd, false, true, false, false);
	size_t ioff = isOdd ? (VOLUME+RAND)/2 : 0;

	for (size_t i_eosub = 0; i_eosub < VOLUME/2; i_eosub+=1) {
		size_t i_eo = i_eosub + ioff;
		size_t i_lexic = g_eo2lexic[i_eo];
		int t = g_coord[i_lexic][0];
		int x = g_coord[i_lexic][1];
		int y = g_coord[i_lexic][2];
		int z = g_coord[i_lexic][3];
		spinor_array64 *sourcespinor = (spinor_array64*)&sourcefield[i_eosub];
		assert(bgq_local2isOdd(t,x,y,z)==isOdd);

		ucoord ih = bgq_local2halfvolume(t,x,y,z);
		ucoord ic = bgq_local2collapsed(t,x,y,z);
		ucoord k = bgq_local2k(t,x,y,z);
		bgq_spinorsite *targetspinor = &targetfield->sec_fullspinor[ic];
		bgq_spinorveck_expect(*targetspinor, k, t,x,y,z);
		for (int v = 0; v < 4; v+=1) {
			for (int c = 0; c < 3; c+=1) {
				targetspinor->s[v][c][k] = sourcespinor->v[v].c[c];
			}
		}
		bgq_spinorveck_written(targetspinor, k, t,x,y,z);
	}
}


bgq_spinor bgq_spinorfield_getspinor(bgq_weylfield_controlblock *field, ucoord t, ucoord x, ucoord y, ucoord z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);

	ucoord ic = bgq_local2collapsed(t,x,y,z);
	ucoord k = bgq_local2k(t,x,y,z);
	if (field->hasFullspinorData) {
		bgq_spinorsite spinor = field->sec_fullspinor[ic];
		return bgq_spinor_fromvec(spinor,k);
	} else if (field->hasWeylfieldData) {
		bgq_su3_spinor_decl(spinor);
		bgq_HoppingMatrix_loadWeyllayout(spinor,&field->sec_collapsed[ic], bgq_t2t(t,0), bgq_t2t(t,1), x,y,z);
		return bgq_spinor_fromqpx(spinor,k);
	} else {
		assert(!"Field contains no data");
		UNREACHABLE
	}
}
