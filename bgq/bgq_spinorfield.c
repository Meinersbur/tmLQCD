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
	bgq_spinorfield_setup(bgqfield, isOdd, readFulllayout, false, !readFulllayout, false);
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
									master_print("Coordinate (%zu,%zu,%zu,%zu)(%zu,%zu): ref=(%f + %fi) != bgb=(%f + %fi) off by %f\n", t, x, y, z, v, c, creal(refvalue), cimag(refvalue), creal(bgqvalue), cimag(bgqvalue), diff);
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
	assert(weyl.s[1][2] == 0);
}






static inline void bgq_HoppingMatrix_worker_datamove_recvxyz(bgq_weylfield_controlblock *spinorfield, bool isOdd, size_t beginj, size_t endj, bool noprefetchstream) {
	if (!noprefetchstream) {
		bgq_prefetch_forward(&spinorfield->sec_recv[XUP][beginj]);
		bgq_prefetch_forward(&spinorfield->destptrFromRecv[beginj]);
	}

	for (size_t j = beginj; j < endj; j+=1) {
#ifndef NDEBUG
		size_t offset = (uint8_t*)&g_bgq_sec_recv[XUP][j] - g_bgq_sec_comm + bgq_weyl_section_offset(sec_comm);
		ucoord index = bgq_offset2index(offset);
		ucoord ic = g_bgq_index2collapsed[isOdd][index];
		ucoord t1 = bgq_collapsed2t(isOdd, ic, 0);
		ucoord t2 = bgq_collapsed2t(isOdd, ic, 1);
		ucoord x = bgq_collapsed2x(isOdd, ic);
		ucoord y = bgq_collapsed2y(isOdd, ic);
		ucoord z = bgq_collapsed2z(isOdd, ic);
		bgq_direction d = bgq_offset2ddst(offset);
#endif

		//TODO: Check strength reduction
		//TODO: Prefetch
		//TODO: Inline assembler
		bgq_weyl_vec *weyladdr_src= &spinorfield->sec_recv[XUP][j]; // Note: overlaps into following sections
		bgq_weyl_vec *weyladdr_dst = spinorfield->destptrFromRecv[j];
		assert((bgq_weyl_vec *)spinorfield->sec_weyl <= weyladdr_dst && weyladdr_dst <  (bgq_weyl_vec *)spinorfield->sec_end);

		bgq_prefetch(&spinorfield->destptrFromRecv[j+1]);
		bgq_su3_weyl_prefetch_double(&spinorfield->sec_recv[XUP][j+1]);

		bgq_su3_weyl_decl(weyl);
		bgq_su3_weyl_load_double(weyl, weyladdr_src);
		bgq_weylqpx_expect(weyl, t1, t2, x, y, z, d, false);

		//bgq_weylvec_expect(*weyladdr_dst, t1,t2,x,y,z,d,false);
		bgq_su3_weyl_store_double(weyladdr_dst, weyl);
		//bgq_weylvec_written(weyladdr_dst, t1, t2, x, y, z, d, false);
	}
}


static inline void bgq_HoppingMatrix_worker_datamove_recvtup(bgq_weylfield_controlblock *spinorfield, bool isOdd, size_t beginj, size_t endj, bool noprefetchstream) {
	assert(COMM_T);
	if (!noprefetchstream) {
		bgq_prefetch_forward(&g_bgq_sec_send[TDOWN][beginj]);
		bgq_prefetch_forward(&g_bgq_sec_recv[TUP][beginj]);
		bgq_prefetch_forward(&spinorfield->consptr_recvtup[beginj]);
	}

	for (size_t j = beginj; j < endj; j+=1) {
#ifndef NDEBUG
		size_t offset_left =  (uint8_t*)&g_bgq_sec_send[TDOWN][j] - g_bgq_sec_comm + bgq_weyl_section_offset(sec_comm);
		ucoord index_left = bgq_offset2index(offset_left);
		size_t offset_right =(uint8_t*)&g_bgq_sec_recv[TUP][j] - g_bgq_sec_comm + bgq_weyl_section_offset(sec_comm);
		ucoord index_right = bgq_offset2index(offset_right);
		ucoord ic_left = g_bgq_index2collapsed[isOdd][index_left];
		ucoord ic_right = g_bgq_index2collapsed[isOdd][index_right];
		assert(ic_left == ic_right);
		ucoord ic = ic_left;
		ucoord t1 = bgq_collapsed2t(isOdd, ic_left, 0);
		ucoord t2 = bgq_collapsed2t(isOdd, ic_right, 1);
		ucoord x = bgq_collapsed2x(isOdd, ic);
		ucoord y = bgq_collapsed2y(isOdd, ic);
		ucoord z = bgq_collapsed2z(isOdd, ic);
		bgq_direction d1 = bgq_offset2ddst(offset_left);
		bgq_direction d2 = bgq_offset2ddst(offset_right);
		assert(d1==d2);
		bgq_direction d = d1;
#endif

		//TODO: Check strength reduction
		//TODO: Inline assembler
		//TODO: Is reading just the 16 used bytes faster?
		bgq_weyl_vec *weyladdr_left = &g_bgq_sec_send[TDOWN][j];
		bgq_weyl_vec *weyladdr_right = &g_bgq_sec_recv[TUP][j];
		bgq_weyl_vec *weyladdr_dst = spinorfield->consptr_recvtup[j];

		bgq_prefetch(&spinorfield->consptr_recvtup[j+1]);
		bgq_su3_weyl_prefetch_double(&g_bgq_sec_send[TDOWN][j+1]);
		bgq_su3_weyl_prefetch_double(&g_bgq_sec_recv[TUP][j+1]);

		bgq_su3_weyl_decl(weyl_left);
		bgq_su3_weyl_load_double(weyl_left, weyladdr_left);
		assert(bgq_cmplxval2(weyl_left_v0_c0)!=0); // for valgrind
		bgq_weylqpxk_expect(weyl_left, 1, t1, x, y, z, d1, false);

		bgq_su3_weyl_decl(weyl_right);
		bgq_su3_weyl_load_double(weyl_right, weyladdr_right);
		assert(bgq_cmplxval1(weyl_left_v0_c0)!=0); // for valgrind
		bgq_weylqpxk_expect(weyl_right, 0, t2, x, y, z, d2, false);

		bgq_su3_weyl_decl(weyl);
		bgq_su3_weyl_merge2(weyl, weyl_left, weyl_right);
		bgq_weylqpx_expect(weyl, t1, t2, x, y, z, d, false);

#ifndef NDEBUG
		//bgq_weylvec_expect(*weyladdr_dst, t1,t2,x,y,z,d,false);
		size_t offset = bgq_pointer2offset(spinorfield, weyladdr_dst);
		ucoord index = bgq_offset2index(offset);
		bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
		assert(sec==sec_body || sec==sec_surface);
		ucoord ic_check = g_bgq_index2collapsed[isOdd][index];
		assert(ic == ic_check);
#endif
		bgq_su3_weyl_store_double(weyladdr_dst, weyl);
		//bgq_weylvec_written(weyladdr_dst, t1, t2, x, y, z, d, false);
	}
}

static inline void bgq_HoppingMatrix_worker_datamove_recvtdown(bgq_weylfield_controlblock *spinorfield, bool isOdd, size_t beginj, size_t endj, bool noprefetchstream) {
	assert(COMM_T);
	if (!noprefetchstream) {
		bgq_prefetch_forward(&g_bgq_sec_recv[TDOWN][beginj]);
		bgq_prefetch_forward(&g_bgq_sec_send[TUP][beginj]);
		bgq_prefetch_forward(&spinorfield->consptr_recvtup[beginj]);
	}

	for (size_t j = beginj; j < endj; j+=1) {
#ifndef NDEBUG
		size_t offset_left = (uint8_t*)&g_bgq_sec_recv[TDOWN][j] - g_bgq_sec_comm + bgq_weyl_section_offset(sec_comm);
		ucoord index_left = bgq_offset2index(offset_left);
		size_t offset_right = (uint8_t*)&g_bgq_sec_send[TUP][j] - g_bgq_sec_comm + bgq_weyl_section_offset(sec_comm);
		ucoord index_right = bgq_offset2index(offset_right);
		ucoord ic_left = g_bgq_index2collapsed[isOdd][index_left];
		ucoord ic_right = g_bgq_index2collapsed[isOdd][index_right];
		assert(ic_left == ic_right);
		ucoord ic = ic_left;
		ucoord t1 = bgq_collapsed2t(isOdd, ic_left, 0);
		ucoord t2 = bgq_collapsed2t(isOdd, ic_right, 1);
		ucoord x = bgq_collapsed2x(isOdd, ic);
		ucoord y = bgq_collapsed2y(isOdd, ic);
		ucoord z = bgq_collapsed2z(isOdd, ic);
		bgq_direction d1 = bgq_offset2ddst(offset_left);
		bgq_direction d2 = bgq_offset2ddst(offset_right);
		assert(d1==d2);
		bgq_direction d = d1;
#endif

		//TODO: Check strength reduction
		//TODO: Prefetch
		//TODO: Inline assembler
		bgq_weyl_vec *weyladdr_left = &g_bgq_sec_recv[TDOWN][j];
		bgq_weyl_vec *weyladdr_right = &g_bgq_sec_send[TUP][j];
		bgq_weyl_vec *weyladdr_dst = spinorfield->consptr_recvtdown[j];

		bgq_prefetch(&spinorfield->consptr_recvtdown[j+1]);
		bgq_su3_weyl_prefetch_double(&g_bgq_sec_recv[TDOWN][j+1]);
		bgq_su3_weyl_prefetch_double(&g_bgq_sec_send[TUP][j+1]);

		bgq_su3_weyl_decl(weyl_left);
		bgq_su3_weyl_load_double(weyl_left, weyladdr_left);
		bgq_weylqpxk_expect(weyl_left, 1, t1, x, y, z, d1, false);
		assert(bgq_elem2(weyl_left_v0_c0)!=0); // for valgrind

		bgq_su3_weyl_decl(weyl_right);
		bgq_su3_weyl_load_double(weyl_right, weyladdr_right);
		bgq_weylqpxk_expect(weyl_right, 0, t2, x, y, z, d2, false);
		assert(bgq_elem0(weyl_right_v0_c0)!=0);// for valgrind

		bgq_su3_weyl_decl(weyl);
		bgq_su3_weyl_merge2(weyl, weyl_left, weyl_right);
		bgq_weylqpx_expect(weyl, t1, t2, x, y, z, d, false);

#ifndef NDEBUG
		//bgq_weylvec_expect(*weyladdr_dst, t1,t2,x,y,z,d,false);
		size_t offset = bgq_pointer2offset(spinorfield, weyladdr_dst);
		ucoord index = bgq_offset2index(offset);
		bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
		assert(sec==sec_body || sec==sec_surface);
		ucoord ic_check = g_bgq_index2collapsed[isOdd][index];
		assert(ic == ic_check);
#endif
		bgq_su3_weyl_store_double(weyladdr_dst, weyl);
		//bgq_weylvec_written(weyladdr_dst, t1, t2, x, y, z, d, false);
	}
}


static inline void bgq_HoppingMatrix_worker_datamovet_recvtup(bgq_weylfield_controlblock *spinorfield, bool isOdd, size_t beginj, size_t endj, bool noprefetchstream) {
	assert(!COMM_T);
	if (!noprefetchstream) {
		bgq_prefetch_forward(&g_bgq_sec_send[TDOWN][beginj]);
		bgq_prefetch_forward(&g_bgq_sec_recv[TUP][beginj]);
		bgq_prefetch_forward(&spinorfield->consptr_recvtup[beginj]);
	}

	for (size_t j = beginj; j < endj; j+=1) {
#ifndef NDEBUG
		size_t offset =(uint8_t*)&g_bgq_sec_send[TDOWN][j] - g_bgq_sec_comm + bgq_weyl_section_offset(sec_comm);
		ucoord index = bgq_offset2index(offset);
		ucoord ic = g_bgq_index2collapsed[isOdd][index];
		ucoord t1 = bgq_collapsed2t(isOdd, ic, 0);
		ucoord t2 = bgq_collapsed2t(isOdd, ic, 1);
		ucoord x = bgq_collapsed2x(isOdd, ic);
		ucoord y = bgq_collapsed2y(isOdd, ic);
		ucoord z = bgq_collapsed2z(isOdd, ic);
		bgq_direction d = bgq_offset2ddst(offset);
#endif

		//TODO: Check strength reduction
		//TODO: Inline assembler
		bgq_weyl_vec *weyladdr = &g_bgq_sec_send[TDOWN][j];
		bgq_weyl_vec *weyladdr_dst = spinorfield->consptr_recvtup[j];

		bgq_prefetch(&spinorfield->consptr_recvtup[j+1]);
		bgq_su3_weyl_prefetch_double(&g_bgq_sec_send[TDOWN][j+1]);

		bgq_su3_weyl_decl(weyl_before);
		bgq_su3_weyl_load_double(weyl_before, weyladdr);
		bgq_weylqpxk_expect(weyl_before, 0, t2, x, y, z, d, false);
		bgq_weylqpxk_expect(weyl_before, 1, t1, x, y, z, d, false);

		bgq_su3_weyl_decl(weyl_after);
		bgq_su3_weyl_merge2(weyl_after, weyl_before, weyl_before); // Just switch first and second complex
		bgq_weylqpx_expect(weyl_after, t1, t2, x, y, z, d, false);
		bgq_su3_weyl_store_double(weyladdr_dst, weyl_after);
	}
}


static inline void bgq_HoppingMatrix_worker_datamovet_recvtdown(bgq_weylfield_controlblock *spinorfield, bool isOdd, size_t beginj, size_t endj, bool noprefetchstream) {
	assert(!COMM_T);

	if (!noprefetchstream) {
		bgq_prefetch_forward(&g_bgq_sec_recv[TDOWN][beginj]);
		bgq_prefetch_forward(&g_bgq_sec_send[TUP][beginj]);
		bgq_prefetch_forward(&spinorfield->consptr_recvtup[beginj]);
	}

		for (size_t j = beginj; j < endj; j+=1) {
#ifndef NDEBUG
			size_t offset = (uint8_t*)&g_bgq_sec_send[TUP][j] - g_bgq_sec_comm + bgq_weyl_section_offset(sec_comm);
			ucoord index = bgq_offset2index(offset);
			ucoord ic = g_bgq_index2collapsed[isOdd][index];
			ucoord t1 = bgq_collapsed2t(isOdd, ic, 0);
			ucoord t2 = bgq_collapsed2t(isOdd, ic, 1);
			ucoord x = bgq_collapsed2x(isOdd, ic);
			ucoord y = bgq_collapsed2y(isOdd, ic);
			ucoord z = bgq_collapsed2z(isOdd, ic);
			bgq_direction d = bgq_offset2ddst(offset);
#endif
			//TODO: Check strength reduction
			//TODO: Prefetch
			//TODO: Inline assembler
			bgq_weyl_vec *weyladdr = &g_bgq_sec_send[TUP][j];
			bgq_weyl_vec *weyladdr_dst = spinorfield->consptr_recvtdown[j];

			bgq_prefetch(&spinorfield->consptr_recvtdown[j+1]);
			bgq_su3_weyl_prefetch_double(&g_bgq_sec_send[TUP][j+1]);

			bgq_su3_weyl_decl(weyl_before);
			bgq_su3_weyl_load_double(weyl_before, weyladdr);
			bgq_weylqpxk_expect(weyl_before, 0, t2, x, y, z, d2, false);
			bgq_weylqpxk_expect(weyl_before, 1, t1, x, y, z, d1, false);

			bgq_su3_weyl_decl(weyl_after);
			bgq_su3_weyl_merge2(weyl_after, weyl_before, weyl_before);
			bgq_weylqpx_expect(weyl_after, t1, t2, x, y, z, d, false);

			bgq_su3_weyl_store_double(weyladdr_dst, weyl_after);
		}
}




void bgq_HoppingMatrix_datamovet_worker(void *arg_untyped, size_t tid, size_t threads) {
	bgq_work_datamove *arg = arg_untyped;
	bgq_weylfield_controlblock *spinorfield = arg->spinorfield;
	bgq_hmflags opts = arg->opts;
	bool isOdd = spinorfield->isOdd;
	bool noprefetchstream = opts & hm_noprefetchstream;
	assert(!COMM_T);

	const size_t workload_recv_tup = LOCAL_HALO_T/PHYSICAL_LP;
	const size_t workload_recv_tdown = workload_recv_tup;
	const size_t workload = workload_recv_tup + workload_recv_tdown;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);
	for (size_t i = begin; i<end; ) {
		WORKLOAD_DECL(i,workload);

		if (WORKLOAD_SPLIT(workload_recv_tup)) {
			size_t beginj = WORKLOAD_PARAM(workload_recv_tup);
			size_t endj = min_sizet(workload_recv_tup,beginj+threadload);
			bgq_HoppingMatrix_worker_datamovet_recvtup(spinorfield,isOdd,beginj,endj,noprefetchstream);
			i += (endj - beginj);
		} else if (WORKLOAD_SPLIT(workload_recv_tdown)) {
			size_t beginj = WORKLOAD_PARAM(workload_recv_tdown);
			size_t endj = min_sizet(workload_recv_tdown,beginj+threadload);
			bgq_HoppingMatrix_worker_datamovet_recvtdown(spinorfield,isOdd,beginj,endj,noprefetchstream);
			i += (endj - beginj);
		} else {
			UNREACHABLE
		}

		WORKLOAD_CHECK
	}
}


static void bgq_HoppingMatrix_worker_datamove(void *arg_untyped, size_t tid, size_t threads) {
	bgq_work_datamove *arg = arg_untyped;
	bgq_weylfield_controlblock *spinorfield = arg->spinorfield;
	bgq_hmflags opts = arg->opts;
	bool isOdd = spinorfield->isOdd;

	bool noprefetchstream = opts & hm_noprefetchstream;

	//const size_t workload_recvt = COMM_T ? 2*LOCAL_HALO_T/PHYSICAL_LP : LOCAL_HALO_T/PHYSICAL_LP;
	const size_t workload_recv_tup = COMM_T * LOCAL_HALO_T/PHYSICAL_LP;
	const size_t workload_recv_tdown = workload_recv_tup;
	const size_t workload_recv = 2*PHYSICAL_HALO_X + 2*PHYSICAL_HALO_Y + 2*PHYSICAL_HALO_Z;
	const size_t workload = workload_recv + 2*workload_recv_tup + 2*workload_recv_tdown;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);
	for (size_t i = begin; i<end; ) {
		WORKLOAD_DECL(i,workload);

		if (!COMM_T || WORKLOAD_SPLIT(workload_recv)) {
			// Do other dimensions
			size_t beginj = WORKLOAD_PARAM(workload_recv);
			size_t endj = min_sizet(workload_recv,beginj+threadload);
			bgq_HoppingMatrix_worker_datamove_recvxyz(spinorfield,isOdd,beginj,endj,noprefetchstream);
			i += (endj - beginj);
		} else if (COMM_T && WORKLOAD_SPLIT(2*workload_recv_tup)) {
			// Count an T-iteration twice; better have few underloaded threads (so remaining SMT-threads have some more resources) then few overloaded threads (so the master thread has to wait for them)
			size_t twobeginj = WORKLOAD_PARAM(2*workload_recv_tup);
			size_t twoendj = min_sizet(2*workload_recv_tup,twobeginj+threadload);
			size_t beginj = twobeginj / 2;
			size_t endj = twoendj / 2;
			bgq_HoppingMatrix_worker_datamove_recvtup(spinorfield,isOdd,beginj,endj,noprefetchstream);
			i += (twoendj - twobeginj);
		} else if (COMM_T && WORKLOAD_SPLIT(2*workload_recv_tdown)) {
			// Count an T-iteration twice; better have few underloaded threads (so remaining SMT-threads have some more resources) then few overloaded threads (so the master thread has to wait for them)
			size_t twobeginj = WORKLOAD_PARAM(2*workload_recv_tdown);
			size_t twoendj = min_sizet(2*workload_recv_tdown,twobeginj+threadload);
			size_t beginj = twobeginj / 2;
			size_t endj = twoendj / 2;
			bgq_HoppingMatrix_worker_datamove_recvtdown(spinorfield,isOdd,beginj,endj,noprefetchstream);
			i += (twoendj - twobeginj);
		} else {
			UNREACHABLE
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


size_t bgq_collapsed2consecutiveoffset(ucoord ic, bgq_direction d) {
	return bgq_weyl_section_offset(sec_collapsed) + ic*sizeof(bgq_weylsite) + d*sizeof(bgq_weyl_vec);
}


static size_t bgq_weylfield_bufferoffset2consecutiveoffset(bool isOdd, size_t offset, size_t k) {
	assert(offset);
	assert(offset % sizeof(bgq_weyl_vec) == 0);
	bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
	size_t index = bgq_offset2index(offset);

	//bgq_direction d = bgq_section2direction(sec);
//assert(d == bgq_offset2ddst(offset));
	bgq_direction d = bgq_offset2ddst(offset);
	ucoord ic = g_bgq_index2collapsed[isOdd][index];
	return bgq_collapsed2consecutiveoffset(ic, d);
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

		if (field->pendingDatamove) {
			actionDatamove = true;
		}
		if (field->waitingForRecv) {
			assert(actionDatamove);
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
		bgq_comm_wait(nospi);
		field->waitingForRecv = false;
	}
	if (actionDatamove) {
		// 5. Move received to correct location
		bgq_master_sync();
		static bgq_work_datamove work;
		work.spinorfield = field;
		work.opts = field->hmflags;
		bgq_master_call(&bgq_HoppingMatrix_worker_datamove, &work);
		field->pendingDatamove = false;
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

		field->sendptr = malloc_aligned(PHYSICAL_VOLUME * sizeof(*field->sendptr), BGQ_ALIGNMENT_L2);
		field->destptrFromHalfvolume = malloc(PHYSICAL_VOLUME * sizeof(*field->destptrFromHalfvolume));
		field->destptrFromSurface = malloc(PHYSICAL_SURFACE * sizeof(*field->destptrFromSurface));
		field->destptrFromBody = malloc(PHYSICAL_BODY * sizeof(*field->destptrFromBody));


		size_t offset_begin = bgq_weyl_section_offset(sec_recv_xup);
		size_t offset_end = bgq_weyl_section_offset(sec_recv_zdown + 1);
		size_t weylCount = (offset_end - offset_begin) / sizeof(bgq_weyl_vec);
		field->destptrFromRecv = malloc_aligned(weylCount * sizeof(bgq_weyl_vec*), BGQ_ALIGNMENT_L2);

		weylCount = bgq_section_size(sec_send_tup) / sizeof(bgq_weyl_vec);
		field->consptr_recvtdown = malloc_aligned(weylCount * sizeof(*field->consptr_recvtdown), BGQ_ALIGNMENT_L2);
		field->consptr_recvtup = malloc_aligned(weylCount * sizeof(*field->consptr_recvtup), BGQ_ALIGNMENT_L2);
	}

	if (actionInitWeylPtrs) {
		uint8_t *weylbase = field->sec_weyl;
		assert(weylbase);

		// For 1st phase (distribute)
		for (ucoord ic_src=0; ic_src < PHYSICAL_VOLUME; ic_src+=1) {
			for (size_t d_src = 0; d_src < PHYSICAL_LD; d_src += 1) {
				ucoord index = g_bgq_collapsed2indexsend[isOdd][ic_src].d[d_src];
				bgq_weyl_vec *ptr =  bgq_index2pointer(weylbase, index);
				assert(bgq_pointer2offset(field,ptr)==bgq_index2offset(index));
				if (bgq_fieldpointer2offset(ptr)!=bgq_index2offset(index)) {
					bgq_fieldpointer2offset(ptr);
				}
				assert(bgq_fieldpointer2offset(ptr)==bgq_index2offset(index));
				field->sendptr[ic_src].d[d_src] = ptr;
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
		 {
			 // actually, it makes no difference
			 bgq_weylfield_section secleft_tup = COMM_T ? sec_recv_tup : sec_send_tdown;
			 bgq_weylfield_section secright_tdown = COMM_T ? sec_recv_tdown : sec_send_tup;

			ucoord endj = bgq_offset2index(bgq_section_size(sec_send_tdown));
			for (int j = 0; j < endj; j+=1) {
				size_t offset_left = bgq_weyl_section_offset(sec_send_tdown) + j*sizeof(bgq_weyl_vec);
				ucoord index_left = bgq_offset2index(offset_left);
				bgq_weylfield_section sec_left = bgq_sectionOfOffset(offset_left);
				size_t offset_right = bgq_weyl_section_offset(secleft_tup) + j*sizeof(bgq_weyl_vec);
				ucoord index_right = bgq_offset2index(offset_right);
				bgq_weylfield_section sec_right = bgq_sectionOfOffset(offset_right);

				size_t consecutive_left = bgq_weylfield_bufferoffset2consecutiveoffset(isOdd, offset_left, 1);
				ucoord index_consecutive_left = bgq_offset2index(consecutive_left);
				size_t consecutive_right = bgq_weylfield_bufferoffset2consecutiveoffset(isOdd, offset_right, 0);
				ucoord index_consecutive_right = bgq_offset2index(consecutive_right);
				assert(consecutive_left == consecutive_right);

				bgq_weyl_vec *ptr = bgq_offset2pointer(weylbase, consecutive_left);
				assert((uintptr_t)ptr % 32 == 0);
				field->consptr_recvtup[j] = ptr;
			}


			endj =  bgq_offset2index(bgq_section_size(sec_send_tup));
			for (int j = 0; j < endj; j+=1) {
				size_t offset_left = bgq_weyl_section_offset(secright_tdown) + j*sizeof(bgq_weyl_vec);
				ucoord index_left = bgq_offset2index(offset_left);
				bgq_weylfield_section sec_left = bgq_sectionOfOffset(offset_left);
				size_t offset_right = bgq_weyl_section_offset(sec_send_tup) + j*sizeof(bgq_weyl_vec);
				ucoord index_right = bgq_offset2index(offset_right);
				bgq_weylfield_section sec_right = bgq_sectionOfOffset(offset_right);

				size_t consecutive_left = bgq_weylfield_bufferoffset2consecutiveoffset(isOdd, offset_left, 1);
				ucoord index_consecutive_left = bgq_offset2index(consecutive_left);
				size_t consecutive_right = bgq_weylfield_bufferoffset2consecutiveoffset(isOdd, offset_right, 0);
				ucoord index_consecutive_right = bgq_offset2index(consecutive_right);
				assert(consecutive_left == consecutive_right);

				bgq_weyl_vec *ptr = bgq_offset2pointer(weylbase, consecutive_left);
				assert((uintptr_t)ptr % 32 == 0);
				field->consptr_recvtdown[j] = ptr;
			}
		}

		size_t offset_begin = bgq_weyl_section_offset(sec_recv_xup);
		size_t offset_end = bgq_weyl_section_offset(sec_recv_zdown + 1);
		size_t weylCount = (offset_end - offset_begin) / sizeof(bgq_weyl_vec);
		for (size_t j = 0; j < weylCount; j+=1) {
			size_t offset = offset_begin + j*sizeof(bgq_weyl_vec); // Overlaps into following recv sections
			bgq_weylfield_section sec = bgq_sectionOfOffset(offset);

			size_t offset_consecutive = bgq_weylfield_bufferoffset2consecutiveoffset(isOdd, offset, -1/*doesn't matter*/);
			field->destptrFromRecv[j] = bgq_offset2pointer(weylbase, offset_consecutive);
		}
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


void bgq_spinorfield_transfer(bool isOdd, bgq_weylfield_controlblock *targetfield, spinor *sourcefield) {
	bgq_spinorfield_setup(targetfield, isOdd, false, true, false, false);
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
		bgq_spinorsite *targetspinor = &targetfield->sec_fullspinor[ic];
		//bgq_spinorveck_expect(*targetspinor, k, t,x,y,z);
		for (int v = 0; v < 4; v+=1) {
			for (int c = 0; c < 3; c+=1) {
				targetspinor->s[v][c][k] = sourcespinor->v[v].c[c];
			}
		}
		bgq_spinorveck_written(targetspinor, k, t,x,y,z);
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

	assert(bgq_local2isOdd(t,x,y,z)==field->isOdd);
	ucoord ic = bgq_local2collapsed(t,x,y,z);
	ucoord k = bgq_local2k(t,x,y,z);
	if (field->hasFullspinorData) {
		bgq_spinorfield_setup(field, field->isOdd, true, false, false, false);
		bgq_spinorsite spinor = field->sec_fullspinor[ic];
		return bgq_spinor_fromvec(spinor,k);
	} else if (field->hasWeylfieldData) {
		bgq_spinorfield_setup(field, field->isOdd, false, false, true, false);
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
}





char *(g_idxdesc[BGQREF_count]);
complexdouble *g_bgqvalue = NULL;
complexdouble *g_refvalue = NULL;


void bgq_initbgqref_impl() {
	int datasize = sizeof(complexdouble) * VOLUME * lengthof(g_idxdesc);
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
void bgq_setbgqvalue_impl(int t, int x, int y, int z, bgqref idx, complexdouble val) {
	if (idx==BGQREF_TUP && t==0 && x==0 && y==0 && z==0) {
		int a = 0;
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
}



size_t bgq_fieldpointer2offset(void *ptr) {
	size_t commsize = bgq_weyl_section_offset(sec_comm_end) - bgq_weyl_section_offset(sec_comm);
	size_t collapsedsize = bgq_weyl_section_offset(sec_collapsed_end) - bgq_weyl_section_offset(sec_collapsed);
	size_t result = -1;
	if (g_bgq_sec_comm <= (uint8_t*)ptr && (uint8_t*)ptr < g_bgq_sec_comm + commsize) {
		result = (uint8_t*)ptr - g_bgq_sec_comm + bgq_weyl_section_offset(sec_comm);
	} else {
		for (size_t i = 0; i < g_bgq_spinorfields_count; i+=1) {
			bgq_weylfield_controlblock *field = &g_bgq_spinorfields[i];
			if (!field->isInitinialized)
				continue;
			if (!field->sec_weyl)
				continue;

			if ((uint8_t*)field->sec_collapsed <= (uint8_t*)ptr && (uint8_t*)ptr < (uint8_t*)field->sec_collapsed + collapsedsize) {
				result = (uint8_t*)ptr - (uint8_t*)field->sec_collapsed + bgq_weyl_section_offset(sec_collapsed);
			}
		}
	}

	assert(result%sizeof(bgq_weyl_vec)==0);
assert(result!=-1);
	return result;
}
