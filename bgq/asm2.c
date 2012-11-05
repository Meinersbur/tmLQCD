
#include "bgq_HoppingMatrix.h"
#include "bgq_comm.h"

static inline bgq_HoppingMatrix_worker_datamove_recvxyz(bgq_weylfield_controlblock *spinorfield, bool isOdd, size_t beginj, size_t endj) {
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

		bgq_su3_weyl_decl(weyl);
		bgq_su3_weyl_load_double(weyl, weyladdr_src);
		bgq_weylqpx_expect(weyl, t1, t2, x, y, z, d, false);

		//bgq_weylvec_expect(*weyladdr_dst, t1,t2,x,y,z,d,false);
		bgq_su3_weyl_store_double(weyladdr_dst, weyl);
		bgq_weylvec_written(weyladdr_dst, t1, t2, x, y, z, d, false);
	}
}


static inline bgq_HoppingMatrix_worker_datamove_recvtup(bgq_weylfield_controlblock *spinorfield, bool isOdd, size_t beginj, size_t endj) {
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
		//TODO: Prefetch
		//TODO: Inline assembler
		bgq_weyl_vec *weyladdr_left = &g_bgq_sec_send[TDOWN][j];
		bgq_weyl_vec *weyladdr_right = &g_bgq_sec_recv[TUP][j];
		bgq_weyl_vec *weyladdr_dst = spinorfield->consptr_recvtup[j];

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

		//bgq_weylvec_expect(*weyladdr_dst, t1,t2,x,y,z,d,false);
		size_t offset = bgq_pointer2offset(spinorfield, weyladdr_dst);
		ucoord index = bgq_offset2index(offset);
		bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
		assert(sec==sec_body || sec==sec_surface);
		ucoord ic_check = g_bgq_index2collapsed[isOdd][index];
		assert(ic == ic_check);
		bgq_su3_weyl_store_double(weyladdr_dst, weyl);
		bgq_weylvec_written(weyladdr_dst, t1, t2, x, y, z, d, false);
	}
}

static inline bgq_HoppingMatrix_worker_datamove_recvtdown(bgq_weylfield_controlblock *spinorfield, bool isOdd, size_t beginj, size_t endj) {
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

		//bgq_weylvec_expect(*weyladdr_dst, t1,t2,x,y,z,d,false);
		size_t offset = bgq_pointer2offset(spinorfield, weyladdr_dst);
		ucoord index = bgq_offset2index(offset);
		bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
		assert(sec==sec_body || sec==sec_surface);
		ucoord ic_check = g_bgq_index2collapsed[isOdd][index];
		assert(ic == ic_check);
		bgq_su3_weyl_store_double(weyladdr_dst, weyl);
		bgq_weylvec_written(weyladdr_dst, t1, t2, x, y, z, d, false);
	}
}

void bgq_HoppingMatrix_worker_datamove(void *argptr, size_t tid, size_t threads) {
	bgq_weylfield_controlblock *spinorfield = argptr;
	bool isOdd = spinorfield->isOdd;

	//const size_t workload_recvt = COMM_T ? 2*LOCAL_HALO_T/PHYSICAL_LP : LOCAL_HALO_T/PHYSICAL_LP;
	const size_t workload_recv_tup = LOCAL_HALO_T/PHYSICAL_LP;
	const size_t workload_recv_tdown = workload_recv_tup;
	const size_t workload_recv = 2*PHYSICAL_HALO_X + 2*PHYSICAL_HALO_Y + 2*PHYSICAL_HALO_Z;
	const size_t workload = workload_recv + 2*workload_recv_tup + 2*workload_recv_tdown;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min(workload, begin+threadload);
	for (size_t i = begin; i<end; ) {
		WORKLOAD_DECL(i,workload);

		if (WORKLOAD_SPLIT(workload_recv)) {
			// Do other dimensions
			size_t beginj = WORKLOAD_PARAM(workload_recv);
			size_t endj = min_sizet(workload_recv,beginj+threadload);
			bgq_HoppingMatrix_worker_datamove_recvxyz(spinorfield,isOdd,beginj,endj);
			i += (endj - beginj);
		} else if (WORKLOAD_SPLIT(2*workload_recv_tup)) {
			// Count an T-iteration twice; better have few underloaded threads (so remaining SMT-threads have some more resources) then few overloaded threads (so the master thread has to wait for them)
			size_t twobeginj = WORKLOAD_PARAM(2*workload_recv_tup);
			size_t twoendj = min_sizet(2*workload_recv_tup,twobeginj+threadload);
			size_t beginj = twobeginj / 2;
			size_t endj = twoendj / 2;
			bgq_HoppingMatrix_worker_datamove_recvtup(spinorfield,isOdd,beginj,endj);
			i += (twoendj - twobeginj);
		} else if (WORKLOAD_SPLIT(2*workload_recv_tdown)) {
			// Count an T-iteration twice; better have few underloaded threads (so remaining SMT-threads have some more resources) then few overloaded threads (so the master thread has to wait for them)
			size_t twobeginj = WORKLOAD_PARAM(2*workload_recv_tdown);
			size_t twoendj = min_sizet(2*workload_recv_tdown,twobeginj+threadload);
			size_t beginj = twobeginj/2;
			size_t endj = twoendj / 2;
			bgq_HoppingMatrix_worker_datamove_recvtdown(spinorfield,isOdd,beginj,endj);
			i += (twoendj - twobeginj);
		} else {
		UNREACHABLE
		}

		WORKLOAD_CHECK
	}
}


