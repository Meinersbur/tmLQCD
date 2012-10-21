/*
 * bgq_field.c
 *
 *  Created on: Oct 13, 2012
 *      Author: meinersbur
 */


#define BGQ_FIELD_C_
#include "bgq_field.h"

#include "bgq_HoppingMatrix.h"
#include "bgq_dispatch.h"
#include "bgq_qpx.h"


static bgq_weylfield_section bgq_section_commbuftran(bgq_weylfield_section sec, bool isSendSide) {
	switch (sec) {
	case sec_send_tup:
	case sec_recv_tdown:
		return isSendSide ? sec_send_tup : sec_recv_tdown;
	case sec_send_tdown:
	case sec_recv_tup:
		return isSendSide ? sec_send_tdown : sec_recv_tup;
	case sec_send_xup:
	case sec_recv_xdown:
		return isSendSide ? sec_send_xup : sec_recv_xdown;
	case sec_send_xdown:
	case sec_recv_xup:
		return isSendSide ? sec_send_xdown : sec_recv_xup;
	case sec_send_yup:
	case sec_recv_ydown:
		return isSendSide ? sec_send_yup : sec_recv_ydown;
	case sec_send_ydown:
	case sec_recv_yup:
		return isSendSide ? sec_send_ydown : sec_recv_yup;
	case sec_send_zup:
	case sec_recv_zdown:
		return isSendSide ? sec_send_zup : sec_recv_zdown;
	case sec_send_zdown:
	case sec_recv_zup:
		return isSendSide ? sec_send_zdown : sec_recv_zup;
	default:
		// No communication to other node, stay in buffer
		return sec;
	}
}



static bgq_weylfield_section bgq_HoppingMatrix_init_source_sectionof_local(size_t t, size_t x, size_t y, size_t z, bgq_direction d) {
	if (COMM_T && (t==LOCAL_LT-1) && (d==TUP)) {
		return sec_recv_tup;
	}
	if (COMM_T && (t==0) && (d==TDOWN)) {
		return sec_recv_tdown;
	}
	if (COMM_X && (x==LOCAL_LX-1) && (d==XUP)) {
		return sec_recv_xup;
	}
	if (COMM_X && (x==0) && (d==XDOWN)) {
		return sec_recv_xdown;
	}
	if (COMM_Y && (y==LOCAL_LY-1) && (d==YUP)) {
		return sec_recv_yup;
	}
	if (COMM_Y && (y==0) && (d==YDOWN)) {
		return sec_recv_ydown;
	}
	if (COMM_Z && (z==LOCAL_LZ-1) && (d==ZUP)) {
		return sec_recv_zup;
	}
	if (COMM_Z && (z==0) && (d==ZDOWN)) {
		return sec_recv_zdown;
	}

	bool isSurface = bgq_local2isSurface(t, x, y, z);
	return isSurface ? sec_surface : sec_body;
}


static bgq_weylfield_section bgq_HoppingMatrix_init_source_sectionof_physical(bool isOdd, size_t tv, size_t x, size_t y, size_t z, bgq_physical_direction pd) {
	bgq_direction d = bgq_physical2direction(pd);
	switch (pd) {
	case P_TDOWN1:
	case P_TUP1:{
		size_t t1 = bgq_physical2t1(isOdd, tv, x,y,z);
		return bgq_HoppingMatrix_init_source_sectionof_local(t1,x,y,z,d);
	}
	case P_TDOWN2:
	case P_TUP2:{
		size_t t2 = bgq_physical2t2(isOdd, tv, x,y,z);
		return bgq_HoppingMatrix_init_source_sectionof_local(t2,x,y,z,d);
	}
	default:{
		size_t t1 = bgq_physical2t1(isOdd, tv, x,y,z);
		bgq_weylfield_section sec1 = bgq_HoppingMatrix_init_source_sectionof_local(t1,x,y,z,d);
		size_t t2 = bgq_physical2t2(isOdd, tv, x,y,z);
		bgq_weylfield_section sec2 = bgq_HoppingMatrix_init_source_sectionof_local(t2,x,y,z,d);
		assert(sec1 == sec2);
		return sec1;
	}
	}
}


void bgq_indices_init() {
	if (g_bgq_indices_initialized)
		return;
	g_bgq_indices_initialized = true; // Take care for uses within this function itself


	size_t body_volume;
	if ((PHYSICAL_LTV<=1) || (PHYSICAL_LX<=2) || (PHYSICAL_LY<=2) || (PHYSICAL_LZ<=2)) {
		body_volume = 0;
	} else {
		size_t body_tv = COMM_T ? (PHYSICAL_LTV-1) : PHYSICAL_LTV;
		size_t body_x = COMM_X ? (PHYSICAL_LX-2) : PHYSICAL_LX;
		size_t body_y = COMM_Y ? (PHYSICAL_LY-2) : PHYSICAL_LY;
		size_t body_z = COMM_Z ? (PHYSICAL_LZ-2) : PHYSICAL_LZ;
		body_volume = body_tv * body_x * body_y * body_z;
	}
	PHYSICAL_BODY = body_volume;


	for (size_t isOdd = false; isOdd <= true; isOdd+=1) {
		g_bgq_index_surface2halfvolume[isOdd] = malloc(PHYSICAL_SURFACE * sizeof(*g_bgq_index_surface2halfvolume[isOdd]));
		g_bgq_index_body2halfvolume[isOdd] = malloc(PHYSICAL_BODY * sizeof(*g_bgq_index_body2halfvolume[isOdd]));
		g_bgq_index_halfvolume2surface[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_index_halfvolume2surface[isOdd]));
		g_bgq_index_halfvolume2body[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_index_halfvolume2body[isOdd]));
		g_bgq_index_halfvolume2surfacebody[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_index_halfvolume2surfacebody[isOdd]));

		g_bgq_offset_fromHalfvolume[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_offset_fromHalfvolume[isOdd]));
		g_bgq_offset_fromSurface[isOdd] = malloc(PHYSICAL_SURFACE * sizeof(*g_bgq_offset_fromSurface[isOdd]));
		g_bgq_offset_fromBody[isOdd] = malloc(PHYSICAL_BODY * sizeof(*g_bgq_offset_fromBody[isOdd]));

		g_bgq_destoffset_fromHalfvolume[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_destoffset_fromHalfvolume[isOdd]));
		g_bgq_destoffset_fromSurface[isOdd] = malloc(PHYSICAL_SURFACE * sizeof(*g_bgq_destoffset_fromSurface[isOdd]));
		g_bgq_destoffset_fromBody[isOdd] = malloc(PHYSICAL_BODY * sizeof(*g_bgq_destoffset_fromBody[isOdd]));
	}


	// Setup array for index translation
	// It would be difficult to find an explicit expression for indices that are only on surface/body
	// Thus, we simply iterate over all locations and allocate locations in order
	for (size_t isOdd = false; isOdd <= true; isOdd+=1) {
		size_t nextIndexSurface = 0;
		size_t nextIndexBody = 0;

		// Iteration order is effectively chosen here
		for (size_t ih = 0; ih < PHYSICAL_VOLUME; ih += 1) {
			size_t t1 = bgq_halfvolume2t1(isOdd, ih);
			size_t t2 = bgq_halfvolume2t2(isOdd, ih);
			size_t x = bgq_halfvolume2x(ih);
			size_t y = bgq_halfvolume2y(ih);
			size_t z = bgq_halfvolume2z(ih);
			assert(ih == bgq_local2halfvolume(t1,x,y,z));
			assert(ih == bgq_local2halfvolume(t2,x,y,z));
			assert(isOdd == bgq_local2isOdd(t1,x,y,z));
			assert(isOdd == bgq_local2isOdd(t2,x,y,z));

			bool isSurface = bgq_local2isSurface(t1, x, y, z) || bgq_local2isSurface(t2, x, y, z);
			if (isSurface) {
				size_t is = nextIndexSurface;
				assert(0 <= is && is < PHYSICAL_SURFACE);
				g_bgq_index_surface2halfvolume[isOdd][is] = ih;
				g_bgq_index_halfvolume2surfacebody[isOdd][ih] = is;
				g_bgq_index_halfvolume2surface[isOdd][ih] = is;
				g_bgq_index_halfvolume2body[isOdd][ih] = (size_t)-1;
				nextIndexSurface += 1;
			} else {
				size_t ib = nextIndexBody;
				assert(0 <= ib && ib < PHYSICAL_BODY);
				g_bgq_index_body2halfvolume[isOdd][ib] = ih;
				g_bgq_index_halfvolume2surfacebody[isOdd][ih] = ib;
				g_bgq_index_halfvolume2surface[isOdd][ih] = (size_t)-1;
				g_bgq_index_halfvolume2body[isOdd][ih] = ib;
				nextIndexBody += 1;
			}
		}
	}

#ifndef NDEBUG
	// Check consistency
	for (size_t isOdd = false; isOdd <= true; isOdd+=1) {
		for (size_t is2 = 0; is2 < PHYSICAL_SURFACE; is2 += 1) {
			size_t ih = bgq_surface2halfvolume(isOdd,is2);
			assert(bgq_halfvolume2isSurface(isOdd,ih));

			assert(is2 == bgq_halfvolume2surface(isOdd,ih));
			assert(ih == bgq_surface2halfvolume(isOdd,is2));
		}
	}
#endif


	// Setup mapping of weyl to some memory offset
	// (where to read a datum for hoppingmatrix)
	for (size_t isOdd = false; isOdd <= true; isOdd+=1) {
		size_t nextoffset[sec_end];
		for (bgq_weylfield_section sec = 0; sec < sec_end; sec+=1) {
			nextoffset[sec] = bgq_weyl_section_offset(sec);
		}

		for (size_t ih = 0; ih < PHYSICAL_VOLUME; ih+=1) {
			size_t tv = bgq_halfvolume2tv(ih);
			size_t t1 = bgq_halfvolume2t1(isOdd, ih);
			size_t t2 = bgq_halfvolume2t2(isOdd, ih);
			size_t x = bgq_halfvolume2x(ih);
			size_t y = bgq_halfvolume2y(ih);
			size_t z = bgq_halfvolume2z(ih);
			bool isSurface = bgq_local2isSurface(t1,x,y,z) || bgq_local2isSurface(t2,x,y,z);
			assert(isSurface == (g_bgq_index_halfvolume2surface[isOdd][ih]!=(size_t)-1));
			assert(isSurface == (g_bgq_index_halfvolume2body[isOdd][ih]==(size_t)-1));
			bgq_weylfield_section mainsec = isSurface ? sec_surface : sec_body;

			if (tv == PHYSICAL_LTV-1) {
				int a = 0;
			}

			for (size_t pd_cnt = P_TUP1; pd_cnt <= P_ZDOWN; pd_cnt+=1) {
				bgq_physical_direction pd = pd_cnt;
				bgq_direction d = bgq_physical2direction(pd);
				bgq_dimension dim = bgq_direction2dimension(d);
				bgq_weylfield_section sec = bgq_HoppingMatrix_init_source_sectionof_physical(isOdd, tv, x, y, z, pd);
				size_t eltsize = (dim==DIM_T) ? sizeof(bgq_weyl_nonvec) : sizeof(bgq_weyl_vec);

				// Reserve some offset in surface/body to ensure consecutive layout
				size_t thisOffset = nextoffset[mainsec];
				assert(bgq_sectionOfOffset(thisOffset) == mainsec);
				//assert(thisOffset == bgq_weyl_section_offset(sec_surface) + g_bgq_index_halfvolume2surfacebody[isOdd][ih]*sizeof(bgq_weylsite) + bgq_offsetof_weylsite[pd]);
				nextoffset[mainsec] += eltsize;

				if (sec != mainsec) {
					assert((sec!=sec_surface) && (sec!=sec_body));
					// If in one of the mpi buffers, also reserve some space there
					thisOffset = nextoffset[sec];
					nextoffset[sec] += eltsize;
				}
				assert(bgq_sectionOfOffset(thisOffset) == sec);
				g_bgq_offset_fromHalfvolume[isOdd][ih].pd[pd] = bgq_encode_offset(thisOffset);
				if (isSurface) {
					size_t is = bgq_halfvolume2surface(isOdd,ih);
					g_bgq_offset_fromSurface[isOdd][is].pd[pd] = bgq_encode_offset(thisOffset);
				} else {
					size_t ib = bgq_halfvolume2body(isOdd,ih);
					g_bgq_offset_fromBody[isOdd][ib].pd[pd] = bgq_encode_offset(thisOffset);
				}
			}
		}
	}

#ifndef NDEBUG
	// Check consistency
	for (size_t isOdd = false; isOdd <= true; isOdd+=1) {
		for (size_t is = 0; is < PHYSICAL_SURFACE; is+=1) {
			for (size_t pd_cnt = P_TUP1; pd_cnt <= P_ZDOWN; pd_cnt+=1) {
				bgq_physical_direction pd = pd_cnt;
				size_t offset_is = g_bgq_offset_fromSurface[isOdd][is].pd[pd];
				assert(offset_is+1); // For valgrind to check initialization

				size_t ih = bgq_surface2halfvolume(isOdd,is);
				size_t offset_ih = g_bgq_offset_fromHalfvolume[isOdd][ih].pd[pd];
				assert(offset_ih+1);

				assert(offset_is == offset_ih);
			}
		}
	}
#endif


	// Determine dest locations (where to write a weyl before communication)
	for (size_t isOdd_cnt = false; isOdd_cnt <= true; isOdd_cnt+=1) {
		bool isOdd_src = isOdd_cnt;
		bool isOdd_dst = !isOdd_cnt;
		for (size_t ih_src = 0; ih_src < PHYSICAL_VOLUME; ih_src+=1) {
			const size_t tv_src = bgq_halfvolume2tv(ih_src);
			const size_t t1_src = bgq_halfvolume2t1(isOdd_src, ih_src);
			const size_t t2_src = bgq_halfvolume2t2(isOdd_src, ih_src);
			const size_t x_src = bgq_halfvolume2x(ih_src);
			const size_t y_src = bgq_halfvolume2y(ih_src);
			const size_t z_src = bgq_halfvolume2z(ih_src);

			const bool isSurface = bgq_physical2isSurface(isOdd_src,tv_src,x_src,y_src,z_src);
			const bool isSecondInVector_src = bgq_physical2eo(isOdd_src, tv_src,x_src,y_src,z_src);
			assert(!isSecondInVector_src == (t1_src==(tv_src*PHYSICAL_LK*PHYSICAL_LP+2)%LOCAL_LT));
			assert(!isSecondInVector_src == (t2_src==(tv_src*PHYSICAL_LK*PHYSICAL_LP+2+2)%LOCAL_LT));

			for (size_t pd_cnt = P_TUP1; pd_cnt <= P_ZDOWN; pd_cnt+=1) {
				const bgq_physical_direction pd_src = pd_cnt;
				const bgq_direction d_src = bgq_physical2direction(pd_src);
				//const bgq_direction d_dst = bgq_direction_revert(d_src);
				//const bgq_dimension dim = bgq_direction2dimension(d_src);
				const size_t t_src = ((pd_src==P_TUP2) || (pd_src==P_TDOWN2)) ? t2_src : t1_src;

				size_t tv_dst = tv_src;
				//size_t t1_dst = t1_src;
				//size_t t2_dst = t2_src;
				size_t t_dst = t_src;
				size_t x_dst = x_src;
				size_t y_dst = y_src;
				size_t z_dst = z_src;
				bgq_physical_direction pd_dst;

				switch (pd_src) {
				case P_TDOWN1:
					t_dst = (t1_src + LOCAL_LT - 1) % LOCAL_LT;
					tv_dst = bgq_local2tv(t_dst,x_dst,y_dst,z_dst);
					pd_dst = !isSecondInVector_src ? P_TUP2 : P_TUP1;
					assert(tv_dst==bgq_local2tv(t_dst,x_dst,y_dst,z_dst));
					break;
				case P_TUP1:
					t_dst = (t1_src + 1) % LOCAL_LT;
					tv_dst = bgq_local2tv(t_dst,x_dst,y_dst,z_dst);
					pd_dst = !isSecondInVector_src ? P_TDOWN1 : P_TDOWN2;
					break;
				case P_TDOWN2:
					t_dst = (t2_src + LOCAL_LT - 1) % LOCAL_LT;
					tv_dst = bgq_local2tv(t_dst,x_dst,y_dst,z_dst);
					pd_dst = !isSecondInVector_src ? P_TUP1 : P_TUP2;
					assert(tv_dst==bgq_local2tv(t_dst,x_dst,y_dst,z_dst));
					break;
				case P_TUP2:
					t_dst = (t2_src + 1) % LOCAL_LT;
					tv_dst = bgq_local2tv(t_dst,x_dst,y_dst,z_dst);
					pd_dst = !isSecondInVector_src ? P_TDOWN2 : P_TDOWN1;
					assert(tv_dst==bgq_local2tv(t_dst,x_dst,y_dst,z_dst));
					break;
				case P_XDOWN:
					x_dst = (x_src + LOCAL_LX - 1) % LOCAL_LX;
					pd_dst = P_XUP;
					break;
				case P_XUP:
					x_dst = (x_src + 1) % LOCAL_LX;
					pd_dst = P_XDOWN;
					break;
				case P_YDOWN:
					y_dst = (y_src + LOCAL_LY - 1) % LOCAL_LY;
					pd_dst = P_YUP;
					break;
				case P_YUP:
					y_dst = (y_src + 1) % LOCAL_LY;
					pd_dst = P_YDOWN;
					break;
				case P_ZDOWN:
					z_dst = (z_src + LOCAL_LZ - 1) % LOCAL_LZ;
					pd_dst = P_ZUP;
					break;
				case P_ZUP:
					z_dst = (z_src + 1) % LOCAL_LZ;
					pd_dst = P_ZDOWN;
					break;
				}

				bgq_weylfield_section sec_dst = bgq_HoppingMatrix_init_source_sectionof_physical(isOdd_dst, tv_dst, x_dst, y_dst, z_dst, pd_dst); // dst is going to read the weyl
				bgq_weylfield_section sec_src = bgq_section_commbuftran(sec_dst, true); // src finds the correct location in preparation for dst

				// Find the offset at which the other node expects the datum
				size_t ih_dst = bgq_physical2halfvolume(tv_dst, x_dst, y_dst, z_dst);
				size_t offset_dst = bgq_decode_offset(g_bgq_offset_fromHalfvolume[isOdd_dst][ih_dst].pd[pd_dst]);
				assert(bgq_sectionOfOffset(offset_dst) == sec_dst);

				size_t offset_src;
				if ((sec_dst==sec_body) || (sec_dst==sec_surface)) {
					// No communication, do everything locally
					assert((sec_src==sec_body) || (sec_src==sec_surface));
					offset_src = offset_dst;
				} else {
					// Communication with MPI buffers, translate recv- to send-buffer
					assert((sec_src!=sec_body) && (sec_src!=sec_surface));
					size_t relOffset = offset_dst - bgq_weyl_section_offset(sec_dst)/*recv*/;
					offset_src = bgq_weyl_section_offset(sec_src)/*send*/ + relOffset;
				}

				assert(0 <= offset_src && offset_src < bgq_weyl_section_offset(sec_end));
				g_bgq_destoffset_fromHalfvolume[isOdd_src][ih_src].pd[pd_src] = bgq_encode_offset(offset_src);
				if (isSurface) {
					size_t is_src = bgq_halfvolume2surface(isOdd_src, ih_src);
					g_bgq_destoffset_fromSurface[isOdd_src][is_src].pd[pd_src] = bgq_encode_offset(offset_src);
				} else {
					size_t ib_src = bgq_halfvolume2body(isOdd_src, ih_src);
					g_bgq_destoffset_fromBody[isOdd_src][ib_src].pd[pd_src] = bgq_encode_offset(offset_src);
				}
			}
		}
	}
}


uint8_t *g_bgq_spinorfields_data = NULL;

void bgq_spinorfields_init(size_t std_count, size_t chi_count) {
	size_t tot_count = std_count + chi_count;

	size_t field_datasize = bgq_weyl_section_offset(sec_end);
	field_datasize = (field_datasize + 127) & ~(size_t)127;
	field_datasize += PHYSICAL_VOLUME * sizeof(bgq_spinorsite); // Fullspinor field
	field_datasize = (field_datasize + 127) & ~(size_t)127;
	field_datasize += 2*PHYSICAL_VOLUME * sizeof(size_t); // destptrs
	field_datasize += PHYSICAL_HALO * sizeof(void*); // destptr from halo
	field_datasize = (field_datasize + 127) & ~(size_t)127; // For padding the next field

	g_bgq_spinorfields = malloc(tot_count * sizeof(bgq_weylsite));
	//memset(g_bgq_spinorfields, -1, tot_count * sizeof(bgq_weylsite));
	g_bgq_spinorfields_data = malloc_aligned(tot_count * field_datasize, 128);
	//memset(g_bgq_spinorfields_data, -1, tot_count * field_datasize);

	for (size_t i = 0; i < tot_count; i += 1) {
		g_bgq_spinorfields[i].isInitinialized = false;
		uint8_t *weylbase = &g_bgq_spinorfields_data[i * field_datasize];
		g_bgq_spinorfields[i].sec_weyl = weylbase;
		g_bgq_spinorfields[i].sec_send_tup = (bgq_weyl_nonvec*)(weylbase + bgq_weyl_section_offset(sec_send_tup));
		g_bgq_spinorfields[i].sec_send_tdown = (bgq_weyl_nonvec*)(weylbase + bgq_weyl_section_offset(sec_send_tdown));
		g_bgq_spinorfields[i].sec_send_xup = (bgq_weyl_vec*)(weylbase + bgq_weyl_section_offset(sec_send_xup));
		g_bgq_spinorfields[i].sec_send_xdown = (bgq_weyl_vec*)(weylbase + bgq_weyl_section_offset(sec_send_xdown));
		g_bgq_spinorfields[i].sec_send_yup = (bgq_weyl_vec*)(weylbase + bgq_weyl_section_offset(sec_send_yup));
		g_bgq_spinorfields[i].sec_send_ydown = (bgq_weyl_vec*)(weylbase + bgq_weyl_section_offset(sec_send_ydown));
		g_bgq_spinorfields[i].sec_send_zup = (bgq_weyl_vec*)(weylbase + bgq_weyl_section_offset(sec_send_zup));
		g_bgq_spinorfields[i].sec_send_zdown = (bgq_weyl_vec*)(weylbase + bgq_weyl_section_offset(sec_send_zdown));
		g_bgq_spinorfields[i].sec_recv_tup = (bgq_weyl_nonvec*)(weylbase + bgq_weyl_section_offset(sec_recv_tup));
		g_bgq_spinorfields[i].sec_recv_tdown = (bgq_weyl_nonvec*)(weylbase + bgq_weyl_section_offset(sec_recv_tdown));
		g_bgq_spinorfields[i].sec_recv_xup = (bgq_weyl_vec*)(weylbase + bgq_weyl_section_offset(sec_recv_xup));
		g_bgq_spinorfields[i].sec_recv_xdown = (bgq_weyl_vec*)(weylbase + bgq_weyl_section_offset(sec_recv_xdown));
		g_bgq_spinorfields[i].sec_recv_yup = (bgq_weyl_vec*)(weylbase + bgq_weyl_section_offset(sec_recv_yup));
		g_bgq_spinorfields[i].sec_recv_ydown = (bgq_weyl_vec*)(weylbase + bgq_weyl_section_offset(sec_recv_ydown));
		g_bgq_spinorfields[i].sec_recv_zup = (bgq_weyl_vec*)(weylbase + bgq_weyl_section_offset(sec_recv_zup));
		g_bgq_spinorfields[i].sec_recv_zdown = (bgq_weyl_vec*)(weylbase + bgq_weyl_section_offset(sec_recv_zdown));
		g_bgq_spinorfields[i].sec_surface = (bgq_weylsite*)(weylbase + bgq_weyl_section_offset(sec_surface));
		g_bgq_spinorfields[i].sec_body = (bgq_weylsite*)(weylbase + bgq_weyl_section_offset(sec_body));
		g_bgq_spinorfields[i].sec_end = weylbase + bgq_weyl_section_offset(sec_end);

		uint8_t *fullbase = weylbase + bgq_weyl_section_offset(sec_end);
		fullbase = ADD_PADDING(fullbase,BGQ_ALIGNMENT_L2);
		g_bgq_spinorfields[i].sec_fullspinor = (bgq_spinorsite*)fullbase;

		uint8_t *ptrBase = fullbase + PHYSICAL_VOLUME * sizeof(bgq_spinorsite);
		ptrBase = ADD_PADDING(ptrBase,BGQ_ALIGNMENT_L2);
		g_bgq_spinorfields[i].destptrFromHalfvolume = (bgq_weyl_ptr_t*)ptrBase;


		uint8_t *ptrSurfaceBase = ptrBase + PHYSICAL_VOLUME * sizeof(*g_bgq_spinorfields[i].destptrFromHalfvolume);
		g_bgq_spinorfields[i].destptrFromSurface = (bgq_weyl_ptr_t*)ptrSurfaceBase;


		uint8_t *ptrBodyBase = ptrSurfaceBase + PHYSICAL_SURFACE * sizeof(*g_bgq_spinorfields[i].destptrFromSurface);
		g_bgq_spinorfields[i].destptrFromBody = (bgq_weyl_ptr_t*)ptrBodyBase;

		uint8_t *ptrHaloBase = ptrBodyBase + PHYSICAL_BODY * sizeof(bgq_weyl_ptr_t);
		uint8_t *ptrHalo = ptrHaloBase;
		if (COMM_T) {
			g_bgq_spinorfields[i].destptrFromRecvTup = (bgq_weyl_nonvec**)ptrHalo;
			ptrHalo += LOCAL_LT * sizeof(bgq_weyl_nonvec*);
			g_bgq_spinorfields[i].destptrFromRecvTdown = (bgq_weyl_nonvec**)ptrHalo;
			ptrHalo += LOCAL_LT * sizeof(bgq_weyl_nonvec*);
		}
		if (COMM_X) {
			g_bgq_spinorfields[i].destptrFromRecvXup = (bgq_weyl_vec**)ptrHalo;
			ptrHalo += LOCAL_LX * sizeof(bgq_weyl_vec*);
			g_bgq_spinorfields[i].destptrFromRecvXdown = (bgq_weyl_vec**)ptrHalo;
			ptrHalo += LOCAL_LX * sizeof(bgq_weyl_vec*);
		}
		if (COMM_Y) {
			g_bgq_spinorfields[i].destptrFromRecvYup = (bgq_weyl_vec**)ptrHalo;
			ptrHalo += LOCAL_LX * sizeof(bgq_weyl_vec*);
			g_bgq_spinorfields[i].destptrFromRecvYdown = (bgq_weyl_vec**)ptrHalo;
			ptrHalo += LOCAL_LX * sizeof(bgq_weyl_vec*);
		}
		if (COMM_Z) {
			g_bgq_spinorfields[i].destptrFromRecvZup = (bgq_weyl_vec**)ptrHalo;
			ptrHalo += LOCAL_LX * sizeof(bgq_weyl_vec*);
			g_bgq_spinorfields[i].destptrFromRecvZdown = (bgq_weyl_vec**)ptrHalo;
			ptrHalo += LOCAL_LX * sizeof(bgq_weyl_vec*);
		}
	}
}


void bgq_gaugefield_init() {
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		g_bgq_gaugefield_fromHalfvolume[isOdd] = malloc_aligned(PHYSICAL_VOLUME * sizeof(*g_bgq_gaugefield_fromHalfvolume), BGQ_ALIGNMENT_L2);
		g_bgq_gaugefield_fromSurface[isOdd] = malloc_aligned(PHYSICAL_VOLUME * sizeof(*g_bgq_gaugefield_fromSurface), BGQ_ALIGNMENT_L2);
		g_bgq_gaugefield_fromBody[isOdd] = malloc_aligned(PHYSICAL_VOLUME * sizeof(*g_bgq_gaugefield_fromBody), BGQ_ALIGNMENT_L2);
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

	bgq_conversion_args args = {
			.field = field
	};
	bgq_master_call(&bgq_spinorfield_fullspinor2weyl_worker, &args); field->hasWeylfieldData = true;
}


void bgq_spinorfield_setup(bgq_weylfield_controlblock *field, bool isOdd, bool readFullspinor, bool writeFullspinor, bool readWeyl, bool writeWeyl) {
	assert(field);
	assert(readFullspinor || writeFullspinor || readWeyl || writeWeyl); // Do something

	bool initDestPtrs;
	bool fullspinorAvailable;
	bool weylAvailable;
	if (field->isInitinialized) {
		if (field->isOdd == isOdd) {
			initDestPtrs = false;
			fullspinorAvailable = field->hasFullspinorData;
			weylAvailable = field->hasWeylfieldData;
		} else {
			// Warning: field has already initialized for a different oddness
			// For optimal performance, always reuse fields with same oddness
			// TODO: Even and Odd fields can be made similar enough such that the pointers are the same? (e.g. Mirror T-Dimension)
			master_print("PERFORMANCE WARNING: Performance loss by reuse of spinorfield with different oddness\n");
			initDestPtrs = true;
			fullspinorAvailable = false;
			weylAvailable = false;
		}
	} else {
		initDestPtrs = true;
		fullspinorAvailable = false;
		weylAvailable = false;
	}

	if (initDestPtrs) {
		// For phase 1. (distribute)
		uint8_t *weylbase = field->sec_weyl;
		for (size_t ih = 0; ih < PHYSICAL_VOLUME; ih += 1) {
			for (size_t pd = 0; pd < P_COUNT; pd += 1) {
				uint8_t *ptr = &weylbase[bgq_decode_offset(g_bgq_destoffset_fromHalfvolume[isOdd][ih].pd[pd])];
				assert(field->sec_weyl <= ptr && ptr < field->sec_end);
				field->destptrFromHalfvolume[ih].pd[pd] = ptr;
			}
		}
		for (size_t is = 0; is < PHYSICAL_SURFACE; is += 1) {
			for (size_t pd = 0; pd < P_COUNT; pd += 1) {
				field->destptrFromSurface[is].pd[pd] = &weylbase[bgq_decode_offset(g_bgq_destoffset_fromSurface[isOdd][is].pd[pd])];
			}
		}
		for (size_t ib = 0; ib < PHYSICAL_BODY; ib += 1) {
			for (size_t pd = 0; pd < P_COUNT; pd += 1) {
				field->destptrFromBody[ib].pd[pd] = &weylbase[bgq_decode_offset(g_bgq_destoffset_fromBody[isOdd][ib].pd[pd])];
			}
		}

		// For 5th phase (datamove)
		memset(field->destptrFromRecvTup, 0, PHYSICAL_HALO_T * sizeof(bgq_weyl_nonvec*));
		for (size_t is = 0; is < PHYSICAL_SURFACE; is += 1) {
			for (size_t pd_cnt = 0; pd_cnt < P_COUNT; pd_cnt += 1) {
				bgq_physical_direction pd = pd_cnt;
				bgq_direction d = bgq_physical2direction(pd);
				bgq_dimension dim = bgq_direction2dimension(d);
				size_t tv = bgq_surface2tv(isOdd, is);
				if (tv == PHYSICAL_LTV-1) {
					int a = 0;
				}

				size_t offset = bgq_decode_offset(g_bgq_offset_fromSurface[isOdd][is].pd[pd]); // Where does the data come from
				assert(offset+1); // For valgrind
				void *ptr_consecutive = bgq_weylsite_getdirection(&field->sec_surface[is], pd);
				bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
				if ((sec == sec_surface) || (sec == sec_body)) {
					assert(ptr_consecutive == &field->sec_weyl[offset]);
					continue; // If the data does not come from a communication buffer, no move required
				}

				size_t relOffset = offset - bgq_weyl_section_offset(sec);
				size_t eltSize = (dim == DIM_T) ? sizeof(bgq_weyl_nonvec) : sizeof(bgq_weyl_vec);
				assert(relOffset % eltSize == 0);
				size_t index = relOffset / eltSize;

				switch (sec) {
				case sec_recv_tup:
					assert(0 <= index && index < PHYSICAL_HALO_T);
					assert(field->destptrFromRecvTup[index] == 0);
					field->destptrFromRecvTup[index] = ptr_consecutive;
					break;
				case sec_recv_tdown:
					field->destptrFromRecvTdown[index] = ptr_consecutive;
					break;
				case sec_recv_xup:
					field->destptrFromRecvXup[index] = ptr_consecutive;
					break;
				case sec_recv_xdown:
					field->destptrFromRecvXdown[index] = ptr_consecutive;
					break;
				case sec_recv_yup:
					field->destptrFromRecvYup[index] = ptr_consecutive;
					break;
				case sec_recv_ydown:
					field->destptrFromRecvYdown[index] = ptr_consecutive;
					break;
				case sec_recv_zup:
					field->destptrFromRecvZup[index] = ptr_consecutive;
					break;
				case sec_recv_zdown:
					field->destptrFromRecvZdown[index] = ptr_consecutive;
					break;
				default:
					assert(!"Not necessary to move somewhere else");
					break;
				}
			}
		}

#ifndef NDEBUG
		// Consistency checks
		if (COMM_T) {
			for (size_t i = 0; i < LOCAL_HALO_T; i += 1) {
				bgq_weyl_nonvec *ptr = field->destptrFromRecvTup[i];
				assert(ptr+1); // For valgrind
				ptrdiff_t offset = (uint8_t*)ptr - field->sec_weyl;
				bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
				assert(sec == sec_surface || sec==sec_body);
			}
		}
#endif

		field->isInitinialized = true;
		field->isOdd = isOdd;
		field->isSloppy = false;
	}

	if (readFullspinor) {
		if (fullspinorAvailable) {
			// Everything's fine
		} else {
			master_error(1, "Try to read uninitialized spinorfield in weyl format");
		}
	}
	if (readWeyl) {
		if (weylAvailable) {
			// Everything's fine
		} else if (fullspinorAvailable) {
			// Do a conversion
			master_print("PERFORMANCE WARNING: Converting fullspinor layout to weyl layout; For best performance, this should have been done in the first place\n");
			bgq_spinorfield_fullspinor2weyl(field);
		} else {
			master_error(1, "Try to read uninitialized spinorfield in weyl format");
		}
	}

	if (writeWeyl || writeFullspinor) {
		field->hasWeylfieldData = writeWeyl;
		field->hasFullspinorData = writeFullspinor;
	}
}

typedef struct {
	double _Complex c[3];
} su3_vector_array64;

typedef struct {
	su3_vector_array64 v[4];
} spinor_array64;

void bgq_spinorfield_transfer(bool isOdd, bgq_weylfield_controlblock *targetfield, spinor* sourcefield) {
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

		size_t ih = bgq_local2halfvolume(t,x,y,z);
		size_t k = bgq_local2k(t,x,y,z);
		bgq_spinorsite *targetspinor = &targetfield->sec_fullspinor[ih];
		for (int v = 0; v < 4; v+=1) {
			for (int c = 0; c < 3; c+=1) {
				targetspinor->s[v][c][k] = sourcespinor->v[v].c[c];
			}
		}
	}
}


