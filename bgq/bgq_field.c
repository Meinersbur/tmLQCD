/*
 * bgq_field.c
 *
 *  Created on: Oct 13, 2012
 *      Author: meinersbur
 */


#define BGQ_FIELD_C_
#include "bgq_field.h"

#include "bgq_qpx.h"


static bgq_weylfield_section bgq_HoppingMatrix_init_source_sectionof(bool isOdd, size_t t, size_t x, size_t y, size_t z, bgq_direction d, bool isSend) {
	if (COMM_T && (t==LOCAL_LT-1) && (d==TUP)) {
		return isSend ? sec_send_tup : sec_recv_tup;
	}
	if (COMM_T && (t==0) && (d==TDOWN)) {
		return isSend ? sec_send_tdown : sec_recv_tdown;
	}
	if (COMM_X && (x==LOCAL_LX-1) && (d==XUP)) {
		return isSend ? sec_send_xup : sec_recv_xup;
	}
	if (COMM_X && (x==0) && (d==XDOWN)) {
		return isSend ? sec_send_xdown : sec_recv_xdown;
	}
	if (COMM_Y && (y==LOCAL_LY-1) && (d==YUP)) {
		return isSend ? sec_send_yup: sec_recv_yup;
	}
	if (COMM_Y && (y==0) && (d==YDOWN)) {
		return isSend ? sec_send_ydown : sec_recv_ydown;
	}
	if (COMM_Z && (z==LOCAL_LZ-1) && (d==ZUP)) {
		return isSend ? sec_send_zup : sec_recv_zup;
	}
	if (COMM_Z && (z==0) && (d==ZDOWN)) {
		return isSend ? sec_send_zdown : sec_recv_zdown;
	}

	bool isSurface = bgq_local2isSurface(t,x,y,z);
	return isSurface ? sec_surface : sec_body;
}





void bgq_indices_init() {
	if (g_bgq_indices_initialized)
		return;


	size_t body_volume;
	if ((PHYSICAL_LTV<=2) || (PHYSICAL_LX<=2) || (PHYSICAL_LY<=2) || (PHYSICAL_LZ<=2)) {
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


	// Setup mapping of weyl to some memory offset
	for (size_t isOdd = false; isOdd <= true; isOdd+=1) {
		size_t nextoffset[sec_end];
		for (bgq_weylfield_section sec = 0; sec < sec_end; sec+=1) {
			nextoffset[sec] = bgq_weyl_section_offset(sec);
		}

		for (size_t ih = 0; ih < PHYSICAL_VOLUME; ih+=1) {
			size_t t1 = bgq_halfvolume2t1(isOdd, ih);
			size_t t2 = bgq_halfvolume2t2(isOdd, ih);
			size_t x = bgq_halfvolume2x(ih);
			size_t y = bgq_halfvolume2y(ih);
			size_t z = bgq_halfvolume2z(ih);
			bool isSurface = bgq_local2isSurface(t1,x,y,z) || bgq_local2isSurface(t2,x,y,z);
			assert(isSurface == (g_bgq_index_halfvolume2surface[isOdd][ih]!=(size_t)-1));
			assert(isSurface == (g_bgq_index_halfvolume2body[isOdd][ih]==(size_t)-1));
			bgq_weylfield_section mainsec = isSurface ? sec_surface : sec_body;

			for (size_t pd = P_TUP1; pd <= P_ZDOWN; pd+=1) {
				bgq_direction d = bgq_physical2direction(d);
				bgq_dimension dim = bgq_direction2dimension(d);
				bgq_weylfield_section sec = bgq_HoppingMatrix_init_source_sectionof(isOdd, ((pd==P_TUP2) || (pd==P_TDOWN2)) ? t2 : t1, x, y, z, d, false);
				size_t eltsize = (dim==DIM_T) ? sizeof(bgq_weyl_nonvec) : sizeof(bgq_weyl_vec);

				// Reserve some offset in surface/body to ensure consecutive layout
				size_t thisOffset = nextoffset[mainsec];
				//assert(thisOffset == bgq_weyl_section_offset(sec_surface) + g_bgq_index_halfvolume2surfacebody[isOdd][ih]*sizeof(bgq_weylsite) + bgq_offsetof_weylsite[pd]);
				nextoffset[mainsec] += eltsize;

				if (sec != mainsec) {
					// If in one of the mpi buffers, also reserve some space there
					thisOffset = nextoffset[sec];
					nextoffset[sec] += eltsize;
				}
				g_bgq_offset_fromHalfvolume[isOdd][ih].pd[pd] = bgq_encode_offset(thisOffset);
				if (isSurface) {
					size_t is =  bgq_halfvolume2surface(isOdd,ih);
					g_bgq_offset_fromSurface[isOdd][is].pd[pd] = bgq_encode_offset(thisOffset);
				} else {
					size_t ib = bgq_halfvolume2body(isOdd,ih);
					g_bgq_offset_fromBody[isOdd][ib].pd[pd] = bgq_encode_offset(thisOffset);
				}
			}
		}
	}


	// Determine dest locations (where to write a weyl when calculated)
	for (size_t isOdd = false; isOdd <= true; isOdd+=1) {
		bool isOdd_src = !isOdd;
		bool isOdd_dst = isOdd;
		for (size_t ih_src = 0; ih_src < PHYSICAL_VOLUME; ih_src+=1) {
			const size_t tv_src = bgq_halfvolume2tv(ih_src);
			const size_t t1_src = bgq_halfvolume2t1(isOdd, ih_src);
			const size_t t2_src = bgq_halfvolume2t2(isOdd, ih_src);
			const size_t x_src = bgq_halfvolume2x(ih_src);
			const size_t y_src = bgq_halfvolume2y(ih_src);
			const size_t z_src = bgq_halfvolume2z(ih_src);

			const bool isSurface = bgq_local2isSurface(t1_src,x_src,y_src,z_src) || bgq_local2isSurface(t2_src,x_src,y_src,z_src);
			const bool isSecondInVector_src = bgq_physical2eo(isOdd_src, tv_src,x_src,y_src,z_src);
			assert(!isSecondInVector_src == (t1_src==(tv_src*PHYSICAL_LK*PHYSICAL_LP+2)));
			assert(!isSecondInVector_src == (t2_src==(tv_src*PHYSICAL_LK*PHYSICAL_LP+2+2)));

			for (size_t pd_cnt = P_TUP1; pd_cnt <= P_ZDOWN; pd_cnt+=1) {
				const bgq_physical_direction pd_src = pd_cnt;
				const bgq_direction d = bgq_physical2direction(pd_src);
				const bgq_dimension dim = bgq_direction2dimension(d);
				size_t t_src = ((pd_src==P_TUP2) || (pd_src==P_TDOWN2)) ? t2_src : t1_src;

				size_t tv_dst = tv_src;
				size_t t1_dst = t1_src;
				size_t t2_dst = t2_src;
				size_t t_dst = t_src;
				size_t x_dst = x_src;
				size_t y_dst = y_src;
				size_t z_dst = z_src;
				bgq_physical_direction pd_dst = pd_src;

				switch (pd_src) {
				case P_TUP1:
					t_dst = (t1_dst + LOCAL_LT - 1) % LOCAL_LT;
					//t2_src += 1;
					tv_dst = !isSecondInVector_src ? tv_src-1 : tv_src;
					pd_dst = !isSecondInVector_src ? P_TUP2 : P_TUP1;
					assert(tv_dst==bgq_local2tv(t1_dst,x_dst,y_dst,z_dst));
					break;
				case P_TDOWN1:
					t_dst = (t1_dst + 1) % LOCAL_LT;
					//t2_src -= 1;
					tv_dst = !isSecondInVector_src ? tv_src : tv_src;
					pd_dst = !isSecondInVector_src ? P_TDOWN1 : P_TDOWN2;
					assert(tv_dst==bgq_local2tv(t1_dst,x_dst,y_dst,z_dst));
					break;
				case P_TUP2:
					//t1_src += 1;
					t_dst = (t2_dst + LOCAL_LT - 1) % LOCAL_LT;
					tv_dst = !isSecondInVector_src ? tv_src : tv_src;
					pd_dst = !isSecondInVector_src ? P_TUP1 : P_TUP2;
					assert(tv_dst==bgq_local2tv(t1_dst,x_dst,y_dst,z_dst));
					break;
				case P_TDOWN2:
					//t1_src -= 1;
					t_dst = (t2_dst + 1) % LOCAL_LT;
					tv_dst = !isSecondInVector_src ? tv_src : tv_src+1;
					pd_dst = !isSecondInVector_src ? P_TDOWN2 : P_TDOWN1;
					assert(tv_dst==bgq_local2tv(t1_dst,x_dst,y_dst,z_dst));
					break;
				case P_XUP:
					x_dst -= 1;
					//pd_src = P_XDOWN;
					break;
				case P_XDOWN:
					x_dst += 1;
					//pd_dst = P_XUP;
					break;
				case P_YUP:
					y_dst -= 1;
					//pd_src = P_YDOWN;
					break;
				case P_YDOWN:
					y_dst += 1;
					//pd_src = P_YUP;
					break;
				case P_ZUP:
					z_dst -= 1;
					//pd_dst = P_ZDOWN;
					break;
				case P_ZDOWN:
					z_dst += 1;
					//pd_src = P_ZUP;
					break;
				}
				bgq_weylfield_section sec_src = bgq_HoppingMatrix_init_source_sectionof(isOdd_dst, t_src, x_src, y_src, z_src, d, false);
				bgq_weylfield_section sec_dst = bgq_HoppingMatrix_init_source_sectionof(isOdd_dst, t_dst, x_dst, y_dst, z_dst, bgq_direction_revert(d), true);

				// Find the offset at which the other node expects the datum
				size_t offset_src = g_bgq_offset_fromHalfvolume[isOdd_src][ih_src].pd[pd_src];
				size_t offset_src_write;
				size_t ih_dst = bgq_physical2halfvolume(tv_dst, x_dst, y_dst, z_dst);
				if ((sec_src==sec_body) || (sec_src==sec_surface)) {
					// No communication, do everything locally
					assert(sec_src==sec_dst);
					offset_src_write = offset_src;
				} else {
					// Communication with MPI buffers, translate recv- to send-buffer
					size_t relOffset = bgq_weyl_section_offset(sec_src)/*recv*/ - offset_src;
					offset_src_write = bgq_weyl_section_offset(sec_dst) + relOffset;
				}

				g_bgq_destoffset_fromHalfvolume[isOdd][ih_dst].pd[pd_dst] = bgq_encode_offset(offset_src_write);
				if (isSurface) {
					size_t is_dst = bgq_halfvolume2surface(isOdd_dst,ih_dst);
					g_bgq_destoffset_fromSurface[isOdd][is_dst].pd[pd_dst] = bgq_encode_offset(offset_src_write);
				} else {
					size_t ib_dst = bgq_halfvolume2body(isOdd_dst,ih_dst);
					g_bgq_destoffset_fromBody[isOdd][ib_dst].pd[pd_dst] = bgq_encode_offset(offset_src_write);
				}
			}
		}
	}

	g_bgq_indices_initialized = true;
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
	memset(g_bgq_spinorfields, 0, tot_count * sizeof(bgq_weylsite));
	g_bgq_spinorfields_data = malloc_aligned(tot_count * field_datasize, 128);

	for (size_t i = 0; i < tot_count; i += 1) {
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

static void bgq_spinorfield_initdstptr_halo(bool isOdd, bgq_weylfield_controlblock *field, bgq_physical_direction pd) {
	bgq_direction d = bgq_physical2direction(pd);
	bgq_dimension dim = bgq_direction2dimension(d);

	size_t begin_tv = 0;
	size_t begin_x = 0;
	size_t begin_y = 0;
	size_t begin_z = 0;
	size_t end_tv = PHYSICAL_LTV;
	size_t end_x = PHYSICAL_LX;
	size_t end_y = PHYSICAL_LY;
	size_t end_z = PHYSICAL_LZ;
	int shift_t = 0;
	int shift_x = 0;
	int shift_y = 0;
	int shift_z = 0;
	bgq_weylfield_section sec_neighbor;
	bgq_weylfield_section sec;

	switch (pd) {
	case P_TUP1:
		case P_TUP2:
		begin_tv = PHYSICAL_LTV - 1;
		end_tv = PHYSICAL_LTV - 1;
		shift_t = +1;
		sec_neighbor = sec_send_tdown;
		sec = sec_recv_tup;
		break;
	case P_TDOWN1:
		case P_TDOWN2:
		begin_tv = PHYSICAL_LTV - 1;
		end_tv = PHYSICAL_LTV - 1;
		shift_t = -1;
		sec_neighbor = sec_send_tup;
		sec = sec_recv_tdown;
		break;
	case P_XUP:
		begin_tv = PHYSICAL_LX - 1;
		end_tv = PHYSICAL_LX - 1;
		shift_x = +1;
		sec_neighbor = sec_send_xdown;
		sec = sec_recv_xup;
		break;
	case P_XDOWN:
		begin_tv = 0;
		end_tv = 0;
		shift_x = -1;
		sec_neighbor = sec_send_xup;
		sec = sec_recv_xdown;
		break;
	case P_YUP:
		begin_tv = PHYSICAL_LY - 1;
		end_tv = PHYSICAL_LY - 1;
		shift_y = +1;
		sec_neighbor = sec_send_ydown;
		sec = sec_recv_yup;
		break;
	case P_YDOWN:
		begin_tv = 0;
		end_tv = 0;
		shift_y = -1;
		sec_neighbor = sec_send_yup;
		sec = sec_recv_ydown;
		break;
	case P_ZUP:
		begin_tv = PHYSICAL_LZ - 1;
		end_tv = PHYSICAL_LZ - 1;
		shift_x = +1;
		sec_neighbor = sec_send_zdown;
		sec = sec_recv_zup;
		break;
	case P_ZDOWN:
		begin_tv = 0;
		end_tv = 0;
		shift_z = -1;
		sec_neighbor = sec_send_zup;
		sec = sec_recv_zdown;
		break;
	}

	for (size_t z = begin_z; z < end_z; z += 1) {
		for (size_t y = begin_y; y < end_y; y += 1) {
			for (size_t x = begin_x; x < end_x; x += 1) {
				for (size_t tv = begin_tv; tv < end_tv; tv += 1) {
					size_t t1 = bgq_physical2t1(isOdd, tv, x, y, z);
					size_t t2 = bgq_physical2t2(isOdd, tv, x, y, z);
					size_t is = bgq_physical2surface(isOdd, tv, x, y, z);

					bool isOdd_neighbor = !isOdd;
					size_t tv_neighbor = tv; // For completeness
					size_t t1_neighbor = (t1+LOCAL_LT+shift_t)%LOCAL_LT;
					size_t t2_neighbor = (t2+LOCAL_LT+shift_t)%LOCAL_LT;
					size_t x_neighbor = (x+LOCAL_LX+shift_x)%LOCAL_LX;
					size_t y_neighbor = (y+LOCAL_LY+shift_y)%LOCAL_LY;
					size_t z_neighbor = (z+LOCAL_LZ+shift_z)%LOCAL_LZ;
					size_t is_neighbor = bgq_physical2surface(isOdd_neighbor, tv_neighbor, x_neighbor, y_neighbor, z_neighbor);
					assert(isOdd_neighbor == bgq_local2isOdd(t1_neighbor, x_neighbor, y_neighbor, z_neighbor));
					assert(isOdd_neighbor == bgq_local2isOdd(t2_neighbor, x_neighbor, y_neighbor, z_neighbor));

					switch (pd) {
					case P_TUP1:
						break;
					case P_TUP2:
						break;
					case P_TDOWN1:
						break;
					case P_TDOWN2:
						break;
					case P_XUP:
						break;
					case P_XDOWN:
						break;
					case P_YUP:
						break;
					case P_YDOWN:
						break;
					case P_ZUP:
						break;
					case P_ZDOWN:
						break;
					}

					// Determine where the neighbor node would write this value
					size_t offset = bgq_decode_offset(g_bgq_destoffset_fromSurface[isOdd_neighbor][is_neighbor].pd[sec_neighbor]);
					assert(bgq_weyl_section_offset(sec_neighbor) <= offset && offset < bgq_weyl_section_offset(sec_neighbor + 1));
					size_t relativeOffset = offset - bgq_weyl_section_offset(sec_neighbor);

					// determine which index it is
					size_t srcidx = relativeOffset / (dim==DIM_T) ? sizeof(bgq_weyl_nonvec) : sizeof(bgq_weyl_vec);

					// Where does it need to be written?
					void *destptr = bgq_weylsite_getdirection(&field->sec_surface[is], pd);

					field->destptrFromRecvZup[srcidx] = destptr;
				}
			}
		}
	}
}


void bgq_spinorfield_reset(bgq_weylfield_controlblock *field, bool isOdd, bool activateWeyl, bool activateFull) {
	bool changedOddness;
	if (field->isInitinialized) {
		changedOddness = (field->isOdd != isOdd);
	} else {
		changedOddness = true;
	}

	if (changedOddness) {
		uint8_t *weylbase = field->sec_weyl;
		for (size_t ih = 0; ih < PHYSICAL_VOLUME; ih += 1) {
			for (size_t pd = 0; pd < P_COUNT; pd += 1) {
				field->destptrFromHalfvolume[ih].pd[pd] = weylbase + bgq_decode_offset(g_bgq_destoffset_fromHalfvolume[isOdd][ih].pd[pd]);
			}
		}
		for (size_t is = 0; is < PHYSICAL_SURFACE; is += 1) {
			for (size_t pd = 0; pd < P_COUNT; pd += 1) {
				field->destptrFromSurface[is].pd[pd] = weylbase + bgq_decode_offset(g_bgq_destoffset_fromSurface[isOdd][is].pd[pd]);
			}
		}
		for (size_t ib = 0; ib < PHYSICAL_BODY; ib += 1) {
			for (size_t pd = 0; pd < P_COUNT; pd += 1) {
				field->destptrFromBody[ib].pd[pd] = weylbase + bgq_decode_offset(g_bgq_destoffset_fromBody[isOdd][ib].pd[pd]);
			}
		}

		if (COMM_T) {
			bgq_spinorfield_initdstptr_halo(isOdd, field, P_TUP1);
			bgq_spinorfield_initdstptr_halo(isOdd, field, P_TUP2);
			bgq_spinorfield_initdstptr_halo(isOdd, field, P_TDOWN1);
			bgq_spinorfield_initdstptr_halo(isOdd, field, P_TDOWN2);
		}
		if (COMM_X) {
			bgq_spinorfield_initdstptr_halo(isOdd, field, P_XUP);
			bgq_spinorfield_initdstptr_halo(isOdd, field, P_XDOWN);
		}
		if (COMM_Y) {
			bgq_spinorfield_initdstptr_halo(isOdd, field, P_YUP);
			bgq_spinorfield_initdstptr_halo(isOdd, field, P_YDOWN);
		}
		if (COMM_Z) {
			bgq_spinorfield_initdstptr_halo(isOdd, field, P_ZUP);
			bgq_spinorfield_initdstptr_halo(isOdd, field, P_ZDOWN);
		}
	}

	field->isInitinialized = true;
	field->isOdd = isOdd;
	field->isSloppy = false;
	field->hasWeylfieldData = activateWeyl;
	field->hasFullspinorData = activateFull;
}

