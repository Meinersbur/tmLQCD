/*
 * bgq_field.c
 *
 *  Created on: Oct 13, 2012
 *      Author: meinersbur
 */


#define BGQ_FIELD_C_
#include "bgq_field.h"

#include "bgq_qpx.h"


static bgq_weylfield_section bgq_HoppingMatrix_init_source_sectionof(bool isOdd, size_t t, size_t x, size_t y, size_t z, bgq_direction d) {
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
		return sec_recv_yup;
	}
	if (COMM_Z && (z==0) && (d==ZDOWN)) {
		return sec_recv_ydown;
	}

	bool isSurface = bgq_local2isSurface(t,x,y,z);
	return isSurface ? sec_surface : sec_body;
}





void bgq_indices_init() {
	if (g_bgq_indices_initialized)
		return;


	size_t body_volume;
	if (PHYSICAL_LTV<=2 || PHYSICAL_LX<=2 || PHYSICAL_LY<=2 || PHYSICAL_LZ<=2) {
		body_volume = 0;
	} else {
		size_t body_tv = COMM_T ? (PHYSICAL_LTV-2) : PHYSICAL_LTV;
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

		g_bgq_offset_fromHalfvolume[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(g_bgq_offset_fromHalfvolume[isOdd]));
		g_bgq_offset_fromSurface[isOdd] = malloc(PHYSICAL_SURFACE * sizeof(g_bgq_offset_fromSurface[isOdd]));
		g_bgq_offset_fromBody[isOdd] = malloc(PHYSICAL_BODY * sizeof(g_bgq_offset_fromBody[isOdd]));

		g_bgq_destoffset_fromHalfvolume[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(g_bgq_destoffset_fromHalfvolume[isOdd]));
		g_bgq_destoffset_fromSurface[isOdd] = malloc(PHYSICAL_SURFACE * sizeof(g_bgq_destoffset_fromSurface[isOdd]));
		g_bgq_destoffset_fromBody[isOdd] = malloc(PHYSICAL_BODY * sizeof(g_bgq_destoffset_fromBody[isOdd]));
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
			assert(isSurface == bgq_local2isSurface(t1, x, y, z));
			assert(isSurface == bgq_local2isSurface(t2, x, y, z));

			if (isSurface) {
				size_t is = nextIndexSurface;
				g_bgq_index_surface2halfvolume[isOdd][is] = ih;
				g_bgq_index_halfvolume2surfacebody[isOdd][ih] = is;
				g_bgq_index_halfvolume2surface[isOdd][ih] = is;
				g_bgq_index_halfvolume2body[isOdd][ih] = -1;
				nextIndexSurface += 1;
			} else {
				size_t ib = nextIndexBody;
				g_bgq_index_body2halfvolume[isOdd][ib] = ih;
				g_bgq_index_halfvolume2surfacebody[isOdd][ih] = ib;
				g_bgq_index_halfvolume2surface[isOdd][ih] = -1;
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
			assert(isSurface == (g_bgq_index_halfvolume2surface[isOdd][ih]!=-1));
			bgq_weylfield_section mainsec = isSurface ? sec_surface : sec_body;

			for (size_t pd = P_TUP1; pd <= P_ZDOWN; pd+=1) {
				bgq_direction d = bgq_physical2direction(d);
				bgq_dimension dim = bgq_direction2dimension(d);
				bgq_weylfield_section sec = bgq_HoppingMatrix_init_source_sectionof(isOdd, ((pd==P_TUP2) || (pd==P_TDOWN2)) ? t2 : t1, x, y, z, d);
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
					size_t is = g_bgq_index_halfvolume2surface[isOdd][ih];
					assert(is!=-1);
					g_bgq_offset_fromSurface[isOdd][is].pd[pd] = bgq_encode_offset(thisOffset);
				} else {
					size_t ib = g_bgq_index_halfvolume2body[isOdd][ih];
					assert(ib!=-1);
					g_bgq_offset_fromBody[isOdd][ib].pd[pd] = bgq_encode_offset(thisOffset);
				}
			}
		}
	}


	// Determine dest locations
	for (size_t isOdd = false; isOdd <= true; isOdd+=1) {
		for (size_t ih = 0; ih < PHYSICAL_VOLUME; ih+=1) {
			size_t tv = PHYSICAL_LEXICAL2TV(isOdd, ih);
			size_t t1 = PHYSICAL_LEXICAL2T1(isOdd, ih);
			size_t t2 = PHYSICAL_LEXICAL2T2(isOdd, ih);
			size_t x = PHYSICAL_LEXICAL2X(isOdd, ih);
			size_t y = PHYSICAL_LEXICAL2Y(isOdd, ih);
			size_t z = PHYSICAL_LEXICAL2Z(isOdd, ih);

			bool isSurface =  bgq_local2isSurface(t1,x,y,z) || bgq_local2isSurface(t2,x,y,z);
			bool isFirstInVector = bgq_physical2eo(isOdd, tv,x,y,z);
			assert(isFirstInVector == (t1==tv*PHYSICAL_LK*PHYSICAL_LP));
			assert(isFirstInVector == (t2==(tv*PHYSICAL_LK*PHYSICAL_LP+2)));

			for (size_t pd = P_TUP1; pd <= P_ZDOWN; pd+=1) {
				bgq_direction d = bgq_physical2direction(d);

				size_t tv_src = tv;
				size_t t1_src = t1;
				size_t t2_src = t2;
				size_t x_src = x;
				size_t y_src = y;
				size_t z_src = z;
				bgq_physical_direction pd_src;
				switch (pd) {
				case P_TUP1:
					t1_src += 1;
					//t2_src += 1;
					tv_src = tv;
					pd_src = isFirstInVector ? P_TDOWN1 : P_TDOWN2;
					assert(tv_src==t1_src/(PHYSICAL_LP*PHYSICAL_LK));
					break;
				case P_TDOWN1:
					t1_src -= 1;
					//t2_src -= 1;
					tv_src = isFirstInVector ? tv-1 : tv;
					pd_src = isFirstInVector ? P_TUP2 : P_TUP1;
					assert(tv_src==t1_src/(PHYSICAL_LP*PHYSICAL_LK));
					break;
				case P_TUP2:
					//t1_src += 1;
					t2_src += 1;
					tv_src = isFirstInVector ? tv : tv+1;
					pd_src = isFirstInVector ? P_TDOWN2 : P_TDOWN1;
					assert(tv_src==t2_src/(PHYSICAL_LP*PHYSICAL_LK));
					break;
				case P_TDOWN2:
					//t1_src -= 1;
					t2_src -= 1;
					tv_src = tv;
					pd_src = isFirstInVector ? P_TUP1 : P_TUP2;
					assert(tv_src==t2_src/(PHYSICAL_LP*PHYSICAL_LK));
					break;
				case P_XUP:
					x_src += 1;
					pd_src = P_XDOWN;
					break;
				case P_XDOWN:
					x_src -= 1;
					pd_src = P_XUP;
					break;
				case P_YUP:
					y_src += 1;
					pd_src = P_YDOWN;
					break;
				case P_YDOWN:
					y_src -= 1;
					pd_src = P_YUP;
					break;
				case P_ZUP:
					z_src += 1;
					pd_src = P_ZDOWN;
					break;
				case P_ZDOWN:
					z_src -= 1;
					pd_src = P_ZUP;
					break;
				}

				size_t ih_src = bgq_physical2halfvolume(tv_src,x_src,y_src,z_src);
				uint32_t encodedOffset = g_bgq_offset_fromHalfvolume[!isOdd][ih_src].pd[pd_src];
				g_bgq_destoffset_fromHalfvolume[isOdd][ih].pd[pd] = encodedOffset;
				if (isSurface) {
					size_t is = g_bgq_index_halfvolume2surface[isOdd][ih];
					assert(is!=-1);
					g_bgq_destoffset_fromSurface[isOdd][is].pd[pd] = encodedOffset;
				} else {
					size_t ib = g_bgq_index_halfvolume2body[isOdd][ih];
					assert(ib!=-1);
					g_bgq_destoffset_fromBody[isOdd][ib].pd[pd] = encodedOffset;
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

