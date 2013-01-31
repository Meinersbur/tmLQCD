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
#include "bgq_comm.h"




static bgq_weylfield_section bgq_HoppingMatrix_init_source_sectionof_local(size_t t, size_t x, size_t y, size_t z, bgq_direction d) {
	// Check whether moving out of the local lattice
	if (HALO_T && (t==LOCAL_LT-1) && (d==TUP)) {
		if (BGQ_UNVECTORIZE || !COMM_T)
			return sec_temp_tup;
		else
			return sec_send_tup;
	}
	if (HALO_T && (t==0) && (d==TDOWN)) {
		if (BGQ_UNVECTORIZE || !COMM_T)
			return sec_temp_tdown;
		else
			return sec_send_tdown;
	}
	if (HALO_X && (x==LOCAL_LX-1) && (d==XUP)) {
		return sec_send_xup;
	}
	if (HALO_X && (x==0) && (d==XDOWN)) {
		return sec_send_xdown;
	}
	if (HALO_Y && (y==LOCAL_LY-1) && (d==YUP)) {
		return sec_send_yup;
	}
	if (HALO_Y && (y==0) && (d==YDOWN)) {
		return sec_send_ydown;
	}
	if (HALO_Z && (z==LOCAL_LZ-1) && (d==ZUP)) {
		return sec_send_zup;
	}
	if (HALO_Z && (z==0) && (d==ZDOWN)) {
		return sec_send_zdown;
	}

	// Stay within the local lattice, determine whether the target is on the surface
	bgq_direction_move_local(&t, &x, &y, &z, d);
	bool isSurface = bgq_local2isSurface(t, x, y, z);
	return isSurface ? sec_surface : sec_body;
}


static bgq_weylfield_section bgq_HoppingMatrix_init_source_sectionof_physical(bool isOdd, size_t tv, size_t x, size_t y, size_t z, bgq_direction d) {
	size_t t1 = bgq_physical2t1(isOdd, tv, x, y, z);
	bgq_weylfield_section sec1 = bgq_HoppingMatrix_init_source_sectionof_local(t1,x,y,z,d);
	size_t t2 = bgq_physical2t2(isOdd, tv, x, y, z);
	bgq_weylfield_section sec2 = bgq_HoppingMatrix_init_source_sectionof_local(t2,x,y,z,d);

	if (sec1 == sec2) {
		return sec1;
	} else if (((sec1 == sec_surface) || (sec1 == sec_body)) && ((sec2 == sec_send_tup) || (sec2 == sec_temp_tup))) {
		// Happens because t1 is not strictly at the border
		// We give t2 priority such that we save both into the buffer
		return sec2;
	} else if (((sec2 == sec_surface) || (sec2 == sec_body)) && ((sec1 == sec_send_tdown) || (sec1 == sec_temp_tdown))) {
		return sec1;
	} else {
		assert(!"Unknown case");
		UNREACHABLE
		return -1;
	}
}


bgq_direction bgq_offset2ddst(size_t offset) {
	bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
	switch (sec) {
	case sec_send_tup:
	case sec_recv_tdown:
	case sec_temp_tup:
		return TDOWN;
	case sec_send_tdown:
	case sec_recv_tup:
	case sec_temp_tdown:
		return TUP;
	case sec_send_xup:
	case sec_recv_xdown:
		return XDOWN;
	case sec_send_xdown:
	case sec_recv_xup:
		return XUP;
	case sec_send_yup:
	case sec_recv_ydown:
		return YDOWN;
	case sec_send_ydown:
	case sec_recv_yup:
		return YUP;
	case sec_send_zup:
	case sec_recv_zdown:
		return ZDOWN;
	case sec_send_zdown:
	case sec_recv_zup:
		return ZUP;
	case sec_surface:
	case sec_body: {
		size_t index = (offset - bgq_weyl_section_offset(sec_surface)) / sizeof(bgq_weyl_vec_double);
		return index % PHYSICAL_LD;
	}
	default:
		UNREACHABLE
		break;
	}
	return -1;
}


bgq_direction bgq_offset2dsrc(size_t offset) {
	bgq_direction d_dst = bgq_offset2ddst(offset);
	return bgq_direction_revert(d_dst);
}


size_t bgq_src2ih_dst(size_t t_src, size_t x_src, size_t y_src, size_t z_src, bgq_direction d_src) {
	size_t t_dst = t_src;
	size_t x_dst = x_src;
	size_t y_dst = y_src;
	size_t z_dst = z_src;

	switch (d_src) {
	case TDOWN:
		t_dst = (t_src + LOCAL_LT - 1) % LOCAL_LT;
		break;
	case TUP:
		t_dst = (t_src + 1) % LOCAL_LT;
		break;
	case XDOWN:
		x_dst = (x_src + LOCAL_LX - 1) % LOCAL_LX;
		break;
	case XUP:
		x_dst = (x_src + 1) % LOCAL_LX;
		break;
	case YDOWN:
		y_dst = (y_src + LOCAL_LY - 1) % LOCAL_LY;
		break;
	case YUP:
		y_dst = (y_src + 1) % LOCAL_LY;
		break;
	case ZDOWN:
		z_dst = (z_src + LOCAL_LZ - 1) % LOCAL_LZ;
		break;
	case ZUP:
		z_dst = (z_src + 1) % LOCAL_LZ;
		break;
	}

	size_t ih_dst = bgq_local2halfvolume(t_dst,x_dst,y_dst,z_dst);
	return ih_dst;
}


size_t bgq_src2k_dst(size_t t_src, size_t x_src, size_t y_src, size_t z_src, bgq_direction d_src) {
	size_t t_dst = t_src;
	size_t x_dst = x_src;
	size_t y_dst = y_src;
	size_t z_dst = z_src;

	switch (d_src) {
	case TDOWN:
		t_dst = (t_src + LOCAL_LT - 1) % LOCAL_LT;
		break;
	case TUP:
		t_dst = (t_src + 1) % LOCAL_LT;
		break;
	case XDOWN:
		x_dst = (x_src + LOCAL_LX - 1) % LOCAL_LX;
		break;
	case XUP:
		x_dst = (x_src + 1) % LOCAL_LX;
		break;
	case YDOWN:
		y_dst = (y_src + LOCAL_LY - 1) % LOCAL_LY;
		break;
	case YUP:
		y_dst = (y_src + 1) % LOCAL_LY;
		break;
	case ZDOWN:
		z_dst = (z_src + LOCAL_LZ - 1) % LOCAL_LZ;
		break;
	case ZUP:
		z_dst = (z_src + 1) % LOCAL_LZ;
		break;
	}

	size_t k_dst = bgq_local2k(t_dst, x_dst, y_dst, z_dst);
	return k_dst;
}


void bgq_indices_init() {
	if (g_bgq_indices_initialized)
		return;
	g_bgq_indices_initialized = true; // Take care for uses within this function itself

	g_comm_t = (g_nb_t_dn != g_proc_id);
	g_comm_x = (g_nb_x_dn != g_proc_id);
	g_comm_y = (g_nb_y_dn != g_proc_id);
	g_comm_z = (g_nb_z_dn != g_proc_id);
	assert(g_comm_t == (g_nb_t_up != g_proc_id));
	assert(g_comm_x == (g_nb_x_up != g_proc_id));
	assert(g_comm_y == (g_nb_y_up != g_proc_id));
	assert(g_comm_z == (g_nb_z_up != g_proc_id));
	g_comm_t = true;
	g_comm_x = true;
	g_comm_y = true;
	g_comm_z = true;
	g_bgq_dimension_isDistributed[DIM_T] = g_comm_t;
	g_bgq_dimension_isDistributed[DIM_X] = g_comm_x;
	g_bgq_dimension_isDistributed[DIM_Y] = g_comm_y;
	g_bgq_dimension_isDistributed[DIM_Z] = g_comm_z;
	g_bgq_dimension_hasHalo[DIM_T] = true;
	g_bgq_dimension_hasHalo[DIM_X] = g_comm_x;
	g_bgq_dimension_hasHalo[DIM_Y] = g_comm_y;
	g_bgq_dimension_hasHalo[DIM_Z] = g_comm_z;
	master_print("BGQ SPI/MPI communication enabled for: %s%s%s%s\n", COMM_T?"T,":"", COMM_X?"X,":"", COMM_Y?"Y,":"", COMM_Z?"Z":"");

	assert(PHYSICAL_LTV>=2);
	PHYSICAL_VOLUME = (PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LY*PHYSICAL_LZ);
	size_t physical_body;
	if ((PHYSICAL_LTV <= 1) || (PHYSICAL_LX <= 2) || (PHYSICAL_LY <= 2) || (PHYSICAL_LZ <= 2)) {
		physical_body = 0;
	} else {
		size_t body_tv = COMM_T ? (PHYSICAL_LTV - 1) : PHYSICAL_LTV;
		size_t body_x = COMM_X ? (PHYSICAL_LX - 2) : PHYSICAL_LX;
		size_t body_y = COMM_Y ? (PHYSICAL_LY - 2) : PHYSICAL_LY;
		size_t body_z = COMM_Z ? (PHYSICAL_LZ - 2) : PHYSICAL_LZ;
		physical_body = body_tv * body_x * body_y * body_z;
	}
	PHYSICAL_SURFACE = PHYSICAL_VOLUME - physical_body;
	PHYSICAL_INNER = (PHYSICAL_LTV - 1)*(PHYSICAL_LX - 2)*(PHYSICAL_LY - 2)*(PHYSICAL_LZ - 2);

	// Setup array for spinor index translation
	// It would be difficult to find an explicit expression for indices that are only on surface/body
	// Thus, we simply iterate over all locations and allocate locations in order
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		g_bgq_collapsed2halfvolume[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_collapsed2halfvolume[isOdd]));
		g_bgq_halfvolume2collapsed[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_halfvolume2collapsed[isOdd]));
		size_t nextIndexSurface = 0;
		size_t nextIndexBody = PHYSICAL_SURFACE;
		//size_t nextIndexOuterBody = PHYSICAL_SURFACE;
		//size_t nextIndexInnerBody = PHYSICAL_SURFACE+PHYSICAL_OUTERBODY;

		for (ucoord tv = 0; tv < PHYSICAL_LTV; tv+=1) {
			for (ucoord x = 0; x < PHYSICAL_LX; x+=1) {
				for (ucoord y = 0; y < PHYSICAL_LY; y+=1) {
					for (ucoord z = 0; z < PHYSICAL_LZ; z+=1) { // Make inner z-lines consecutive if possible
						ucoord ih = bgq_physical2halfvolume(tv,x,y,z);
						size_t t1 = bgq_physical2t1(isOdd,tv,x,y,z);
						size_t t2 = bgq_physical2t2(isOdd,tv,x,y,z);
						assert(isOdd == bgq_local2isOdd(t1,x,y,z));
						assert(isOdd == bgq_local2isOdd(t2,x,y,z));

						bool isSurface = bgq_physical2isSurface(isOdd, tv, x, y, z);
						//bool isInner = bgq_physical2isInner(isOdd, tv, x, y, z);

						size_t ic;
						if (isSurface) {
							ic = nextIndexSurface;
							nextIndexSurface += 1;
							assert(bgq_collapsed2isSurface(ic));
						} else {
							ic = nextIndexBody;
							nextIndexBody += 1;
							assert(bgq_collapsed2isBody(ic));
						}
						g_bgq_collapsed2halfvolume[isOdd][ic] = ih;
						g_bgq_halfvolume2collapsed[isOdd][ih] = ic;
					}
				}
			}
		}
	}


	// Setup indices for communication spinors
	for (size_t isOdd = false; isOdd <= true; isOdd+=1) {
		g_bgq_collapsed2recvidx[isOdd] = malloc_aligned(PHYSICAL_VOLUME * sizeof(*g_bgq_collapsed2recvidx[isOdd]), BGQ_ALIGNMENT_L2);
		memset(g_bgq_collapsed2recvidx[isOdd], 0xFE, PHYSICAL_VOLUME * sizeof(*g_bgq_collapsed2recvidx[isOdd]));
		g_bgq_collapsed2sendidx[isOdd] = malloc_aligned(PHYSICAL_VOLUME * sizeof(*g_bgq_collapsed2sendidx[isOdd]), BGQ_ALIGNMENT_L2);
		memset(g_bgq_collapsed2sendidx[isOdd], 0xFE, PHYSICAL_VOLUME * sizeof(*g_bgq_collapsed2sendidx[isOdd]));

		for (ucoord d = TUP; d <= ZDOWN; d+=1) {
			bgq_dimension dim = bgq_direction2dimension(d);
			size_t length = bgq_dimension_physicalLengthOfHalo(dim, false);
			g_bgq_recvidx2collapsed[isOdd][d] = malloc(length * sizeof(*g_bgq_recvidx2collapsed[isOdd][d]));
			memset(g_bgq_recvidx2collapsed[isOdd][d], 0xFF, length * sizeof(*g_bgq_recvidx2collapsed[isOdd][d]));
			g_bgq_sendidx2collapsed[isOdd][d] = malloc(length * sizeof(*g_bgq_sendidx2collapsed[isOdd][d]));
			memset(g_bgq_sendidx2collapsed[isOdd][d], 0xFF, length * sizeof(*g_bgq_sendidx2collapsed[isOdd][d]));
		}
	}
	// Setup receive buffer indices
	for (size_t isOdd = false; isOdd <= true; isOdd+=1) {
		bool isOdd_src = isOdd;
		bool isOdd_dst = !isOdd;

		for (size_t d_dst = TUP; d_dst <= ZDOWN; d_dst+=1) {
			size_t nextRecvIdx = 0;
			for (ucoord ic_dst = 0; ic_dst < PHYSICAL_VOLUME; ic_dst+=1) {
				bgq_direction d_src = bgq_direction_revert(d_dst);
				size_t idx = -1;
				if (bgq_direction_isCrossingBorder_collapsed(isOdd_dst, ic_dst, d_dst)) {
					idx = nextRecvIdx;
					nextRecvIdx += 1;
				}
				assert(g_bgq_collapsed2recvidx[isOdd_dst][ic_dst][d_dst] == 0xFEFEFEFEFEFEFEFEull);
				g_bgq_collapsed2recvidx[isOdd_dst][ic_dst][d_dst] = idx;

				if (idx!=-1) {
					assert(g_bgq_recvidx2collapsed[isOdd_dst][d_dst][idx] == -1);
					g_bgq_recvidx2collapsed[isOdd_dst][d_dst][idx] = ic_dst;
				}
			}
		}
	}
	// Setup send buffer indices
	for (size_t isOdd = false; isOdd <= true; isOdd+=1) {
		bool isOdd_src = isOdd;
		bool isOdd_dst = !isOdd;

		for (size_t d_dst = TUP; d_dst <= ZDOWN; d_dst+=1) {
			bgq_direction d_src = bgq_direction_revert(d_dst);
			for (ucoord ic_src = 0; ic_src < PHYSICAL_VOLUME; ic_src+=1) {
				ucoord t1_src = bgq_collapsed2t1(isOdd_src, ic_src);
				ucoord t2_src = bgq_collapsed2t2(isOdd_src, ic_src);
				ucoord x_src = bgq_collapsed2x(isOdd_src, ic_src);
				ucoord y_src = bgq_collapsed2y(isOdd_src, ic_src);
				ucoord z_src = bgq_collapsed2z(isOdd_src, ic_src);

				size_t idx = -1;
				if (bgq_direction_isCrossingBorder_collapsed(isOdd_src, ic_src, d_src)) {
					// Find the index at which the receiving node expect the data
					ucoord ic_dst = bgq_direction_move_collapsed(isOdd_src, ic_src, d_src);
					idx = g_bgq_collapsed2recvidx[isOdd_dst][ic_dst][d_dst];
					assert(idx < 0xFEFEFEFEFEFEFEFEull);
				}
				assert(g_bgq_collapsed2sendidx[isOdd_src][ic_src][d_dst] == 0xFEFEFEFEFEFEFEFEULL);
				g_bgq_collapsed2sendidx[isOdd_src][ic_src][d_dst] = idx;

				if (idx!=-1) {
					assert(g_bgq_sendidx2collapsed[isOdd_src][d_dst][idx] == -1);
					g_bgq_sendidx2collapsed[isOdd_src][d_dst][idx] = ic_src;
				}
			}
		}
	}

#ifndef NDEBUG
	for (size_t isOdd = false; isOdd <= true; isOdd+=1) {
		bool isOdd_src = isOdd;
		bool isOdd_dst = !isOdd;

		for (ucoord ic_src = 0; ic_src < PHYSICAL_VOLUME; ic_src+=1) {
			for (ucoord d = TUP; d <= ZDOWN; d+=1) {
				bgq_direction d_dst = d;
				bgq_direction d_src = bgq_direction_revert(d);
				ucoord ic_dst = bgq_direction_move_collapsed(isOdd_src, ic_src, d_src);

				assert(g_bgq_collapsed2sendidx[isOdd_src][ic_src][d_dst] == g_bgq_collapsed2recvidx[isOdd_dst][ic_dst][d_dst]);
			}
		}

		for (ucoord d = TUP; d <= ZDOWN; d+=1) {
			bgq_dimension dim = bgq_direction2dimension(d);
			bgq_direction d_dst = d;
			bgq_direction d_src = bgq_direction_revert(d);

			size_t length = bgq_dimension_physicalLengthOfHalo(dim, false);
			for (ucoord idx = 0; idx < length; idx+=1) {
				ucoord ic_src = g_bgq_sendidx2collapsed[isOdd_src][d_dst][idx];
				ucoord ic_dst = g_bgq_recvidx2collapsed[isOdd_dst][d_dst][idx];
				assert(bgq_direction_move_collapsed(isOdd_src, ic_src, d_src) == ic_dst);
				assert(bgq_direction_move_collapsed(isOdd_dst, ic_dst, d_dst) == ic_src);
			}
		}
	}
#endif


	// Setup mapping of weyl to some memory offset
	// (where to read a datum for hoppingmatrix)
	// Note: although we arrange data ordered as in ih_dst, the field contains data ih_src
	//TODO: Use previous index allocation for buffers
	assert(bgq_weyl_section_offset(sec_end) % sizeof(bgq_weyl_vec_double) == 0);
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		size_t indices = bgq_weyl_section_offset(sec_end) / sizeof(bgq_weyl_vec_double);
		g_bgq_index2collapsed[isOdd] = malloc(indices * sizeof(*g_bgq_index2collapsed[isOdd]));
		g_bgq_collapsed2indexsend[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_collapsed2indexsend[isOdd]));
	}
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		bool isOdd_src = !isOdd;
		bool isOdd_dst = isOdd;

		size_t nextoffset[sec_end];
		for (bgq_weylfield_section sec = 0; sec < sec_end; sec += 1) {
			nextoffset[sec] = bgq_weyl_section_offset(sec);
		}

		for (ucoord ic_dst = 0; ic_dst < PHYSICAL_VOLUME; ic_dst += 1) {
			//ucoord ic_dst = bgq_halfvolume2collapsed(isOdd_dst, ih_dst);
			ucoord ih_dst = bgq_collapsed2halfvolume(isOdd_dst, ic_dst);
			size_t tv_dst = bgq_halfvolume2tv(ih_dst);
			size_t x_dst = bgq_halfvolume2x(ih_dst);
			size_t y_dst = bgq_halfvolume2y(ih_dst);
			size_t z_dst = bgq_halfvolume2z(ih_dst);
			bool isSurface_dst = bgq_collapsed2isSurface(ic_dst);
			bgq_weylfield_section mainsec_dst = isSurface_dst ? sec_surface : sec_body;

			for (size_t d_dst = 0; d_dst < PHYSICAL_LD; d_dst += 1) {
				bgq_direction d_src = bgq_direction_revert(d_dst);
				bgq_dimension dim = bgq_direction2dimension(d_src);
				ucoord ic_src = bgq_collapsed_dst2src(isOdd_dst, ic_dst, d_dst);
				ucoord ih_src = bgq_collapsed2halfvolume(isOdd_src, ic_src);
				size_t tv_src = bgq_halfvolume2tv(ih_src);
				size_t x_src = bgq_halfvolume2x(ih_src);
				size_t y_src = bgq_halfvolume2y(ih_src);
				size_t z_src = bgq_halfvolume2z(ih_src);
				bool isSurface_src = bgq_collapsed2isSurface(ic_src);

				//bgq_weylfield_section mainsec_src = isSurface_src ? sec_surface : sec_body;
				bgq_weylfield_section sec_write = bgq_HoppingMatrix_init_source_sectionof_physical(isOdd_src, tv_src, x_src, y_src, z_src, d_src);

				// Reserve some offset in surface/body to ensure consecutive layout
				size_t offset_main = nextoffset[mainsec_dst];
				nextoffset[mainsec_dst] += sizeof(bgq_weyl_vec_double);
				assert(bgq_collapsed2consecutiveoffset(ic_dst, d_dst) == offset_main);

				ucoord index_main = bgq_offset2index(offset_main);
				if (ic_dst == 60) {
					int a = 0;
				}
				g_bgq_index2collapsed[isOdd_dst][index_main] = ic_dst;
				assert(bgq_sectionOfOffset(offset_main) == mainsec_dst);

				if (sec_write != mainsec_dst) {
					// If in one of the comm send buffers, also reserve some space there
					assert((sec_write!=sec_surface) && (sec_write!=sec_body));

					size_t offset_write = nextoffset[sec_write];
					nextoffset[sec_write] += sizeof(bgq_weyl_vec_double);

					ucoord index_write = bgq_offset2index(offset_write);
					g_bgq_collapsed2indexsend[isOdd_dst][ic_src/*!!!*/].d[d_dst] = index_write;
					g_bgq_index2collapsed[isOdd_dst][index_write] = ic_dst;
				} else {
					g_bgq_collapsed2indexsend[isOdd_dst][ic_src/*!!!*/].d[d_dst] = index_main;
				}

			}
		}
	}
}



