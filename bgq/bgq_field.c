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

#pragma GCC diagnostic ignored "-Wunused-variable"
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


static bgq_weylfield_section bgq_HoppingMatrix_init_source_sectionof_physical(bool isOdd, size_t tv, size_t x, size_t y, size_t z, bgq_direction d) {
	size_t t1 = bgq_physical2t1(isOdd, tv, x,y,z);
	bgq_weylfield_section sec1 = bgq_HoppingMatrix_init_source_sectionof_local(t1,x,y,z,d);
	size_t t2 = bgq_physical2t2(isOdd, tv, x,y,z);
	bgq_weylfield_section sec2 = bgq_HoppingMatrix_init_source_sectionof_local(t2,x,y,z,d);

	if (sec1 == sec2) {
		return sec1;
	} else if ((sec1 == sec_surface) && (sec2 == sec_recv_tup)) {
		// Happens because t1 is not strictly at the border
		// We give t2 priority such that we save both into the buffer
		return sec_recv_tup;
	} else if ((sec2 == sec_surface) && (sec1 == sec_recv_tdown)) {
		return sec_recv_tdown;
	} else {
		assert(!"Unknown case");
		exit(1);
	}
}


static size_t bgq_src2ih_dst(size_t t_src, size_t x_src, size_t y_src, size_t z_src, bgq_direction d_src) {
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

static size_t bgq_src2k_dst(size_t t_src, size_t x_src, size_t y_src, size_t z_src, bgq_direction d_src) {
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


static size_t bgq_offset_send2recv(size_t offset_send) {
	bgq_weylfield_section sec_send = bgq_sectionOfOffset(offset_send);
	if ((sec_send==sec_surface)||(sec_send==sec_body))
		return offset_send;
	bgq_weylfield_section sec_recv = bgq_section_commbuftran(sec_send, false);
	assert(sec_recv!=sec_send);

	size_t offset_recv = offset_send - bgq_weyl_section_offset(sec_send) + bgq_weyl_section_offset(sec_recv);
	assert(offset_recv);
	return offset_recv;
}


static size_t bgq_offset_recv2send(size_t offset_recv) {
	bgq_weylfield_section sec_recv = bgq_sectionOfOffset(offset_recv);
	if ((sec_recv==sec_surface)||(sec_recv==sec_body))
		return offset_recv;
	bgq_weylfield_section sec_send = bgq_section_commbuftran(sec_recv, true);
	assert(sec_recv!=sec_send); // Happens if offset_recv is send sendbuf

	size_t offset_send = offset_recv - bgq_weyl_section_offset(sec_recv) + bgq_weyl_section_offset(sec_send);
	assert(offset_send);
	return offset_send;
}



// Find the offset in the weylfield where a specific value is expected for the HoppingMatrix kernel
// Used for 1st phase (distribute)
// i.e. body -> body,surface
//      surface -> body,surface or communication buffer to be sent to other node
static size_t bgq_weylfield_destoffsetForWeyl(bool isOdd, size_t ih_src, bgq_direction d_src) {
	bool isOdd_src = isOdd;
	bool isOdd_dst = !isOdd;
	size_t tv_src = bgq_halfvolume2tv(ih_src);
	size_t t1_src = bgq_halfvolume2t1(isOdd_src, ih_src);
	size_t t2_src = bgq_halfvolume2t2(isOdd_src, ih_src);
	size_t x_src = bgq_halfvolume2x(ih_src);
	size_t y_src = bgq_halfvolume2y(ih_src);
	size_t z_src = bgq_halfvolume2z(ih_src);
	//bool isSurface = bgq_physical2isSurface(isOdd_src, tv_src, x_src, y_src, z_src);
	//bgq_dimension dim = bgq_direction2dimension(d_src);

	size_t tv_dst = tv_src;
	size_t t1_dst = t1_src;
	size_t t2_dst = t2_src;
	size_t x_dst = x_src;
	size_t y_dst = y_src;
	size_t z_dst = z_src;
	bgq_direction d_dst = bgq_direction_revert(d_src);

	// Find the coordinate that consumes this value for its calculation
	// Change the coordinate to go into that direction
	//TODO: Use bgq_src2dst function
	switch (d_src) {
	case TDOWN:
		t1_dst = (t1_src + LOCAL_LT - 1) % LOCAL_LT;
		t2_dst = (t2_src + LOCAL_LT - 1) % LOCAL_LT;
		tv_dst = !bgq_physical2eo(isOdd_src, tv_src, x_src, y_src, z_src) ? (tv_src + PHYSICAL_LTV - 1) % PHYSICAL_LTV : tv_src;
		break;
	case TUP:
		t1_dst = (t1_src + 1) % LOCAL_LT;
		t2_dst = (t2_src + 1) % LOCAL_LT;
		tv_dst = !bgq_physical2eo(isOdd_src, tv_src, x_src, y_src, z_src) ? tv_src : (tv_src + 1) % PHYSICAL_LTV;
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
	assert(bgq_local2tv(t1_dst,x_dst,y_dst,z_dst) == tv_dst);
	assert(bgq_local2tv(t2_dst,x_dst,y_dst,z_dst) == tv_dst);
	bgq_weylfield_section sec_read = bgq_HoppingMatrix_init_source_sectionof_physical(isOdd_dst, tv_dst, x_dst, y_dst, z_dst, d_dst); // dst is going to read the weyl
	assert(!bgq_section2isSend(sec_read));
	// Not that the section might logically belong to a different MPI node
	bgq_weylfield_section sec_write = bgq_section_commbuftran(sec_read, true); // src finds the correct location in preparation for dst
	assert(!bgq_section2isRecv(sec_write));
	// Find the offset at which the src node (can be this very same node) expects the datum
	size_t ih_dst = bgq_physical2halfvolume(tv_dst, x_dst, y_dst, z_dst);
	size_t offset_read = bgq_decode_offset(g_bgq_ih_dst2offset[isOdd][ih_dst].d[d_dst]);
	assert(bgq_sectionOfOffset(offset_read) == sec_read);

	size_t offset_write;
	if ((sec_read == sec_body) || (sec_read == sec_surface)) {
		// No communication, do everything locally
		assert((sec_write==sec_body) || (sec_write==sec_surface));
		offset_write = offset_read;
	} else {
		// Communication with MPI buffers, translate recv- to send-buffer
		// i.e. dst expects the datum in the recv buffer, so we have to write it into the send buffer in order to be written there
		assert((sec_write!=sec_body) && (sec_write!=sec_surface));
		offset_write = bgq_offset_recv2send(offset_read);
	}
	assert(offset_write);
	assert(0 <= offset_write && offset_write < bgq_weyl_section_offset(sec_end));
	assert(!bgq_section2isRecv(bgq_sectionOfOffset(offset_write)));
	return offset_write;
}


static size_t bgq_weyllayout_halfvolume2consecutiveoffset(bool isOdd_dst, size_t ih_dst, bgq_direction d_dst) {
	bool isSurface = bgq_halfvolume2isSurface(isOdd_dst, ih_dst);
	if (isSurface) {
		size_t is_dst =bgq_halfvolume2surface(isOdd_dst, ih_dst);
		return bgq_weyl_section_offset(sec_surface) + is_dst*sizeof(bgq_weylsite) + d_dst * sizeof(bgq_weyl_vec);
	} else {
		size_t ib_dst = bgq_halfvolume2body(isOdd_dst, ih_dst);
		return bgq_weyl_section_offset(sec_body) + ib_dst*sizeof(bgq_weylsite) + d_dst * sizeof(bgq_weyl_vec);
	}
}


static bgq_direction bgq_offset2ddst(size_t offset) {
	bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
	switch (sec) {
	case sec_send_tup:
	case sec_recv_tdown:
		return TDOWN;
	case sec_send_tdown:
	case sec_recv_tup:
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
		size_t index = (offset - bgq_weyl_section_offset(sec_surface)) / sizeof(bgq_weyl_vec);
		return index % PHYSICAL_LD;
	}
	default:
		UNREACHABLE
	}
}


static bgq_direction bgq_offset2dsrc(size_t offset) {
	bgq_direction d_dst = bgq_offset2ddst(offset);
	return bgq_direction_revert(d_dst);
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


void bgq_indices_init() {
	if (g_bgq_indices_initialized)
		return;
	g_bgq_indices_initialized = true; // Take care for uses within this function itself

	assert(PHYSICAL_LTV>=2);
	if ((PHYSICAL_LTV <= 1) || (PHYSICAL_LX <= 2) || (PHYSICAL_LY <= 2) || (PHYSICAL_LZ <= 2)) {
		PHYSICAL_BODY = 0;
	} else {
		size_t body_tv = COMM_T ? (PHYSICAL_LTV - 1) : PHYSICAL_LTV;
		size_t body_x = COMM_X ? (PHYSICAL_LX - 2) : PHYSICAL_LX;
		size_t body_y = COMM_Y ? (PHYSICAL_LY - 2) : PHYSICAL_LY;
		size_t body_z = COMM_Z ? (PHYSICAL_LZ - 2) : PHYSICAL_LZ;
		PHYSICAL_BODY = body_tv * body_x * body_y * body_z;
	}

	// Setup array for index translation
	// It would be difficult to find an explicit expression for indices that are only on surface/body
	// Thus, we simply iterate over all locations and allocate locations in order
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		g_bgq_index_surface2halfvolume[isOdd] = malloc(PHYSICAL_SURFACE * sizeof(*g_bgq_index_surface2halfvolume[isOdd]));
		g_bgq_index_body2halfvolume[isOdd] = malloc(PHYSICAL_BODY * sizeof(*g_bgq_index_body2halfvolume[isOdd]));
		g_bgq_index_halfvolume2surface[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_index_halfvolume2surface[isOdd]));
		g_bgq_index_halfvolume2body[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_index_halfvolume2body[isOdd]));
		g_bgq_index_halfvolume2surfacebody[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_index_halfvolume2surfacebody[isOdd]));
		size_t nextIndexSurface = 0;
		size_t nextIndexBody = 0;

		// Iteration order is effectively chosen here
		for (size_t ih = 0; ih < PHYSICAL_VOLUME; ih += 1) {
			size_t tv = bgq_halfvolume2tv(ih);
			size_t t1 = bgq_halfvolume2t1(isOdd, ih);
			size_t t2 = bgq_halfvolume2t2(isOdd, ih);
			size_t x = bgq_halfvolume2x(ih);
			size_t y = bgq_halfvolume2y(ih);
			size_t z = bgq_halfvolume2z(ih);
			assert(ih == bgq_local2halfvolume(t1,x,y,z));
			assert(ih == bgq_local2halfvolume(t2,x,y,z));
			assert(isOdd == bgq_local2isOdd(t1,x,y,z));
			assert(isOdd == bgq_local2isOdd(t2,x,y,z));

			bool isSurface = bgq_physical2isSurface(isOdd, tv, x, y, z);
			if (isSurface) {
				size_t is = nextIndexSurface;
				assert(0 <= is && is < PHYSICAL_SURFACE);
				g_bgq_index_surface2halfvolume[isOdd][is] = ih;
				g_bgq_index_halfvolume2surfacebody[isOdd][ih] = is;
				g_bgq_index_halfvolume2surface[isOdd][ih] = is;
				g_bgq_index_halfvolume2body[isOdd][ih] = (size_t) -1;
				nextIndexSurface += 1;
			} else {
				size_t ib = nextIndexBody;
				assert(0 <= ib && ib < PHYSICAL_BODY);
				g_bgq_index_body2halfvolume[isOdd][ib] = ih;
				g_bgq_index_halfvolume2surfacebody[isOdd][ih] = ib;
				g_bgq_index_halfvolume2surface[isOdd][ih] = (size_t) -1;
				g_bgq_index_halfvolume2body[isOdd][ih] = ib;
				nextIndexBody += 1;
			}
		}
	}

#ifndef NDEBUG
	// Check consistency
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		for (size_t is2 = 0; is2 < PHYSICAL_SURFACE; is2 += 1) {
			size_t ih = bgq_surface2halfvolume(isOdd, is2);
			assert(bgq_halfvolume2isSurface(isOdd,ih));

			assert(is2 == bgq_halfvolume2surface(isOdd,ih));
			assert(ih == bgq_surface2halfvolume(isOdd,is2));
		}
	}
#endif

	// Setup mapping of weyl to some memory offset
	// (where to read a datum for hoppingmatrix)
	// Note: although we arrange data ordered as in ih_dst, the field contains data ih_src
	assert(bgq_weyl_section_offset(sec_end) % sizeof(bgq_weyl_vec) == 0);
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		size_t indices = bgq_weyl_section_offset(sec_end) / sizeof(bgq_weyl_vec);
		g_bgq_index2ih_dst[isOdd] = malloc(indices * sizeof(*g_bgq_index2ih_dst[isOdd]));
		g_bgq_index2d_dst[isOdd] = malloc(indices * sizeof(*g_bgq_index2d_dst[isOdd]));

		g_bgq_ih_dst2offset[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_ih_dst2offset[isOdd]));
		g_bgq_is_dst2offset[isOdd] = malloc(PHYSICAL_SURFACE * sizeof(*g_bgq_is_dst2offset[isOdd]));
		g_bgq_ib_dst2offset[isOdd] = malloc(PHYSICAL_BODY * sizeof(*g_bgq_ib_dst2offset[isOdd]));
	}
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		bool isOdd_src = isOdd;
		bool isOdd_dst = !isOdd;

		size_t nextoffset[sec_end];
		for (bgq_weylfield_section sec = 0; sec < sec_end; sec += 1) {
			nextoffset[sec] = bgq_weyl_section_offset(sec);
		}

		for (size_t ih_dst = 0; ih_dst < PHYSICAL_VOLUME; ih_dst += 1) {
			size_t tv_dst = bgq_halfvolume2tv(ih_dst);
			size_t x_dst = bgq_halfvolume2x(ih_dst);
			size_t y_dst = bgq_halfvolume2y(ih_dst);
			size_t z_dst = bgq_halfvolume2z(ih_dst);
			bool isSurface = bgq_physical2isSurface(isOdd_dst, tv_dst, x_dst, y_dst, z_dst);
			assert(isSurface == (g_bgq_index_halfvolume2surface[isOdd_dst][ih_dst]!=(size_t)-1));
			assert(isSurface == (g_bgq_index_halfvolume2body[isOdd_dst][ih_dst]==(size_t)-1));
			bgq_weylfield_section mainsec = isSurface ? sec_surface : sec_body;

			if (ih_dst==126) {
				int g = 0;
			}

			for (size_t d_dst = 0; d_dst < PHYSICAL_LD; d_dst += 1) {
				bgq_dimension dim = bgq_direction2dimension(d_dst);
				bgq_weylfield_section sec = bgq_HoppingMatrix_init_source_sectionof_physical(isOdd_dst, tv_dst, x_dst, y_dst, z_dst, d_dst);

				// Reserve some offset in surface/body to ensure consecutive layout
				size_t thisOffset = nextoffset[mainsec];
				if (thisOffset/sizeof(bgq_weyl_vec)==32) {
					int d = 0;
				}
				assert(bgq_sectionOfOffset(thisOffset) == mainsec);
				assert(thisOffset == bgq_weyllayout_halfvolume2consecutiveoffset(isOdd_dst, ih_dst, d_dst));
				g_bgq_index2ih_dst[isOdd][thisOffset / sizeof(bgq_weyl_vec)] = ih_dst;
				g_bgq_index2d_dst[isOdd][thisOffset / sizeof(bgq_weyl_vec)] = d_dst;
				nextoffset[mainsec] += sizeof(bgq_weyl_vec);

				if (sec != mainsec) {
					// If in one of the mpi buffers, also reserve some space there
					assert((sec!=sec_surface) && (sec!=sec_body));

					thisOffset = nextoffset[sec];
					if (thisOffset/sizeof(bgq_weyl_vec)==32) {
						int d = 0;
					}
					g_bgq_index2ih_dst[isOdd][thisOffset / sizeof(bgq_weyl_vec)] = ih_dst;
					g_bgq_index2d_dst[isOdd][thisOffset / sizeof(bgq_weyl_vec)] = d_dst;
					nextoffset[sec] += sizeof(bgq_weyl_vec);

					// Sender part (this is on the inter-node sending the data to this node)
					size_t offset_sender = bgq_offset_recv2send(thisOffset);
					assert(bgq_sectionOfOffset(offset_sender) == bgq_section_commbuftran(sec,true));
					if (offset_sender/sizeof(bgq_weyl_vec)==32) {
						int d = 0;
					}
					g_bgq_index2ih_dst[isOdd][offset_sender / sizeof(bgq_weyl_vec)] = ih_dst;
					g_bgq_index2d_dst[isOdd][offset_sender / sizeof(bgq_weyl_vec)] = d_dst;
				}

				assert(bgq_sectionOfOffset(thisOffset) == sec);
				g_bgq_ih_dst2offset[isOdd][ih_dst].d[d_dst] = bgq_encode_offset(thisOffset);
				if (isSurface) {
					size_t is_dst = bgq_halfvolume2surface(isOdd_dst, ih_dst);
					g_bgq_is_dst2offset[isOdd][is_dst].d[d_dst] = bgq_encode_offset(thisOffset);
				} else {
					size_t ib_dst = bgq_halfvolume2body(isOdd_dst, ih_dst);
					g_bgq_ib_dst2offset[isOdd][ib_dst].d[d_dst] = bgq_encode_offset(thisOffset);
				}
			}
		}
	}

#ifndef NDEBUG
	// Check consistency
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		bool isOdd_dst = !isOdd;
		bool isOdd_src = isOdd;

		for (size_t is_dst = 0; is_dst < PHYSICAL_SURFACE; is_dst += 1) {
			for (size_t d_dst = 0; d_dst < PHYSICAL_LD; d_dst += 1) {
				size_t offset_is = g_bgq_is_dst2offset[isOdd][is_dst].d[d_dst];
				assert(offset_is);
				// For valgrind to check initialization

				size_t ih_dst = bgq_surface2halfvolume(isOdd_dst, is_dst);
				size_t offset_ih = g_bgq_ih_dst2offset[isOdd][ih_dst].d[d_dst];
				assert(offset_ih);

				assert(offset_is == offset_ih);
			}
		}
	}
#endif


	// Determine dest locations (where to write a weyl before communication)
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		g_bgq_ihsrc2offsetwrite[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_ihsrc2offsetwrite[isOdd]));
		g_bgq_issrc2offsetwrite[isOdd] = malloc(PHYSICAL_SURFACE * sizeof(*g_bgq_issrc2offsetwrite[isOdd]));
		g_bgq_ibsrc2offsetwrite[isOdd] = malloc(PHYSICAL_BODY * sizeof(*g_bgq_ibsrc2offsetwrite[isOdd]));

		size_t indices = bgq_weyl_section_offset(sec_end) / sizeof(bgq_weyl_vec);
		g_bgq_index2ih_src[isOdd] = malloc(indices * sizeof(*g_bgq_index2ih_src[isOdd]));
		g_bgq_index2d_src[isOdd] = malloc(indices * sizeof(*g_bgq_index2d_src[isOdd]));
	}
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		bool isOdd_src = isOdd;
		bool isOdd_dst = !isOdd;
		for (size_t ih_src = 0; ih_src < PHYSICAL_VOLUME; ih_src += 1) {
			for (size_t d_src = 0; d_src < PHYSICAL_LD; d_src += 1) {
				if (ih_src==6 && d_src==XUP) {
					int k = 0;
				}
				size_t offset_send = bgq_weylfield_destoffsetForWeyl(isOdd, ih_src, d_src);
				size_t index_send = offset_send / sizeof(bgq_weyl_vec);
				if (43200 == offset_send) {
					int b = 0;
				}

				assert(!bgq_section2isRecv(bgq_sectionOfOffset(offset_send)));

				size_t encodedOffset = bgq_encode_offset(offset_send);
				g_bgq_ihsrc2offsetwrite[isOdd][ih_src].d[d_src] = encodedOffset;
				bool isSurface = bgq_halfvolume2isSurface(isOdd_src, ih_src);
				if (isSurface) {
					size_t is_src = bgq_halfvolume2surface(isOdd_src, ih_src);
					g_bgq_issrc2offsetwrite[isOdd][is_src].d[d_src] = encodedOffset;
				} else {
					size_t ib_src = bgq_halfvolume2body(isOdd, ih_src);
					g_bgq_ibsrc2offsetwrite[isOdd][ib_src].d[d_src] = encodedOffset;
				}

				// Reverse lookup
				// 3 different offset calculations: send,recv,consecutive; may overlap some regions
				size_t index_write = offset_send / sizeof(bgq_weyl_vec);
				if (index_write==516) {
					int e = 0;
				}
				g_bgq_index2ih_src[isOdd][index_write] = ih_src;
				g_bgq_index2d_src[isOdd][index_write] = d_src;
				assert(d_src == bgq_offset2dsrc(offset_send));

				size_t offset_receiver = bgq_offset_send2recv(offset_send);
				if (18624 == offset_receiver) {
					int a = 0;
				}
				assert(!bgq_section2isRecv(offset_receiver));
				assert(bgq_offset_recv2send(offset_receiver) == offset_send);
				size_t index_receiver = offset_receiver / sizeof(bgq_weyl_vec);
				if (index_receiver==516) {
					int e = 0;
				}
				g_bgq_index2ih_src[isOdd][index_receiver] = ih_src;
				g_bgq_index2d_src[isOdd][index_receiver] = d_src;
				assert(d_src == bgq_offset2dsrc(offset_receiver));

				size_t tv_src = bgq_halfvolume2tv(ih_src);
				size_t t1_src = bgq_halfvolume2t1(isOdd_dst, ih_src);
				size_t t2_src = bgq_halfvolume2t2(isOdd_dst, ih_src);
				size_t x_src = bgq_halfvolume2x(ih_src);
				size_t y_src = bgq_halfvolume2y(ih_src);
				size_t z_src = bgq_halfvolume2z(ih_src);
				bgq_direction d_dst = bgq_direction_revert(d_src);
				size_t ih_dst = bgq_src2ih_dst(t1_src, x_src, y_src, z_src, d_src);
				assert(ih_dst == bgq_src2ih_dst(t2_src, x_src, y_src, z_src, d_src));
				size_t offset_consecutive = bgq_weyllayout_halfvolume2consecutiveoffset(isOdd_dst,ih_dst,d_dst);
				size_t index_consecutive = offset_consecutive / sizeof(bgq_weyl_vec);
				g_bgq_index2ih_src[isOdd][index_consecutive] = ih_src;
				g_bgq_index2d_src[isOdd][index_consecutive] = d_src;
			}
		}
	}

#ifndef NDEBUG
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		for (size_t offset = bgq_weyl_section_offset(0); offset < bgq_weyl_section_offset(sec_end); offset += sizeof(bgq_weyl_vec)) {
			size_t index = offset / sizeof(bgq_weyl_vec);
			if (index == 516) {
				int f = 0;
			}
			if (6336 == offset) {
				int b= 0;
			}
			if (18624 == offset) {
				int a = 0;
			}
			bgq_weylfield_section sec = bgq_sectionOfOffset(offset);
			if (bgq_section2isRecv(sec)) {
				size_t offsend = bgq_offset_recv2send(offset);
			}


			size_t ih_dst = g_bgq_index2ih_dst[isOdd][index];
			assert(0 <= ih_dst && ih_dst < PHYSICAL_VOLUME);
			bgq_direction d_dst = g_bgq_index2d_dst[isOdd][index];
			assert(d_dst == bgq_offset2ddst(offset));

			size_t ih_src = g_bgq_index2ih_src[isOdd][index];
			assert(0 <= ih_src && ih_src < PHYSICAL_VOLUME);
			bgq_direction d_src = g_bgq_index2d_src[isOdd][index];
			assert(d_src == bgq_offset2dsrc(offset));

			assert(d_dst == bgq_direction_revert(d_src));
			//assert(bgq_src2ih_dst());
		}
	}
#endif

}


//uint8_t *g_bgq_spinorfields_data = NULL;

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
	//g_bgq_spinorfields_data = malloc_aligned(tot_count * field_datasize, 128);
	//memset(g_bgq_spinorfields_data, -1, tot_count * field_datasize);

	for (size_t i = 0; i < tot_count; i += 1) { // Defer initialization to when it is used the first time
		g_bgq_spinorfields[i].isInitinialized = false;
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

	bgq_spinorfield_setup(field, field->isOdd, true, false, false, true);
	bgq_conversion_args args = {
			.field = field
	};
	bgq_master_call(&bgq_spinorfield_fullspinor2weyl_worker, &args); field->hasWeylfieldData = true;
}


static bgq_weyl_vec *bgq_offset2pointer(uint8_t *weylbase, size_t offset) {
	assert(weylbase);
	assert(offset);
	assert(bgq_weyl_section_offset(0) <= offset && offset < bgq_weyl_section_offset(sec_end));
	assert(offset % sizeof(bgq_weyl_vec) == 0);
	return (bgq_weyl_vec *)(weylbase + offset);
}


static bgq_weyl_vec *bgq_encodedoffset2pointer(uint8_t *weylbase, size_t code) {
	size_t offset = bgq_decode_offset(code);
	return bgq_offset2pointer(weylbase,offset);
}


void bgq_spinorfield_setup(bgq_weylfield_controlblock *field, bool isOdd, bool readFullspinor, bool writeFullspinor, bool readWeyl, bool writeWeyl) {
	assert(field);
	assert(readFullspinor || writeFullspinor || readWeyl || writeWeyl);
	// Do something

	//bool initDestPtrs;
	bool fullspinorAvailable;
	bool weylAvailable;
	if (field->isInitinialized) {
		if (field->isOdd == isOdd) {
			//initDestPtrs = false;
			fullspinorAvailable = field->hasFullspinorData;
			weylAvailable = field->hasWeylfieldData;
		} else {
			// Warning: field has already initialized for a different oddness
			// For optimal performance, always reuse fields with same oddness
			// TODO: Even and Odd fields can be made similar enough such that the pointers are the same? (e.g. Mirror T-Dimension)
			master_print("PERFORMANCE WARNING: Performance loss by reuse of spinorfield with different oddness\n");
			//initDestPtrs = true;
			fullspinorAvailable = false;
			weylAvailable = false;
		}
	} else {
		//initDestPtrs = true;
		fullspinorAvailable = false;
		weylAvailable = false;
		field->sec_weyl = NULL;
		field->sec_fullspinor = NULL;
	}

	if ((writeWeyl && !field->sec_weyl) || (field->sec_weyl && field->isOdd != isOdd)) {
		bool allocate = field->sec_weyl==NULL;

		size_t weylfieldsize = bgq_weyl_section_offset(sec_end);
		if (allocate) {
			size_t displacement = bgq_weyl_section_offset(0);
			field->sec_weyl = (uint8_t*)malloc_aligned(weylfieldsize - displacement, BGQ_ALIGNMENT_L2) + displacement/*so that offset points to invalid memory*/;
			//TODO: Also store a proper pointer to tell Valgrind we haven't forgotten this pointer
		}

		uint8_t *weylbase = field->sec_weyl;
		for (size_t d_cnt = 0; d_cnt < PHYSICAL_LD; d_cnt += 1) {
			bgq_direction d = d_cnt;
			field->sec_send[d] = (bgq_weyl_vec*) (weylbase + bgq_weyl_section_offset(bgq_direction2section(d,true)));
			field->sec_recv[d] = (bgq_weyl_vec*) (weylbase + bgq_weyl_section_offset(bgq_direction2section(d,false)));
		}
		field->sec_surface = (bgq_weylsite*) (weylbase + bgq_weyl_section_offset(sec_surface));
		field->sec_body = (bgq_weylsite*) (weylbase + bgq_weyl_section_offset(sec_body));
		field->sec_end = weylbase + bgq_weyl_section_offset(sec_end);

		// For 1st phase (distribute)
		if (allocate)
			field->destptrFromHalfvolume = malloc(PHYSICAL_VOLUME * sizeof(*field->destptrFromHalfvolume));
		for (size_t ih_src = 0; ih_src < PHYSICAL_VOLUME; ih_src += 1) {
			for (size_t d_src = 0; d_src < PHYSICAL_LD; d_src += 1) {
				field->destptrFromHalfvolume[ih_src].d[d_src] = bgq_encodedoffset2pointer(weylbase, g_bgq_ihsrc2offsetwrite[isOdd][ih_src].d[d_src]);
			}
		}

		if (allocate)
			field->destptrFromSurface = malloc(PHYSICAL_SURFACE * sizeof(*field->destptrFromSurface));
		for (size_t is_src = 0; is_src < PHYSICAL_SURFACE; is_src += 1) {
			for (size_t d_src = 0; d_src < PHYSICAL_LD; d_src += 1) {
				field->destptrFromSurface[is_src].d[d_src] = bgq_encodedoffset2pointer(weylbase, g_bgq_issrc2offsetwrite[isOdd][is_src].d[d_src]);
				assert(field->destptrFromSurface[is_src].d[d_src] == bgq_encodedoffset2pointer(weylbase, g_bgq_ihsrc2offsetwrite[isOdd][bgq_surface2halfvolume(isOdd, is_src)].d[d_src]));
			}
		}

		if (allocate)
			field->destptrFromBody = malloc(PHYSICAL_SURFACE * sizeof(*field->destptrFromBody));
		for (size_t ib_src = 0; ib_src < PHYSICAL_BODY; ib_src += 1) {
			for (size_t d_src = 0; d_src < PHYSICAL_LD; d_src += 1) {
				field->destptrFromBody[ib_src].d[d_src] = bgq_encodedoffset2pointer(weylbase, g_bgq_ibsrc2offsetwrite[isOdd][ib_src].d[d_src]);
				assert(field->destptrFromBody[ib_src].d[d_src] == bgq_encodedoffset2pointer(weylbase, g_bgq_ihsrc2offsetwrite[isOdd][bgq_body2halfvolume(isOdd, ib_src)].d[d_src]));
			}
		}

		// For 5th phase (datamove)
		if (COMM_T) {
			if (allocate)
				field->destptrFromTRecv = malloc((LOCAL_HALO_T/PHYSICAL_LP + 0 + LOCAL_HALO_T/PHYSICAL_LP) * sizeof(bgq_weyl_vec*));

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
		if (allocate)
			field->destptrFromRecv = malloc(weylCount * sizeof(bgq_weyl_vec*));
		for (size_t j = 0; j < weylCount; j+=1) {
			size_t offset = bgq_weyl_section_offset(sec_recv_xup) + j*sizeof(bgq_weyl_vec); // Overlaps into following recv sections
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

	}

	field->isInitinialized = true;
	field->isOdd = isOdd;
	field->isSloppy = false;

	if (writeFullspinor && !field->sec_fullspinor) {
		size_t fullfieldsize = PHYSICAL_VOLUME * sizeof(bgq_spinorsite);
		field->sec_fullspinor = malloc_aligned(fullfieldsize,BGQ_ALIGNMENT_L2);
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
			bgq_spinorfield_fullspinor2weyl(field); // Calls bgq_spinorfield_setup recursively
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


