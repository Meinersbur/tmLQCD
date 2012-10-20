/*
 * bgq_field.h
 *
 *  Created on: Oct 13, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_FIELD_H_
#define BGQ_FIELD_H_

#include "bgq_utils.h"

#include "../global.h"

#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>


#ifndef BGQ_FIELD_C_
#define EXTERN_INLINE EXTERN_INLINE_DECLARATION
#define EXTERN_FIELD extern
#define EXTERN_INIT(val)
#else
#define EXTERN_INLINE EXTERN_INLINE_DEFINITION
#define EXTERN_FIELD
#define EXTERN_INIT(val) = (val)
#endif



/* global logical coordinates */

#define GLOBAL_LT ((size_t)T_global)
#define GLOBAL_LX ((size_t)LX*(size_t)N_PROC_X)
#define GLOBAL_LY ((size_t)LY*(size_t)N_PROC_Y)
#define GLOBAL_LZ ((size_t)LZ*(size_t)N_PROC_Z)
#define GLOBAL_VOLUME (GLOBAL_LT*GLOBAL_LX*GLOBAL_LY*GLOBAL_LZ)


/* local logical coordinates */

#define LOCAL_LT ((size_t)T)
#define LOCAL_LX ((size_t)LX)
#define LOCAL_LY ((size_t)LY)
#define LOCAL_LZ ((size_t)LZ)
#define LOCAL_VOLUME ((size_t)VOLUME)
#define LOCAL_LD 8 /* Number of directions */

#define LOCAL_HALO_T (LOCAL_LX*LOCAL_LY*LOCAL_LZ)
#define LOCAL_HALO_X (LOCAL_LT*LOCAL_LY*LOCAL_LZ)
#define LOCAL_HALO_Y (LOCAL_LT*LOCAL_LX*LOCAL_LZ)
#define LOCAL_HALO_Z (LOCAL_LT*LOCAL_LX*LOCAL_LZ)


/* local physical coordinates */

#define PHYSICAL_LP 2 /* Even/Odd */
#define PHYSICAL_LTV (LOCAL_LT/(PHYSICAL_LP*PHYSICAL_LK))
#define PHYSICAL_LX LOCAL_LX
#define PHYSICAL_LY LOCAL_LY
#define PHYSICAL_LZ LOCAL_LZ
#define PHYSICAL_LK 2 /* Vector unit width (2 complex = 4 reals) */
#define PHYSICAL_LD 8/*10*/ /* Number of directions */
#define PHYSICAL_VOLUME (PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LY*PHYSICAL_LZ) /* (VOLUME/(PHYSICAL_LP*PHYSICAL_LK)) */
EXTERN_FIELD size_t PHYSICAL_BODY;
#define PHYSICAL_SURFACE (PHYSICAL_VOLUME-PHYSICAL_BODY)
#define PHYSICAL_HALO (COMM_T*2*LOCAL_HALO_T + COMM_X*2*LOCAL_HALO_X + COMM_Y*2*LOCAL_HALO_Y + COMM_Z*2*LOCAL_HALO_Z)

#define PHYSICAL_INDEX_LEXICAL(isOdd, tv, x, y, z) (tv + PHYSICAL_LTV*(x + PHYSICAL_LX*(y + PHYSICAL_LY*(z))))
#define PHYSICAL_LEXICAL2TV(isOdd,ih) ((ih)%PHYSICAL_LTV)
#define PHYSICAL_LEXICAL2X(isOdd,ih) (((ih)/PHYSICAL_LTV)%PHYSICAL_LX)
#define PHYSICAL_LEXICAL2Y(isOdd,ih) (((ih)/(PHYSICAL_LTV*PHYSICAL_LX))%PHYSICAL_LY)
#define PHYSICAL_LEXICAL2Z(isOdd,ih) ((ih)/(PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LY))
#define PHYSICAL_LEXICAL2T1(isOdd,ih) (4*PHYSICAL_LEXICAL2TV(isOdd,ih) + ((PHYSICAL_LEXICAL2X(isOdd,ih)+PHYSICAL_LEXICAL2Y(isOdd,ih)+PHYSICAL_LEXICAL2Z(isOdd,ih)+isOdd)&1))
#define PHYSICAL_LEXICAL2T2(isOdd,ih) (PHYSICAL_LEXICAL2T1(isOdd,ih) + 2)

#define PHYSICAL_HALO_T (COMM_T*PHYSICAL_LX*PHYSICAL_LY*PHYSICAL_LZ/PHYSICAL_LP)
#define PHYSICAL_HALO_X (COMM_X*PHYSICAL_LTV*PHYSICAL_LY*PHYSICAL_LZ)
#define PHYSICAL_HALO_Y (COMM_Y*PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LZ)
#define PHYSICAL_HALO_Z (COMM_Z*PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LY)



#if defined(PARALLELT)
#define COMM_T 1
#define COMM_X 0
#define COMM_Y 0
#define COMM_Z 0
#elif defined(PARALLELX)
#define COMM_T 0
#define COMM_X 1
#define COMM_Y 0
#define COMM_Z 0
#elif defined(PARALLELXT)
#define COMM_T 1
#define COMM_X 1
#define COMM_Y 0
#define COMM_Z 0
#elif defined(PARALLELXY)
#define COMM_T 0
#define COMM_X 1
#define COMM_Y 1
#define COMM_Z 0
#elif defined(PARALLELXYT)
#define COMM_T 1
#define COMM_X 1
#define COMM_Y 1
#define COMM_Z 0
#elif defined(PARALLELXYZ)
#define COMM_T 0
#define COMM_X 1
#define COMM_Y 1
#define COMM_Z 1
#elif defined(PARALLELXYZT)
#define COMM_T 1
#define COMM_X 1
#define COMM_Y 1
#define COMM_Z 1
#else
#error Unknown inter-node parallelization
#endif


#define COMPLEX_PRECISION complexdouble

typedef enum {
	TUP = 0,
	TDOWN = 1,
	XUP = 2,
	XDOWN = 3,
	YUP = 4,
	YDOWN = 5,
	ZUP = 6,
	ZDOWN = 7,

	DIR_UP = 0,
	DIR_DOWN = 1
} bgq_direction;

typedef enum {
	P_TUP1,
	P_TUP2,
	P_TDOWN1,
	P_TDOWN2,
	P_XUP,
	P_XDOWN,
	P_YUP,
	P_YDOWN,
	P_ZUP,
	P_ZDOWN
} bgq_physical_direction;
#define P_COUNT 10

typedef struct {
	uint32_t pd[P_COUNT];
} bgq_weyl_offsets_t;

typedef struct {
	void *pd[P_COUNT];
} bgq_weyl_ptr_t;

typedef enum {
	DIM_T,
	DIM_X,
	DIM_Y,
	DIM_Z
} bgq_dimension;




EXTERN_INLINE bgq_direction bgq_physical2direction(bgq_physical_direction pd) {
	switch (pd) {
	case P_TUP1:
	case P_TUP2:
		return TUP;
	case P_TDOWN1:
	case P_TDOWN2:
		return TDOWN;
	case P_XUP:
		return XUP;
	case P_XDOWN:
		return XDOWN;
	case P_YUP:
		return YUP;
	case P_YDOWN:
		return YDOWN;
	case P_ZUP:
		return ZUP;
	case P_ZDOWN:
		return ZDOWN;
	}
	assert(!"Unreachable");
	return -1;
}

void bgq_indices_init();
void bgq_spinorfields_init(size_t std_count, size_t chi_count);
void bgq_gaugefield_init();



EXTERN_INLINE bgq_dimension bgq_direction2dimension(bgq_direction d) {
	return d/2;
}


typedef struct {
	COMPLEX_PRECISION c[3][3][PHYSICAL_LK]; /* 3*3*2*sizeof(COMPLEX_PRECISION) = 288;144 bytes (4.5;2.25 L1 cache lines) */
	//COMPLEX_PRECISION padding[6];
} bgq_gaugesu3;
typedef struct {
	bgq_gaugesu3 su3[8];
} bgq_gaugesite;
typedef struct {
	bgq_gaugesu3 *(eodir[PHYSICAL_LP][PHYSICAL_LD]);
} bgq_gaugeeodir;
typedef bgq_gaugeeodir (*bgq_gaugefield);

typedef struct {
	COMPLEX_PRECISION s[4][3][PHYSICAL_LK]; /* 4*3*2*sizeof(COMPLEX_PRECISION) = 384;192 bytes (6;3 L1 cache lines) */
} bgq_spinorsite;
typedef bgq_spinorsite (*bgq_spinorfield);

typedef struct {
	COMPLEX_PRECISION s[2][3][PHYSICAL_LK]; // 96 byte (3 L1 cache lines)
} bgq_weyl_vec;

typedef struct {
	COMPLEX_PRECISION s[2][3]; // 48 byte
	COMPLEX_PRECISION _padding[2]; // 32 byte
} bgq_weyl_nonvec; // 64 byte (One L1 cache line)

typedef struct {
	bgq_weyl_nonvec tup1;
	bgq_weyl_nonvec tup2;
	bgq_weyl_nonvec tdown1;
	bgq_weyl_nonvec tdown2;
	bgq_weyl_vec xup;
	bgq_weyl_vec xdown;
	bgq_weyl_vec yup;
	bgq_weyl_vec ydown;
	bgq_weyl_vec zup;
	bgq_weyl_vec zdown;
} bgq_weylsite;
#if 0
EXTERN_FIELD size_t bgq_offsetof_weylsite[P_COUNT]
EXTERN_INIT((
		{offsetof(bgq_weylsite,tup1),offsetof(bgq_weylsite,tup2),offsetof(bgq_weylsite,tdown1),offsetof(bgq_weylsite,tdown2),
		offsetof(bgq_weylsite,xup),offsetof(bgq_weylsite,xdown),offsetof(bgq_weylsite,yup),offsetof(bgq_weylsite,ydown),offsetof(bgq_weylsite,zup),offsetof(bgq_weylsite,zdown)}
	));
#endif

EXTERN_INLINE void *bgq_weylsite_getdirection(bgq_weylsite *site, bgq_physical_direction d) {
	switch (d) {
	case P_TUP1:
		return &site->tup1;
	case P_TUP2:
		return &site->tup2;
	case P_TDOWN1:
		return &site->tdown1;
	case P_TDOWN2:
		return &site->tdown2;
	case P_XUP:
		return &site->xup;
	case P_XDOWN:
		return &site->xdown;
	case P_YUP:
		return &site->yup;
	case P_YDOWN :
		return &site->ydown;
	case P_ZUP:
		return &site->zup;
	case P_ZDOWN:
		return &site->zdown;
	}
	assert(!"Unreachable");
	exit(1);
}

typedef struct {
	COMPLEX_PRECISION s[2];
} bgq_weyl;

typedef struct {
	bool isSloppy;
	bool hasWeylfieldData;
	bool hasFullspinorData;
} bgq_weylfield_status;

typedef enum {
	//sec_status,
	sec_send_tup,
	sec_send_tdown,
	sec_send_xup,
	sec_send_xdown,
	sec_send_yup,
	sec_send_ydown,
	sec_send_zup,
	sec_send_zdown,
	sec_recv_tup,
	sec_recv_tdown,
	sec_recv_xup,
	sec_recv_xdown,
	sec_recv_yup,
	sec_recv_ydown,
	sec_recv_zup,
	sec_recv_zdown,
	sec_surface,
	sec_body,
	//sec_fullspinor,
	sec_end
} bgq_weylfield_section;

EXTERN_INLINE bgq_direction bgq_section2direction(bgq_weylfield_section sec) {
	switch (sec) {
	case sec_send_tup:
	case sec_recv_tup:
		return TUP;
	case sec_send_tdown:
	case sec_recv_tdown:
		return TDOWN;
	case sec_send_xup:
	case sec_recv_xup:
		return XUP;
	case sec_send_xdown:
	case sec_recv_xdown:
		return XDOWN;
	case sec_send_yup:
	case sec_recv_yup:
		return YUP;
	case sec_send_ydown:
	case sec_recv_ydown:
		return YDOWN;
	case sec_send_zup:
	case sec_recv_zup:
		return ZUP;
	case sec_send_zdown:
	case sec_recv_zdown:
		return ZDOWN;
	default:
		assert(!"This section has no direction");
		exit(1);
	}
}

typedef struct {
	bool isInitinialized;
	bool isOdd;
	bool isSloppy; // To be implemented
	bool hasWeylfieldData;
	bool hasFullspinorData;

	uint8_t *sec_weyl; // corresponds to offset 0 for the following fields
	bgq_weyl_nonvec *sec_send_tup;
	bgq_weyl_nonvec *sec_send_tdown;
	bgq_weyl_vec *sec_send_xup;
	bgq_weyl_vec *sec_send_xdown;
	bgq_weyl_vec *sec_send_yup;
	bgq_weyl_vec *sec_send_ydown;
	bgq_weyl_vec *sec_send_zup;
	bgq_weyl_vec *sec_send_zdown;
	bgq_weyl_nonvec *sec_recv_tup;
	bgq_weyl_nonvec *sec_recv_tdown;
	bgq_weyl_vec *sec_recv_xup;
	bgq_weyl_vec *sec_recv_xdown;
	bgq_weyl_vec *sec_recv_yup;
	bgq_weyl_vec *sec_recv_ydown;
	bgq_weyl_vec *sec_recv_zup;
	bgq_weyl_vec *sec_recv_zdown;
	bgq_weylsite *sec_surface;
	bgq_weylsite *sec_body;

	bgq_spinorsite *sec_fullspinor;

	//TODO: We may even interleave these with the data itself, but may cause alignment issues
	// Idea: sizeof(bgq_weyl_ptr_t)==10*8==80, so one bgq_weyl_ptr_t every 2(5;10) spinors solves the issue
	// In case we write as fullspinor layout, the are not needed
	bgq_weyl_ptr_t *destptrFromHalfvolume;
	bgq_weyl_ptr_t *destptrFromSurface;
	bgq_weyl_ptr_t *destptrFromBody;

	bgq_weyl_nonvec **destptrFromRecvTup;
	bgq_weyl_nonvec **destptrFromRecvTdown;
	bgq_weyl_vec **destptrFromRecvXup;
	bgq_weyl_vec **destptrFromRecvXdown;
	bgq_weyl_vec **destptrFromRecvYup;
	bgq_weyl_vec **destptrFromRecvYdown;
	bgq_weyl_vec **destptrFromRecvZup;
	bgq_weyl_vec **destptrFromRecvZdown;
} bgq_weylfield_controlblock;

// Index translations
EXTERN_FIELD size_t *g_bgq_index_surface2halfvolume[PHYSICAL_LP];
EXTERN_FIELD size_t *g_bgq_index_body2halfvolume[PHYSICAL_LP];
EXTERN_FIELD size_t *g_bgq_index_halfvolume2surface[PHYSICAL_LP]; // -1 if not surface
EXTERN_FIELD size_t *g_bgq_index_halfvolume2body[PHYSICAL_LP]; // -1 if not body
EXTERN_FIELD size_t *g_bgq_index_halfvolume2surfacebody[PHYSICAL_LP];

// Mapping of weyls to memory offsets
EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_offset_fromHalfvolume[PHYSICAL_LP];
EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_offset_fromSurface[PHYSICAL_LP];
EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_offset_fromBody[PHYSICAL_LP];

// Offsets of weyl when going into the corresponding direction
EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_destoffset_fromHalfvolume[PHYSICAL_LP];
EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_destoffset_fromSurface[PHYSICAL_LP];
EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_destoffset_fromBody[PHYSICAL_LP];

// The gaugefield as GAUGE_COPY
EXTERN_FIELD bgq_weylfield_controlblock *g_bgq_spinorfields EXTERN_INIT(NULL);
EXTERN_FIELD bgq_gaugesite *g_bgq_gaugefield_fromHalfvolume[PHYSICAL_LP];
EXTERN_FIELD bgq_gaugesite *g_bgq_gaugefield_fromSurface[PHYSICAL_LP];
EXTERN_FIELD bgq_gaugesite *g_bgq_gaugefield_fromBody[PHYSICAL_LP];


EXTERN_FIELD bool g_bgq_indices_initialized EXTERN_INIT(false);
void bgq_indices_init();
void bgq_spinorfields_init(size_t std_count, size_t chi_count);
void bgq_spinorfield_reset(bgq_weylfield_controlblock *field, bool isOdd, bool activateWeyl, bool activateFull);
void bgq_gaugefield_init();

EXTERN_INLINE size_t bgq_local2global_t(size_t t) {
	assert(0 <= t && t < LOCAL_LT);
	return LOCAL_LT*g_proc_coords[0] + t;
}
EXTERN_INLINE size_t bgq_local2global_x(size_t x) {
	assert(0 <= x && x < LOCAL_LX);
	return LOCAL_LX*g_proc_coords[1] + x;
}
EXTERN_INLINE size_t bgq_local2global_y(size_t y) {
	assert(0 <= y && y < LOCAL_LY);
	return LOCAL_LY*g_proc_coords[2] + y;
}
EXTERN_INLINE size_t bgq_local2global_z(size_t z) {
	assert(0 <= z && z < LOCAL_LZ);
	return LOCAL_LZ*g_proc_coords[3] + z;
}


EXTERN_INLINE size_t bgq_local2tv(size_t t, size_t x, size_t y, size_t z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	size_t tv = (t + LOCAL_LT - 2) / (PHYSICAL_LP*PHYSICAL_LK);
	assert(0 <= tv && tv < PHYSICAL_LTV);
	return tv;
}
EXTERN_INLINE bool bgq_local2k(size_t t, size_t x, size_t y, size_t z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	size_t teo = (t + LOCAL_LT - 2) / PHYSICAL_LP;
	size_t k = teo % PHYSICAL_LK;
	assert((k==0) || (k==1));
	return k;
}
EXTERN_INLINE bool bgq_local2eo(size_t t, size_t x, size_t y, size_t z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	return t&1;
}
EXTERN_INLINE bool bgq_local2isOdd(size_t t, size_t x, size_t y, size_t z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	return (t+x+y+z)%PHYSICAL_LP;
}
EXTERN_INLINE bool bgq_local2isSurface(size_t t, size_t x, size_t y, size_t z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);

	if (COMM_T && (t==0))
		return true;
	if (COMM_T && (t==LOCAL_LT-1))
		return true;
	if (COMM_X && (x==0))
		return true;
	if (COMM_X && (x==LOCAL_LX-1))
		return true;
	if (COMM_Y && (y==0))
		return true;
	if (COMM_Y && (y==LOCAL_LY-1))
		return true;
	if (COMM_Z && (z==0))
		return true;
	if (COMM_Z && (z==LOCAL_LZ-1))
		return true;
	return false;
}

EXTERN_INLINE bool bgq_physical2eo(bool isOdd, size_t tv, size_t x, size_t y, size_t z) {
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	return (x+y+z+isOdd)&1;
}

EXTERN_INLINE size_t bgq_physical2t1(bool isOdd, size_t tv, size_t x, size_t y, size_t z) {
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	size_t t1 = PHYSICAL_LP*PHYSICAL_LK*tv + bgq_physical2eo(isOdd,tv,x,y,z);
	assert(0 <= t1 && t1 < LOCAL_LT);
	return t1;
}
EXTERN_INLINE size_t bgq_physical2t2(bool isOdd, size_t tv, size_t x, size_t y, size_t z) {
	size_t t1 = bgq_physical2t1(isOdd, tv, x, y, z);
	return (t1+2)%LOCAL_LT;
}
EXTERN_INLINE size_t bgq_physical2halfvolume(size_t tv, size_t x, size_t y, size_t z) {
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= x && x < PHYSICAL_LX);
	assert(0 <= y && y < PHYSICAL_LY);
	assert(0 <= z && z < PHYSICAL_LZ);
	size_t ih = tv+PHYSICAL_LTV*(x+PHYSICAL_LX*(y+PHYSICAL_LY*(z)));
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	return ih;
}
EXTERN_INLINE size_t bgq_local2halfvolume(size_t t, size_t x, size_t y, size_t z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	size_t tv = t / (PHYSICAL_LP*PHYSICAL_LK);
	return bgq_physical2halfvolume(tv,x,y,z);
}
EXTERN_INLINE size_t bgq_halfvolume2tv(size_t ih) {
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	return ih%PHYSICAL_LTV;
}
EXTERN_INLINE size_t bgq_halfvolume2x(size_t ih) {
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	return (ih/PHYSICAL_LTV)%PHYSICAL_LX;
}
EXTERN_INLINE size_t bgq_halfvolume2y(size_t ih) {
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	return (ih/(PHYSICAL_LTV*PHYSICAL_LX))%PHYSICAL_LY;
}
EXTERN_INLINE size_t bgq_halfvolume2z(size_t ih) {
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	size_t result = ih/(PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LY);
	assert(result < PHYSICAL_LZ);
	return result;
}
EXTERN_INLINE size_t bgq_halfvolume2t1(bool isOdd, size_t ih) {
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	size_t tv = bgq_halfvolume2tv(ih);
	size_t x = bgq_halfvolume2x(ih);
	size_t y = bgq_halfvolume2y(ih);
	size_t z = bgq_halfvolume2z(ih);
	size_t t1 = bgq_physical2t1(isOdd,tv,x,y,z);
	return t1;
}
EXTERN_INLINE size_t bgq_halfvolume2t2(bool isOdd, size_t ih) {
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	size_t tv = bgq_halfvolume2tv(ih);
	size_t x = bgq_halfvolume2x(ih);
	size_t y = bgq_halfvolume2y(ih);
	size_t z = bgq_halfvolume2z(ih);
	size_t t2 = bgq_physical2t2(isOdd,tv,x,y,z);
	return t2;
}
EXTERN_INLINE bool bgq_halfvolume2isSurface(bool isOdd, size_t ih) {
	assert(g_bgq_indices_initialized);
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	return g_bgq_index_halfvolume2surface[isOdd][ih] != -1;
}
EXTERN_INLINE bool bgq_halfvolume2surface(bool isOdd, size_t ih) {
	assert(g_bgq_indices_initialized);
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	assert(g_bgq_index_halfvolume2surface[isOdd][ih] != -1);
	return g_bgq_index_halfvolume2surface[isOdd][ih];
}
EXTERN_INLINE size_t bgq_surface2halfvolume(bool isOdd, size_t is) {
	assert(g_bgq_indices_initialized);
	assert(0 <= is && is < PHYSICAL_SURFACE);
	size_t ih = g_bgq_index_surface2halfvolume[isOdd][is];
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	return ih;
}
EXTERN_INLINE size_t bgq_body2halfvolume(bool isOdd, size_t ib) {
	assert(g_bgq_indices_initialized);
	assert(0 <= ib && ib < PHYSICAL_BODY);
	size_t ih = g_bgq_index_body2halfvolume[isOdd][ib];
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	return ih;
}
EXTERN_INLINE size_t bgq_physical2surface(bool isOdd, size_t tv, size_t x, size_t y, size_t z) {
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= x && x < PHYSICAL_LX);
	assert(0 <= y && y < PHYSICAL_LY);
	assert(0 <= z && z < PHYSICAL_LZ);
	size_t ih = bgq_physical2halfvolume(tv,x,y,z);
	size_t is = bgq_halfvolume2surface(isOdd, ih);
	assert(0 <= is && is < PHYSICAL_SURFACE);
	return is;
}




EXTERN_INLINE size_t bgq_weyl_section_offset(bgq_weylfield_section section) {
	size_t result = 0;

	// TODO: loop-switch, more stable of order changes
	if (section > sec_send_tup) {
		result += PHYSICAL_HALO_T * sizeof(bgq_weyl_nonvec);
		result = (result + 127) & ~127;
	}
	if (section > sec_send_tdown) {
		result = PHYSICAL_HALO_T * sizeof(bgq_weyl_nonvec);
		result = (result + 127) & ~127;
	}
	if (section > sec_send_xup) {
		result += PHYSICAL_HALO_X * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_send_xdown) {
		result += PHYSICAL_HALO_X * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_send_yup) {
		result += PHYSICAL_HALO_Y * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_send_ydown) {
		result += PHYSICAL_HALO_Y * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_send_zup) {
		result += PHYSICAL_HALO_Z * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_send_zdown) {
		result += PHYSICAL_HALO_Z * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_recv_tup) {
		result += PHYSICAL_HALO_T * sizeof(bgq_weyl_nonvec);
		result = (result + 127) & ~127;
	}
	if (section > sec_recv_tdown) {
		result = PHYSICAL_HALO_T * sizeof(bgq_weyl_nonvec);
		result = (result + 127) & ~127;
	}
	if (section > sec_recv_xup) {
		result += PHYSICAL_HALO_X * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_recv_xdown) {
		result += PHYSICAL_HALO_X * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_recv_yup) {
		result += PHYSICAL_HALO_Y * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_recv_ydown) {
		result += PHYSICAL_HALO_Y * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_recv_zup) {
		result += PHYSICAL_HALO_Z * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_recv_zdown) {
		result += PHYSICAL_HALO_Z * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_surface) {
		result += PHYSICAL_SURFACE * sizeof(bgq_weylsite);
		result = (result + 127) & ~127;
	}
	if (section > sec_body) {
		result += PHYSICAL_BODY * sizeof(bgq_weylsite);
		result = (result + 127) & ~127;
	}
	return result;
}


static inline bool bgq_isOdd(size_t t, size_t x, size_t y, size_t z) {
	return (t+x+y+z)&1;
}


EXTERN_INLINE bool bgq_isCommSurface(size_t t, size_t x, size_t y, size_t z) {
#ifdef COMM_T
if ((t==0) || (t==LOCAL_LT-1))
	return true;
#endif
#ifdef COMM_X
if ((x==0) || (x==LOCAL_LX-1))
	return true;
#endif
#ifdef COMM_Y
if ((y==0) || (y==LOCAL_LY-1))
	return true;
#endif
#ifdef COMM_Z
if ((z==0) || (z==LOCAL_LZ-1))
	return true;
#endif
return false;
}


EXTERN_INLINE bool bgq_isCommBody(size_t t, size_t x, size_t y, size_t z) {
	return !bgq_isCommSurface(t,x,y,z);
}


// We compress offsets so they fit into an 32 bit integer
EXTERN_INLINE uint32_t bgq_encode_offset(size_t index) {
	assert((index & (32-1)) == 0);
	return index >> 5; // Always 32-bit aligned
}

EXTERN_INLINE size_t bgq_decode_offset(uint32_t code) {
	return (size_t)code << 5;
}



#undef EXTERN_INLINE
#undef EXTERN_FIELD
#undef EXTERN_INIT

#endif /* BGQ_FIELD_H_ */


