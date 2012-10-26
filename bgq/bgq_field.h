/*
 * bgq_field.h
 *
 *  Created on: Oct 13, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_FIELD_H_
#define BGQ_FIELD_H_
#include "../global.h"

#include "bgq_utils.h"
#include "bgq_qpx.h"

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
/*
A point in the lattice

Coordinate: (t,x,y,z) where
assert(0 <= t && t < GLOBAL_LT);
assert(0 <= x && x < GLOBAL_LX);
assert(0 <= y && y < GLOBAL_LY);
assert(0 <= z && z < GLOBAL_LZ);
Logically, it is a torus, therefore (t+i*GLOBAL_LT)==t (mod GLOBAL_LT) for every dimension
*/

#define GLOBAL_LT ((size_t)T_global)
#define GLOBAL_LX ((size_t)LX*(size_t)N_PROC_X)
#define GLOBAL_LY ((size_t)LY*(size_t)N_PROC_Y)
#define GLOBAL_LZ ((size_t)LZ*(size_t)N_PROC_Z)
#define GLOBAL_VOLUME (GLOBAL_LT*GLOBAL_LX*GLOBAL_LY*GLOBAL_LZ)


/* local logical coordinates */
/*
A point in the lattice stored on this MPI rank, relative to the top left 4D-rectangular coordinate on this rank

Spinor coordinate: (g_proc_coords[0],g_proc_coords[1],g_proc_coords[2],g_proc_coords[3]) (t,x,y,z) where
assert(0 <= t && t < LOCAL_LT);
assert(0 <= x && x < LOCAL_LX);
assert(0 <= y && y < LOCAL_LY);
assert(0 <= z && z < LOCAL_LZ);
Coordinates out of this rectangle may be stores on other MPI ranks
Note the the global lattice is a torus and not every dimension is needed to span multiple nodes

Weyl coordinate: (g_proc_coords[0],g_proc_coords[1],g_proc_coords[2],g_proc_coords[3]) (t,x,y,z,d) where
assert(0 <= t && t < LOCAL_LT);
assert(0 <= x && x < LOCAL_LX);
assert(0 <= y && y < LOCAL_LY);
assert(0 <= z && z < LOCAL_LZ);
assert(TUP <= d && z <= ZDOWN);
Specifies the weyl component when the spinor at (t,x,y,z) is decompositioned into its directions
*/

#define LOCAL_LT ((size_t)T)
#define LOCAL_LX ((size_t)LX)
#define LOCAL_LY ((size_t)LY)
#define LOCAL_LZ ((size_t)LZ)
#define LOCAL_LD 8 /* Number of directions */
#define LOCAL_VOLUME ((size_t)VOLUME)

#define LOCAL_HALO_T (LOCAL_LX*LOCAL_LY*LOCAL_LZ)
#define LOCAL_HALO_X (LOCAL_LT*LOCAL_LY*LOCAL_LZ)
#define LOCAL_HALO_Y (LOCAL_LT*LOCAL_LX*LOCAL_LZ)
#define LOCAL_HALO_Z (LOCAL_LT*LOCAL_LX*LOCAL_LZ)


/* local physical coordinates */
/*
Used to compute the memory location where specific lattice datum is stored

Spinor coordinate: (g_proc_coords[0],g_proc_coords[1],g_proc_coords[2],g_proc_coords[3]) (isOdd,tv,x,y,z,k) where
assert(false <= isOdd && isOdd <= true); assert(0 <= isOdd && isOdd < PHYSICAL_LP);
assert(0 <= tv && tv < PHYSICAL_LTV);
assert(0 <= x && x < PHYSICAL_LX);
assert(0 <= y && y < PHYSICAL_LY);
assert(0 <= z && z < PHYSICAL_LZ);
assert(0 <= k && k < PHYSICAL_LK);
Even and odd locations are stored independently (P stands for for parity)
Every memory location stores two logical lattice location, two locations in T-dimension are always processed together (V stands for Vector) These are usually called t1 and t2, for k==0 and k==1 respectively
Note that tv is shifted such that tv==PHYSICAL_LTV-1 are the two logical locations on the surface (i.e. t2==0 and t1==LOCAL_LT-1)
The other dimensions are handled like local coordinates


Weyl coordinate: (g_proc_coords[0],g_proc_coords[1],g_proc_coords[2],g_proc_coords[3]) (isOdd,tv,x,y,z,k,pd)
assert(false <= isOdd && isOdd <= true); assert(0 <= isOdd && isOdd < PHYSICAL_LP);
assert(0 <= tv && tv < PHYSICAL_LTV);
assert(0 <= x && x < PHYSICAL_LX);
assert(0 <= y && y < PHYSICAL_LY);
assert(0 <= z && z < PHYSICAL_LZ);
assert(0 <= k && k < PHYSICAL_LK);
assert(P_TUP1 <= pd && pd <= P_ZDOWN); assert(0 <= pd && pd < P_COUNT);
Denotes the weyl components used assemble the spinor in the next HoppingMatrix iteration (note the difference to the local coordinate system)
The t-directions are not vectorized here because they originate from different locations. Therefore, the k-coordinate is split up into two new "physical" dimension


Halfvolume coordinate: (g_proc_coords[0],g_proc_coords[1],g_proc_coords[2],g_proc_coords[3]) (isOdd,ih) where
assert(false <= isOdd && isOdd <= true); assert(0 <= isOdd && isOdd < PHYSICAL_LP);
assert(0 <= ih && ih < PHYSICAL_VOLUME);
The tv,x,y,z-coordinates linearized into the ih coordinate


Surface/Body coordinate: (g_proc_coords[0],g_proc_coords[1],g_proc_coords[2],g_proc_coords[3]) (isOdd,isSurface,is/ib) where
assert(false <= isSurface && isSurface <= true);
assert(0 <= is && is <= PHYSCIAL_SURFACE);
assert(0 <= ib && ib <= PHYSCIAL_BODY);
The spinor at the local lattice's surface are stored in a different memory block than the locations not at a border

*/

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

//#define PHYSICAL_INDEX_LEXICAL(isOdd, tv, x, y, z) (tv + PHYSICAL_LTV*(x + PHYSICAL_LX*(y + PHYSICAL_LY*(z))))
//#define PHYSICAL_LEXICAL2TV(isOdd,ih) ((ih)%PHYSICAL_LTV)
//#define PHYSICAL_LEXICAL2X(isOdd,ih) (((ih)/PHYSICAL_LTV)%PHYSICAL_LX)
//#define PHYSICAL_LEXICAL2Y(isOdd,ih) (((ih)/(PHYSICAL_LTV*PHYSICAL_LX))%PHYSICAL_LY)
//#define PHYSICAL_LEXICAL2Z(isOdd,ih) ((ih)/(PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LY))
//#define PHYSICAL_LEXICAL2T1(isOdd,ih) (4*PHYSICAL_LEXICAL2TV(isOdd,ih) + ((PHYSICAL_LEXICAL2X(isOdd,ih)+PHYSICAL_LEXICAL2Y(isOdd,ih)+PHYSICAL_LEXICAL2Z(isOdd,ih)+isOdd)&1))
//#define PHYSICAL_LEXICAL2T2(isOdd,ih) (PHYSICAL_LEXICAL2T1(isOdd,ih) + 2)

#define PHYSICAL_HALO_T (COMM_T*PHYSICAL_LX*PHYSICAL_LY*PHYSICAL_LZ/PHYSICAL_LP) //TODO: Remove? Interpretation not clear. With vectorization? without?
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

#define COMM_ORD_T 0
#define COMM_ORD_TSEND (COMM_ORD_T+0)
#define COMM_ORD_TRECV (COMM_ORD_T+1)
#define COMM_ORD_X (COMM_ORD_T+COMM_T)
#define COMM_ORD_XSEND (COMM_ORD_X+0)
#define COMM_ORD_XRECV (COMM_ORD_X+1)
#define COMM_ORD_Y (COMM_ORD_X+COMM_X)
#define COMM_ORD_YSEND (COMM_ORD_Y+0)
#define COMM_ORD_YRECV (COMM_ORD_Y+1)
#define COMM_ORD_Z (COMM_ORD_Y+COMM_Y)
#define COMM_ORD_ZSEND (COMM_ORD_Z+0)
#define COMM_ORD_ZRECV (COMM_ORD_X+1)

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

#if 0
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
#endif



typedef enum {
	DIM_T,
	DIM_X,
	DIM_Y,
	DIM_Z
} bgq_dimension;



#if 0
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
#endif

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
	COMPLEX_PRECISION s[2][3][PHYSICAL_LK]; // 192 byte (3 L1 cache lines)
} bgq_weyl_vec;

#if 0
typedef struct {
	COMPLEX_PRECISION s[2][3]; // 96 byte
	COMPLEX_PRECISION _padding[2]; // 32 byte
} bgq_weyl_nonvec; // 128 byte (2 L1 cache lines)
#endif

typedef struct {
	bgq_weyl_vec d[PHYSICAL_LD];
#if 0
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
#endif
} bgq_weylsite;
#if 0
EXTERN_FIELD size_t bgq_offsetof_weylsite[P_COUNT]
EXTERN_INIT((
		{offsetof(bgq_weylsite,tup1),offsetof(bgq_weylsite,tup2),offsetof(bgq_weylsite,tdown1),offsetof(bgq_weylsite,tdown2),
		offsetof(bgq_weylsite,xup),offsetof(bgq_weylsite,xdown),offsetof(bgq_weylsite,yup),offsetof(bgq_weylsite,ydown),offsetof(bgq_weylsite,zup),offsetof(bgq_weylsite,zdown)}
	));
#endif

typedef struct {
	//uint32_t pd[P_COUNT];
	uint32_t d[PHYSICAL_LD];
} bgq_weyl_offsets_t;

typedef struct {
	//void *pd[P_COUNT];
	bgq_weyl_vec *d[PHYSICAL_LD];
} bgq_weyl_ptr_t;

#if 0
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
#endif

typedef struct {
	COMPLEX_PRECISION s[2];
} bgq_weyl;

typedef struct {
	bool isSloppy;
	bool hasWeylfieldData;
	bool hasFullspinorData;
} bgq_weylfield_status;

typedef enum {
	// Yes, the order is important for COMM_X==1

	// k==0
	/* BEGIN consecutive */
	sec_send_tup,
	sec_recv_tup,
	/* END consecutive */

	// k==1
	/* BEGIN consecutive */
	sec_recv_tdown,
	sec_send_tdown,
	/* END consecutive */

	/* BEGIN consecutive */
	sec_send_xup,
	sec_send_xdown,
	sec_send_yup,
	sec_send_ydown,
	sec_send_zup,
	sec_send_zdown,
	/* END consecutive */

	// Recv sections are unused without inter-node communication
	// However, keep the here in order not to break code
	sec_recv_xup,
	sec_recv_xdown,
	sec_recv_yup,
	sec_recv_ydown,
	sec_recv_zup,
	sec_recv_zdown,

	sec_surface,
	sec_body,

	sec_end
} bgq_weylfield_section;//TODO: rename weyllayout

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
		UNREACHABLE
	}
}

EXTERN_INLINE bool bgq_section2isSend(bgq_weylfield_section sec) {
	switch (sec) {
		case sec_send_tup:
		case sec_send_tdown:
		case sec_send_xup:
		case sec_send_xdown:
		case sec_send_yup:
		case sec_send_ydown:
		case sec_send_zup:
		case sec_send_zdown:
			return true;
		case sec_recv_tup:
		case sec_recv_tdown:
		case sec_recv_xup:
		case sec_recv_xdown:
		case sec_recv_yup:
		case sec_recv_ydown:
		case sec_recv_zup:
		case sec_recv_zdown:
		return false;
		default:
			return false;
		}
}
EXTERN_INLINE bool bgq_section2isRecv(bgq_weylfield_section sec) {
	switch (sec) {
		case sec_send_tup:
		case sec_send_tdown:
		case sec_send_xup:
		case sec_send_xdown:
		case sec_send_yup:
		case sec_send_ydown:
		case sec_send_zup:
		case sec_send_zdown:
			return false;
		case sec_recv_tup:
		case sec_recv_tdown:
		case sec_recv_xup:
		case sec_recv_xdown:
		case sec_recv_yup:
		case sec_recv_ydown:
		case sec_recv_zup:
		case sec_recv_zdown:
		return true;
		default:
			return false;
		}
}


typedef struct {
	bool isInitinialized;
	bool isOdd;
	bool isSloppy; // To be implemented
	bool hasWeylfieldData;
	bool waitingForRecv; /* true==Need to wait for SPI recv and then copy data to consecutive area; false==All data available in sec_surface and sec_body */
	bool hasFullspinorData;

	uint8_t *sec_weyl; //TODO: can be made bgq_weyl_vec // corresponds to offset 0 for the following fields
	bgq_weyl_vec *sec_send[PHYSICAL_LD];
	bgq_weyl_vec *sec_recv[PHYSICAL_LD];
	bgq_weylsite *sec_surface;
	bgq_weylsite *sec_body;
	uint8_t *sec_end;

	bgq_spinorsite *sec_fullspinor;
	bgq_spinorsite *sec_fullspinor_surface;
	bgq_spinorsite *sec_fullspinor_body;

	//TODO: We may even interleave these with the data itself, but may cause alignment issues
	// Idea: sizeof(bgq_weyl_ptr_t)==10*8==80, so one bgq_weyl_ptr_t every 2(5;10) spinors solves the issue
	// In case we write as fullspinor layout, the are not needed
	bgq_weyl_ptr_t *destptrFromHalfvolume;
	bgq_weyl_ptr_t *destptrFromSurface;
	bgq_weyl_ptr_t *destptrFromBody;


	bgq_weyl_vec **destptrFromTRecv;
	bgq_weyl_vec **destptrFromRecv;

	//bgq_weyl_vec **destptrFromRecv[PHYSICAL_LD];
} bgq_weylfield_controlblock;

// Index translations
EXTERN_FIELD size_t *g_bgq_index_surface2halfvolume[PHYSICAL_LP];
EXTERN_FIELD size_t *g_bgq_index_body2halfvolume[PHYSICAL_LP];
EXTERN_FIELD size_t *g_bgq_index_halfvolume2surface[PHYSICAL_LP]; // -1 if not surface
EXTERN_FIELD size_t *g_bgq_index_halfvolume2body[PHYSICAL_LP]; // -1 if not body
EXTERN_FIELD size_t *g_bgq_index_halfvolume2surfacebody[PHYSICAL_LP];

// Mapping of dst weyls to memory offsets
EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_ih_dst2offset[PHYSICAL_LP];
EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_is_dst2offset[PHYSICAL_LP];
EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_ib_dst2offset[PHYSICAL_LP];
EXTERN_FIELD size_t *g_bgq_index2ih_dst[PHYSICAL_LP];
EXTERN_FIELD bgq_direction *g_bgq_index2d_dst[PHYSICAL_LP];

// Offsets of weyl when going into the corresponding direction
EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_ihsrc2offsetwrite[PHYSICAL_LP];
EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_issrc2offsetwrite[PHYSICAL_LP];
EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_ibsrc2offsetwrite[PHYSICAL_LP];
EXTERN_FIELD size_t *g_bgq_index2ih_src[PHYSICAL_LP];
EXTERN_FIELD bgq_direction *g_bgq_index2d_src[PHYSICAL_LP];


// The gaugefield as GAUGE_COPY
EXTERN_FIELD bgq_weylfield_controlblock *g_bgq_spinorfields EXTERN_INIT(NULL);
EXTERN_FIELD bgq_gaugesite *g_bgq_gaugefield_fromHalfvolume[PHYSICAL_LP];
EXTERN_FIELD bgq_gaugesite *g_bgq_gaugefield_fromSurface[PHYSICAL_LP];
EXTERN_FIELD bgq_gaugesite *g_bgq_gaugefield_fromBody[PHYSICAL_LP];


EXTERN_FIELD bool g_bgq_indices_initialized EXTERN_INIT(false);
void bgq_indices_init();
void bgq_spinorfields_init(size_t std_count, size_t chi_count);
void bgq_spinorfield_setup(bgq_weylfield_controlblock *field, bool isOdd, bool readFullspinor, bool writeFullspinor, bool readWeyl, bool writeWeyl);
//void bgq_spinorfield_reset(bgq_weylfield_controlblock *field, bool isOdd, bool activateWeyl, bool activateFull);
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


EXTERN_INLINE bool bgq_local2isOdd(size_t t, size_t x, size_t y, size_t z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	return (t+x+y+z)%PHYSICAL_LP;
}
EXTERN_INLINE bool bgq_physical2eo(bool isOdd, size_t tv, size_t x, size_t y, size_t z) {
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= x && x < PHYSICAL_LX);
	assert(0 <= y && y < PHYSICAL_LY);
	assert(0 <= z && z < PHYSICAL_LZ);
	return (x+y+z+isOdd)%PHYSICAL_LK;
}

EXTERN_INLINE size_t bgq_physical2t(bool isOdd, size_t tv, size_t x, size_t y, size_t z, size_t k) {
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= x && x < PHYSICAL_LX);
	assert(0 <= y && y < PHYSICAL_LY);
	assert(0 <= z && z < PHYSICAL_LZ);
	assert(0 <= k && k < PHYSICAL_LK);
	size_t t = PHYSICAL_LP*tv + k*LOCAL_LT/2 + bgq_physical2eo(isOdd,tv,x,y,z);
	assert(0 <= t && t < LOCAL_LT);
	return t;
}
EXTERN_INLINE size_t bgq_local2tv(size_t t, size_t x, size_t y, size_t z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	size_t t1 = t % (LOCAL_LT/PHYSICAL_LK); // Normalize to left virtual halflattice
	size_t tv = t1 / PHYSICAL_LP; // Remove isOdd information
	assert(0 <= tv && tv < PHYSICAL_LTV);
	return tv;
}
EXTERN_INLINE bool bgq_local2k(size_t t, size_t x, size_t y, size_t z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	size_t k = (t >= (LOCAL_LT/2));
	assert(t==bgq_physical2t(bgq_local2isOdd(t,x,y,z),bgq_local2tv(t,x,y,z),x,y,z,k));
	return k;
}
EXTERN_INLINE bool bgq_local2eo(size_t t, size_t x, size_t y, size_t z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	return t%PHYSICAL_LP;
}


//Do not use! "isSurface" is a property of a physical (vectorized) site
EXTERN_INLINE bool bgq_local2isSurface_raw(size_t t, size_t x, size_t y, size_t z) {
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


EXTERN_INLINE size_t bgq_physical2t1(bool isOdd, size_t tv, size_t x, size_t y, size_t z) {
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= x && x < PHYSICAL_LX);
	assert(0 <= y && y < PHYSICAL_LY);
	assert(0 <= z && z < PHYSICAL_LZ);
	size_t t1 = bgq_physical2t(isOdd,tv,x,y,z,0);
	assert(0 <= t1 && t1 < LOCAL_LT/2);
	return t1;
}
EXTERN_INLINE size_t bgq_physical2t2(bool isOdd, size_t tv, size_t x, size_t y, size_t z) {
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= x && x < PHYSICAL_LX);
	assert(0 <= y && y < PHYSICAL_LY);
	assert(0 <= z && z < PHYSICAL_LZ);
	size_t t2 = bgq_physical2t(isOdd,tv,x,y,z,1);
	assert(LOCAL_LT/2 <= t2 && t2 < LOCAL_LT);
	return t2;
}

EXTERN_INLINE bool bgq_physical2isSurface(bool isOdd, size_t tv, size_t x, size_t y, size_t z) {
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= x && x < PHYSICAL_LX);
	assert(0 <= y && y < PHYSICAL_LY);
	assert(0 <= z && z < PHYSICAL_LZ);
	size_t t1 = bgq_physical2t1(isOdd,tv,x,y,z);
	size_t t2 = bgq_physical2t2(isOdd,tv,x,y,z);
	bool eo = bgq_physical2eo(isOdd,tv,x,y,z);
	bool isSurface = false;
	if (COMM_T)
		isSurface = isSurface || (!eo && (tv==0)) || (eo && (tv==PHYSICAL_LTV-1));
	if (COMM_X)
		isSurface = isSurface || (x==0) || (x==PHYSICAL_LX-1);
	if (COMM_Y)
		isSurface = isSurface || (y==0) || (y==PHYSICAL_LY-1);
	if (COMM_Z)
		isSurface = isSurface || (z==0) || (z==PHYSICAL_LZ-1);
	assert(isSurface == (bgq_local2isSurface_raw(t1,x,y,z) || bgq_local2isSurface_raw(t2,x,y,z)));
	return isSurface;
}
EXTERN_INLINE bool bgq_local2isSurface(size_t t, size_t x, size_t y, size_t z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	bool isOdd = bgq_local2isOdd(t,x,y,z);
	size_t tv = bgq_local2tv(t,x,y,z);
	return bgq_physical2isSurface(isOdd,tv,x,y,z);
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
	size_t tv = bgq_local2tv(t,x,y,z);
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
	size_t tv = bgq_halfvolume2tv(ih);
	size_t x = bgq_halfvolume2x(ih);
	size_t y = bgq_halfvolume2y(ih);
	size_t z = bgq_halfvolume2z(ih);
	size_t t2 = bgq_physical2t2(isOdd,tv,x,y,z);
	return t2;
}
EXTERN_INLINE size_t bgq_halfvolume2t(bool isOdd, size_t ih, size_t k) {
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	assert(0 <= k && k < PHYSICAL_LK);
	size_t tv = bgq_halfvolume2tv(ih);
	size_t x = bgq_halfvolume2x(ih);
	size_t y = bgq_halfvolume2y(ih);
	size_t z = bgq_halfvolume2z(ih);
	size_t t = bgq_physical2t(isOdd, tv, x,y,z,k);
	return t;
}
EXTERN_INLINE bool bgq_halfvolume2isSurface(bool isOdd, size_t ih) {
	assert(g_bgq_indices_initialized);
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	return g_bgq_index_halfvolume2surface[isOdd][ih] != (size_t)-1;
}
EXTERN_INLINE size_t bgq_halfvolume2surface(bool isOdd, size_t ih) {
	assert(g_bgq_indices_initialized);
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	size_t is2 = g_bgq_index_halfvolume2surface[isOdd][ih];
	assert(0 <= is2 && is2 < PHYSICAL_SURFACE);
	return is2;
}
EXTERN_INLINE size_t bgq_halfvolume2body(bool isOdd, size_t ih) {
	assert(g_bgq_indices_initialized);
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	size_t ib = g_bgq_index_halfvolume2body[isOdd][ih];
	assert(0 <= ib && ib < PHYSICAL_BODY);
	return ib;
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
EXTERN_INLINE size_t bgq_surface2tv(bool isOdd, size_t is) {
	assert(0 <= is && is < PHYSICAL_SURFACE);
	size_t ih = bgq_surface2halfvolume(isOdd, is);
	return bgq_halfvolume2tv(ih);
}

EXTERN_INLINE size_t bgq_halfvolume2volume(bool isOdd, size_t ih, size_t k) {
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	assert(0 <= k && k < PHYSICAL_LK);
	return isOdd + PHYSICAL_LP*(k + PHYSICAL_LK*ih);
}
EXTERN_INLINE size_t bgq_volume2halfvolume(size_t iv) {
	assert(0 <= iv && iv < LOCAL_VOLUME);
	return iv / (PHYSICAL_LP*PHYSICAL_LK);
}
EXTERN_INLINE size_t bgq_volume2isOdd(size_t iv) {
	assert(0 <= iv && iv < LOCAL_VOLUME);
	return iv % PHYSICAL_LP;
}
EXTERN_INLINE size_t bgq_volume2k(size_t iv) {
	assert(0 <= iv && iv < LOCAL_VOLUME);
	return (iv / PHYSICAL_LP) % PHYSICAL_LK;
}
EXTERN_INLINE size_t bgq_volume2t(size_t iv) {
	assert(0 <= iv && iv < LOCAL_VOLUME);
	bool isOdd = bgq_volume2isOdd(iv);
	size_t k = bgq_volume2k(iv);
	size_t ih = bgq_volume2halfvolume(iv);
	size_t tv = bgq_halfvolume2tv(ih);
	size_t x = bgq_halfvolume2x(ih);
	size_t y = bgq_halfvolume2y(ih);
	size_t z = bgq_halfvolume2z(ih);
	return bgq_physical2t(isOdd, tv, x, y, z, k);
}
EXTERN_INLINE size_t bgq_volume2x(size_t iv) {
	assert(0 <= iv && iv < LOCAL_VOLUME);
	size_t ih = bgq_volume2halfvolume(iv);
	size_t x = bgq_halfvolume2x(ih);
	return x;
}
EXTERN_INLINE size_t bgq_volume2y(size_t iv) {
	assert(0 <= iv && iv < LOCAL_VOLUME);
	size_t ih = bgq_volume2halfvolume(iv);
	size_t y = bgq_halfvolume2y(ih);
	return y;
}
EXTERN_INLINE size_t bgq_volume2z(size_t iv) {
	assert(0 <= iv && iv < LOCAL_VOLUME);
	size_t ih = bgq_volume2halfvolume(iv);
	size_t z = bgq_halfvolume2z(ih);
	return z;
}
EXTERN_INLINE size_t bgq_local2volume(size_t t, size_t x, size_t y, size_t z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	bool isOdd = bgq_local2isOdd(t,x,y,z);
	size_t k = bgq_local2k(t,x,y,z);
	size_t ih = bgq_local2halfvolume(t,x,y,z);
	return bgq_halfvolume2volume(isOdd, ih, k);
}

EXTERN_INLINE bgq_weylfield_section bgq_direction2section(bgq_direction d, bool isSend) {
	switch (d) {
	case TUP:
		return isSend ? sec_send_tup : sec_recv_tup;
	case TDOWN:
		return isSend ? sec_send_tdown : sec_recv_tdown;
	case XUP:
		return isSend ? sec_send_tup : sec_recv_tup;
	case XDOWN:
		return isSend ? sec_send_tdown : sec_recv_tdown;
	case YUP:
		return isSend ? sec_send_tup : sec_recv_tup;
	case YDOWN:
		return isSend ? sec_send_tdown : sec_recv_tdown;
	case ZUP:
		return isSend ? sec_send_tup : sec_recv_tup;
	case ZDOWN:
		return isSend ? sec_send_tdown : sec_recv_tdown;
	default:
		UNREACHABLE
	}
}


EXTERN_INLINE size_t bgq_weyl_section_offset(bgq_weylfield_section section) {
	size_t result = sizeof(bgq_weyl_vec); // Do not use the first entry; "0" should signify an error

	for (bgq_weylfield_section sec = 0; sec < sec_end; sec += 1) {
		if (section == sec)
			return result;

		size_t secsize;
		switch (sec) {
		case sec_send_tup:
			case sec_send_tdown:
			case sec_recv_tup:
			case sec_recv_tdown:
			secsize = LOCAL_HALO_T * sizeof(bgq_weyl_vec)/PHYSICAL_LP;
			break;
		case sec_send_xup:
			case sec_send_xdown:
			case sec_recv_xup:
			case sec_recv_xdown:
			secsize = PHYSICAL_HALO_X * sizeof(bgq_weyl_vec);
			break;
		case sec_send_yup:
			case sec_send_ydown:
			case sec_recv_yup:
			case sec_recv_ydown:
			secsize = PHYSICAL_HALO_Y * sizeof(bgq_weyl_vec);
			break;
		case sec_send_zup:
			case sec_send_zdown:
			case sec_recv_zup:
			case sec_recv_zdown:
			secsize = PHYSICAL_HALO_Z * sizeof(bgq_weyl_vec);
			break;
		case sec_surface:
			secsize = PHYSICAL_SURFACE * sizeof(bgq_weylsite);
			break;
		case sec_body:
			secsize = PHYSICAL_BODY * sizeof(bgq_weylsite);
			break;
		default:
			assert(!"Never should get here");
			secsize = 0;
			break;
		}

		result += secsize;
		//result = (result + (BGQ_ALIGNMENT_L2-1)) & ~(BGQ_ALIGNMENT_L2-1); // Padding for alignment
	}
	assert(section==sec_end);
	return result;
}


EXTERN_INLINE bgq_weylfield_section bgq_sectionOfOffset(size_t offset) {
	for (bgq_weylfield_section sec = 0; sec < sec_end; sec+=1) {
		if ((bgq_weyl_section_offset(sec) <= offset) && (offset < bgq_weyl_section_offset(sec+1)))
			return sec;
	}
	assert(!"Out of range");
	return sec_end;
}







EXTERN_INLINE bgq_direction bgq_direction_revert(bgq_direction d) {
	return d^1;
}


// We compress offsets so they fit into an 32 bit integer
EXTERN_INLINE uint32_t bgq_encode_offset(size_t index) { assert(index+1);
	assert((index & (32-1)) == 0);
	return index >> 5; // Always 32-bit aligned
}


EXTERN_INLINE size_t bgq_decode_offset(uint32_t code) { assert(code+1);
	return (size_t)code << 5;
}


void bgq_spinorfield_transfer(bool isOdd, bgq_weylfield_controlblock *targetfield, spinor* sourcefield);

#undef EXTERN_INLINE
#undef EXTERN_FIELD
#undef EXTERN_INIT

#endif /* BGQ_FIELD_H_ */


