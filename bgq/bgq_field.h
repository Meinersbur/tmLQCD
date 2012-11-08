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

typedef uint_fast32_t ucoord;
typedef int_fast32_t scoord;
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

#define GLOBAL_LT ((ucoord)T_global)
#define GLOBAL_LX ((ucoord)LX*(ucoord)N_PROC_X)
#define GLOBAL_LY ((ucoord)LY*(ucoord)N_PROC_Y)
#define GLOBAL_LZ ((ucoord)LZ*(ucoord)N_PROC_Z)
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

#define LOCAL_LT ((ucoord)T)
#define LOCAL_LX ((ucoord)LX)
#define LOCAL_LY ((ucoord)LY)
#define LOCAL_LZ ((ucoord)LZ)
#define LOCAL_LD 8 /* Number of directions */
#define LOCAL_VOLUME ((ucoord)VOLUME)

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

#define COMMDIR_COUNT 2*(COMM_T+COMM_X+COMM_Y+COMM_Z)

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
	DIM_T,
	DIM_X,
	DIM_Y,
	DIM_Z
} bgq_dimension;


EXTERN_INLINE bgq_dimension bgq_direction2dimension(bgq_direction d) {
	return d/2;
}


EXTERN_INLINE bool bgq_dimension_isDistributed(bgq_dimension dim) {
	switch (dim) {
	case DIM_T:
		return COMM_T;
	case DIM_X:
		return COMM_X;
	case DIM_Y:
		return COMM_Y;
	case DIM_Z:
		return COMM_Z;
	default:
		UNREACHABLE
	}
}

EXTERN_INLINE bgq_direction bgq_dimenstion2direction_up(bgq_dimension dim) {
	switch (dim) {
	case DIM_T:
		return TUP;
	case DIM_X:
		return XUP;
	case DIM_Y:
		return YUP;
	case DIM_Z:
		return ZUP;
	default:
		UNREACHABLE
	}
}
EXTERN_INLINE bgq_direction bgq_dimenstion2direction_down(bgq_dimension dim) {
	switch (dim) {
	case DIM_T:
		return TDOWN;
	case DIM_X:
		return XDOWN;
	case DIM_Y:
		return YDOWN;
	case DIM_Z:
		return ZDOWN;
	default:
		UNREACHABLE
	}
}


EXTERN_INLINE bool bgq_direction_isDistributed(bgq_direction d) {
	return bgq_dimension_isDistributed(bgq_direction2dimension(d));
}



void bgq_indices_init(void);
void bgq_spinorfields_init(size_t std_count, size_t chi_count);
void bgq_gaugefield_init(void);






typedef struct {
	COMPLEX_PRECISION s[4][3][PHYSICAL_LK]; /* 4*3*2*sizeof(COMPLEX_PRECISION) = 384;192 bytes (6;3 L1 cache lines) */
} bgq_spinorsite;
typedef bgq_spinorsite bgq_spinor_vec;
typedef bgq_spinorsite (*bgq_spinorfield);

typedef struct {
	COMPLEX_PRECISION s[2][3][PHYSICAL_LK]; // 192 byte (3 L1 cache lines)
} bgq_weyl_vec;



typedef struct {
	bgq_weyl_vec d[PHYSICAL_LD];
} bgq_weylsite;

typedef struct {
	uint32_t d[PHYSICAL_LD];
} bgq_weyl_offsets_t;

typedef struct {
	bgq_weyl_vec *d[PHYSICAL_LD];
} bgq_weyl_ptr_t;


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

	sec_surface,
	sec_body,

	// k==0
	/* BEGIN consecutive */
	//sec_send_tup,
	//sec_recv_tup,
	/* END consecutive */

	// k==1
	/* BEGIN consecutive */
	//sec_recv_tdown,
	//sec_send_tdown,
	/* END consecutive */

	sec_send_tup,
	sec_send_tdown,
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
	sec_recv_tup,
	sec_recv_tdown,
	/* BEGIN consecutive */
	sec_recv_xup,
	sec_recv_xdown,
	sec_recv_yup,
	sec_recv_ydown,
	sec_recv_zup,
	sec_recv_zdown,
	/* END consecutive */

	sec_end,

	sec_collapsed=sec_surface,
	sec_collapsed_end=sec_send_tup,
	sec_comm=sec_send_tup,
	sec_comm_end=sec_end,
	sec_send_begin=sec_send_tup,
	sec_send_end=sec_send_zdown+1,
	sec_recv_begin=sec_recv_tup,
	sec_recv_end=sec_recv_zdown+1
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


typedef enum {
	hm_nocom = 1 << 0,
	hm_nooverlap = 1 << 1,
	hm_nokamul = 1 << 2,
	hm_fixedoddness = 1 << 3, // obsolete

	hm_noprefetchexplicit = 1 << 4, // obsolete
	hm_noprefetchlist = 1 << 5,
	hm_noprefetchstream = 1 << 6,

	hm_noweylsend = 1 << 7, // obsolete
	hm_nobody = 1 << 8,
	hm_nosurface = 1 << 9, // obsolete (->hm_nodistribute)

	hm_l1pnonstoprecord = 1 << 10,
	hm_experimental = 1 << 11, // obsolete

	hm_prefetchimplicitdisable = 1 << 12,
	hm_prefetchimplicitoptimistic = 2 << 12,
	hm_prefetchimplicitconfirmed = 3 << 12,

	hm_withcheck = 1 << 14,
	hm_nodistribute = 1 << 15,
	hm_nodatamove = 1 << 16,
	hm_nospi = 1 << 17
} bgq_hmflags;


//TODO: make incomplete type
typedef struct {
	bool isInitinialized;
	bool isOdd;
	bool isSloppy; // To be implemented
	bool hasWeylfieldData;
	bool waitingForRecv; /* true==Need to wait for SPI recv and then copy data to consecutive area; false==All data available in sec_surface and sec_body */
	//bool waitingForRecvNoSPI;
	bgq_hmflags hmflags;
	bool pendingDatamove;
	bool hasFullspinorData;

	uint8_t *sec_weyl;
	bgq_weyl_vec *sec_index; // obsolete
	bgq_weyl_vec *sec_send[PHYSICAL_LD]; // obsolete
	bgq_weyl_vec *sec_recv[PHYSICAL_LD]; // obsolete
	bgq_weylsite *sec_collapsed;
	bgq_weylsite *sec_surface;
	bgq_weylsite *sec_body;
	uint8_t *sec_end;

	bgq_spinorsite *sec_fullspinor;
	bgq_spinorsite *sec_fullspinor_surface;
	bgq_spinorsite *sec_fullspinor_body;

	//TODO: We may even interleave these with the data itself, but may cause alignment issues
	// Idea: sizeof(bgq_weyl_ptr_t)==10*8==80, so one bgq_weyl_ptr_t every 2(5;10) spinors solves the issue
	// In case we write as fullspinor layout, the are not needed
	bgq_weyl_ptr_t *destptrFromHalfvolume; // obsolete
	bgq_weyl_ptr_t *destptrFromSurface; // obsolete
	bgq_weyl_ptr_t *destptrFromBody; // obsolete
	bgq_weyl_ptr_t *sendptr;


	bgq_weyl_vec **consptr_recvtup;
	bgq_weyl_vec **consptr_recvtdown;
	//bgq_weyl_vec **destptrFromTRecv;
	bgq_weyl_vec **destptrFromRecv;

	//bgq_weyl_vec **destptrFromRecv[PHYSICAL_LD];
} bgq_weylfield_controlblock;

// Index translations
EXTERN_FIELD size_t *g_bgq_collapsed2halfvolume[PHYSICAL_LP];
EXTERN_FIELD size_t *g_bgq_halfvolume2collapsed[PHYSICAL_LP];

EXTERN_FIELD size_t *g_bgq_index_surface2halfvolume[PHYSICAL_LP]; // deprecated
EXTERN_FIELD size_t *g_bgq_index_body2halfvolume[PHYSICAL_LP];// deprecated
EXTERN_FIELD size_t *g_bgq_index_halfvolume2surface[PHYSICAL_LP];// deprecated // -1 if not surface
EXTERN_FIELD size_t *g_bgq_index_halfvolume2body[PHYSICAL_LP]; // deprecated// -1 if not body
EXTERN_FIELD size_t *g_bgq_index_halfvolume2surfacebody[PHYSICAL_LP];// deprecated
//EXTERN_FIELD size_t *g_bgq_index_ordered2halfvolume[PHYSICAL_LP];

// Mapping of dst weyls to memory offsets
EXTERN_FIELD size_t *g_bgq_index2collapsed[PHYSICAL_LP];
//EXTERN_FIELD bgq_direction *g_bgq_index2d[PHYSICAL_LP];
EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_collapsed2indexrecv[PHYSICAL_LP];
EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_collapsed2indexsend[PHYSICAL_LP];

EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_ih_dst2offset[PHYSICAL_LP]; // deprecated
EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_is_dst2offset[PHYSICAL_LP]; // deprecated
EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_ib_dst2offset[PHYSICAL_LP]; // deprecated
EXTERN_FIELD size_t *g_bgq_index2ih_dst[PHYSICAL_LP]; // deprecated
EXTERN_FIELD bgq_direction *g_bgq_index2d_dst[PHYSICAL_LP]; // deprecated

// Offsets of weyl when going into the corresponding direction
EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_ihsrc2offsetwrite[PHYSICAL_LP];// deprecated
EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_issrc2offsetwrite[PHYSICAL_LP];// deprecated
EXTERN_FIELD bgq_weyl_offsets_t *g_bgq_ibsrc2offsetwrite[PHYSICAL_LP];// deprecated
EXTERN_FIELD size_t *g_bgq_index2ih_src[PHYSICAL_LP];// deprecated
EXTERN_FIELD bgq_direction *g_bgq_index2d_src[PHYSICAL_LP];// deprecated


// The gaugefield as GAUGE_COPY
EXTERN_FIELD bgq_weylfield_controlblock *g_bgq_spinorfields EXTERN_INIT(NULL);



EXTERN_FIELD bool g_bgq_indices_initialized EXTERN_INIT(false);
void bgq_indices_init(void);
void bgq_spinorfields_init(size_t std_count, size_t chi_count);
//void bgq_spinorfield_reset(bgq_weylfield_controlblock *field, bool isOdd, bool activateWeyl, bool activateFull);
void bgq_gaugefield_init(void);

EXTERN_INLINE ucoord bgq_local2global_t(scoord t) {
	assert(0 <= t && t < LOCAL_LT);
	return LOCAL_LT*g_proc_coords[0] + t;
}
EXTERN_INLINE ucoord bgq_local2global_x(scoord x) {
	assert(0 <= x && x < LOCAL_LX);
	return LOCAL_LX*g_proc_coords[1] + x;
}
EXTERN_INLINE ucoord bgq_local2global_y(scoord y) {
	assert(0 <= y && y < LOCAL_LY);
	return LOCAL_LY*g_proc_coords[2] + y;
}
EXTERN_INLINE ucoord bgq_local2global_z(scoord z) {
	assert(0 <= z && z < LOCAL_LZ);
	return LOCAL_LZ*g_proc_coords[3] + z;
}


EXTERN_INLINE ucoord bgq_t2t(ucoord t, ucoord k) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= k && k < PHYSICAL_LK);
	size_t result = (t % (LOCAL_LT/PHYSICAL_LK)) + k*(LOCAL_LT/PHYSICAL_LK);
	assert(0 <= result && result < LOCAL_LT);
	return result;
}


EXTERN_INLINE bool bgq_local2isOdd(size_t t, size_t x, size_t y, size_t z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	return (t+x+y+z)%PHYSICAL_LP;
}
EXTERN_INLINE bool bgq_physical2eo(bool isOdd, ucoord tv, ucoord x, ucoord y, ucoord z) {
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= x && x < PHYSICAL_LX);
	assert(0 <= y && y < PHYSICAL_LY);
	assert(0 <= z && z < PHYSICAL_LZ);
	return (x+y+z+isOdd)%PHYSICAL_LK;
}

EXTERN_INLINE ucoord bgq_physical2t(bool isOdd, ucoord tv, ucoord x, ucoord y, ucoord z, ucoord k) {
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= x && x < PHYSICAL_LX);
	assert(0 <= y && y < PHYSICAL_LY);
	assert(0 <= z && z < PHYSICAL_LZ);
	assert(0 <= k && k < PHYSICAL_LK);
	ucoord t = PHYSICAL_LP*tv + k*LOCAL_LT/2 + bgq_physical2eo(isOdd,tv,x,y,z);
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
EXTERN_INLINE ucoord bgq_halfvolume2tv(ucoord ih) {
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	return ih%PHYSICAL_LTV;
}
EXTERN_INLINE ucoord bgq_halfvolume2x(ucoord ih) {
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	return (ih/PHYSICAL_LTV)%PHYSICAL_LX;
}
EXTERN_INLINE ucoord bgq_halfvolume2y(ucoord ih) {
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	return (ih/(PHYSICAL_LTV*PHYSICAL_LX))%PHYSICAL_LY;
}
EXTERN_INLINE ucoord bgq_halfvolume2z(ucoord ih) {
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	ucoord result = ih/(PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LY);
	assert(0 <= result && result < PHYSICAL_LZ);
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
EXTERN_INLINE ucoord bgq_halfvolume2t(bool isOdd, ucoord ih, ucoord k) {
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	assert(0 <= k && k < PHYSICAL_LK);
	ucoord tv = bgq_halfvolume2tv(ih);
	ucoord x = bgq_halfvolume2x(ih);
	ucoord y = bgq_halfvolume2y(ih);
	ucoord z = bgq_halfvolume2z(ih);
	ucoord t = bgq_physical2t(isOdd, tv, x,y,z,k);
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

EXTERN_INLINE size_t bgq_halfvolume2collapsed(bool isOdd, size_t ih) {
	assert(g_bgq_indices_initialized);
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	ucoord ic = g_bgq_halfvolume2collapsed[isOdd][ih];
	assert(0 <= ic && ic < PHYSICAL_VOLUME);
	return ic;
}

EXTERN_INLINE bool bgq_collapsed2isSurface(size_t ic) {
	assert(0 <= ic && ic < PHYSICAL_VOLUME);
	return ic < PHYSICAL_SURFACE;
}

EXTERN_INLINE size_t bgq_collapsed2halfvolume(bool isOdd, size_t ic) {
	assert(g_bgq_indices_initialized);
	assert(0 <= ic && ic < PHYSICAL_VOLUME);
	ucoord ih = g_bgq_collapsed2halfvolume[isOdd][ic];
	assert(0 <= ih && ih < PHYSICAL_VOLUME);
	return ih;
}

EXTERN_INLINE size_t bgq_surface2collapsed(size_t is) {
	assert(0 <= is && is < PHYSICAL_SURFACE);
	return is;
}
EXTERN_INLINE size_t bgq_surface2halfvolume(bool isOdd, size_t is) {
	assert(0 <= is && is < PHYSICAL_SURFACE);
	//size_t ih = g_bgq_index_surface2halfvolume[isOdd][is];
	//assert(0 <= ih && ih < PHYSICAL_VOLUME);
	ucoord ic = bgq_surface2collapsed(is);
	ucoord ih = bgq_collapsed2halfvolume(isOdd, ic) ;
	return ih;
}
EXTERN_INLINE size_t bgq_surface2tv(bool isOdd, size_t is) {
	assert(0 <= is && is < PHYSICAL_SURFACE);
	size_t ih = bgq_surface2halfvolume(isOdd, is);
	return bgq_halfvolume2tv(ih);
}

EXTERN_INLINE size_t bgq_collapsed2surface(size_t ic) {
	assert(0 <= ic && ic < PHYSICAL_VOLUME);
	assert(bgq_collapsed2isSurface(ic));
	return ic;
}

EXTERN_INLINE size_t bgq_body2collapsed(size_t ib) {
	assert(0 <= ib && ib < PHYSICAL_BODY);
	size_t ic = PHYSICAL_SURFACE + ib;
	assert(0 <= ic && ic < PHYSICAL_VOLUME);
	return ic;
}
EXTERN_INLINE size_t bgq_collapsed2body(size_t ic) {
	assert(0 <= ic && ic < PHYSICAL_VOLUME);
	assert(!bgq_collapsed2isSurface(ic));
	return ic - PHYSICAL_SURFACE;
}

EXTERN_INLINE size_t bgq_local2collapsed(size_t t, size_t x, size_t y, size_t z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	bool isOdd = bgq_local2isOdd(t,x,y,z);
	size_t ih = bgq_local2halfvolume(t,x,y,z);
	return bgq_halfvolume2collapsed(isOdd, ih);
}

EXTERN_INLINE ucoord bgq_collapsed2t(bool isOdd, ucoord ic, ucoord k) {
	assert(0 <= ic && ic < PHYSICAL_VOLUME);
	assert(0 <= k && k < PHYSICAL_LK);
	ucoord ih = bgq_collapsed2halfvolume(isOdd, ic);
	return bgq_halfvolume2t(isOdd, ih, k);
}
EXTERN_INLINE size_t bgq_collapsed2x(bool isOdd, size_t ic) {
	assert(0 <= ic && ic < PHYSICAL_VOLUME);
	size_t ih = bgq_collapsed2halfvolume(isOdd, ic);
	return bgq_halfvolume2x(ih);
}
EXTERN_INLINE size_t bgq_collapsed2y(bool isOdd, size_t ic) {
	assert(0 <= ic && ic < PHYSICAL_VOLUME);
	size_t ih = bgq_collapsed2halfvolume(isOdd, ic);
	return bgq_halfvolume2y(ih);
}
EXTERN_INLINE size_t bgq_collapsed2z(bool isOdd, size_t ic) {
	assert(0 <= ic && ic < PHYSICAL_VOLUME);
	size_t ih = bgq_collapsed2halfvolume(isOdd, ic);
	return bgq_halfvolume2z(ih);
}


EXTERN_INLINE size_t bgq_local2halfvolume_neighbor(size_t t, size_t x, size_t y, size_t z, bgq_direction d) {
	switch (d) {
	case TUP:
		t = (t + 1) % LOCAL_LT;
		break;
	case TDOWN:
		t = (t + LOCAL_LT - 1) % LOCAL_LT;
		break;
	case XUP:
		x = (x + 1) % LOCAL_LX;
		break;
	case XDOWN:
		x = (x + LOCAL_LX - 1) % LOCAL_LX;
		break;
	case YUP:
		y = (y + 1) % LOCAL_LY;
		break;
	case YDOWN:
		y = (y + LOCAL_LY - 1) % LOCAL_LY;
		break;
	case ZUP:
		z = (z + 1) % LOCAL_LZ;
		break;
	case ZDOWN:
		z = (z + LOCAL_LZ - 1) % LOCAL_LZ;
		break;
	}

	size_t ih_src = bgq_local2halfvolume(t,x,y,z);
	return ih_src;
}

EXTERN_INLINE size_t bgq_localdst2ksrc(size_t t, size_t x, size_t y, size_t z, bgq_direction d) {
	switch (d) {
	case TUP:
		t = (t + 1) % LOCAL_LT;
		break;
	case TDOWN:
		t = (t + LOCAL_LT - 1) % LOCAL_LT;
		break;
	case XUP:
		x = (x + 1) % LOCAL_LX;
		break;
	case XDOWN:
		x = (x + LOCAL_LX - 1) % LOCAL_LX;
		break;
	case YUP:
		y = (y + 1) % LOCAL_LY;
		break;
	case YDOWN:
		y = (y + LOCAL_LY - 1) % LOCAL_LY;
		break;
	case ZUP:
		z = (z + 1) % LOCAL_LZ;
		break;
	case ZDOWN:
		z = (z + LOCAL_LZ - 1) % LOCAL_LZ;
		break;
	}

	size_t k_src = bgq_local2k(t,x,y,z);
	return k_src;
}


EXTERN_INLINE bgq_weylfield_section bgq_direction2section(bgq_direction d, bool isSend) {
	switch (d) {
	case TUP:
		return isSend ? sec_send_tup : sec_recv_tup;
	case TDOWN:
		return isSend ? sec_send_tdown : sec_recv_tdown;
	case XUP:
		return isSend ? sec_send_xup : sec_recv_xup;
	case XDOWN:
		return isSend ? sec_send_xdown : sec_recv_xdown;
	case YUP:
		return isSend ? sec_send_yup : sec_recv_yup;
	case YDOWN:
		return isSend ? sec_send_ydown : sec_recv_ydown;
	case ZUP:
		return isSend ? sec_send_zup : sec_recv_zup;
	case ZDOWN:
		return isSend ? sec_send_zdown : sec_recv_zdown;
	default:
		UNREACHABLE
	}
}


EXTERN_INLINE size_t bgq_weyl_section_offset(bgq_weylfield_section section) {
	size_t result = 0;

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
EXTERN_INLINE size_t bgq_spinorfield_indexOfSection(bgq_weylfield_section sec) {
	return bgq_weyl_section_offset(sec)/sizeof(bgq_weyl_vec);
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

EXTERN_INLINE size_t bgq_weyllayout_collapsed2consecutiveoffset(bool isOdd, ucoord ic, bgq_direction d) {
	return bgq_weyl_section_offset(sec_collapsed) + ic*sizeof(bgq_weylsite) + d*sizeof(bgq_weyl_vec);
}
size_t bgq_weyllayout_halfvolume2consecutiveoffset(bool isOdd_dst, size_t ih_dst, bgq_direction d_dst);
bgq_direction bgq_offset2ddst(size_t offset);
bgq_direction bgq_offset2dsrc(size_t offset);
size_t bgq_src2ih_dst(size_t t_src, size_t x_src, size_t y_src, size_t z_src, bgq_direction d_src);
size_t bgq_src2k_dst(size_t t_src, size_t x_src, size_t y_src, size_t z_src, bgq_direction d_src);


EXTERN_INLINE void bgq_direction_move_local(ucoord *t, ucoord *x, ucoord *y, ucoord *z, bgq_direction d) {
	switch (d) {
	case TDOWN:
		*t = (*t + LOCAL_LT - 1) % LOCAL_LT;
		break;
	case TUP:
		*t = (*t + 1) % LOCAL_LT;
		break;
	case XDOWN:
		*x = (*x + LOCAL_LX - 1) % LOCAL_LX;
		break;
	case XUP:
		*x = (*x + 1) % LOCAL_LX;
		break;
	case YDOWN:
		*y = (*y + LOCAL_LY - 1) % LOCAL_LY;
		break;
	case YUP:
		*y = (*y + 1) % LOCAL_LY;
		break;
	case ZDOWN:
		*z = (*z + LOCAL_LZ - 1) % LOCAL_LZ;
		break;
	case ZUP:
		*z = (*z + 1) % LOCAL_LZ;
		break;
	}
}


EXTERN_INLINE void bgq_direction_move_physical(bool isOdd, ucoord *tv, ucoord *x, ucoord *y, ucoord *z, bgq_direction d) {
	switch (d) {
	case TDOWN:
		*tv = bgq_physical2eo(isOdd, *tv, *x, *y, *z) ? (*tv) : ((*tv + PHYSICAL_LTV - 1) % PHYSICAL_LTV);
		break;
	case TUP:
		*tv = bgq_physical2eo(isOdd, *tv, *x, *y, *z) ? ((*tv + 1) % PHYSICAL_LTV) : (*tv);
		break;
	case XDOWN:
		*x = (*x + LOCAL_LX - 1) % LOCAL_LX;
		break;
	case XUP:
		*x = (*x + 1) % LOCAL_LX;
		break;
	case YDOWN:
		*y = (*y + LOCAL_LY - 1) % LOCAL_LY;
		break;
	case YUP:
		*y = (*y + 1) % LOCAL_LY;
		break;
	case ZDOWN:
		*z = (*z + LOCAL_LZ - 1) % LOCAL_LZ;
		break;
	case ZUP:
		*z = (*z + 1) % LOCAL_LZ;
		break;
	}
}


EXTERN_INLINE void bgq_direction_move_global(ucoord *t, ucoord *x, ucoord *y, ucoord *z, bgq_direction d) {
	switch (d) {
	case TDOWN:
		*t = (*t + GLOBAL_LT - 1) % GLOBAL_LT;
		break;
	case TUP:
		*t = (*t + 1) % GLOBAL_LT;
		break;
	case XDOWN:
		*x = (*x + GLOBAL_LX - 1) % GLOBAL_LX;
		break;
	case XUP:
		*x = (*x + 1) % GLOBAL_LX;
		break;
	case YDOWN:
		*y = (*y + GLOBAL_LY - 1) % GLOBAL_LY;
		break;
	case YUP:
		*y = (*y + 1) % GLOBAL_LY;
		break;
	case ZDOWN:
		*z = (*z + GLOBAL_LZ - 1) % GLOBAL_LZ;
		break;
	case ZUP:
		*z = (*z + 1) % GLOBAL_LZ;
		break;
	}
}


EXTERN_INLINE ucoord bgq_offset2index(size_t offset) {
	assert(bgq_weyl_section_offset(0) <= offset && offset <= bgq_weyl_section_offset(sec_end));
	assert(offset % sizeof(bgq_weyl_vec) == 0);
	return offset / sizeof(bgq_weyl_vec);
}


EXTERN_INLINE size_t bgq_index2offset(ucoord index) {
	size_t offset = index * sizeof(bgq_weyl_vec);
	assert(bgq_weyl_section_offset(0) <= offset && offset < bgq_weyl_section_offset(sec_end));
	return offset;
}


EXTERN_INLINE ucoord bgq_halfvolume_src2dst(bool isOdd_src, ucoord ih_src, bgq_direction d_src) {
	assert(0 <= ih_src && ih_src < PHYSICAL_VOLUME);
	ucoord tv_src = bgq_halfvolume2tv(ih_src);
	ucoord x_src = bgq_halfvolume2x(ih_src);
	ucoord y_src = bgq_halfvolume2y(ih_src);
	ucoord z_src = bgq_halfvolume2z(ih_src);

	ucoord tv_dst = tv_src;
	ucoord x_dst = x_src;
	ucoord y_dst = y_src;
	ucoord z_dst = z_src;
	bgq_direction_move_physical(isOdd_src, &tv_dst, &x_dst, &y_dst, &z_dst, d_src);
	ucoord ih_dst = bgq_physical2halfvolume(tv_dst, x_dst, y_dst, z_dst);
	return ih_dst;
}


EXTERN_INLINE ucoord bgq_collapsed_src2dst(bool isOdd_src, ucoord ic_src, bgq_direction d_src) {
	assert(0 <= ic_src && ic_src < PHYSICAL_VOLUME);
	bool isOdd_dst = !isOdd_src;
	ucoord ih_src = bgq_collapsed2halfvolume(isOdd_src, ic_src);
	ucoord ih_dst = bgq_halfvolume_src2dst(isOdd_src, ih_src, d_src);
	return bgq_halfvolume2collapsed(isOdd_dst, ih_dst);
}



size_t bgq_pointer2offset(bgq_weylfield_controlblock *field, void *ptr) ;


EXTERN_FIELD size_t g_bgq_spinorfields_count;


EXTERN_INLINE size_t bgq_section_size(bgq_weylfield_section sec) {
	return bgq_weyl_section_offset(sec+1) - bgq_weyl_section_offset(sec);
}

EXTERN_INLINE size_t bgq_sectionrange_size(bgq_weylfield_section sec_begin, bgq_weylfield_section sec_stop/*inclusive*/) {
	return bgq_weyl_section_offset(sec_stop+1) - bgq_weyl_section_offset(sec_stop);
}


#undef EXTERN_INLINE
#undef EXTERN_FIELD
#undef EXTERN_INIT

#endif /* BGQ_FIELD_H_ */


