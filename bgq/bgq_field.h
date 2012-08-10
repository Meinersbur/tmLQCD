#ifndef BGQ_H_INCLUDED
#define BGQ_H_INCLUDED

#include "bgq_utils.h"
#include "complex_c99.h"
#include "../global.h"
#include "../read_input.h"
#include <stdbool.h>
#include <assert.h>


#ifndef BGQ_FIELD_C_H_
#define EXTERN_INLINE inline
#define EXTERN_FIELD extern
#else
#define EXTERN_INLINE extern inline
#define EXTERN_FIELD
#endif

#define VECTOR_WIDTH 2 /* 2 complex values = 4 reals */

#define LOCAL_LP 2 /* Even/Odd */
#define LOCAL_LT T
#define LOCAL_LX LX
#define LOCAL_LY LY
#define LOCAL_LZ LZ

#define PHYSICAL_LP 2 /* Even/Odd */
#define PHYSICAL_LTV (LOCAL_LT/(PHYSICAL_LP*PHYSICAL_LK))
#define PHYSICAL_LX LOCAL_LX
#define PHYSICAL_LXV (LOCAL_LX/(PHYSICAL_LP*PHYSICAL_LK))
#define PHYSICAL_LY LOCAL_LY
#define PHYSICAL_LZ LOCAL_LZ
#define PHYSICAL_LK 2 /* Vector unit width (2 complex = 4 reals) */
#define PHYSICAL_LD 5 /* No of directions (for gauge field) = 4 up-directions + 1 for shifted t-line (note that it is one vector-element longer for wraparound values) */



// faces, 2 of each
#define SURFACE_T_SITES ((1)*(PHYSICAL_LX-2)*(PHYSICAL_LY-2)*(PHYSICAL_LZ))
#define SURFACE_X_SITES ((PHYSICAL_LTV-2)*(1)*(PHYSICAL_LY-2)*(PHYSICAL_LZ))
#define SURFACE_Y_SITES ((PHYSICAL_LTV-2)*(PHYSICAL_LX-2)*(1)*(PHYSICAL_LZ))
#define SURFACE_FACE_SITES (2*SURFACE_T_SITES + 2*SURFACE_X_SITES + 2*SURFACE_Y_SITES)

#define SURFACE_T_ZLINES ((1)*(PHYSICAL_LX-2)*(PHYSICAL_LY-2))
#define SURFACE_X_ZLINES ((PHYSICAL_LTV-2)*(1)*(PHYSICAL_LY-2))
#define SURFACE_Y_ZLINES ((PHYSICAL_LTV-2)*(PHYSICAL_LX-2)*(1))
#define SURFACE_FACE_ZLINES (2*SURFACE_T_ZLINES + 2*SURFACE_X_ZLINES + 2*SURFACE_Y_ZLINES)

// edges, 4 of each
#define SURFACE_TX_SITES ((1)*(1)*(PHYSICAL_LY-2)*(PHYSICAL_LZ))
#define SURFACE_TY_SITES ((1)*(PHYSICAL_LX-2)*(1)*(PHYSICAL_LZ))
#define SURFACE_XY_SITES ((PHYSICAL_LTV-2)*(1)*(1)*(PHYSICAL_LZ))
#define SURFACE_EDGE_SITES (4*SURFACE_TX_SITES + 4*SURFACE_TY_SITES + 4*SURFACE_XY_SITES)

#define SURFACE_TX_ZLINES ((1)*(1)*(PHYSICAL_LY-2))
#define SURFACE_TY_ZLINES ((1)*(PHYSICAL_LX-2)*(1))
#define SURFACE_XY_ZLINES ((PHYSICAL_LTV-2)*(1)*(1))
#define SURFACE_EDGE_ZLINES (4*SURFACE_TX_ZLINES + 4*SURFACE_TY_ZLINES + 4*SURFACE_XY_ZLINES)

// vertices, 8
#define SURFACE_TXY_SITES ((1)*(1)*(1)*(PHYSICAL_LZ))
#define SURFACE_VERTICE_SITES (8*SURFACE_TXY_SITES)

#define SURFACE_TXY_ZLINES ((1)*(1)*(1))
#define SURFACE_VERTICE_ZLINES (8*SURFACE_TXY_ZLINES)

// total
#define SURFACE_SITES (SURFACE_FACE_SITES+SURFACE_FACE_SITES+SURFACE_EDGE_SITES)
#define BODY_SITES ((PHYSICAL_LTV-2)*(PHYSICAL_LX-2)*(PHYSICAL_LY-2)*(PHYSICAL_LZ))
#define VOLUME_SITES (VOLUME/(PHYSICAL_LP*PHYSICAL_LK))

#define SURFACE_ZLINES (SURFACE_FACE_ZLINES+SURFACE_FACE_ZLINES+SURFACE_EDGE_ZLINES)
#define BODY_ZLINES ((PHYSICAL_LTV-2)*(PHYSICAL_LX-2)*(PHYSICAL_LY-2))


#define GAUGE_VOLUME ((LOCAL_LT+1)*(LOCAL_LX+1)*(LOCAL_LY+1)*(LOCAL_LZ)) /* LOCAL volume */
#define GAUGE_EOVOLUME ((PHYSICAL_LTV+1)*(PHYSICAL_LX+1)*(PHYSICAL_LY+1)*(PHYSICAL_LZ+1/*for wraparound*/)) /* This is not hole-free; we could also define an accessor function per direction*/


typedef struct {
	double _Complex s[4][3][PHYSICAL_LK];
} bgq_spinorsite_double;
typedef struct {
	double _Complex s[2][3][PHYSICAL_LK];
} bgq_weylsite_double;
typedef bgq_spinorsite_double (*bgq_spinorfield_double);
typedef bgq_weylsite_double (*bgq_weylfield_double);

typedef struct {
	double _Complex c[3][3][PHYSICAL_LK];
} bgq_gaugesite_double;
typedef struct {
	bgq_gaugesite_double *(eodir[PHYSICAL_LP][PHYSICAL_LD]);
} bgq_gaugeeodir_double;
typedef bgq_gaugeeodir_double (*bgq_gaugefield_double);


typedef enum {
	TUP = 0,
	TDOWN = 1,
	XUP = 2,
	XDOWN = 3,
	YUP = 4,
	YDOWN = 5,
	ZUP = 6,
	ZDOWN = 7,

	TUP_SHIFT = 8,
	TDOWN_SHIFT = 9,

	DIR_UP = 0,
	DIR_DOWN = 1
} direction;


void *malloc_aligned(size_t size, size_t alignment);



////////////////////////////////////////////////////////////////////////////////
// Spinorfields

#ifndef NDEBUG
EXTERN_FIELD bgq_spinorsite_double *g_spinorfields_doubledata_coords;
#endif
EXTERN_FIELD bgq_spinorsite_double *g_spinorfields_doubledata;
EXTERN_FIELD bgq_spinorfield_double *g_spinorfields_double;

void bgq_init_spinorfields(int count);
void bgq_free_spinofields();

EXTERN_INLINE bgq_spinorfield_double bgq_translate_spinorfield(spinor * const field) {
	const int V = even_odd_flag ? VOLUMEPLUSRAND / 2 : VOLUMEPLUSRAND;
	const int fieldsize  = V * sizeof(*field);
// This computes the original index address of the passed field; be aware that its correctness depends on the implementation of init_spinorfield
	long long offset = (char*)field - (char*)g_spinor_field[0];
	assert(offset >= 0);
	assert(offset % fieldsize == 0);
	int result = offset  / fieldsize;
	assert(g_spinorfields_double[result] >= g_spinorfields_doubledata);
	return g_spinorfields_double[result];
}

EXTERN_FIELD int g_num_spinorfields;
EXTERN_FIELD int *g_spinorfield_isOdd;


void bgq_transfer_spinorfield(bool isOdd, bgq_spinorfield_double targetfield, spinor *sourcefield);
void bgq_spinorfield_resetcoord(bgq_spinorfield_double spinorfield, bool isOdd, int expected_reads_min, int expected_reads_max, int expected_writes_min, int expected_writes_max);





bool assert_spinorfield_coord(bgq_spinorfield_double spinorfield, bool isOdd, int t, int x, int y, int z, int tv, int k, int v, int c, bool isRead, bool isWrite);
bool assert_spinorcoord(bgq_spinorfield_double spinorfield, bool isOdd, int t, int x, int y, int z, int tv, int k, bool isRead, bool isWrite);

#define BGQ_SPINORSITE_ACCESS(spinorfield, isOdd, tv, x, y, z) \
	(assert(tv <= 0 && tv < PHYSICAL_LTV),      \
	 assert(x <= 0 && x < PHYSICAL_LX),         \
	 assert(y <= 0 && y < PHYSICAL_LY),         \
	 assert(z <= 0 && z < PHYSICAL_LZ),         \
	 &spinorfield[(((tv)*PHYSICAL_LX + (x))*PHYSICAL_LY + (y))*PHYSICAL_LZ + (z)])

#define BGQ_SPINORSITE(spinorfield, isOdd, tv, x, y, z, t1, t2, isWrite)                         \
		(assert(assert_spinorcoord(spinorfield, isOdd, tv, x, y, z, t1, 0, !isWrite, isWrite)), \
		 assert(assert_spinorcoord(spinorfield, isOdd, tv, x, y, z, t2, 1, !isWrite, isWrite)), \
		 BGQ_SPINORSITE_ACCESS(spinorfield, isOdd, tv, x, y, z))

#define BGQ_SPINORSITE_LEFT(spinorfield, isOdd, tv, x, y, z, t1, t2, isWrite)                         \
		(assert(assert_spinorcoord(spinorfield, isOdd, tv, x, y, z, t1, 0, !isWrite, isWrite)), \
		 BGQ_SPINORSITE_ACCESS(spinorfield, isOdd, tv, x, y, z))

#define BGQ_SPINORSITE_RIGHT(spinorfield, isOdd, tv, x, y, z, t1, t2, isWrite)                         \
		(assert(assert_spinorcoord(spinorfield, isOdd, tv, x, y, z, t2, 1, !isWrite, isWrite)), \
		 BGQ_SPINORSITE_ACCESS(spinorfield, isOdd, tv, x, y, z))

#define BGQ_SPINORVAL(spinorfield,isOdd,t,x,y,z,tv,k, v,c, isRead,isWrite) \
	(assert(assert_spinorfield_coord(spinorfield,isOdd,t,x,y,z,tv,k,v,c,isRead,isWrite)), \
	 &(BGQ_SPINORSITE_ACCESS(spinorfield,isOdd,tv,x,y,z)->s[v][c][k]))

//TODO: move to .c
EXTERN_INLINE _Complex double *bgq_spinorfield_double_ref(bgq_spinorfield_double spinorfield, bool isOdd, int t, int x, int y, int z, int v, int c, bool isRead, bool isWrite) {
	const int teo = t / PHYSICAL_LP;
	const int tv = teo/PHYSICAL_LK;
	const int k = mod(teo,PHYSICAL_LK);

	_Complex double *val = BGQ_SPINORVAL(spinorfield, isOdd, t, x, y, z, tv, k, v, c, isRead, isWrite);
	return val;
}


EXTERN_INLINE _Complex double bgq_spinorfield_double_get(bgq_spinorfield_double spinorfield, bool isOdd, int t, int x, int y, int z, int v, int c) {
	_Complex double *ptr = bgq_spinorfield_double_ref(spinorfield,isOdd,t,x,y,z,v,c,true,false);
	return *ptr;
}


EXTERN_INLINE void bgq_spinorfield_double_set(bgq_spinorfield_double spinorfield, bool isOdd, int t, int x, int y, int z, int v, int c, double _Complex value) {
	_Complex double *ptr = bgq_spinorfield_double_ref(spinorfield,isOdd,t,x,y,z,v,c,true,false);
	*ptr = value;
}


////////////////////////////////////////////////////////////////////////////////
// Gaugefield

#ifndef NDEBUG
//TODO: privatize
EXTERN_FIELD bgq_gaugeeodir_double g_gaugefield_doubledata_debug;
#endif
EXTERN_FIELD bgq_gaugeeodir_double g_gaugefield_doubledata;
EXTERN_FIELD bgq_gaugefield_double g_gaugefield_double;


void bgq_init_gaugefield();
void bgq_free_gaugefield();

void bgq_transfer_gaugefield(bgq_gaugefield_double targetfield, su3 **sourcefield);

void bgq_gaugefield_resetcoord(bgq_gaugefield_double gaugefield, int expected_reads, int expected_writes);


bool assert_gaugeval(bgq_gaugefield_double gaugefield, bool isOdd, int t, int x, int y, int z, int tv, int k, direction dir, int i, int l, bool isRead, bool isWrite);
bool assert_gaugesite(bgq_gaugefield_double gaugefield, bool isOdd, int t, int x, int y, int z, int tv, int k, direction dir, bool isRead, bool isWrite);

#define BGQ_GAUGESITE_ACCESS(gaugefield,isOdd,tv,x,y,z,dir) \
	(assert(0 <= tv && tv <= PHYSICAL_LTV),        \
	 assert(-1 <= x && x < PHYSICAL_LX),        \
	 assert(-1 <= y && y < PHYSICAL_LY),        \
	 assert(-1 <= z && z < PHYSICAL_LZ),        \
	 &gaugefield->eodir[(isOdd)][(dir)/2][(((tv)*(PHYSICAL_LX+1) + ((x)+1))*(PHYSICAL_LY+1) + ((y)+1))*(PHYSICAL_LZ+1) + ((z)+1)])

#define BGQ_GAUGESITE(gaugefield,isOdd,tv,x,y,z,dir,t1,t2,isRead,isWrite)         \
	(assert(assert_gaugesite(gaugefield,isOdd,t1,x,y,z,tv,0,dir,isRead,isWrite)),  \
	 assert(assert_gaugesite(gaugefield,isOdd,t2,x,y,z,tv,1,dir,isRead,isWrite)),  \
	 BGQ_GAUGESITE_ACCESS(gaugefield,isOdd,tv,x,y,z,dir))

#define BGQ_GAUGESITE_LEFT(gaugefield,isOdd,tv,x,y,z,dir,t1,t2,isRead,isWrite)         \
	(assert(assert_gaugesite(gaugefield,isOdd,t1,x,y,z,tv,0,dir,isRead,isWrite)),  \
	 BGQ_GAUGESITE_ACCESS(gaugefield,isOdd,tv,x,y,z,dir))

#define BGQ_GAUGESITE_RIGHT(gaugefield,isOdd,tv,x,y,z,dir,t1,t2,isRead,isWrite)      \
	(assert(assert_gaugesite(gaugefield,isOdd,t2,x,y,z,tv,1,dir,isRead,isWrite)),  \
	 BGQ_GAUGESITE_ACCESS(gaugefield,isOdd,tv,x,y,z,dir))

#define BGQ_GAUGEVAL(gaugefield,isOdd,t,x,y,z,tv,k,dir,i,l, isRead,isWrite) \
	(assert(assert_gaugeval(gaugefield,isOdd,t,x,y,z,tv,k,dir,i,l,isRead,isWrite)), \
	 &(BGQ_GAUGESITE_ACCESS(gaugefield,isOdd,tv,x,y,z,dir)->c[i][l][k]))


//TODO: move to .c
EXTERN_INLINE _Complex double *bgq_gaugefield_double_ref(bgq_gaugefield_double gaugefield, bool isOdd, int t, int x, int y, int z, direction d, int i, int l, bool isRead, bool isWrite) {
	assert(gaugefield);
	assert(false <= isOdd && isOdd <= true);
	assert(isOdd == ((t+x+y+z)&1));
	assert(-1 <= t && t < LOCAL_LT);
	assert(-1 <= x && x < LOCAL_LX);
	assert(-1 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	assert(TUP <= d && d <= ZDOWN);
	assert(0 <= i && i < 3);
	assert(0 <= l && l < 3);
	assert((d&1)==0); // Must be up-direction

	const int teo = (t+1) / PHYSICAL_LP;
	const int tv = teo / PHYSICAL_LK;
	const int k = mod(teo, PHYSICAL_LK);


	assert(assert_gaugeval(gaugefield,isOdd,t,x,y,z,tv,k,d,i,l,isRead,isWrite));
	bgq_gaugesite_double *site = BGQ_GAUGESITE_ACCESS(gaugefield, isOdd, tv,x,y,z,d);
	_Complex double *val = &site->c[i][l][k];

	return val;
}

EXTERN_INLINE _Complex double bgq_gaugefield_double_get(bgq_gaugefield_double gaugefield, bool isOdd, int t, int x, int y, int z, direction d, int i, int l) {
	bool adjoint = false;
	switch (d) {
	case TDOWN:
	case TDOWN_SHIFT:
		t -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = TUP;
		break;
	case XDOWN:
		x -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = XUP;
		break;
	case YDOWN:
		y -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = YUP;
		break;
	case ZDOWN:
		z -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = ZUP;
		break;
	case TUP_SHIFT:
		d = TUP;
		break;
	default:
		break;
	}

	if (adjoint) {
		const int tmp = i;
		i = l;
		l = tmp;
	}

	_Complex double result = *bgq_gaugefield_double_ref(gaugefield, isOdd, t,x,y,z,d,i,l, true, false);
	if (adjoint) {
		result = conj(result);
	}
	return result;
}



EXTERN_INLINE void bgq_gaugefield_double_set(bgq_gaugefield_double gaugefield, bool isOdd, int t, int x, int y, int z, direction d, int i, int l, double _Complex value) {
	assert(gaugefield);
	assert(false <= isOdd && isOdd <= true);
	assert(isOdd == ((t+x+y+z)&1));
	assert(-1 <= t && t < LOCAL_LT);
	assert(-1 <= x && x < LOCAL_LX);
	assert(-1 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	assert(TUP <= d && d <= ZDOWN);
	assert(0 <= i && i < 3);
	assert(0 <= l && l < 3);

	bool adjoint = false;
	switch (d) {
	case TDOWN:
	case TDOWN_SHIFT:
		t -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = TUP;
		break;
	case XDOWN:
		x -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = XUP;
		break;
	case YDOWN:
		y -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = YUP;
		break;
	case ZDOWN:
		z -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = ZUP;
		break;
	case TUP_SHIFT:
		d = TUP;
		break;
	default:
		break;
	}

	const int teo = (t+1) / PHYSICAL_LP;
	const int tv = teo / PHYSICAL_LK;
	const int k = mod(teo, PHYSICAL_LK);

	if (adjoint) {
		int tmp = i;
		i = l;
		l = tmp;
		crealf(value);
		value = conj(value);
	}

	{
	_Complex double *valptr = BGQ_GAUGEVAL(gaugefield, isOdd, t,x,y,z, tv,k, d, i,l, false,true);
	*valptr = value;
	}

	if (d == ZUP) {
		// For z-direction...
		if (z==LOCAL_LZ-1) {
			// wraparound, also store as z=-1 coordinate
			_Complex double *wrapptr = BGQ_GAUGEVAL(gaugefield, isOdd, t,x,y,-1, tv,k, d, i,l, false,true);
			*wrapptr = value;
		}
	}

	if (d == TUP) {
			// For t-direction...
			// also store as shifted

			// Move one to the right
		    const int teo_shift = mod(teo+1, 1+LOCAL_LT/PHYSICAL_LP);
			const int tv_shift = teo_shift / PHYSICAL_LK;
			const int k_shift = mod(teo_shift, PHYSICAL_LK);

			if (t == LOCAL_LT-1) {
				assert(tv_shift == 0);
				assert(k_shift == 0);
			} else {	assert((tv_shift==tv && k_shift==1) || (tv_shift==tv+1 && k_shift==0));
			}
			_Complex double *shiftptr = BGQ_GAUGEVAL(gaugefield, isOdd, t,x,y,z, tv_shift,k_shift, TUP_SHIFT, i,l, false,true);
			*shiftptr = value;
		}
}


////////////////////////////////////////////////////////////////////////////////
// Weylfields


bool assert_weylfield_t(bgq_weylfield_double weylfield, bool isOdd, int t, int x, int y, int z, int xv, int k);
#define BGQ_WEYLSITE_T(weylfield, isOdd, t, xv, y, z, x1, x2)  \
	(assert(assert_weylfield_t(weylfield,isOdd,t,x1,y,z,xv,0)), \
	 assert(assert_weylfield_t(weylfield,isOdd,t,x2,y,z,xv,1)), \
	 &weylfield[((xv)*PHYSICAL_LY + (y))*PHYSICAL_LZ + (z)])

bool assert_weylfield_x(bgq_weylfield_double weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k);
#define BGQ_WEYLSITE_X(weylfield, isOdd, tv, x, y, z, t1, t2)  \
	(assert(assert_weylfield_t(weylfield,isOdd,t1,x,y,z,tv,0)), \
	 assert(assert_weylfield_t(weylfield,isOdd,t2,x,y,z,tv,1)), \
	 &weylfield[((tv)*PHYSICAL_LY + (y))*PHYSICAL_LZ + (z)])

bool assert_weylfield_y(bgq_weylfield_double weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k);
#define BGQ_WEYLSITE_Y(weylfield, isOdd, tv, x, y, z, t1, t2)  \
	(assert(assert_weylfield_t(weylfield,isOdd,t1,x,y,z,tv,0)), \
	 assert(assert_weylfield_t(weylfield,isOdd,t2,x,y,z,tv,1)), \
	 &weylfield[((tv)*PHYSICAL_LX + (x))*PHYSICAL_LZ + (z)])








#undef EXTERN_INLINE
#undef EXTERN_FIELD
#endif // BGQ_H_INCLUDED
