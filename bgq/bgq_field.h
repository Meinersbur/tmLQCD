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
#define PHYSICAL_LTV (LOCAL_LZ/(PHYSICAL_LP*PHYSICAL_LK))
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
#define VOLUME_SITES (VOLUME)

#define SURFACE_ZLINES (SURFACE_FACE_ZLINES+SURFACE_FACE_ZLINES+SURFACE_EDGE_ZLINES)
#define BODY_ZLINES ((PHYSICAL_LTV-2)*(PHYSICAL_LX-2)*(PHYSICAL_LY-2))


#define GAUGE_VOLUME ((LOCAL_LT+1)*(LOCAL_LX+1)*(LOCAL_LY+1)*(LOCAL_LZ))


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
	T_UP = 0,
	T_DOWN = 1,
	X_UP = 2,
	X_DOWN = 3,
	Y_UP = 4,
	Y_DOWN = 5,
	Z_UP = 6,
	Z_DOWN = 7,

	T_UP_SHIFT = 8,
	T_DOWN_SHIFT = 9,

	DIR_UP = 0,
	DIR_DOWN = 1
} direction;


void *malloc_aligned(size_t size, size_t alignment);




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

EXTERN_FIELD bgq_gaugeeodir_double g_gaugefield_doubledata;
EXTERN_FIELD bgq_gaugefield_double g_gaugefield_double;

void bgq_init_gaugefield();
void bgq_free_gaugefield();

void bgq_transfer_gaugefield(bgq_gaugefield_double targetfield, su3 **sourcefield);













bool assert_spinorcoord(bgq_spinorfield_double spinorfield, bool isOdd, int t, int x, int y, int z, int tv, int k);
#define BGQ_SPINORSITE(spinorfield, isOdd, tv, x, y, z, t1, t2)                        \
		(assert(assert_spinorcoord(spinorfield, isOdd, tv, x, y, z, t1, 0)),			\
		 assert(assert_spinorcoord(spinorfield, isOdd, tv, x, y, z, t2, 1)),		  	\
		 &spinorfield[(((tv)*PHYSICAL_LX + (x))*PHYSICAL_LY + (y))*PHYSICAL_LZ + (z)])

EXTERN_INLINE void bgq_spinorfield_double_set(bgq_spinorfield_double spinorfield, bool isOdd, int t, int x, int y, int z, int v, int c, double _Complex value) {
	const int teo = t / PHYSICAL_LP;
	const int tv = teo/PHYSICAL_LK;
	const int k = mod(teo,PHYSICAL_LK);

	bgq_spinorsite_double *site = BGQ_SPINORSITE(spinorfield, isOdd, t, x, y, z, tv, k);
	site->s[v][c][k] = value;
}


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



bool assert_gaugesite(bgq_gaugefield_double gaugefield, bool isOdd, int t, int x, int y, int z, int tv, int k, direction dir);
#define BGQ_GAUGESITE(gaugefield,isOdd,tv,x,y,z,dir,t1,t2)         \
	(assert(assert_gaugesite(gaugefield,isOdd,t1,x,y,z,tv,0,dir)),  \
	 assert(assert_gaugesite(gaugefield,isOdd,t2,x,y,z,tv,1,dir)),  \
	 &gaugefield->eodir[(isOdd)][(dir)/2][(((tv)*PHYSICAL_LX + ((x)+1))*PHYSICAL_LY + ((y)+1))*PHYSICAL_LZ + (z)])

EXTERN_INLINE void bgq_gaugefield_double_set(bgq_gaugefield_double gaugefield, bool isOdd, int t, int x, int y, int z, direction d, int i, int l, double _Complex value) {
	assert(gaugefield);
	assert(false <= isOdd && isOdd <= true);
	assert(isOdd == (8+t+x+y+z)%2);
	assert(-1 <= t && t < LOCAL_LT);
	assert(-1 <= x && x < LOCAL_LX);
	assert(-1 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	assert(T_UP <= d && d <= Z_DOWN);
	assert(0 <= i && i < 3);
	assert(0 <= l && l < 3);

	bool adjoint = false;
	switch (d) {
	case T_DOWN:
	case T_DOWN_SHIFT:
		t -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = T_UP;
		break;
	case X_DOWN:
		x -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = X_UP;
		break;
	case Y_DOWN:
		y -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = Y_UP;
		break;
	case Z_DOWN:
		z -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = Z_UP;
		break;
	case T_UP_SHIFT:
		d = T_UP;
		break;
	default:
		break;
	}

	const int teo = t / PHYSICAL_LP;
	const int tv = (PHYSICAL_LK + teo) / PHYSICAL_LK;
	const int k = mod(PHYSICAL_LK + teo, PHYSICAL_LK);

	if (adjoint) {
		int tmp = i;
		i = l;
		l = tmp;
		crealf(value);
		value = conj(value);
	}

	bgq_gaugesite_double *site = BGQ_GAUGESITE(gaugefield, isOdd, tv,x,y,z,d, t-k,t+(1-k));
	site->c[i][l][k] = value;

	if (d == T_UP) {
		// Do some special things for Z_UP, values exist in multiple places

		// Move one to the right
		const int tv_shift = (teo+1) / PHYSICAL_LK;
		const int k_shift = mod(teo+1, PHYSICAL_LK);

		bgq_gaugesite_double *raggedsite = BGQ_GAUGESITE(gaugefield, isOdd, tv_shift,x,y,z,T_UP_SHIFT, t-k_shift,t+(1-k_shift));
		raggedsite->c[i][l][k_shift] = value;

		if (t == LOCAL_LT-1) {
			// wraparound
			bgq_gaugesite_double *wraparoundsite = BGQ_GAUGESITE(gaugefield, isOdd, 0,x,y,z,T_UP_SHIFT, 0,1);
			wraparoundsite->c[i][l][0] = value;
		}

		if (t == 0) {
			// wraparound
			bgq_gaugesite_double *wraparoundsite = BGQ_GAUGESITE(gaugefield, isOdd, PHYSICAL_LTV,x,y,z,T_UP_SHIFT, LOCAL_LT-2,LOCAL_LT-1);
			wraparoundsite->c[i][l][1] = value;
		}
	}
}






#undef EXTERN_INLINE
#undef EXTERN_FIELD
#endif // BGQ_H_INCLUDED
