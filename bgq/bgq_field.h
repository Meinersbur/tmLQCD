#ifndef BGQ_H_INCLUDED
#define BGQ_H_INCLUDED

#include "../global.h"
#include <stdbool.h>
#include <assert.h>
#include "complex_c99.h"

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
#define PHYSICAL_LT LOCAL_LT
#define PHYSICAL_LX LOCAL_LX
#define PHYSICAL_LY LOCAL_LY
#define PHYSICAL_LZV (LOCAL_LZ/(PHYSICAL_LP*PHYSICAL_LK))
#define PHYSICAL_LK 2 /* Vector unit width (2 complex = 4 reals) */
#define PHYSICAL_LD 5 /* No of directions (for gauge field) = 4 up-directions + 1 for shifted z-line (note that it is one vector-element longer for wraparound values) */

// volume without borders
#define BODY_ZLINES ((PHYSICAL_LT-2)*(PHYSICAL_LX-2)*(PHYSICAL_LY-2))

// faces, 2 of each
#define BORDER_ZLINES_T ((1)*(PHYSICAL_LX-2)*(PHYSICAL_LY-2))
#define BORDER_ZLINES_X ((PHYSICAL_LT-2)*(1)*(PHYSICAL_LY-2))
#define BORDER_ZLINES_Y ((PHYSICAL_LT-2)*(PHYSICAL_LX-2)*(1))
#define BORDER_ZLINES_ALLFACES (2*BORDER_ZLINES_T+2*BORDER_ZLINES_X+2*BORDER_ZLINES_Y)

// edges, 4 of each
#define BORDER_ZLINES_TX ((1)*(1)*(PHYSICAL_LY-2))
#define BORDER_ZLINES_TY ((1)*(PHYSICAL_LX-2)*(1))
#define BORDER_ZLINES_XY ((PHYSICAL_LT-2)*(1)*(1))
#define BORDER_ZLINES_ALLEDGES (4*BORDER_ZLINES_TX+4*BORDER_ZLINES_TY+4*BORDER_ZLINES_XY)

// vertices, 8
#define BORDER_ZLINES_TXY ((1)*(1)*(1))
#define BORDER_ZLINES_ALLVERTICES (8*BORDER_ZLINES_TXY)

#define BORDER_ZLINES_TOTAL (BORDER_ZLINES_ALLFACES+BORDER_ZLINES_ALLEDGES+BORDER_ZLINES_ALLVERTICES)
#define VOLUME_ZLINES ((PHYSICAL_LT)*(PHYSICAL_LX)*(PHYSICAL_LY))

#define SURFACE_ZLINES_T ((1)*(PHYSICAL_LX)*(PHYSICAL_LY))
#define SURFACE_ZLINES_X ((PHYSICAL_LT)*(1)*(PHYSICAL_LY))
#define SURFACE_ZLINES_Y ((PHYSICAL_LT)*(PHYSICAL_LX)*(1))
#define SURFACE_ZLINES_TOTAL (2*(SURFACE_ZLINES_T+SURFACE_ZLINES_X+SURFACE_ZLINES_Y))

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

	Z_UP_SHIFT = 8,
	Z_DOWN_SHIFT = 9,

	DIR_UP = 0,
	DIR_DOWN = 1
} direction;

void *malloc_aligned(size_t size, size_t alignment);

#define BGQ_SPINORSITE(spinorfield, isOdd, t, x, y, zv) \
		(&spinorfield[((t*PHYSICAL_LX + x)*PHYSICAL_LY + y)*PHYSICAL_LZV + zv]);

EXTERN_INLINE double _Complex *bgq_spinorfield_double_local_to_physical(bgq_spinorfield_double spinorfield, bool isOdd, int t, int x, int y, int z, int s, int c) {
	assert(spinorfield);
	assert(false <= isOdd && isOdd <= true);
	assert(isOdd == (x+y+z+t)%2);
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	assert(0 <= s && s < 4);
	assert(0 <= c && c < 3);

	int zeo = z / PHYSICAL_LP;
	int zv = zeo / PHYSICAL_LK;
	int k = zeo % PHYSICAL_LK;
	bgq_spinorsite_double *site = BGQ_SPINORSITE(spinorfield, isOdd, t, x, y, zv);
	return &site->s[s][c][k];
}


#define BGQ_WEYLSITE_T(weylfield, isOdd, t, x, y, zv) \
	(&weylfield[(x*PHYSICAL_LY + y)*PHYSICAL_LZV + zv]);

#define BGQ_WEYLSITE_X(weylfield, isOdd, t, x, y, zv) \
	(&weylfield[(t*PHYSICAL_LY + y)*PHYSICAL_LZV + zv]);

#define BGQ_WEYLSITE_Y(weylfield, isOdd, t, x, y, zv) \
	(&weylfield[(t*PHYSICAL_LX + x)*PHYSICAL_LZV + zv]);

#define BGQ_GAUGESITE(gaugefield,isOdd,t,x,y,zv,direction) \
		(&gaugefield->eodir[(isOdd)][(direction)/2][((((t)+1)*PHYSICAL_LX + ((x)+1))*PHYSICAL_LY + ((y)+1))*PHYSICAL_LZV + (zv)]);

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
		case Z_DOWN_SHIFT:
		z -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = Z_UP;
		break;
	case Z_UP_SHIFT:
		d = Z_UP;
		break;
	default:
		break;
	}

	int zeo = z / PHYSICAL_LP;
	int zv = (PHYSICAL_LK + zeo) / PHYSICAL_LK;
	int k = (PHYSICAL_LK + zeo) % PHYSICAL_LK;

	if (adjoint) {
		int tmp = i;
		i = l;
		l = tmp;
		crealf(value);
		value = conj(value);
	}

	bgq_gaugesite_double *site = BGQ_GAUGESITE(gaugefield, isOdd, t,x,y,zv,d);
	site->c[i][l][k] = value;

	if (d == Z_UP) {
		// Do some special things for Z_UP, values exist in multiple places

		// Move one to the right
		int zv_shift = (zeo + 1) / PHYSICAL_LK;
		int k_shift = (zeo + 1) % PHYSICAL_LK;

		bgq_gaugesite_double *raggedsite = BGQ_GAUGESITE(gaugefield, isOdd, t,x,y,zv_shift, Z_UP_SHIFT);
		raggedsite->c[i][l][k_shift] = value;

		if (z == LOCAL_LZ - 1) {
			// wraparound
			bgq_gaugesite_double *site = BGQ_GAUGESITE(gaugefield, isOdd, t,x,y,0, Z_UP_SHIFT);
			site->c[i][l][0] = value;
		}

		if (z == 0) {
			// wraparound
			bgq_gaugesite_double *site = BGQ_GAUGESITE(gaugefield, isOdd, t,x,y,PHYSICAL_LZV, Z_UP_SHIFT);
			site->c[i][l][1] = value;
		}
	}
}

EXTERN_FIELD bgq_spinorsite_double *g_spinorfields_doubledata;
EXTERN_FIELD bgq_spinorfield_double *g_spinorfields_double;

void bgq_init_spinorfields(int count);
void bgq_free_spinofields();

void bgq_transfer_spinorfield(bool isOdd, bgq_spinorfield_double targetfield, spinor *sourcefield);

EXTERN_FIELD bgq_gaugeeodir_double g_gaugefield_doubledata;
EXTERN_FIELD bgq_gaugefield_double g_gaugefield_double;

void bgq_init_gaugefield();
void bgq_free_gaugefield();

void bgq_transfer_gaugefield(bgq_gaugefield_double targetfield, su3 **sourcefield);

#undef EXTERN_INLINE
#undef EXTERN_FIELD
#endif // BGQ_H_INCLUDED
