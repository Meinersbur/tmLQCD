#ifndef BGQ_H_INCLUDED
#define BGQ_H_INCLUDED

#include "../global.h"
#include <stdbool.h>
#include <complex.h>
#include <assert.h>

#ifndef EXTERN
#define EXTERN extern
#endif

#define VECTOR_WIDTH 2 /* 2 complex values = 4 reals */

#define LOCAL_LP 2 /* Even/Odd */
#define LOCAL_LX LX
#define LOCAL_LY LY
#define LOCAL_LZ LZ
#define LOCAL_LT T
#define LOCAL_LK 2 /* Vector unit width (2 complex = 4 reals) */

#define PHYSICAL_LP 2 /* Even/Odd */
#define PHYSICAL_LX LOCAL_LX
#define PHYSICAL_LY LOCAL_LY
#define PHYSICAL_LZ LOCAL_LZ
#define PHYSICAL_LTV (LOCAL_LT/(PHYSICAL_LP*PHYSICAL_LK))
#define PHYSICAL_LK 2 /* Vector unit width (2 complex = 4 reals) */
#define PHYSICAL_LD 5 /* No of directions (for gauge field) */


// 2 of each
#define COUNT_FACE_XY ((PHYSICAL_LX-2)*(PHYSICAL_LY-2)*(1))
#define COUNT_FACE_XZ ((PHYSICAL_LX-2)*(1)*(PHYSICAL_LZ-2))
#define COUNT_FACE_YZ ((1)*(PHYSICAL_LY-2)*(PHYSICAL_LZ-2))
#define COUNT_FACES (2*(COUNT_FACE_XY + COUNT_FACE_XZ + COUNT_FACE_YZ))

// 4 of each
#define COUNT_EDGE_X ((PHYSICAL_LX-2)*(1)*(1))
#define COUNT_EDGE_Y ((1)*(PHYSICAL_LY-2)*(1))
#define COUNT_EDGE_Z ((1)*(1)*(PHYSICAL_LZ-2))
#define COUNT_EDGES  (4*(COUNT_EDGE_X + COUNT_EDGE_Y + COUNT_EDGE_Z))

// 8
#define COUNT_VERTEX ((1)*(1)*(1))
#define COUNT_VERTICES (8*COUNT_VERTEX)

#define TOTAL_BORDER (COUNT_FACES + COUNT_EDGES + COUNT_VERTICES)



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
	bgq_gaugesite_double (*eodir)[PHYSICAL_LP][PHYSICAL_LD];
} bgq_gaugeeodir_double;
typedef bgq_gaugeeodir_double (*bgq_gaugefield_double);


typedef enum {
	DIR_DOWN = 0,
	DIR_UP = 1,

	X_DOWN = 0,
	X_UP = 1,
	Y_DOWN = 2,
	Y_UP = 3,
	Z_DOWN = 4,
	Z_UP = 5,
	T_DOWN = 6,
	T_UP = 7,

	T_RAGGED_DOWN = 8,
	T_RAGGED_UP = 9
} direction;


#define BGQ_SPINORSITE(spinorfield, isOdd, x, y, z, tv) \
		(&spinorfield[((x*PHYSICAL_LY + y)*PHYSICAL_LZ + z)*PHYSICAL_LTV + tv]);

EXTERN inline double _Complex *bgq_spinorfield_double_local_to_physical(bgq_spinorfield_double spinorfield, bool isOdd, int x, int y, int z, int t, int s, int c) {
	assert(spinorfield);
	assert(0 <= isOdd && isOdd < LOCAL_LP);
	assert(isOdd == (x+y+z+t)%2);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= s && s < 4);
	assert(0 <= c && c < 3);

	int teo = t / PHYSICAL_LP;
	int tv = teo / PHYSICAL_LK;
	int k = teo % PHYSICAL_LK;
	bgq_spinorsite_double *site = BGQ_SPINORSITE(spinorfield, isOdd, x, y, z, tv);
	return &site->s[s][c][k];
}

#define BGQ_WEYLSITE_X(weylfield, isOdd, x, y, z, tv) \
	(&weylfield[(y*PHYSICAL_LZ + z)*PHYSICAL_LTV + tv]);

#define BGQ_WEYLSITE_Y(weylfield, isOdd, x, y, z, tv) \
	(&weylfield[(x*PHYSICAL_LZ + z)*PHYSICAL_LTV + tv]);

#define BGQ_WEYLSITE_Z(weylfield, isOdd, x, y, z, tv) \
	(&weylfield[(x*PHYSICAL_LY + y)*PHYSICAL_LTV + tv]);




#define BGQ_GAUGESITE(gaugefield,isOdd,x,y,z,tv,direction) \
		(&gaugefield->eodir[isOdd][direction][((x*PHYSICAL_LY + y)*PHYSICAL_LZ + z)*PHYSICAL_LTV + tv]);

EXTERN inline double _Complex *bgq_gaugefield_double_local_to_physical(bgq_gaugefield_double gaugefield, bool isOdd, int x, int y, int z, int t, direction d, int i, int l) {
	assert(gaugefield);
	assert(false <= isOdd && isOdd <= true);
	assert(isOdd == (x+y+z+t)%2);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	assert(0 <= t && t < LOCAL_LT);
	assert(X_DOWN <= d && d <= T_UP);
	assert(0 <= i && i < 3);
	assert(0 <= l && l < 3);

	int teo = t / LOCAL_LP;
	int tv = teo / LOCAL_LK;
	int k = teo % LOCAL_LK;

	switch (d) {
	case X_UP:
		x -= 1;
		break;
	case Y_UP:
		y -= 1;
		break;
	case Z_UP:
		z -= 1;
		break;
	case T_UP:
		t -= 1;
		break;
	default:
		break;
	}

	bgq_gaugesite_double *site = BGQ_GAUGESITE(gaugefield,isOdd,x,y,z,tv,d/2);
	return &site->c[i][l][k];
}


#undef EXTERN
#endif // BGQ_H_INCLUDED


