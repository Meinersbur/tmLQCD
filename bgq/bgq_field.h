#ifndef BGQ_H_INCLUDED
#define BGQ_H_INCLUDED

#include "../global.h"
#include <stdbool.h>
#include <complex.h>

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
#define PHYSICAL_LD 4 /* No of directions (for gauge field) */

typedef struct {
	double _Complex s[4][3][PHYSICAL_LK];
} bgq_spinorsite_double;

typedef struct {
	double _Complex c[3][3][PHYSICAL_LK];
} bgq_gaugesite_double;

typedef struct {
	double _Complex s[2][3][PHYSICAL_LK];
} bgq_weylsite_double;

typedef bgq_spinorsite_double (*bgq_spinorfield_double);
typedef bgq_weylsite_double (*bgq_weylfield_double);
typedef bgq_gaugesite_double (*bgq_gaugefield_double)[PHYSICAL_LD];
typedef  {
	bgq_gaugefield_double eo[PHYSICAL_LP];
} bgq_gaugefieldeo_double;






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
} direction;


EXTERN inline bgq_spinorsite_double *bgq_spinorsite_double_physical_pointer(bgq_spinorfield_double spinorfield, bool isOdd, int x, int y, int z, int tv) {
	assert(spinorfield);
	assert(0 <= isOdd && isOdd < PHYSICAL_LP);
	assert(0 <= x && x < PHYSICAL_LX);
	assert(0 <= y && y < PHYSICAL_LY);
	assert(0 <= z && z < PHYSICAL_LZ);
	assert(0 <= tv && tv < PHYSICAL_LTV);

	return &spinorfield[((x*LOCAL_LY + y)*LOCAL_LZ + z)*PHYSICAL_LTV + tv];
}


EXTERN inline complex *bgq_spinorfield_double_local_to_physical(bgq_spinorfield_double spinorfield, bool isOdd, int x, int y, int z, int t, int s, int c) {
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
	bgq_spinorfield_double site = bgq_spinorsite_double_physical_pointer(spinorfield, isOdd, x, y, z, tv);
	return &(*site)[s][c][k];
}

EXTERN inline bgq_gaugesite_double *bgq_gaugesiteeo_double_physical_pointer(bgq_gaugefieldeo_double gaugefield, bool isOdd, int x, int y, int z, int tv, direction d) {
	assert(gaugefield);
	assert(0 <= isOdd && isOdd <= PHYSICAL_LP);
	assert(0 <= x && x < PHYSICAL_LX);
	assert(0 <= y && y < PHYSICAL_LY);
	assert(0 <= z && z < PHYSICAL_LZ);
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(X_DOWN <= direction && direction <= T_UP);
	assert(direction % 2 == DIR_UP && "Only up-directions are stored");

	return &gaugefield[isOdd][((x*LOCAL_LY + y)*LOCAL_LZ + z)*PHYSICAL_LTV + tv][direction/2];
}


EXTERN inline bgq_gaugesite_double *bgq_gaugesite_double_physical_pointer(bgq_gaugefield_double gaugefield, bool isOdd, int x, int y, int z, int tv, direction d) {
	assert(gaugefield);
	assert(0 <= isOdd && isOdd <= PHYSICAL_LP);
	assert(0 <= x && x < PHYSICAL_LX);
	assert(0 <= y && y < PHYSICAL_LY);
	assert(0 <= z && z < PHYSICAL_LZ);
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(X_DOWN <= direction && direction <= T_UP);
	assert(direction % 2 == DIR_UP && "Only up-directions are stored");

	return &gaugefield[((x*LOCAL_LY + y)*LOCAL_LZ + z)*PHYSICAL_LTV + tv][direction/2];
}

EXTERN inline complex bgq_gaugefield_double_local_to_physical(bgq_gaugefield_double gaugefield, bool isOdd, int x, int y, int z, int t, direction d, int i, int l) {
	assert(gaugefield);
	assert(false <= isOdd && isOdd <= true);
	assert(isOdd == (x+y+z+t)%2);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	assert(0 <= t && t < LOCAL_LT);
	assert(X_DOWN <= direction && direction <= T_UP);
	assert(0 <= i && i < 3);
	assert(0 <= l && l < 3);

	int teo = t / LOCAL_LP;
	int tv = t / LOCAL_LK;
	int k = t % LOCAL_LK;

	switch (direction) {
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
	}

	bgq_gaugesite_double *site = bgq_gaugesite_double_physical_pointer(gaugefield, isOdd, x, y, z, tv, (direction/2)*2);
	return &(*site)[i][l][k];
}


#undef EXTERN
#endif // BGQ_H_INCLUDED


