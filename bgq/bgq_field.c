
#define BGQ_FIELD_C_H_

#include "bgq_field.h"


#include "bgq.h"
#include "../global.h"
#include <string.h>



void *malloc_aligned(size_t size, size_t alignment) {
	void *result=NULL;
	int errcode = posix_memalign(&result, alignment, size) ;
	if (errcode != 0) {
		fprintf(stderr ,"malloc returned %d\n", errcode);
		exit(10);
	}
	memset(result, 0, size);
	return result;
}


void bgq_init_gaugefield() {
	for (bool isOdd = false; isOdd <= true; isOdd+=1) {
		for (int dir = T_UP; dir <= Z_UP_SHIFT; dir+=2) {
			int zlinelength = PHYSICAL_LZV;
			if (dir == Z_UP_SHIFT)
				zlinelength +=1;

			g_gaugefield_doubledata.eodir[isOdd][dir/2] = malloc_aligned((PHYSICAL_LT+1)*(PHYSICAL_LX+1)*(PHYSICAL_LY+1) * zlinelength * sizeof(bgq_gaugesite_double), 128 /*L2 cache line size*/);
		}
	}

	g_gaugefield_double = &g_gaugefield_doubledata;
}

void bgq_free_gaugefield() {
	for (bool isOdd = false; isOdd <= true; isOdd+=1) {
			for (int dir = T_UP; dir <= Z_UP_SHIFT; dir+=2) {
				free(g_gaugefield_doubledata.eodir[isOdd][dir/2]);
				g_gaugefield_doubledata.eodir[isOdd][dir/2] = NULL;
			}
	}
	g_gaugefield_double = NULL;
}


void bgq_transfer_gaugefield(bgq_gaugefield_double targetfield, su3 **sourcefield) {

#pragma omp parallel for schedule(static)
	for (int xyz = 0; xyz < GAUGE_VOLUME; xyz+=1) {
		WORKLOAD_DECL(GAUGE_VOLUME);
		const int z = WORKLOAD_PARAM(LOCAL_LZ+1)-1;
		const int y = WORKLOAD_PARAM(LOCAL_LY+1)-1;
		const int x = WORKLOAD_PARAM(LOCAL_LX+1)-1;
		const int t = WORKLOAD_PARAM(LOCAL_LT+1)-1;
		const bool isOdd = (t+x+y+z)%2;

		const int ix = g_ipt[t][x][y][z];

		for (int d = X_UP; d <= T_UP; d+=1) {
			su3 *m = &g_gauge_field[ix][d/2];
			for (int i = 0; i < 3; i+=1) {
				for (int l = 0; l < 3; l+=1) {
					complex *c = ((complex*)m)+i*3+l;
					bgq_gaugefield_double_set(targetfield, isOdd, x, y, z, t, d, i,l, cs2c99(*c));
				}
			}
		}
	}

#if 0
#pragma omp parallel for schedule(static)
	for (int xyz = 0; xyz < VOLUME; xyz+=1) {
		WORKLOAD_DECL(VOLUME);

		int i;
		bool isOdd;
		if (WORKLOAD_SPLIT(VOLUME/2)) {
			// Even sites
			i = WORKLOAD_PARAM(VOLUME/2);
			isOdd = false;
		} else {
			// Odd sites
			assert (xyz_total == VOLUME/2);
			i = WORKLOAD_PARAM(VOLUME/2)+(VOLUME+RAND)/2;
			isOdd = true;
		}
		assert(xyz == 0);
		assert(xyz_total == 1);

		const int x = g_x[i];
		const int y = g_y[i];
		const int z = g_y[i];
		const int t = g_t[i];

		for (int d = T_UP; d <= Z_UP; d+=2) {
			int kb = g_idn[g_eo2lexic[i]][d/2];
			su3 *m = &sourcefield[kb][d/2];
			bgq_gaugefield_double_set(targetfield, isOdd, x, y, z, t, d, 0,0, cs2c99(m->c00));
			bgq_gaugefield_double_set(targetfield, isOdd, x, y, z, t, d, 0,1, cs2c99(m->c01));
			bgq_gaugefield_double_set(targetfield, isOdd, x, y, z, t, d, 0,2, cs2c99(m->c02));
			bgq_gaugefield_double_set(targetfield, isOdd, x, y, z, t, d, 1,0, cs2c99(m->c10));
			bgq_gaugefield_double_set(targetfield, isOdd, x, y, z, t, d, 1,1, cs2c99(m->c11));
			bgq_gaugefield_double_set(targetfield, isOdd, x, y, z, t, d, 1,2, cs2c99(m->c12));
			bgq_gaugefield_double_set(targetfield, isOdd, x, y, z, t, d, 2,0, cs2c99(m->c20));
			bgq_gaugefield_double_set(targetfield, isOdd, x, y, z, t, d, 2,1, cs2c99(m->c21));
			bgq_gaugefield_double_set(targetfield, isOdd, x, y, z, t, d, 2,2, cs2c99(m->c22));
		}
	}
#endif
}





