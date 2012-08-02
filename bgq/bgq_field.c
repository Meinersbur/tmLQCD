#define BGQ_FIELD_C_H_

#include "bgq_field.h"

#include "bgq.h"
#include "../geometry_eo.h"
#include "../global.h"

#include <string.h>

void *malloc_aligned(size_t size, size_t alignment) {
	void *result = NULL;
	int errcode = posix_memalign(&result, alignment, size);
	if (errcode != 0) {
		fprintf(stderr, "malloc returned %d\n", errcode);
		exit(10);
	}
	memset(result, 0, size);
	return result;
}

void bgq_init_gaugefield() {
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (direction dir = T_UP; dir <= Z_UP_SHIFT; dir += 2) {
			int zlinelength = PHYSICAL_LZV;
			if (dir == Z_UP_SHIFT)
				zlinelength += 1;

			g_gaugefield_doubledata.eodir[isOdd][dir / 2] = malloc_aligned((PHYSICAL_LT + 1) * (PHYSICAL_LX + 1) * (PHYSICAL_LY + 1) * zlinelength * sizeof(bgq_gaugesite_double), 128 /*L2 cache line size*/);
		}
	}

	g_gaugefield_double = &g_gaugefield_doubledata;
}

void bgq_free_gaugefield() {
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (int dir = T_UP; dir <= Z_UP_SHIFT; dir += 2) {
			free(g_gaugefield_doubledata.eodir[isOdd][dir / 2]);
			g_gaugefield_doubledata.eodir[isOdd][dir / 2] = NULL;
		}
	}
	g_gaugefield_double = NULL;
}
typedef struct {
	double _Complex c[3][3];
} su3_array64;

void bgq_transfer_gaugefield(bgq_gaugefield_double targetfield, su3 **sourcefield) {

#pragma omp parallel for schedule(static)
	for (int tzy = 0; tzy < GAUGE_VOLUME; tzy += 1) {
		WORKLOAD_DECL(tzy, GAUGE_VOLUME);
		const int z = WORKLOAD_PARAM(LOCAL_LZ);
		const int y = WORKLOAD_PARAM(LOCAL_LY+1) - 1;
		const int x = WORKLOAD_PARAM(LOCAL_LX+1) - 1;
		const int t = WORKLOAD_PARAM(LOCAL_LT+1) - 1;
		WORKLOAD_CHECK
		;

		const bool isOdd = (8 + t + x + y + z) % 2;
		const int ix = Index(t, x, y, z);

		for (direction d = T_UP; d <= Z_UP; d += 2) {
			su3_array64 *m = (su3_array64*) &g_gauge_field[ix][d / 2];
			for (int i = 0; i < 3; i += 1) {
				for (int l = 0; l < 3; l += 1) {
					//complex *c = ((complex*)m)+i*3+l;
					bgq_gaugefield_double_set(targetfield, isOdd, t, x, y, z, d, i, l, m->c[i][l]);
				}
			}
		}
	}
}
