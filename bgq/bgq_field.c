

#include "bgq.h"
#include "../global.h"

#define EXTERN
#include "bgq_field.h"


void bgq_transfer_gaugefield(bgq_gaugefield_double targetfield, spinor *sourcefield) {

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

		for (int d = X_UP; d <= T_UP; d+=2) {
			int kb = g_idn[g_eo2lexic[i]][d/2];
			su3 *m = &g_gauge_field[kb][d/2];
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
}





