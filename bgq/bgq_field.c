

#include "bgq.h"
#include "../global.h"

#define EXTERN
#include "bgq_field.h"


void bgq_transfer_gaugefield(bgq_gaugefield_double targetfield, spinor *sourcefield) {

#pragma parallel for schedule(static)
	for (int xyz = 0; xyz < VOLUME; xyz+=1) {
		WORKLOAD_DECL(VOLUME);

		int i;
		if (WORKLOAD_SPLIT(VOLUME/2)) {
			// Even sites
			i = WORKLOAD_PARAM(VOLUME/2);
		} else {
			// Odd sites
			assert (xyz_total == VOLUME/2);
			i = WORKLOAD_PARAM(VOLUME/2)+(VOLUME+RAND)/2;
		}
		assert(xyz == 0);
		assert(xyz_total == 1);

		for (int d = X_UP; d <= T_UP; d+=2) {
			int kb = g_idn[g_eo2lexic[i]][d/2];
			su3 *m = &g_gauge_field[kb][d/2];
			bgq_gaugefield_double_set(targetfield, true, g_x[i], g_y[i], g_z[i], g_t[i], d, 0,0,m->c00);
			bgq_gaugefield_double_set(targetfield, true, g_x[i], g_y[i], g_z[i], g_t[i], d, 0,1,m->c01);
			bgq_gaugefield_double_set(targetfield, true, g_x[i], g_y[i], g_z[i], g_t[i], d, 0,2,m->c02);
			bgq_gaugefield_double_set(targetfield, true, g_x[i], g_y[i], g_z[i], g_t[i], d, 1,0,m->c10);
			bgq_gaugefield_double_set(targetfield, true, g_x[i], g_y[i], g_z[i], g_t[i], d, 1,1,m->c11);
			bgq_gaugefield_double_set(targetfield, true, g_x[i], g_y[i], g_z[i], g_t[i], d, 1,2,m->c12);
			bgq_gaugefield_double_set(targetfield, true, g_x[i], g_y[i], g_z[i], g_t[i], d, 2,0,m->c20);
			bgq_gaugefield_double_set(targetfield, true, g_x[i], g_y[i], g_z[i], g_t[i], d, 2,1,m->c21);
			bgq_gaugefield_double_set(targetfield, true, g_x[i], g_y[i], g_z[i], g_t[i], d, 2,2,m->c22);
		}
	}
}





