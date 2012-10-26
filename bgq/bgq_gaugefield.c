/*
 * bgq_gaugefield.c
 *
 *  Created on: Oct 26, 2012
 *      Author: meinersbur
 */

#define BGQ_GAUGEFIELD_C_
#include "bgq_gaugefield.h"

#include "bgq_dispatch.h"

#include "geometry_eo.h"


void bgq_gaugefield_init() {
	bgq_indices_init();

	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		g_bgq_gaugefield_fromHalfvolume[isOdd] = malloc_aligned(PHYSICAL_VOLUME * sizeof(*g_bgq_gaugefield_fromHalfvolume[isOdd]), BGQ_ALIGNMENT_L2);
		g_bgq_gaugefield_fromSurface[isOdd] = malloc_aligned(PHYSICAL_SURFACE * sizeof(*g_bgq_gaugefield_fromSurface[isOdd]), BGQ_ALIGNMENT_L2);
		g_bgq_gaugefield_fromBody[isOdd] = malloc_aligned(PHYSICAL_BODY * sizeof(*g_bgq_gaugefield_fromBody[isOdd]), BGQ_ALIGNMENT_L2);
	}
}



typedef struct {
	double _Complex c[3][3];
} su3_array64;

static void bgq_gaugefield_worker_transferfrom(void *arg_untyped, size_t tid, size_t threads) {
	su3 **sourcefield = (su3**) arg_untyped;

	const size_t workload = PHYSICAL_VOLUME * PHYSICAL_LP;
	const size_t threadload = (workload + threads - 1) / threads;
	const size_t begin = tid * threadload;
	const size_t end = min(workload, begin + threadload);
	for (size_t it = begin; it < end; it += 1) {
		WORKLOAD_DECL(it, workload);
		bool isOdd = WORKLOAD_CHUNK(2);
		size_t ih_src = WORKLOAD_PARAM(PHYSICAL_VOLUME);
		WORKLOAD_CHECK

		bool isOdd_src = isOdd;
		bool isOdd_dst = !isOdd;

		for (size_t k_src = 0; k_src < PHYSICAL_LK; k_src += 1) {
			size_t t_src = bgq_halfvolume2t(isOdd_src, ih_src, k_src);
			size_t x_src = bgq_halfvolume2x(ih_src);
			size_t y_src = bgq_halfvolume2y(ih_src);
			size_t z_src = bgq_halfvolume2y(ih_src);
			bool isSurface = bgq_halfvolume2isSurface(isOdd_src, ih_src);

			for (size_t d_src = 0; d_src < PHYSICAL_LD; d_src += 1) {
				bgq_dimension dim = bgq_direction2dimension(d_src);
				bgq_direction d_dst = bgq_direction_revert(d_src);

				size_t t = t_src;
				size_t x = x_src;
				size_t y = y_src;
				size_t z = z_src;
				switch (d_src) {
				case TDOWN:
					t -= 1;
					break;
				case XDOWN:
					x -= 1;
					break;
				case YDOWN:
					y -= 1;
					break;
				case ZDOWN:
					z -= 1;
					break;
				}

				size_t ix = Index(t, x, y, z);/* lexic coordinate; g_ipt[t][x][y][z] is not defined for -1 coordinates */
				su3_array64 *m = (su3_array64*) &sourcefield[ix][dim];
				for (size_t i = 0; i < 3; i += 1) {
					for (size_t l = 0; l < 3; l += 1) {
						g_bgq_gaugefield_fromHalfvolume[isOdd][ih_src].su3[d_src].c[i][l][k_src] = m->c[i][l];
						if (isSurface) {
							size_t is_src = bgq_halfvolume2surface(isOdd_src, ih_src);
							g_bgq_gaugefield_fromSurface[isOdd][is_src].su3[d_src].c[i][l][k_src] = m->c[i][l];
						} else {
							size_t ib_src = bgq_halfvolume2body(isOdd_src, ih_src);
							g_bgq_gaugefield_fromBody[isOdd][ib_src].su3[d_src].c[i][l][k_src] = m->c[i][l];
						}
					}
				}
			}
		}
	}
}
void bgq_gaugefield_transferfrom(su3 **sourcefield) {
	bgq_master_call(&bgq_gaugefield_worker_transferfrom, sourcefield);
}
