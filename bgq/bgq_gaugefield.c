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


static bool g_bgq_gaugefield_initialized = false;

void bgq_gaugefield_init() {
	if (g_bgq_gaugefield_initialized)
		return;
	g_bgq_gaugefield_initialized = true;

	bgq_indices_init();
	for (size_t isOdd = false; isOdd <= true; isOdd += 1) {
		g_bgq_gaugefield_fromCollapsed_src[isOdd] = malloc_aligned(PHYSICAL_VOLUME * sizeof(*g_bgq_gaugefield_fromCollapsed_src[isOdd]), BGQ_ALIGNMENT_L2);
		g_bgq_gaugefield_fromCollapsed_dst[isOdd] = malloc_aligned(PHYSICAL_VOLUME * sizeof(*g_bgq_gaugefield_fromCollapsed_dst[isOdd]), BGQ_ALIGNMENT_L2);
	}
}


static bgq_su3matrix bgq_gauge_coord_encode(scoord t, scoord x, scoord y, scoord z, bgq_direction d, bool isSrc) {
	size_t t_global = bgq_local2global_t(t);
	size_t x_global = bgq_local2global_x(x);
	size_t y_global = bgq_local2global_y(y);
	size_t z_global = bgq_local2global_z(z);
	if (isSrc) {
		bgq_direction_move_global(&t_global,&x_global,&y_global,&z_global,d);
		d = bgq_direction_revert(d);
	}

	bgq_su3matrix gauge={{{0}}};
	gauge.c[0][0] = t_global;
	gauge.c[0][1] = x_global;
	gauge.c[0][2] = y_global;
	gauge.c[1][0] = z_global;
	gauge.c[1][1] = d;

	return gauge;
}


static void bgq_gaugeveck_store(bgq_gaugesu3 *target, ucoord k, bgq_su3matrix data) {
	for (size_t i = 0; i < 3; i += 1) {
		for (size_t l = 0; l < 3; l += 1) {
			target->c[i][l][k] = data.c[i][l];
		}
	}
}


static void bgq_gaugeveck_written(bgq_gaugesu3 *target, ucoord k, ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc) {
#ifdef BGQ_COORDCHECK
	bgq_su3matrix coord = bgq_gauge_coord_encode(t,x,y,z,d,isSrc);
	bgq_gaugeveck_store(target, k, coord);
#endif
}



typedef struct {
	double _Complex c[3][3];
} su3_array64;

static void bgq_gaugefield_worker_transferfrom(void *arg_untyped, size_t tid, size_t threads) {
	su3 **sourcefield = (su3**)arg_untyped;

	const size_t workload = PHYSICAL_VOLUME * PHYSICAL_LP;
	const size_t threadload = (workload + threads - 1) / threads;
	const size_t begin = tid * threadload;
	const size_t end = min_sizet(workload, begin + threadload);
	for (size_t it = begin; it < end; it += 1) {
		WORKLOAD_DECL(it, workload);
		bool isOdd = WORKLOAD_CHUNK(2);
		size_t ih_src = WORKLOAD_PARAM(PHYSICAL_VOLUME);
		WORKLOAD_CHECK

		bool isOdd_src = isOdd;
		bool isOdd_dst = !isOdd;
		size_t ic_src = bgq_halfvolume2collapsed(isOdd_src, ih_src);

		for (size_t k_src = 0; k_src < PHYSICAL_LK; k_src += 1) {
			size_t t_src = bgq_halfvolume2t(isOdd_src, ih_src, k_src);
			size_t x_src = bgq_halfvolume2x(ih_src);
			size_t y_src = bgq_halfvolume2y(ih_src);
			size_t z_src = bgq_halfvolume2z(ih_src);

			for (size_t d_src = 0; d_src < PHYSICAL_LD; d_src += 1) {
				//bgq_dimension dim = bgq_direction2dimension(d_src);
				bgq_direction d_dst = bgq_direction_revert(d_src);

				ucoord t_dst = t_src;
				ucoord x_dst = x_src;
				ucoord y_dst = y_src;
				ucoord z_dst = z_src;
				bgq_direction_move_local(&t_dst, &x_dst, &y_dst, &z_dst, d_src);
				ucoord ic_dst = bgq_local2collapsed(t_dst, x_dst, y_dst, z_dst);
				ucoord k_dst = bgq_local2k(t_dst, x_dst, y_dst, z_dst);
				ucoord tv_dst = bgq_local2tv(t_dst, x_dst, y_dst, z_dst);

				if (ic_dst==0 && !isOdd_dst && d_dst==TDOWN) {
					int a = 0;
				}

				scoord t = t_src;
				scoord x = x_src;
				scoord y = y_src;
				scoord z = z_src;
				bgq_dimension dim;
				switch (d_dst) {
				case TUP:
					t -= 1;
					dim = DIM_T;
					break;
				case TDOWN:
					//t += 1;
					dim = DIM_T;
					break;
				case XUP:
					x -= 1;
					dim = DIM_X;
					break;
				case XDOWN:
					//x += 1;
					dim = DIM_X;
					break;
				case YUP:
					y -= 1;
					dim = DIM_Y;
					break;
				case YDOWN:
					//y += 1;
					dim = DIM_Y;
					break;
				case ZUP:
					z -= 1;
					dim = DIM_Z;
					break;
				case ZDOWN:
					//z += 1;
					dim = DIM_Z;
					break;
				default:
					UNREACHABLE
					break;
				}

				size_t ix = Index(t, x, y, z);/* lexic coordinate; g_ipt[t][x][y][z] is not defined for -1 coordinates */
				su3_array64 *m = (su3_array64*) &sourcefield[ix][dim];
				for (size_t i = 0; i < 3; i += 1) {
					for (size_t l = 0; l < 3; l += 1) {
						COMPLEX_PRECISION val = m->c[i][l];
						//TODO: We may multiply with ka0,ka1,ka2,ka3 here instead of inside HoppingMatrix
						g_bgq_gaugefield_fromCollapsed_src[isOdd_src][ic_src].su3[d_dst].c[i][l][k_src] = val;
						g_bgq_gaugefield_fromCollapsed_dst[isOdd_dst][ic_dst].su3[d_dst].c[i][l][k_dst] = val;
					}
				}
				bgq_gaugeveck_written(&g_bgq_gaugefield_fromCollapsed_src[isOdd_src][ic_src].su3[d_dst], k_src, t_src, x_src, y_src, z_src, d_src, true);
				bgq_gaugeveck_written(&g_bgq_gaugefield_fromCollapsed_dst[isOdd_dst][ic_dst].su3[d_dst], k_dst, t_src, x_src, y_src, z_src, d_src, true);
			}
		}
	}
}
void bgq_gaugefield_transferfrom(su3 **sourcefield) {
	bgq_master_call(&bgq_gaugefield_worker_transferfrom, sourcefield);
}


void bgq_gauge_expect(bgq_su3matrix gauge, ucoord t, ucoord x, ucoord y, ucoord z, bgq_direction d, bool isSrc) {
#ifdef BGQ_COORDCHECK
	size_t t_global = bgq_local2global_t(t);
	size_t x_global = bgq_local2global_x(x);
	size_t y_global = bgq_local2global_y(y);
	size_t z_global = bgq_local2global_z(z);
	if (isSrc) {
		bgq_direction_move_global(&t_global, &x_global, &y_global, &z_global, d);
		d = bgq_direction_revert(d);
	}

	assert(gauge.c[0][0] == t_global);
	assert(gauge.c[0][1] == x_global);
	assert(gauge.c[0][2] == y_global);
	assert(gauge.c[1][0] == z_global);
	assert(gauge.c[1][1] == d);
	assert(gauge.c[1][2] == 0);
	assert(gauge.c[2][0] == 0);
	assert(gauge.c[2][1] == 0);
	assert(gauge.c[2][2] == 0);
#endif
}
