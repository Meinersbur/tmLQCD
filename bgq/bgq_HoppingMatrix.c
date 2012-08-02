/*
 * bgq_HoppingMatrix.c
 *
 *  Created on: Jul 27, 2012
 *      Author: meinersbur
 */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif

#include "bgq_field.h"
#include "bgq.h"
#include "global.h"
#include "boundary.h"
#include <mpi.h>
#include <omp.h>

#include "bgq_HoppingMatrix.h"


#define SETUP_PERSISTENT_SEND(ISODD, DIM, DIR) \
	MPI_CHECK(MPI_Rsend_init(&sendsurface->dim ## DIM[(ISODD)][(-(DIR)+1)/3], sizeof(sendsurface->dim ## DIM[(ISODD)][(-(DIR)+1)/3]), MPI_BYTE, g_cart_shift[joinDimdir((DIM),(DIR))], g_cart_tag[(ISODD)][joinDimdir((DIM),(DIR))], g_cart, &g_cart_send_requests[(ISODD)][joinDimdir((DIM),(DIR))]))

#define SETUP_PERSISTENT_RECV(ISODD, DIM, DIR) \
	MPI_CHECK(MPI_Recv_init(&recvsurface->dim ## DIM[(ISODD)][(-(DIR)+1)/3], sizeof(recvsurface->dim ## DIM[(ISODD)][(-(DIR)+1)/3]), MPI_BYTE, g_cart_shift[joinDimdir((DIM),(DIR))], g_cart_tag[(ISODD)][joinDimdir((DIM),(DIR))], g_cart, &g_cart_recv_requests[(ISODD)][joinDimdir((DIM),(DIR))]))







bgq_weylfield_double weylxchange_recv_double[6];
bgq_weylfield_double weylxchange_send_double[6];
size_t weylxchange_size_double[3];
int weylexchange_destination[6];

void bgq_init() {
	weylxchange_size_double[T_UP/2] = PHYSICAL_LX*PHYSICAL_LY * PHYSICAL_LZV * sizeof(bgq_weylsite_double);
	weylxchange_size_double[X_UP/2] = PHYSICAL_LT*PHYSICAL_LY * PHYSICAL_LZV * sizeof(bgq_weylsite_double);
	weylxchange_size_double[Y_UP/2] = PHYSICAL_LT*PHYSICAL_LX * PHYSICAL_LZV * sizeof(bgq_weylsite_double);

	for (int d = T_UP; d<=Y_DOWN; d+=1) {
		size_t size = weylxchange_size_double[d/2];
		weylxchange_recv_double[d] = (bgq_weylfield_double)malloc_aligned(size,128);
		weylxchange_send_double[d] = (bgq_weylfield_double)malloc_aligned(size,128);
	}

	weylexchange_destination[T_UP] = g_nb_t_up;
	weylexchange_destination[T_DOWN] = g_nb_t_dn;
	weylexchange_destination[X_UP] = g_nb_x_up;
	weylexchange_destination[X_DOWN] = g_nb_x_dn;
	weylexchange_destination[Y_UP] = g_nb_y_up;
	weylexchange_destination[Y_DOWN] = g_nb_y_dn;
}

void bgq_destroy() {
	for (int d = T_UP; d<=Y_DOWN; d+=1) {
		free(weylxchange_recv_double[d]);
		weylxchange_recv_double[d] = NULL;
		free(weylxchange_send_double[d]);
		weylxchange_send_double[d] = NULL;
	}
}




#define BGQ_HM_NOFUNC 1
#define BGQ_HM_ZLINE_NOFUNC 1
#define BGQ_HM_SITE_NOFUNC 1
#define BGQ_HM_DIR_NOFUNC 1

#define BGQ_HM_BORDER_NOFUNC 1
#define BGQ_HM_BORDERDIST_NOFUNC 1

#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
//#pragma GCC diagnostic ignored "-Wunused-variable"

void bgq_HoppingMatrix_borderzline_double(bool isOdd, bgq_spinorfield_double spinorfield, bgq_spinorfield_double targetfield, bgq_gaugefield_double gaugefield, int t, int x, int y) {
#define BGQ_HM_ZLINE_ZLINEINDENT -1
#define BGQ_HM_ZUP_WEYLREAD -1
#define BGQ_HM_ZDOWN_WEYLREAD -1
#define BGQ_HM_XUP_WEYLREAD -1
#define BGQ_HM_XDOWN_WEYLREAD -1
#define BGQ_HM_YUP_WEYLREAD -1
#define BGQ_HM_YDOWN_WEYLREAD -1
#include "bgq_HoppingMatrix_zline.inc.c"
}



void bgq_HoppingMatrix_double(bool isOdd, bgq_spinorfield_double spinorfield, bgq_spinorfield_double targetfield, bgq_gaugefield_double gaugefield) {

	MPI_Request request_recv[6];
	for (int d = T_UP; d<=Y_DOWN; d+=1) {
		MPI_CHECK(MPI_Irecv(weylxchange_recv_double[d], weylxchange_size_double[d/2], MPI_BYTE, weylexchange_destination[d], d, MPI_COMM_WORLD, &request_recv[d]));
	}

	MPI_CHECK(MPI_Barrier(g_cart_grid)); // To ensure that all ranks started the receive requests (necessary? how expensive is this?)


#pragma omp parallel for schedule(static)
	for (int xyz = 0; xyz < SURFACE_ZLINES_TOTAL*PHYSICAL_LZV; xyz+=1) {
		WORKLOAD_DECL(xyz,SURFACE_ZLINES_TOTAL*PHYSICAL_LZV);
		const int dir = WORKLOAD_CHUNK(6);
		const int zv = WORKLOAD_PARAM(PHYSICAL_LZV);
		switch (dir) {
		case T_UP:{
			const int y = WORKLOAD_PARAM(PHYSICAL_LY);
			const int x = WORKLOAD_PARAM(PHYSICAL_LX);
			const int t = PHYSICAL_LT-1;

#define BGQ_HM_TUP_WEYL_SEND 1
#include "bgq_HoppingMatrix_tup.inc.c"

		}break;
		case T_DOWN:{
			const int y = WORKLOAD_PARAM(PHYSICAL_LY);
						const int x = WORKLOAD_PARAM(PHYSICAL_LX);
						const int t = PHYSICAL_LT-1;

#define BGQ_HM_TDOWN_WEYL_SEND 1
#include "bgq_HoppingMatrix_tdown.inc.c"

		}break;
		case X_UP:{
			const int y = WORKLOAD_PARAM(PHYSICAL_LY);
						const int x = PHYSICAL_LX-1;
						const int t = WORKLOAD_PARAM(PHYSICAL_LT);

#define BGQ_HM_XUP_WEYL_SEND 1
#include "bgq_HoppingMatrix_xup.inc.c"

		}break;
		case X_DOWN:{
			const int y = WORKLOAD_PARAM(PHYSICAL_LY);
						const int x = PHYSICAL_LX-1;
						const int t = WORKLOAD_PARAM(PHYSICAL_LT);

#define BGQ_HM_XDOWN_WEYL_SEND 1
#include "bgq_HoppingMatrix_xdown.inc.c"

		}break;
		case Y_UP:{
			const int y = PHYSICAL_LY-1;
								const int x = WORKLOAD_PARAM(PHYSICAL_LX);
								const int t = WORKLOAD_PARAM(PHYSICAL_LT);

#define BGQ_HM_YUP_WEYL_SEND 1
#include "bgq_HoppingMatrix_yup.inc.c"

		}break;
		case Y_DOWN:{
			const int y = 0;
									const int x = WORKLOAD_PARAM(PHYSICAL_LX);
									const int t = WORKLOAD_PARAM(PHYSICAL_LT);

#define BGQ_HM_XDOWN_WEYL_SEND 1
#include "bgq_HoppingMatrix_ydown.inc.c"

		}break;
		default:
			assert(false);
			break;
		}

		assert(xyz==0);
		assert(xyz_total==1);
	}

	MPI_Request request_send[6];
	for (int d = T_UP; d<=Y_DOWN; d+=1) {
		MPI_CHECK(MPI_Isend(weylxchange_send_double[d], weylxchange_size_double[d/2], MPI_BYTE, weylexchange_destination[d], d, MPI_COMM_WORLD, &request_send[d]));
	}




	// Body volume kernel, flush lines
#pragma omp parallel for schedule(static,1)
	for (int xyz = 0; xyz < BODY_ZLINES/2; xyz += 1) {
		WORKLOAD_DECL(xyz,BODY_ZLINES/2);
		const int y = WORKLOAD_PARAM((PHYSICAL_LY-2)/2)*2+isOdd+1;
		const int x = WORKLOAD_PARAM(PHYSICAL_LX-2)+1;
		const int t = WORKLOAD_PARAM(PHYSICAL_LT-2)+1;

#define BGQ_HM_ZLINE_ZLINEINDENT 0
		#include "bgq_HoppingMatrix_zline.inc.c"
	}

	// Body volume kernel, ragged lines
#pragma omp parallel for schedule(static,1)
	for (int xyz = 0; xyz < BODY_ZLINES/2; xyz += 1) {
		WORKLOAD_DECL(xyz,BODY_ZLINES/2);
		const int y = WORKLOAD_PARAM((PHYSICAL_LY-2)/2)*2+!isOdd+1;
		const int x = WORKLOAD_PARAM(PHYSICAL_LX-2)+1;
		const int t = WORKLOAD_PARAM(PHYSICAL_LT-2)+1;

#define BGQ_HM_ZLINE_ZLINEINDENT 1
		#include "bgq_HoppingMatrix_zline.inc.c"
	}

	MPI_Status weylxchange_recv_status[6];
#if 1
	MPI_CHECK(MPI_Waitall(6, request_recv, weylxchange_recv_status));
#ifndef NDEBUG
	for (int d = T_UP; d<=Y_DOWN; d+=1) {
		assert(get_MPI_count(&weylxchange_recv_status[d]) == weylxchange_size_double[d]);
	}
#endif
#else
	for (int d = T_UP; d<=Y_DOWN; d+=1) {
		MPI_CHECK(MPI_Wait(&request_recv[d], &weylxchange_recv_status[d]));
		assert(get_MPI_count(&weylxchange_recv_status[d]) == weylxchange_size_double[d]);
	}
#endif

#pragma omp parallel for schedule(static,1)
	for (int xyz = 0; xyz < BORDER_ZLINES_TOTAL; xyz += 1) {
		WORKLOAD_DECL(xyz,BORDER_ZLINES_TOTAL);

		int t;
		int x;
		int y;
		if (WORKLOAD_SPLIT(BORDER_ZLINES_ALLFACES)) {
			const bool p = WORKLOAD_CHUNK(2);
			if (WORKLOAD_SPLIT(2*BORDER_ZLINES_T)) {
				y = WORKLOAD_PARAM((PHYSICAL_LY-2)/2)*2+p+1;
				x = WORKLOAD_PARAM(PHYSICAL_LX-2)+1;
				t = WORKLOAD_CHUNK(2)*(PHYSICAL_LT-1);
			} else if (WORKLOAD_SPLIT(2*BORDER_ZLINES_X)) {
				y = WORKLOAD_PARAM((PHYSICAL_LY-2)/2)*2+p+1;
				t = WORKLOAD_PARAM(PHYSICAL_LT-2)+1;
				x = WORKLOAD_CHUNK(2)*(PHYSICAL_LX-1);
			} else {
				WORKLOAD_SPLIT(2*BORDER_ZLINES_Y);

				x = WORKLOAD_PARAM((PHYSICAL_LX-2)/2)*2+p+1;
				t = WORKLOAD_PARAM(PHYSICAL_LT-2)+1;
				y = WORKLOAD_CHUNK(2)*(PHYSICAL_LX-1);
			}
		} else if (WORKLOAD_SPLIT(BORDER_ZLINES_ALLEDGES)) {
			const bool p = WORKLOAD_CHUNK(2);
			if (WORKLOAD_SPLIT(4*BORDER_ZLINES_TX)) {
				const bool p = WORKLOAD_TILE(2);
				y = WORKLOAD_PARAM((PHYSICAL_LY-2)/2)*2+p+1;
				x = WORKLOAD_CHUNK(2)*(PHYSICAL_LX-1);
				t = WORKLOAD_CHUNK(2)*(PHYSICAL_LT-1);
			} else if (WORKLOAD_SPLIT(4*BORDER_ZLINES_TY)) {
				const bool p = WORKLOAD_TILE(2);
				x = WORKLOAD_PARAM((PHYSICAL_LX-2)/2)*2+p+1;
				y = WORKLOAD_CHUNK(2)*(PHYSICAL_LY-1);
				t = WORKLOAD_CHUNK(2)*(PHYSICAL_LT-1);
			} else {
				WORKLOAD_SPLIT(4*BORDER_ZLINES_XY);

				t = WORKLOAD_PARAM((PHYSICAL_LT-2)/2)*2+p+1;
				y = WORKLOAD_CHUNK(2)*(PHYSICAL_LY-1);
				x = WORKLOAD_CHUNK(2)*(PHYSICAL_LX-1);
			}
		} else {
			// vertices
			int bits = WORKLOAD_PARAM(8);

			y = (bits & 0x4)*(PHYSICAL_LY-1);
			x = (bits & 0x2)*(PHYSICAL_LX-1);
			t = (bits & 0x1)*(PHYSICAL_LT-1);
		}
		assert(xyz_total == 1);
		assert(xyz == 0);

		bgq_HoppingMatrix_borderzline_double(isOdd,spinorfield,targetfield,gaugefield,t,x,y);
	}
}

#undef BGQ_HM_NOFUNC
#undef BGQ_HM_ZLINE_NOFUNC
#undef BGQ_HM_SITE_NOFUNC
#undef BGQ_HM_DIR_NOFUNC

