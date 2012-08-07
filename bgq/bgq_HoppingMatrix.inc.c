
#ifndef BGQ_HOPPINGMATRIX_C_
#include "bgq_field.h"
#include "bgq.h"
#include "global.h"
#include "boundary.h"
#include <mpi.h>
#include <omp.h>
#include <stdbool.h>

#include "bgq_HoppingMatrix.h"

#define BGQ_HM_SUFFIX inc
#endif

#ifndef BGQ_HM_SUFFIX
#error Define a name suffix for HoppingMatrix
#define BGQ_HM_SUFFIX
#endif

#ifndef BGQ_HM_PRECISION
#define BGQ_HM_PRECISION double
#endif

#ifndef BGQ_HM_NOCOM
#define BGQ_HM_NOCOM 0
#endif

#define BGQ_HM_NOFUNC 1
#define BGQ_HM_ZLINE_NOFUNC 1
#define BGQ_HM_SITE_NOFUNC 1
#define BGQ_HM_DIR_NOFUNC 1

#define BGQ_HM_NAME CONCAT(bgq_HoppingMatrix_ ,BGQ_HM_SUFFIX)



// isOdd refers to the oddness of targetfield; spinorfield will have the opposite oddness
void BGQ_HM_NAME(bool isOdd, bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield) {

#if !BGQ_HM_NOCOM && defined(MPI)
	//if (g_proc_id == 0)
	//	fprintf(stderr, "MK HM Irecv\n");
	MPI_Request request_recv[6];
	for (int d = T_UP; d <= Y_DOWN; d += 1) {
		MPI_CHECK(MPI_Irecv(weylxchange_recv_double[d], weylxchange_size_double[d/2], MPI_BYTE, weylexchange_destination[d], d^1, MPI_COMM_WORLD, &request_recv[d]));
	}

	MPI_CHECK(MPI_Barrier(g_cart_grid)); // To ensure that all ranks started the receive requests (necessary? how expensive is this?)
#endif

	// Load constants
	bgq_vector4double_decl(qka0); // t
	bgq_cconst(qka0, ka0.re, ka0.im);
	bgq_vector4double_decl(qka1); // x
	bgq_cconst(qka1, ka1.re, ka1.im);
	bgq_vector4double_decl(qka2); // y
	bgq_cconst(qka2, ka2.re, ka2.im);
	bgq_vector4double_decl(qka3); // z
	bgq_cconst(qka3, ka3.re, ka3.im);

#pragma omp parallel for schedule(static)
	for (int xyz = 0; xyz < SURFACE_ZLINES_TOTAL*PHYSICAL_LZV; xyz+=1) {
		WORKLOAD_DECL(xyz, SURFACE_ZLINES_TOTAL*PHYSICAL_LZV);

		if (WORKLOAD_SPLIT(2*SURFACE_ZLINES_T*PHYSICAL_LZV)) {
			const int zv = WORKLOAD_PARAM(PHYSICAL_LZV);

			if (WORKLOAD_SPLIT(SURFACE_ZLINES_T)) {
				const int y = WORKLOAD_PARAM(PHYSICAL_LY);
				const int x = WORKLOAD_PARAM(PHYSICAL_LX);
				const int t = PHYSICAL_LT;
				const int z1 = zv*4 + ((t+x+y+isOdd)&1);
				const int z2 = zv*4 + ((t+x+y+isOdd)&1) + 2;

				#define BGQ_HM_TDOWN_WEYL_SEND 1
				#include "bgq_HoppingMatrix_tdown.inc.c"
			} else {
				const int y = WORKLOAD_PARAM(PHYSICAL_LY);
				const int x = PHYSICAL_LX;
				const int t = WORKLOAD_PARAM(PHYSICAL_LT);
				const int z1 = zv*4 + ((t+x+y+isOdd)&1);
				const int z2 = zv*4 + ((t+x+y+isOdd)&1) + 2;

				#define BGQ_HM_XDOWN_WEYL_SEND 1
				#include "bgq_HoppingMatrix_xdown.inc.c"
			}
		} else if (WORKLOAD_SPLIT(2*SURFACE_ZLINES_X*PHYSICAL_LZV)) {
			const int zv = WORKLOAD_PARAM(PHYSICAL_LZV);

			if (WORKLOAD_SPLIT(SURFACE_ZLINES_X)) {
				const int y = WORKLOAD_PARAM(PHYSICAL_LY);
				const int x = PHYSICAL_LX;
				const int t = WORKLOAD_PARAM(PHYSICAL_LT);
				const int z1 = zv*4 + ((t+x+y+isOdd)&1);
				const int z2 = zv*4 + ((t+x+y+isOdd)&1) + 2;

				#define BGQ_HM_XDOWN_WEYL_SEND 1
				#include "bgq_HoppingMatrix_xdown.inc.c"
			} else {
				const int y = WORKLOAD_PARAM(PHYSICAL_LY);
				const int x = -1;
				const int t = WORKLOAD_PARAM(PHYSICAL_LT);
				const int z1 = zv*4 + ((t+x+y+isOdd)&1);
				const int z2 = zv*4 + ((t+x+y+isOdd)&1) + 2;

				#define BGQ_HM_XUP_WEYL_SEND 1
				#include "bgq_HoppingMatrix_xup.inc.c"
			}
		} else {
			WORKLOAD_SPLIT(2*SURFACE_ZLINES_Y*PHYSICAL_LZV);
			const int zv = WORKLOAD_PARAM(PHYSICAL_LZV);

			if (WORKLOAD_SPLIT(SURFACE_ZLINES_Y)) {
				const int y = PHYSICAL_LY;
				const int x = WORKLOAD_PARAM(PHYSICAL_LX);
				const int t = WORKLOAD_PARAM(PHYSICAL_LT);
				const int z1 = zv*4 + ((t+x+y+isOdd)&1);
				const int z2 = zv*4 + ((t+x+y+isOdd)&1) + 2;

				#define BGQ_HM_YDOWN_WEYL_SEND 1
				#include "bgq_HoppingMatrix_ydown.inc.c"
			} else {
				const int y = -1;
				const int x = WORKLOAD_PARAM(PHYSICAL_LX);
				const int t = WORKLOAD_PARAM(PHYSICAL_LT);
				const int z1 = zv*4 + ((t+x+y+isOdd)&1);
				const int z2 = zv*4 + ((t+x+y+isOdd)&1) + 2;

				#define BGQ_HM_YUP_WEYL_SEND 1
				#include "bgq_HoppingMatrix_yup.inc.c"
			}
		}
		WORKLOAD_CHECK
	}


#if !BGQ_HM_NOCOM && defined(MPI)
	//if (g_proc_id == 0)
	//	fprintf(stderr, "MK HM Isend\n");
	MPI_Request request_send[6];
	for (int d = T_UP; d <= Y_DOWN; d += 1) {
		MPI_CHECK(MPI_Isend(weylxchange_send_double[d], weylxchange_size_double[d/2], MPI_BYTE, weylexchange_destination[d], d, MPI_COMM_WORLD, &request_send[d]));
	}
#endif

	// Body volume kernel, flush lines, wavefronting
#pragma omp parallel for schedule(static,1)
	for (int xyz = 0; xyz < BODY_ZLINES/2; xyz += 1) {
		WORKLOAD_DECL(xyz, BODY_ZLINES/2);
		const int t = WORKLOAD_CHUNK(PHYSICAL_LT-2) + 1;
		const int x = WORKLOAD_CHUNK(PHYSICAL_LX-2) + 1;
		const int y = WORKLOAD_CHUNK((PHYSICAL_LY-2)/2)*2 + ((t+x+!isOdd)&1) + 1;
		WORKLOAD_CHECK

#define BGQ_HM_ZLINE_ZLINEINDENT 0
#include "bgq_HoppingMatrix_zline.inc.c"
	}

	// Body volume kernel, ragged lines, wavefronting
#pragma omp parallel for schedule(static,1)
	for (int xyz = 0; xyz < BODY_ZLINES / 2; xyz += 1) {
		WORKLOAD_DECL(xyz, BODY_ZLINES/2);
		const int t = WORKLOAD_CHUNK(PHYSICAL_LT-2) + 1;
		const int x = WORKLOAD_CHUNK(PHYSICAL_LX-2) + 1;
		const int y = WORKLOAD_CHUNK((PHYSICAL_LY-2)/2)*2 + ((t+x+isOdd)&1) + 1;
		WORKLOAD_CHECK

#define BGQ_HM_ZLINE_ZLINEINDENT 1
#include "bgq_HoppingMatrix_zline.inc.c"
	}

#if !BGQ_HM_NOCOM && defined(MPI)
	//if (g_proc_id == 0)
	//	fprintf(stderr, "MK HM Waitall\n");
	MPI_Status weylxchange_recv_status[6];
	MPI_CHECK(MPI_Waitall(6, request_recv, weylxchange_recv_status));
#ifndef NDEBUG
	for (int d = T_UP; d <= Y_DOWN; d += 1) {
		//fprintf(stderr, "MK(rank: %d) Waitall got: %d, expected: %d\n", g_proc_id, get_MPI_count(&weylxchange_recv_status[d]), weylxchange_size_double[d]);
		assert(get_MPI_count(&weylxchange_recv_status[d]) == weylxchange_size_double[d/2]);
	}
#endif
#endif


#pragma omp parallel for schedule(static)
	for (int xyz = 0; xyz < BORDER_ZLINES_TOTAL; xyz += 1) {
		WORKLOAD_DECL(xyz, BORDER_ZLINES_TOTAL);

		int t;
		int x;
		int y;
		if (WORKLOAD_SPLIT(BORDER_ZLINES_ALLFACES)) {
			const bool p = WORKLOAD_CHUNK(2);
			if (WORKLOAD_SPLIT(BORDER_ZLINES_T)) {
				y = WORKLOAD_PARAM((PHYSICAL_LY-2)/2) * 2 + p + 1;
				x = WORKLOAD_PARAM(PHYSICAL_LX-2) + 1;
				t = WORKLOAD_CHUNK(2) * (PHYSICAL_LT - 1);
			} else if (WORKLOAD_SPLIT(BORDER_ZLINES_X)) {
				y = WORKLOAD_PARAM((PHYSICAL_LY-2)/2) * 2 + p + 1;
				t = WORKLOAD_PARAM(PHYSICAL_LT-2) + 1;
				x = WORKLOAD_CHUNK(2) * (PHYSICAL_LX - 1);
			} else {
				WORKLOAD_SPLIT(BORDER_ZLINES_Y);

				x = WORKLOAD_PARAM((PHYSICAL_LX-2)/2) * 2 + p + 1;
				t = WORKLOAD_PARAM(PHYSICAL_LT-2) + 1;
				y = WORKLOAD_CHUNK(2) * (PHYSICAL_LX - 1);
			}
		} else if (WORKLOAD_SPLIT(BORDER_ZLINES_ALLEDGES)) {
			const bool p = WORKLOAD_CHUNK(2);
			if (WORKLOAD_SPLIT(4*BORDER_ZLINES_TX/2)) {
				y = WORKLOAD_PARAM((PHYSICAL_LY-2)/2) * 2 + p + 1;
				x = WORKLOAD_CHUNK(2) * (PHYSICAL_LX - 1);
				t = WORKLOAD_CHUNK(2) * (PHYSICAL_LT - 1);
			} else if (WORKLOAD_SPLIT(4*BORDER_ZLINES_TY/2)) {
				x = WORKLOAD_PARAM((PHYSICAL_LX-2)/2) * 2 + p + 1;
				y = WORKLOAD_CHUNK(2) * (PHYSICAL_LY - 1);
				t = WORKLOAD_CHUNK(2) * (PHYSICAL_LT - 1);
			} else {
				WORKLOAD_SPLIT(4*BORDER_ZLINES_XY/2);

				t = WORKLOAD_PARAM((PHYSICAL_LT-2)/2) * 2 + p + 1;
				y = WORKLOAD_CHUNK(2) * (PHYSICAL_LY - 1);
				x = WORKLOAD_CHUNK(2) * (PHYSICAL_LX - 1);
			}
		} else {
			// vertices
			int bits = WORKLOAD_PARAM(8);

			y = ((bits>>2) & 0x1) * (PHYSICAL_LY - 1);
			x = ((bits>>1) & 0x1) * (PHYSICAL_LX - 1);
			t = ((bits>>0) & 0x1) * (PHYSICAL_LT - 1);
		}
		WORKLOAD_CHECK


		#define BGQ_HM_ZLINE_ZLINEINDENT -1
		#define BGQ_HM_TUP_WEYLREAD -1
		#define BGQ_HM_TDOWN_WEYLREAD -1
		#define BGQ_HM_XUP_WEYLREAD -1
		#define BGQ_HM_XDOWN_WEYLREAD -1
		#define BGQ_HM_YUP_WEYLREAD -1
		#define BGQ_HM_YDOWN_WEYLREAD -1
		#include "bgq_HoppingMatrix_zline.inc.c"
	}
}


#undef BGQ_HM_SUFFIX
#undef BGQ_HM_PRECISION
#undef BGQ_HM_NOCOM

#undef BGQ_HM_NAME

#undef BGQ_HM_NOFUNC
#undef BGQ_HM_ZLINE_NOFUNC
#undef BGQ_HM_SITE_NOFUNC
#undef BGQ_HM_DIR_NOFUNC
