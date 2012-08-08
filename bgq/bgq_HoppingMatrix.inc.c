
#ifndef BGQ_HOPPINGMATRIX_C_
#include "bgq_field.h"
#include "bgq.h"
#include "global.h"
#include "boundary.h"
#include <mpi.h>
#include <omp.h>
#include <stdbool.h>

#include "bgq_HoppingMatrix.h"
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



// isOdd refers to the oddness of targetfield; spinorfield will have the opposite oddness
void HoppingMatrix(bool isOdd, bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield) {

#if !BGQ_HM_NOCOM && defined(MPI)
	//TODO: Persistent communication
	//if (g_proc_id == 0)
	//	fprintf(stderr, "MK HM Irecv\n");
	MPI_Request request_recv[6];
	for (int d = T_UP; d <= Y_DOWN; d += 1) {
		MPI_CHECK(MPI_Irecv(weylxchange_recv_double[d], weylxchange_size_double[d/2], MPI_BYTE, weylexchange_destination[d], d^1, MPI_COMM_WORLD, &request_recv[d]));
	}

	//MPI_CHECK(MPI_Barrier(g_cart_grid)); // To ensure that all ranks started the receive requests (necessary? how expensive is this?)
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

// TDOWN
#pragma omp parallel for schedule(static)
	for (int xy = 0; xy < PHYSICAL_LXV*PHYSICAL_LY*PHYSICAL_LZ; xy+=1) {
		WORKLOAD_DECL(xy,PHYSICAL_LXV*PHYSICAL_LY*PHYSICAL_LZ);
		const int t = -1;
		const int y = WORKLOAD_PARAM(PHYSICAL_LY);
		const int z = WORKLOAD_PARAM(PHYSICAL_LZ);
		const int xv = WORKLOAD_PARAM(PHYSICAL_LXV);
		WORKLOAD_CHECK

		const int x1 = ((isOdd+t+y+z)&1)+xv*2;
		bgq_spinorsite_double *spinorsite1_tup = BGQ_SPINORSITE(spinorfield, !isOdd, t+1, x1, y, z, 0, 2);
		bgq_su3_spinor_decl(spinor1_tup);
		bgq_su3_spinor_double_load_left(spinor1_tup, spinorsite1_tup);

		const int x2 = x1+2;
		bgq_spinorsite_double *spinorsite2_tup = BGQ_SPINORSITE(spinorfield, !isOdd, t+1, x2, y, z, 0, 2);
		bgq_su3_spinor_decl(spinor2_tup);
		bgq_su3_spinor_double_load_left(spinor2_tup, spinorsite2_tup);

		bgq_su3_spinor_decl(spinor_tup);
		bgq_su3_spinor_merge(spinor_tup, spinor1_tup, spinor2_tup);

		// Compute its halfspinor
		bgq_su3_weyl_decl(weyl_tup);
		bgq_su3_vadd(weyl_tup_v0, spinor_tup_v0, spinor_tup_v2);
		bgq_su3_vadd(weyl_tup_v1, spinor_tup_v1, spinor_tup_v3);

		// Store the halfspinor to be transfered to the neighbor node
		bgq_weylsite_double *weylsite_tup = BGQ_WEYLSITE_T(weylxchange_send_double[T_DOWN/*!!!*/], isOdd, t, xv, y, z, x1, x2);
		bgq_su3_weyl_double_store(weylsite_tup, weyl_tup);
	}


	// TUP
#pragma omp parallel for schedule(static)
	for (int xyz = 0; xyz < PHYSICAL_LXV*PHYSICAL_LY*PHYSICAL_LZ; xyz+=1) {
		WORKLOAD_DECL(xyz,PHYSICAL_LXV*PHYSICAL_LY*PHYSICAL_LZ);
		const int t = LOCAL_LT;
		const int y = WORKLOAD_PARAM(PHYSICAL_LY);
		const int z = WORKLOAD_PARAM(PHYSICAL_LZ);
		const int xv = WORKLOAD_PARAM(PHYSICAL_LXV);
		WORKLOAD_CHECK

		const int x1 = ((isOdd+t+y+z)&1)+xv*2;
		bgq_spinorsite_double *spinorsite1_tdown = BGQ_SPINORSITE(spinorfield, !isOdd, t-1, x1, y, z, LOCAL_LT-3, LOCAL_LT-1);
		bgq_su3_spinor_decl(spinor1_tdown);
		bgq_su3_spinor_double_load_right(spinor1_tdown, spinorsite1_tdown);

		const int x2 = x1+2;
		bgq_spinorsite_double *spinorsite2_tdown = BGQ_SPINORSITE(spinorfield, !isOdd, t-1, x2, y, z, LOCAL_LT-3, LOCAL_LT-1);
		bgq_su3_spinor_decl(spinor2_tdown);
		bgq_su3_spinor_double_load_right(spinor2_tdown, spinorsite2_tdown);

		bgq_su3_spinor_decl(spinor_tdown);
		bgq_su3_spinor_merge(spinor_tdown, spinor1_tdown, spinor2_tdown);

		// Compute its halfspinor
		bgq_su3_weyl_decl(weyl_tdown);
		bgq_su3_vsub(weyl_tdown_v0, spinor_tdown_v0, spinor_tdown_v2);
		bgq_su3_vsub(weyl_tdown_v1, spinor_tdown_v1, spinor_tdown_v3);

		// Store the halfspinor to be transfered to the neighbor node
		bgq_weylsite_double *weylsite_tdown = BGQ_WEYLSITE_T(weylxchange_send_double[T_UP/*!!!*/], isOdd, t, xv, y, z, x1, x2);
		bgq_su3_weyl_double_store(weylsite_tdown, weyl_tdown);
	}


	// XDOWN
#pragma omp parallel for schedule(static)
	for (int tyz = 0; tyz < PHYSICAL_LTV*PHYSICAL_LY*PHYSICAL_LZ; tyz += 1) {
		WORKLOAD_DECL(tyz,PHYSICAL_LTV*PHYSICAL_LY*PHYSICAL_LZ);
		const int tv = WORKLOAD_PARAM(PHYSICAL_LTV);
		const int x = -1;
		const int y = WORKLOAD_PARAM(PHYSICAL_LY);
		const int z = WORKLOAD_PARAM(PHYSICAL_LZ);
		WORKLOAD_CHECK

		const int t1 = ((isOdd+x+y+z)&1)+tv*2;
		const int t2 = t1 + 2;

#define BGQ_HM_XUP_WEYL_SEND 1
#include "bgq_HoppingMatrix_xup.inc.c"
	}


	// XUP
#pragma omp parallel for schedule(static)
	for (int tyz = 0; tyz < PHYSICAL_LTV*PHYSICAL_LY*PHYSICAL_LZ; tyz += 1) {
		WORKLOAD_DECL(tyz,PHYSICAL_LTV*PHYSICAL_LY*PHYSICAL_LZ);
		const int tv = WORKLOAD_PARAM(PHYSICAL_LTV);
		const int x = LOCAL_LX;
		const int y = WORKLOAD_PARAM(PHYSICAL_LY);
		const int z = WORKLOAD_PARAM(PHYSICAL_LZ);
		WORKLOAD_CHECK

		const int t1 = ((isOdd+x+y+z)&1)+tv*2;
		const int t2 = t1 + 2;

#define BGQ_HM_XDOWN_WEYL_SEND 1
#include "bgq_HoppingMatrix_xdown.inc.c"
	}



	// YDOWN
#pragma omp parallel for schedule(static)
	for (int txz = 0; txz < PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LZ; txz += 1) {
		WORKLOAD_DECL(txz,PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LZ);
		const int tv = WORKLOAD_PARAM(PHYSICAL_LTV);
		const int x = WORKLOAD_PARAM(PHYSICAL_LX);
		const int y = -1;
		const int z = WORKLOAD_PARAM(PHYSICAL_LZ);
		WORKLOAD_CHECK

		const int t1 = ((isOdd+x+y+z)&1)+tv*2;
		const int t2 = t1 + 2;

#define BGQ_HM_YUP_WEYL_SEND 1
#include "bgq_HoppingMatrix_yup.inc.c"
	}


	// YUP
#pragma omp parallel for schedule(static)
	for (int txz = 0; txz < PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LZ; txz += 1) {
		WORKLOAD_DECL(txz,PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LZ);
		const int tv = WORKLOAD_PARAM(PHYSICAL_LTV);
		const int x = WORKLOAD_PARAM(PHYSICAL_LX);
		const int y = LOCAL_LY;
		const int z = WORKLOAD_PARAM(PHYSICAL_LZ);
		WORKLOAD_CHECK

		const int t1 = ((isOdd+x+y+z)&1)+tv*2;
		const int t2 = t1 + 2;

#define BGQ_HM_YDOWN_WEYL_SEND 1
#include "bgq_HoppingMatrix_ydown.inc.c"
	}


#if !BGQ_HM_NOCOM && defined(MPI)
	//if (g_proc_id == 0)
	//	fprintf(stderr, "MK HM Isend\n");
	MPI_Request request_send[6];
	for (int d = T_UP; d <= Y_DOWN; d += 1) {
		MPI_CHECK(MPI_Isend(weylxchange_send_double[d], weylxchange_size_double[d/2], MPI_BYTE, weylexchange_destination[d], d^1, MPI_COMM_WORLD, &request_send[d]));
	}
#endif


#pragma omp parallel for schedule(static,1)
	for (int txy = 0; txy < (PHYSICAL_LTV-2)*(PHYSICAL_LX-2)*(PHYSICAL_LY-2); txy +=1) {
		WORKLOAD_DECL(txy, (PHYSICAL_LTV-2)*(PHYSICAL_LX-2)*(PHYSICAL_LY-2));
		const int tv = WORKLOAD_CHUNK(PHYSICAL_LTV-2) + 1;
		const int x = WORKLOAD_CHUNK(PHYSICAL_LX-2) + 1;
		const int y = WORKLOAD_CHUNK(PHYSICAL_LY-2) + 1;

		int t1 = tv*PHYSICAL_LP*PHYSICAL_LK;
		int t2 = tv*PHYSICAL_LP*PHYSICAL_LK + 2;

#define BGQ_HM_ZLINE_ID KERNEL
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
	for (int xyzt = 0; xyzt < SURFACE_ZLINES; xyzt+=1) {
		WORKLOAD_DECL(xyzt, SURFACE_ZLINES);
		int tv;
		int x;
		int y;
		if (WORKLOAD_SPLIT(SURFACE_FACE_ZLINES)) {
			if (WORKLOAD_SPLIT(2*SURFACE_T_ZLINES)) {
				tv = WORKLOAD_PARAM(2)*(PHYSICAL_LTV-1);
				x = WORKLOAD_PARAM(PHYSICAL_LX-2)+1;
				y = WORKLOAD_PARAM(PHYSICAL_LY-2)+1;
			} else if (WORKLOAD_SPLIT(2*SURFACE_X_SITES)) {
				tv = WORKLOAD_PARAM(PHYSICAL_LTV-2)+1;
				x = WORKLOAD_PARAM(2)*(PHYSICAL_LX-1);
				y = WORKLOAD_PARAM(PHYSICAL_LY-2)+1;
			} else {
				WORKLOAD_SPLIT(2*SURFACE_Y_SITES);

				tv = WORKLOAD_PARAM(PHYSICAL_LTV-2)+1;
				x = WORKLOAD_PARAM(PHYSICAL_LX-2)+1;
				y = WORKLOAD_PARAM(2)*(PHYSICAL_LY-1);
			}
		} else if (WORKLOAD_SPLIT(SURFACE_EDGE_SITES)) {
			if (WORKLOAD_SPLIT(4*SURFACE_TX_SITES)) {
				tv = WORKLOAD_PARAM(2)*(PHYSICAL_LTV-1);
				x = WORKLOAD_PARAM(2)*(PHYSICAL_LX-1);
				y = WORKLOAD_PARAM(PHYSICAL_LY-2)+1;
			} else if (WORKLOAD_SPLIT(4*SURFACE_TY_SITES)) {
				tv = WORKLOAD_PARAM(2)*(PHYSICAL_LTV-1);
				x = WORKLOAD_PARAM(PHYSICAL_LX-2)+1;
				y = WORKLOAD_PARAM(2)*(PHYSICAL_LY-1);
			} else {
				WORKLOAD_SPLIT(4*SURFACE_XY_SITES);

				tv = WORKLOAD_PARAM(PHYSICAL_LTV-2)+1;
				x = WORKLOAD_PARAM(2)*(PHYSICAL_LX-1);
				y = WORKLOAD_PARAM(2)*(PHYSICAL_LY-1);
			}
		} else {
			// vertices
			int bits = WORKLOAD_PARAM(SURFACE_VERTICE_SITES);

			tv = ((bits>>0)&1) * (PHYSICAL_LTV - 1);
			x = ((bits>>1)&1) * (PHYSICAL_LX - 1);
			y = ((bits>>2)&1) * (PHYSICAL_LY - 1);
		}
		WORKLOAD_CHECK

		int t1 = ((isOdd+x+y)&1) + 2*tv;
		int t2 = t1 + 2;

#define BGQ_HM_ZLINE_ID SURFACE
#define BGQ_HM_TUP_WEYLREAD -1
#define BGQ_HM_TDOWN_WEYLREAD -1
#define BGQ_HM_XUP_WEYLREAD -1
#define BGQ_HM_XDOWN_WEYLREAD -1
#define BGQ_HM_YUP_WEYLREAD -1
#define BGQ_HM_YDOWN_WEYLREAD -1
#include "bgq_HoppingMatrix_zline.inc.c"
	}
}


#undef BGQ_HM_PRECISION
#undef BGQ_HM_NOCOM

#undef BGQ_HM_NOFUNC
#undef BGQ_HM_ZLINE_NOFUNC
#undef BGQ_HM_SITE_NOFUNC
#undef BGQ_HM_DIR_NOFUNC