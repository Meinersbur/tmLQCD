
#if 1
#undef BGQ_SPINORSITE
#define BGQ_SPINORSITE(spinorfield, isOdd, tv, x, y, z, t1, t2, isRead, isWrite) (bgq_spinorsite*)(isRead ? /*read*/(tmp=addr,addr+=2*4*3*2,tmp) : (isWrite ? /*write*/(tmp=writeaddr,writeaddr+=2*4*3*2,tmp) : /*prefetch*/(addr) ) )
#undef BGQ_SPINORSITE_LEFT
#define BGQ_SPINORSITE_LEFT BGQ_SPINORSITE
#undef BGQ_SPINORSITE_RIGHT
#define BGQ_SPINORSITE_RIGHT BGQ_SPINORSITE

#undef BGQ_WEYLSITE_T
#define BGQ_WEYLSITE_T(weylfield, isOdd, t, xv, y, z, x1, x2, isRead, isWrite) (bgq_weylsite*)(isRead ? /*read*/(tmp=addr,addr+=2*2*3*2,tmp) : (isWrite ? /*write*/(tmp=writeaddr,writeaddr+=2*2*3*2,tmp) : /*prefetch*/(addr) ) )
#undef BGQ_WEYLSITE_X
#define BGQ_WEYLSITE_X BGQ_WEYLSITE_T
#undef BGQ_WEYLSITE_Y
#define BGQ_WEYLSITE_Y BGQ_WEYLSITE_T
#undef BGQ_WEYLSITE_Z
#define BGQ_WEYLSITE_Z BGQ_WEYLSITE_T
#undef BGQ_WEYLSITE_T_LEFT
#define BGQ_WEYLSITE_T_LEFT BGQ_WEYLSITE_T
#undef BGQ_WEYLSITE_T_RIGHT
#define BGQ_WEYLSITE_T_RIGHT BGQ_WEYLSITE_T

#undef BGQ_GAUGESITE
#define BGQ_GAUGESITE(gaugefield,isOdd,tv,x,y,z,dir,t1,t2,isRead,isWrite) (bgq_gaugesite*)(isRead ? /*read*/(tmp=addr,addr+=2*3*3*2,tmp) : (isWrite ? /*write*/(tmp=writeaddr,writeaddr+=2*3*3*2,tmp) : /*prefetch*/(addr) ) )
#undef BGQ_GAUGESITE_LEFT
#define BGQ_GAUGESITE_LEFT BGQ_GAUGESITE
#undef BGQ_GAUGESITE_RIGHT
#define BGQ_GAUGESITE_RIGHT BGQ_GAUGESITE
#endif

#ifndef BGQ_HM_NOFUNC
#include "bgq_HoppingMatrix.h"
#include "bgq_field_double.h"
#include "bgq_field_float.h"
#include "bgq.h"
#include "global.h"
#include "boundary.h"
#include <mpi.h>
#include <omp.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include "bgq_utils.h"

#define BGQ_HM_ZLINE_NOFUNC 1
#define BGQ_HM_SITE_NOFUNC 1
#define BGQ_HM_DIR_NOFUNC 1
//#define DD1_L1P_Workaround 1


// isOdd refers to the oddness of targetfield; spinorfield will have the opposite oddness
void bgq_HoppingMatrix(bool isOdd, bgq_spinorfield_double targetfield, bgq_spinorfield_double spinorfield, bgq_gaugefield_double gaugefield, bool nocom, bool nooverlapcom, bool nokamul, bool withprefetchlist, bool withprefetchstream, bool withprefetchexplicit, bool coordcheck) {
#else
	{
#endif

	bgq_spinorfield_setOdd(targetfield, isOdd, true);
	bgq_spinorfield_setOdd(spinorfield, !isOdd, false);

	assert(!omp_in_parallel() && "Should not be called while in #pragma omp parallel");
	assert(omp_get_thread_num() == 0 && "Should be in a #pragma omp master");
	//assert(omp_in_parallel() && "Should be called while in #pragma omp parallel");

#if BGQ_FIELD_COORDCHECK
#pragma omp master
	{
		bgq_spinorfield_resetcoord(targetfield, isOdd, -1, -1, -1, -1);
		bgq_spinorfield_resetcoord(spinorfield, !isOdd, -1, -1, -1, -1);
		bgq_gaugefield_resetcoord(gaugefield, -1, -1, -1, -1);

		bgq_weylfield_t_resetcoord(weylxchange_send[TUP], LOCAL_LT - 1, !isOdd, -1, -1, -1, -1);
		bgq_weylfield_t_resetcoord(weylxchange_send[TDOWN], 0, !isOdd, -1, -1, -1, -1);
		bgq_weylfield_t_resetcoord(weylxchange_recv[TUP], LOCAL_LT, !isOdd, -1, -1, -1, -1);
		bgq_weylfield_t_resetcoord(weylxchange_recv[TDOWN], -1, !isOdd, -1, -1, -1, -1);

		bgq_weylfield_x_resetcoord(weylxchange_send[XUP], LOCAL_LX - 1, !isOdd, -1, -1, -1, -1);
		bgq_weylfield_x_resetcoord(weylxchange_send[XDOWN], 0, !isOdd, -1, -1, -1, -1);
		bgq_weylfield_x_resetcoord(weylxchange_recv[XUP], LOCAL_LX, !isOdd, -1, -1, -1, -1);
		bgq_weylfield_x_resetcoord(weylxchange_recv[XDOWN], -1, !isOdd, -1, -1, -1, -1);

		bgq_weylfield_y_resetcoord(weylxchange_send[YUP], LOCAL_LY - 1, !isOdd, -1, -1, -1, -1);
		bgq_weylfield_y_resetcoord(weylxchange_send[YDOWN], 0, !isOdd, -1, -1, -1, -1);
		bgq_weylfield_y_resetcoord(weylxchange_recv[YUP], LOCAL_LY, !isOdd, -1, -1, -1, -1);
		bgq_weylfield_y_resetcoord(weylxchange_recv[YDOWN], -1, !isOdd, -1, -1, -1, -1);
	}
#endif


#pragma omp parallel
	{
		const int threads = omp_get_num_threads();
		const int tid = omp_get_thread_num();

		PRECISION *tmp;
		PRECISION *addr = ((PRECISION*)spinorfield) + (size_t)tid * READTOTLENGTH / threads;
		addr = (PRECISION*)((size_t)addr & ~(32-1)); // alignment

		PRECISION *writeaddr = ((PRECISION*)targetfield) + (size_t)tid * WRITETOTLENGTH / threads;
		writeaddr = (PRECISION*)((size_t)writeaddr & ~(32-1)); // alignment


	// Load kamul constants
	bgq_vector4double_decl(qka0); // t
	bgq_cconst(qka0, ka0.re, ka0.im);
	bgq_vector4double_decl(qka1); // x
	bgq_cconst(qka1, ka1.re, ka1.im);
	bgq_vector4double_decl(qka2); // y
	bgq_cconst(qka2, ka2.re, ka2.im);
	bgq_vector4double_decl(qka3); // z
	bgq_cconst(qka3, ka3.re, ka3.im);


#pragma omp master
		{
#ifdef MPI
			if (!nocom) {
#ifndef NDEBUG
				for (direction d = TUP; d <= YDOWN; d += 1) {
					memset(weylxchange_recv[d], 0xFE, weylxchange_size[d / 2]);
					memset(weylxchange_send[d], 0xFF, weylxchange_size[d / 2]);
				}
#endif
				MPI_CHECK(MPI_Startall(lengthof(weylexchange_request_recv), weylexchange_request_recv));
			} else {
				for (direction d = TUP; d <= YDOWN; d += 1) {
					// To make it deterministic
					memset(weylxchange_recv[d], 0, weylxchange_size[d / 2]);
				}
			}
#endif
		}


#if BGQ_PREFETCH_LIST
		if (!noprefetchlist) {
			int targetindex = bgq_spinorfield_find_index(targetfield);
			int spinorindex = bgq_spinorfield_find_index(spinorfield);
			int totThreads = omp_get_num_threads();
			L1P_Pattern_t **patternhandle = &bgq_listprefetch_handle[targetindex + spinorindex][isOdd][nokamul][Kernel_ProcessorID()][nobody][nosurface][ilog(totThreads)];
			bool l1p_first = false;
			if (!*patternhandle) {
				//const uint64_t prefsize = 10000/*how to choose this?*/;
				const uint64_t prefsize = 12 * VOLUME / totThreads;
				int retval = L1P_GetPattern(patternhandle);
				if (retval == L1P_NOTCONFIGURED) {
					L1P_CHECK(L1P_PatternConfigure(prefsize));
					L1P_CHECK(L1P_GetPattern(patternhandle));
				} else {
					L1P_CHECK(retval);
					L1P_CHECK(L1P_AllocatePattern(prefsize,patternhandle));
				}
				l1p_first = true;
			}
			assert(*patternhandle);

			// Don't have acces to L1P_CFG_PF_USR (segfault)
			//L1P_CHECK(L1P_PatternSetEnable(1));

			L1P_CHECK(L1P_SetPattern(*patternhandle));
			L1P_CHECK(L1P_PatternStart(true));

			// Don't have acces to L1P_CFG_PF_USR (segfault)
			//int pattern_enable = 0;
			//L1P_CHECK(L1P_PatternGetEnable(&pattern_enable));
			//if (!pattern_enable)
			//	printf("error: List prefetcher not enabled\n");
		}
#endif


		if (!noweylsend) {
// TDOWN
#pragma omp for schedule(static)
				for (int xyz = 0; xyz < PHYSICAL_LXV * PHYSICAL_LY * PHYSICAL_LZ; xyz += 1) {
					WORKLOAD_DECL(xyz, PHYSICAL_LXV*PHYSICAL_LY*PHYSICAL_LZ);
					const int t = -1;
					const int y = WORKLOAD_PARAM(PHYSICAL_LY);
					const int z = WORKLOAD_PARAM(PHYSICAL_LZ);
					const int xv = WORKLOAD_PARAM(PHYSICAL_LXV);
					WORKLOAD_CHECK

					const int x1 = ((isOdd + t + y + z) & 1) + xv * PHYSICAL_LP * PHYSICAL_LK;
					bgq_spinorsite *spinorsite1_tup = BGQ_SPINORSITE_LEFT(spinorfield, !isOdd, 0, x1, y, z, t+1, 2, true,false);
					bgq_su3_spinor_decl(spinor1_tup);
					bgq_su3_spinor_load_left(spinor1_tup, spinorsite1_tup);

					const int x2 = x1 + 2;
					bgq_spinorsite *spinorsite2_tup = BGQ_SPINORSITE_LEFT(spinorfield, !isOdd, 0, x2, y, z, t+1, 2, true,false);
					bgq_su3_spinor_decl(spinor2_tup);
					bgq_su3_spinor_load_left(spinor2_tup, spinorsite2_tup);

					bgq_su3_spinor_decl(spinor_tup);
					bgq_su3_spinor_merge(spinor_tup, spinor1_tup, spinor2_tup);

					// Compute its halfspinor
					bgq_su3_weyl_decl(weyl_tup);
					bgq_su3_vadd(weyl_tup_v0, spinor_tup_v0, spinor_tup_v2);
					bgq_su3_vadd(weyl_tup_v1, spinor_tup_v1, spinor_tup_v3);

					// Store the halfspinor to be transfered to the neighbor node
					bgq_weylsite *weylsite_tup = BGQ_WEYLSITE_T(weylxchange_send[TDOWN/*!!!*/], !isOdd, t+1, xv, y, z, x1, x2, false, true);
					bgq_su3_weyl_zeroload(weylsite_tup);
					bgq_su3_weyl_store(weylsite_tup, weyl_tup);
					bgq_su3_weyl_flush(weylsite_tup);
					int a = 0;
				}

			// TUP
#pragma omp for schedule(static) nowait
			for (int xyz = 0; xyz < PHYSICAL_LXV * PHYSICAL_LY * PHYSICAL_LZ; xyz += 1) {
				WORKLOAD_DECL(xyz, PHYSICAL_LXV*PHYSICAL_LY*PHYSICAL_LZ);
				const int t = LOCAL_LT;
				const int y = WORKLOAD_PARAM(PHYSICAL_LY);
				const int z = WORKLOAD_PARAM(PHYSICAL_LZ);
				const int xv = WORKLOAD_PARAM(PHYSICAL_LXV);
				WORKLOAD_CHECK

				const int x1 = ((isOdd + t + y + z) & 1) + xv * PHYSICAL_LP * PHYSICAL_LK;
				bgq_spinorsite *spinorsite1_tdown = BGQ_SPINORSITE_RIGHT(spinorfield, !isOdd, PHYSICAL_LTV-1, x1, y, z, LOCAL_LT-3, t-1, true, false);
				bgq_su3_spinor_decl(spinor1_tdown);
				bgq_su3_spinor_load_right(spinor1_tdown, spinorsite1_tdown);

				const int x2 = x1 + 2;
				bgq_spinorsite *spinorsite2_tdown = BGQ_SPINORSITE_RIGHT(spinorfield, !isOdd, PHYSICAL_LTV-1, x2, y, z, LOCAL_LT-3, t-1, true, false);
				bgq_su3_spinor_decl(spinor2_tdown);
				bgq_su3_spinor_load_right(spinor2_tdown, spinorsite2_tdown);

				bgq_su3_spinor_decl(spinor_tdown);
				bgq_su3_spinor_merge(spinor_tdown, spinor1_tdown, spinor2_tdown);

				// Compute its halfspinor
				bgq_su3_weyl_decl(weyl_tdown);
				bgq_su3_vsub(weyl_tdown_v0, spinor_tdown_v0, spinor_tdown_v2);
				bgq_su3_vsub(weyl_tdown_v1, spinor_tdown_v1, spinor_tdown_v3);

				// Store the halfspinor to be transfered to the neighbor node
				bgq_weylsite *weylsite_tdown = BGQ_WEYLSITE_T(weylxchange_send[TUP/*!!!*/], !isOdd, t-1, xv, y, z, x1, x2, false, true);
				bgq_su3_weyl_zeroload(weylsite_tdown);
				bgq_su3_weyl_store(weylsite_tdown, weyl_tdown);
				bgq_su3_weyl_flush(weylsite_tdown);
				int a = 0;
			}

			// XDOWN
#pragma omp for schedule(static) nowait
			for (int tyz = 0; tyz < PHYSICAL_LTV * PHYSICAL_LY * PHYSICAL_LZ; tyz += 1) {
				WORKLOAD_DECL(tyz, PHYSICAL_LTV*PHYSICAL_LY*PHYSICAL_LZ);
				const int tv = WORKLOAD_PARAM(PHYSICAL_LTV);
				const int x = -1;
				const int y = WORKLOAD_PARAM(PHYSICAL_LY);
				const int z = WORKLOAD_PARAM(PHYSICAL_LZ);
				WORKLOAD_CHECK

				const int t1 = ((isOdd + x + y + z) & 1) + tv * PHYSICAL_LP * PHYSICAL_LK;
				const int t2 = t1 + 2;

#define BGQ_HM_XUP_WEYLREAD 0
#define BGQ_HM_XUP_WEYL_SEND 1
#include "bgq_HoppingMatrix_xup.inc.c"
#undef BGQ_HM_XUP_WEYLREAD
			}

			// XUP
#pragma omp for schedule(static) nowait
			for (int tyz = 0; tyz < PHYSICAL_LTV * PHYSICAL_LY * PHYSICAL_LZ; tyz += 1) {
				WORKLOAD_DECL(tyz, PHYSICAL_LTV*PHYSICAL_LY*PHYSICAL_LZ);
				const int tv = WORKLOAD_PARAM(PHYSICAL_LTV);
				const int x = LOCAL_LX;
				const int y = WORKLOAD_PARAM(PHYSICAL_LY);
				const int z = WORKLOAD_PARAM(PHYSICAL_LZ);
				WORKLOAD_CHECK

				const int t1 = ((isOdd + x + y + z) & 1) + tv * PHYSICAL_LP * PHYSICAL_LK;
				const int t2 = t1 + 2;

#define BGQ_HM_XDOWN_WEYLREAD 0
#define BGQ_HM_XDOWN_WEYL_SEND 1
#include "bgq_HoppingMatrix_xdown.inc.c"
#undef BGQ_HM_XDOWN_WEYLREAD
			}

			// YDOWN
#pragma omp for schedule(static) nowait
			for (int txz = 0; txz < PHYSICAL_LTV * PHYSICAL_LX * PHYSICAL_LZ; txz += 1) {
				WORKLOAD_DECL(txz, PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LZ);
				const int tv = WORKLOAD_PARAM(PHYSICAL_LTV);
				const int x = WORKLOAD_PARAM(PHYSICAL_LX);
				const int y = -1;
				const int z = WORKLOAD_PARAM(PHYSICAL_LZ);
				WORKLOAD_CHECK

				const int t1 = ((isOdd + x + y + z) & 1) + tv * PHYSICAL_LP * PHYSICAL_LK;
				const int t2 = t1 + 2;

#define BGQ_HM_YUP_WEYLREAD 0
#define BGQ_HM_YUP_WEYL_SEND 1
#include "bgq_HoppingMatrix_yup.inc.c"
#undef BGQ_HM_YUP_WEYLREAD
			}

			// YUP
#pragma omp for schedule(static)
			for (int txz = 0; txz < PHYSICAL_LTV * PHYSICAL_LX * PHYSICAL_LZ; txz += 1) {
				WORKLOAD_DECL(txz, PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LZ);
				const int tv = WORKLOAD_PARAM(PHYSICAL_LTV);
				const int x = WORKLOAD_PARAM(PHYSICAL_LX);
				const int y = LOCAL_LY;
				const int z = WORKLOAD_PARAM(PHYSICAL_LZ);
				WORKLOAD_CHECK

				const int t1 = ((isOdd + x + y + z) & 1) + tv * PHYSICAL_LP * PHYSICAL_LK;
				const int t2 = t1 + 2;

#define BGQ_HM_YDOWN_WEYLREAD 0
#define BGQ_HM_YDOWN_WEYL_SEND 1
#include "bgq_HoppingMatrix_ydown.inc.c"
#undef BGQ_HM_YDOWN_WEYLREAD
			}
		}


#if BGQ_PREFETCH_LIST
		if (!noprefetchlist) {
			L1P_PatternPause();
		}
#endif


#if BGQ_FIELD_COORDCHECK
#pragma omp master
		{
			if (!noweylsend) {
				bgq_weylfield_t_resetcoord(weylxchange_send[TDOWN], 0, !isOdd, 0, 0, 1, 1);
				bgq_weylfield_t_resetcoord(weylxchange_send[TUP], LOCAL_LT - 1, !isOdd, 0, 0, 1, 1);
				bgq_weylfield_x_resetcoord(weylxchange_send[XDOWN], 0, !isOdd, 0, 0, 1, 1);
				bgq_weylfield_x_resetcoord(weylxchange_send[XUP], LOCAL_LX - 1, !isOdd, 0, 0, 1, 1);
				bgq_weylfield_y_resetcoord(weylxchange_send[YDOWN], 0, !isOdd, 0, 0, 1, 1);
				bgq_weylfield_y_resetcoord(weylxchange_send[YUP], LOCAL_LY - 1, !isOdd, 0, 0, 1, 1);
			}
		}
#endif

#ifdef MPI
		if (!nocom && !nooverlap) {
#pragma omp master
			{
				//bgq_weylfield_foreach(weylxchange_send[TDOWN], TDOWN, true, isOdd, &bgq_setbgqval, BGQREF_TDOWN_SENDBUF);
				//bgq_weylfield_foreach(weylxchange_send[TUP], TUP, true, isOdd, &bgq_setbgqval, BGQREF_TUP_SENDBUF);
				//bgq_weylfield_foreach(weylxchange_send[XDOWN], XDOWN, true, isOdd, &bgq_setbgqval, BGQREF_XDOWN_SENDBUF);
				//bgq_weylfield_foreach(weylxchange_send[XUP], XUP, true, isOdd, &bgq_setbgqval, BGQREF_XUP_SENDBUF);

				MPI_CHECK(MPI_Startall(lengthof(weylexchange_request_send), weylexchange_request_send));
			}
		}
#endif


#if BGQ_PREFETCH_LIST
		if (!noprefetchlist) {
			L1P_PatternResume();
		}
#endif


////////////////////////////////////////////////////////////////////////////////
// Body kernel

		if (!nobody) {
#pragma omp for schedule(static,1)
			for (int txy = 0; txy < BODY_ZLINES; txy += 1) {
				WORKLOAD_DECL(txy, BODY_ZLINES);
				const int tv = WORKLOAD_PARAM(PHYSICAL_LTV-2) + 1;
				const int x = WORKLOAD_PARAM(PHYSICAL_LX-2) + 1;
				const int y = WORKLOAD_PARAM(PHYSICAL_LY-2) + 1;
				WORKLOAD_CHECK
				int t1 = tv*PHYSICAL_LP*PHYSICAL_LK + ((isOdd+x+y+0)&1);
				int t2 = t1 + PHYSICAL_LP;

#define BGQ_HM_ZLINE_ID KERNEL
#define BGQ_HM_TUP_WEYLREAD 0
#define BGQ_HM_TDOWN_WEYLREAD 0
#define BGQ_HM_XUP_WEYLREAD 0
#define BGQ_HM_XDOWN_WEYLREAD 0
#define BGQ_HM_YUP_WEYLREAD 0
#define BGQ_HM_YDOWN_WEYLREAD 0
#include "bgq_HoppingMatrix_zline.inc.c"
#undef BGQ_HM_TUP_WEYLREAD
#undef BGQ_HM_TDOWN_WEYLREAD
#undef BGQ_HM_XUP_WEYLREAD
#undef BGQ_HM_XDOWN_WEYLREAD
#undef BGQ_HM_YUP_WEYLREAD
#undef BGQ_HM_YDOWN_WEYLREAD
			}
		}

////////////////////////////////////////////////////////////////////////////////


#if BGQ_PREFETCH_LIST
		if (!noprefetchlist) {
			L1P_PatternPause();
		}
#endif


#ifdef MPI
		if (!nocom && nooverlap) {
#pragma omp master
			{
				//bgq_weylfield_foreach(weylxchange_send[TDOWN], TDOWN, true, isOdd, &bgq_setbgqval, BGQREF_TDOWN_SENDBUF);
				//bgq_weylfield_foreach(weylxchange_send[TUP], TUP, true, isOdd, &bgq_setbgqval, BGQREF_TUP_SENDBUF);
				//bgq_weylfield_foreach(weylxchange_send[XDOWN], XDOWN, true, isOdd, &bgq_setbgqval, BGQREF_XDOWN_SENDBUF);
				//bgq_weylfield_foreach(weylxchange_send[XUP], XUP, true, isOdd, &bgq_setbgqval, BGQREF_XUP_SENDBUF);

				MPI_CHECK(MPI_Startall(lengthof(weylexchange_request_send), weylexchange_request_send));
			}
		}


		if (!nocom) {
#pragma omp master
			{
				MPI_Status weylxchange_recv_statuses[6];
				MPI_CHECK(MPI_Waitall(lengthof(weylexchange_request_recv), weylexchange_request_recv, weylxchange_recv_statuses));
				//master_print("length=%d weylexchange_request_recv=%d weylxchange_recv_statuses=%d\n", (int)lengthof(weylexchange_request_recv), (int)weylexchange_request_recv, (int)weylxchange_recv_statuses);
				//bgq_weylfield_foreach(weylxchange_recv[TUP], TUP, false, isOdd, &bgq_setbgqval, BGQREF_TUP_RECVBUF);
				//bgq_weylfield_foreach(weylxchange_recv[TDOWN], TDOWN, false, isOdd, &bgq_setbgqval, BGQREF_TDOWN_RECVBUF);
				//bgq_weylfield_foreach(weylxchange_recv[XUP], XUP, false, isOdd, &bgq_setbgqval, BGQREF_XUP_RECVBUF);
				//bgq_weylfield_foreach(weylxchange_recv[XDOWN], XDOWN, false, isOdd, &bgq_setbgqval, BGQREF_XDOWN_RECVBUF);

				MPI_Status weylxchange_send_statuses[6];
				MPI_CHECK(MPI_Waitall(lengthof(weylexchange_request_send), weylexchange_request_send, weylxchange_send_statuses));
#ifndef NDEBUG
				for (int d = TUP; d <= YDOWN; d += 1) {
					assert(get_MPI_count(&weylxchange_recv_statuses[d]) == weylxchange_size[d/2]);
				}
#endif
			}
		}
#endif
		// Let all threads wait for the master thread which by itself waits for the MPI communication to finish
		// omp master has no implicit barrier
#pragma omp barrier


#if BGQ_PREFETCH_LIST
		if (!noprefetchlist) {
			L1P_PatternResume();
		}
#endif


		if (!nosurface) {
#pragma omp for schedule(static) nowait
			for (int ixyz = 0; ixyz < SURFACE_ZLINES; ixyz += 1) {
				//master_print("MK HM xyz=%d\n", xyz);

				WORKLOAD_DECL(ixyz, SURFACE_ZLINES);
				int tv;
				int x;
				int y;
				if (WORKLOAD_SPLIT(SURFACE_FACE_ZLINES)) {
					if (WORKLOAD_SPLIT(2*SURFACE_T_ZLINES)) {
						tv = WORKLOAD_PARAM(2) * (PHYSICAL_LTV - 1);
						x = WORKLOAD_PARAM(PHYSICAL_LX-2) + 1;
						y = WORKLOAD_PARAM(PHYSICAL_LY-2) + 1;
					} else if (WORKLOAD_SPLIT(2*SURFACE_X_ZLINES)) {
						tv = WORKLOAD_PARAM(PHYSICAL_LTV-2) + 1;
						x = WORKLOAD_PARAM(2) * (PHYSICAL_LX - 1);
						y = WORKLOAD_PARAM(PHYSICAL_LY-2) + 1;
					} else {
						WORKLOAD_SPLIT(2*SURFACE_Y_ZLINES);

						tv = WORKLOAD_PARAM(PHYSICAL_LTV-2) + 1;
						x = WORKLOAD_PARAM(PHYSICAL_LX-2) + 1;
						y = WORKLOAD_PARAM(2) * (PHYSICAL_LY - 1);
					}
				} else if (WORKLOAD_SPLIT(SURFACE_EDGE_ZLINES)) {
					if (WORKLOAD_SPLIT(4*SURFACE_TX_ZLINES)) {
						tv = WORKLOAD_PARAM(2) * (PHYSICAL_LTV - 1);
						x = WORKLOAD_PARAM(2) * (PHYSICAL_LX - 1);
						y = WORKLOAD_PARAM(PHYSICAL_LY-2) + 1;
					} else if (WORKLOAD_SPLIT(4*SURFACE_TY_ZLINES)) {
						tv = WORKLOAD_PARAM(2) * (PHYSICAL_LTV - 1);
						x = WORKLOAD_PARAM(PHYSICAL_LX-2) + 1;
						y = WORKLOAD_PARAM(2) * (PHYSICAL_LY - 1);
					} else {
						WORKLOAD_SPLIT(4*SURFACE_XY_ZLINES);

						tv = WORKLOAD_PARAM(PHYSICAL_LTV-2) + 1;
						x = WORKLOAD_PARAM(2) * (PHYSICAL_LX - 1);
						y = WORKLOAD_PARAM(2) * (PHYSICAL_LY - 1);
					}
				} else {
					// vertices
					int bits = WORKLOAD_PARAM(SURFACE_VERTICE_ZLINES);

					tv = ((bits >> 0) & 1) * (PHYSICAL_LTV - 1);
					x = ((bits >> 1) & 1) * (PHYSICAL_LX - 1);
					y = ((bits >> 2) & 1) * (PHYSICAL_LY - 1);
				}
				WORKLOAD_CHECK

				int t1 = ((isOdd + x + y) & 1) + tv * PHYSICAL_LP * PHYSICAL_LK;
				int t2 = t1 + 2;

#define BGQ_HM_ZLINE_ID SURFACE
#define BGQ_HM_TUP_WEYLREAD -1
#define BGQ_HM_TDOWN_WEYLREAD -1
#define BGQ_HM_XUP_WEYLREAD -1
#define BGQ_HM_XDOWN_WEYLREAD -1
#define BGQ_HM_YUP_WEYLREAD -1
#define BGQ_HM_YDOWN_WEYLREAD -1
#include "bgq_HoppingMatrix_zline.inc.c"
#undef BGQ_HM_TUP_WEYLREAD
#undef BGQ_HM_TDOWN_WEYLREAD
#undef BGQ_HM_XUP_WEYLREAD
#undef BGQ_HM_XDOWN_WEYLREAD
#undef BGQ_HM_YUP_WEYLREAD
#undef BGQ_HM_YDOWN_WEYLREAD
			}
		}


#if BGQ_PREFETCH_LIST
		if (!noprefetchlist) {
			int tid = omp_get_thread_num();
			L1P_PatternGetCurrentDepth(&fetch_depth[tid], &generate_depth[tid]);
			L1P_PatternStatus(&l1plist_status[tid]);
			L1P_PatternStop();
			//master_print("L1P_LIST: maxed=%d abandoned=%d finished=%d endoflist=%d fetch_depth=%d generate_depth=%d \n", st.s.maximum, st.s.abandoned, st.s.finished, st.s.endoflist, (int)fetch_depth, (int)generate_depth);
		}
#endif
	} // #pragma omp parallel



#if BGQ_FIELD_COORDCHECK
	if (!nobody && !nosurface && !noweylsend) {
		bgq_spinorfield_resetcoord(targetfield, isOdd, 0, 0, 1, 1);
		bgq_spinorfield_resetcoord(spinorfield, !isOdd, 8, 8, 0, 0);
		bgq_gaugefield_resetcoord(gaugefield, 0/*every second -1 coordinate is not read in even/odd iteration*/, 1, 0, 0);
	}

	if (!nosurface) {
		bgq_weylfield_t_resetcoord(weylxchange_send[TUP], LOCAL_LT - 1, !isOdd, 0, 0, 0, 0);
		bgq_weylfield_t_resetcoord(weylxchange_send[TDOWN], 0, !isOdd, 0, 0, 0, 0);
		bgq_weylfield_t_resetcoord(weylxchange_recv[TUP], LOCAL_LT, !isOdd, 2, 2, 0, 0);
		bgq_weylfield_t_resetcoord(weylxchange_recv[TDOWN], -1, !isOdd, 2, 2, 0, 0);

		bgq_weylfield_x_resetcoord(weylxchange_send[XUP], LOCAL_LX - 1, !isOdd, 0, 0, 0, 0);
		bgq_weylfield_x_resetcoord(weylxchange_send[XDOWN], 0, !isOdd, 0, 0, 0, 0);
		bgq_weylfield_x_resetcoord(weylxchange_recv[XUP], LOCAL_LX, !isOdd, 1, 1, 0, 0);
		bgq_weylfield_x_resetcoord(weylxchange_recv[XDOWN], -1, !isOdd, 1, 1, 0, 0);

		bgq_weylfield_y_resetcoord(weylxchange_send[YUP], LOCAL_LY - 1, !isOdd, 0, 0, 0, 0);
		bgq_weylfield_y_resetcoord(weylxchange_send[YDOWN], 0, !isOdd, 0, 0, 0, 0);
		bgq_weylfield_y_resetcoord(weylxchange_recv[YUP], LOCAL_LY, !isOdd, 1, 1, 0, 0);
		bgq_weylfield_y_resetcoord(weylxchange_recv[YDOWN], -1, !isOdd, 1, 1, 0, 0);
	}
#endif

// Just to be sure that if something follows this, it also holds up-to-data
// Only requires if the above is not within a #pragma omp parallel
//#pragma omp barrier
}


#ifdef BGQ_HM_NOFUNC
#undef BGQ_HM_ZLINE_NOFUNC
#undef BGQ_HM_SITE_NOFUNC
#undef BGQ_HM_DIR_NOFUNC
#endif

