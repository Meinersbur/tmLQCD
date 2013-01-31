/*
 * bgq_HoppingMatrixStencil.inc.c
 *
 *  Created on: Jan 23, 2013
 *      Author: meinersbur
 */

#define DIRECTION_LOWERCASE_TUP tup
#define DIRECTION_LOWERCASE_TDOWN tdown
#define DIRECTION_LOWERCASE_XUP xup
#define DIRECTION_LOWERCASE_XDOWN xdown
#define DIRECTION_LOWERCASE_YUP yup
#define DIRECTION_LOWERCASE_YDOWN ydown
#define DIRECTION_LOWERCASE_ZUP zup
#define DIRECTION_LOWERCASE_ZDOWN zdown
#define DIRECTION_LOWERCASE(DIR) NAME2(DIRECTION_LOWERCASE,DIR)

#ifndef BGQ_HOPPINGMATRIXSTENCIL_INC_
#include "bgq_qpx.h"
#include "bgq_spinorfield.h"
#include "bgq_gaugefield.h"
#include "bgq_comm.h"

#include <stdbool.h>

#define PRECISION double
//#define HASWEYL(DIRECTION) false
//#define GETSPINORFIELDPTR(DIR) (NAME2(pspinor,DIRECTION_LOWERCASE(DIR)))
//#define GETWEYLFIELDPTR(DIR) NULL

void bgq_HoppingMatrix_stencil_raw(bool gaugemul, ucoord ic_dst, bool isOdd_dst, bgq_weylfield_controlblock *spinorfield, bgq_spinor_vec *pspinor_tup, bgq_spinor_vec *pspinor_tdown, bgq_spinor_vec *pspinor_xup, bgq_spinor_vec *pspinor_xdown, bgq_spinor_vec *pspinor_yup, bgq_spinor_vec *pspinor_ydown, bgq_spinor_vec *pspinor_zup, bgq_spinor_vec *pspinor_zdown, bgq_gaugesite *gaugesite, bgq_spinor_vec *target, bgq_params(qka0), bgq_params(qka1), bgq_params(qka2), bgq_params(qka3), ucoord t1_dst, ucoord t2_dst, ucoord x_dst, ucoord y_dst, ucoord z_dst, bool kamul)
#endif
{
#ifndef NDEBUG
#endif

	bgq_su3_spinor_decl(result);

	// T+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_spinor_vec *spinorptr_tup = &spinorfield->BGQ_SEC_FULLLAYOUT[bgq_direction_move_collapsed(isOdd_dst, ic_dst, TUP)];
		bgq_su3_spinor_decl(spinor_tup);
		bgq_su3_spinor_load(spinor_tup, spinorptr_tup);
		bgq_su3_weyl_decl(weyl_tup);
		bgq_su3_reduce_weyl_tup(weyl_tup, spinor_tup);

		if (bgq_direction_isCrossingBorder_collapsed(isOdd_dst, ic_dst, TUP)) {
			if (COMM_T) {
						bgq_spinorqpxk_direxpect(spinor_tup, 1, t1_dst, x_dst, y_dst, z_dst, TUP);
						bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_TUP_SOURCE, bgq_cmplxval2(spinor_tup_v0_c0));

				// read data received from neighbor node
				size_t idxrecv_tup = g_bgq_collapsed2recvidx[isOdd_dst][ic_dst][TUP];
				bgq_weyl_nonvec *weylptr_tup = &g_bgq_sec_recv_unvectorized[TUP][idxrecv_tup];

				// Load weyl
				bgq_su3_weyl_decl(weyl_tmp);
				bgq_su3_weyl_load_nonvec(weyl_tmp, weylptr_tup);
#if 0
				bgq_vector4double_decl(weyl_v0c0_v0c1);
				bgq_vector4double_decl(weyl_v0c2_v1c0);
				bgq_vector4double_decl(weyl_v1c1_v1c2);
				void *tmpptr = weylptr_tup;
				bgq_lda(weyl_v0c0_v0c1, 0, tmpptr);
				bgq_qvlfuxa(weyl_v0c2_v1c0, tmpptr, PRECISION_SIZEOF * 4);
				bgq_qvlfuxa(weyl_v1c1_v1c2, tmpptr, PRECISION_SIZEOF * 4);
#endif

				//  merge it; use left part of spinor, throw away right part which has been sent to other node
				bgq_su3_weyl_merge2(weyl_tup, weyl_tup/*q2q3*/, weyl_tmp);
#if 0
				bgq_su3_weyl_decl(weyl_tmp);
				bgq_merge(weyl_tmp_v0_c0, weyl_v0c0_v0c1/*q0q1*/, weyl_tup_v0_c0/*q2q3*/);
				bgq_rmerge(weyl_tmp_v0_c1, weyl_v0c0_v0c1/*q2q3*/, weyl_tup_v0_c1/*q2q3*/);
				bgq_merge(weyl_tmp_v0_c2, weyl_v0c2_v1c0/*q0q1*/, weyl_tup_v0_c1/*q2q3*/);
				bgq_rmerge(weyl_tmp_v1_c0, weyl_v0c2_v1c0/*q2q3*/, weyl_tup_v1_c0/*q2q3*/);
				bgq_merge(weyl_tmp_v1_c1, weyl_v1c1_v1c2/*q0q1*/, weyl_tup_v1_c1/*q2q3*/);
				bgq_rmerge(weyl_tmp_v1_c2, weyl_v1c1_v1c2/*q2q3*/, weyl_tup_v1_c2/*q2q3*/);
#endif
						bgq_weylqpxk_expect(weyl_tmp, 1, t2_dst, x_dst, y_dst, z_dst, TUP, false);
			} else {
						bgq_spinorqpx_direxpect(spinor_tup, t2_dst, t1_dst, x_dst, y_dst, z_dst, TUP);
						bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_TUP_SOURCE, bgq_cmplxval2(spinor_tup_v0_c0));
						bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_TUP_SOURCE, bgq_cmplxval1(spinor_tup_v0_c0));

				// No communication necessary, but need to swap the vector
				bgq_su3_weyl_merge2(weyl_tup, weyl_tup, weyl_tup);
			}
		} else {
			bgq_spinorqpx_direxpect(spinor_tup, t1_dst, t2_dst, x_dst, y_dst, z_dst, TUP);
			bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_TUP_SOURCE, bgq_cmplxval1(spinor_tup_v0_c0));
			bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_TUP_SOURCE, bgq_cmplxval2(spinor_tup_v0_c0));
		}
				bgq_setdesc(BGQREF_TUP, "BGQREF_TUP");
				bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_TUP, bgq_cmplxval1(weyl_tup_v0_c0));
				bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_TUP, bgq_cmplxval2(weyl_tup_v0_c0));

		if (gaugemul) {
			bgq_su3_mdecl(gauge_tup);
			bgq_qvlfduxa(gauge_tup_c00, gaugesite, 32);
			bgq_qvlfduxa(gauge_tup_c01, gaugesite, 32);
			bgq_qvlfduxa(gauge_tup_c02, gaugesite, 32);
			bgq_qvlfduxa(gauge_tup_c10, gaugesite, 32);
			bgq_qvlfduxa(gauge_tup_c11, gaugesite, 32);
			bgq_qvlfduxa(gauge_tup_c12, gaugesite, 32);
			bgq_qvlfduxa(gauge_tup_c20, gaugesite, 32);
			bgq_qvlfduxa(gauge_tup_c21, gaugesite, 32);
			bgq_qvlfduxa(gauge_tup_c22, gaugesite, 32);
					bgq_gaugeqpx_expect(gauge_tup, t1_dst, t2_dst, x_dst, y_dst, z_dst, TUP, false);
					bgq_setdesc(BGQREF_TUP_GAUGE, "BGQREF_TUP_GAUGE");
					bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_TUP_GAUGE, bgq_cmplxval1(gauge_tup_c00));
					bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_TUP_GAUGE, bgq_cmplxval2(gauge_tup_c00));
			bgq_su3_weyl_mvmul(weyl_tup, gauge_tup, weyl_tup);
					bgq_setdesc(BGQREF_TUP_WEYL,"BGQREF_TUP_WEYL");
					bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_TUP_WEYL, bgq_cmplxval1(weyl_tup_v0_c0));
					bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_TUP_WEYL, bgq_cmplxval2(weyl_tup_v0_c0));
		}
		if (kamul) {
			bgq_su3_weyl_cmul(weyl_tup, qka0, weyl_tup);
		}
				bgq_setdesc(BGQREF_TUP_KAMUL,"BGQREF_TUP_KAMUL");
				bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_TUP_KAMUL, bgq_cmplxval1(weyl_tup_v0_c0));
				bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_TUP_KAMUL, bgq_cmplxval2(weyl_tup_v0_c0));

		bgq_su3_expand_weyl_tup(result, weyl_tup);
				bgq_setdesc(BGQREF_TUP_ACCUM, "BGQREF_TUP_ACCUM");
				bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_TUP_ACCUM, bgq_cmplxval1(result_v0_c0));
				bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_TUP_ACCUM, bgq_cmplxval2(result_v0_c0));
	}


	// T- /////////////////////////////////////////////////////////////////////////
	{
		bgq_spinor_vec *spinorptr_tdown = &spinorfield->BGQ_SEC_FULLLAYOUT[bgq_direction_move_collapsed(isOdd_dst, ic_dst, TDOWN)];
		bgq_su3_spinor_decl(spinor_tdown);
		bgq_su3_spinor_load(spinor_tdown, spinorptr_tdown);
		bgq_su3_weyl_decl(weyl_tdown);
		bgq_su3_reduce_weyl_tdown(weyl_tdown, spinor_tdown);

		if (bgq_direction_isCrossingBorder_collapsed(isOdd_dst, ic_dst, TDOWN)) {
			if (COMM_T) {
						bgq_spinorqpxk_direxpect(spinor_tdown, 0, t2_dst, x_dst, y_dst, z_dst, TDOWN);
						bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_TDOWN_SOURCE, bgq_cmplxval1(spinor_tdown_v0_c0));
				size_t idxrecv_tdown = g_bgq_collapsed2recvidx[isOdd_dst][ic_dst][TDOWN];
				bgq_weyl_nonvec *weylptr_tdown = &g_bgq_sec_recv_unvectorized[TDOWN][idxrecv_tdown];
				bgq_su3_weyl_decl(weyl_tmp);
				bgq_su3_weyl_load_nonvec(weyl_tmp, weylptr_tdown);

				bgq_su3_weyl_merge2(weyl_tdown, weyl_tmp, weyl_tdown);
			} else {
						bgq_spinorqpx_direxpect(spinor_tdown, t2_dst, t1_dst, x_dst, y_dst, z_dst, TDOWN);
						bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_TDOWN_SOURCE, bgq_cmplxval2(spinor_tdown_v0_c0));
						bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_TDOWN_SOURCE, bgq_cmplxval1(spinor_tdown_v0_c0));
				bgq_su3_weyl_merge2(weyl_tdown, weyl_tdown, weyl_tdown);
			}
		} else {
			bgq_spinorqpx_direxpect(spinor_tdown, t1_dst, t2_dst, x_dst, y_dst, z_dst, TDOWN);
			bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_TDOWN_SOURCE, bgq_cmplxval1(spinor_tdown_v0_c0));
			bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_TDOWN_SOURCE, bgq_cmplxval2(spinor_tdown_v0_c0));
		}
				bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_TDOWN, bgq_cmplxval1(weyl_tdown_v0_c2));
				bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_TDOWN, bgq_cmplxval2(weyl_tdown_v0_c2));

		if (gaugemul) {
			bgq_su3_mdecl(gauge_tdown);
			bgq_qvlfduxa(gauge_tdown_c00, gaugesite, 32);
			bgq_qvlfduxa(gauge_tdown_c01, gaugesite, 32);
			bgq_qvlfduxa(gauge_tdown_c02, gaugesite, 32);
			bgq_qvlfduxa(gauge_tdown_c10, gaugesite, 32);
			bgq_qvlfduxa(gauge_tdown_c11, gaugesite, 32);
			bgq_qvlfduxa(gauge_tdown_c12, gaugesite, 32);
			bgq_qvlfduxa(gauge_tdown_c20, gaugesite, 32);
			bgq_qvlfduxa(gauge_tdown_c21, gaugesite, 32);
			bgq_qvlfduxa(gauge_tdown_c22, gaugesite, 32);
					bgq_gaugeqpx_expect(gauge_tdown, t1_dst, t2_dst, x_dst, y_dst, z_dst, TDOWN, false);
					bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_TDOWN_GAUGE, bgq_cmplxval1(gauge_tdown_c00));
					bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_TDOWN_GAUGE, bgq_cmplxval2(gauge_tdown_c00));
			bgq_su3_weyl_mvinvmul(weyl_tdown, gauge_tdown, weyl_tdown);
					bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_TDOWN_WEYL, bgq_cmplxval1(weyl_tdown_v0_c0));
					bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_TDOWN_WEYL, bgq_cmplxval2(weyl_tdown_v0_c0));
		}
		if (kamul) {
			bgq_su3_weyl_cmul(weyl_tdown, qka0, weyl_tdown);
		}
				bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_TDOWN_KAMUL, bgq_cmplxval1(weyl_tdown_v0_c0));
				bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_TDOWN_KAMUL, bgq_cmplxval2(weyl_tdown_v0_c0));

		bgq_su3_accum_weyl_tdown(result, weyl_tdown);
				bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_TDOWN_ACCUM, bgq_cmplxval1(result_v0_c0));
				bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_TDOWN_ACCUM, bgq_cmplxval2(result_v0_c0));
	}


	// X+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xup);
		if (COMM_X && bgq_direction_isCrossingBorder_collapsed(isOdd_dst, ic_dst, XUP)) {
			size_t idxrecv_xup = g_bgq_collapsed2recvidx[isOdd_dst][ic_dst][XUP];
			bgq_weyl_vec *weylptr_xup = &g_bgq_sec_recv[XUP][idxrecv_xup];
			bgq_su3_weyl_load(weyl_xup, weylptr_xup);
					bgq_weylqpx_expect(weyl_xup, t1_dst, t2_dst, x_dst, y_dst, z_dst, XUP, false);
		} else {
			bgq_spinor_vec *spinorptr_xup = &spinorfield->BGQ_SEC_FULLLAYOUT[bgq_direction_move_collapsed(isOdd_dst, ic_dst, XUP)];
			bgq_su3_spinor_decl(spinor_xup);
			bgq_su3_spinor_load(spinor_xup, spinorptr_xup);
					bgq_spinorqpx_direxpect(spinor_xup, t1_dst, t2_dst, x_dst, y_dst, z_dst, XUP);
			bgq_su3_reduce_weyl_xup(weyl_xup, spinor_xup);
		}

		if (gaugemul) {
			bgq_su3_mdecl(gauge_xup);
			bgq_qvlfduxa(gauge_xup_c00, gaugesite, 32);
			bgq_qvlfduxa(gauge_xup_c01, gaugesite, 32);
			bgq_qvlfduxa(gauge_xup_c02, gaugesite, 32);
			bgq_qvlfduxa(gauge_xup_c10, gaugesite, 32);
			bgq_qvlfduxa(gauge_xup_c11, gaugesite, 32);
			bgq_qvlfduxa(gauge_xup_c12, gaugesite, 32);
			bgq_qvlfduxa(gauge_xup_c20, gaugesite, 32);
			bgq_qvlfduxa(gauge_xup_c21, gaugesite, 32);
			bgq_qvlfduxa(gauge_xup_c22, gaugesite, 32);
					bgq_gaugeqpx_expect(gauge_xup, t1_dst, t2_dst, x_dst, y_dst, z_dst, XUP, false);
			bgq_su3_weyl_mvmul(weyl_xup, gauge_xup, weyl_xup);
		}
		if (kamul) {
			bgq_su3_weyl_cmul(weyl_xup, qka1, weyl_xup);
		}
				bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_XUP_KAMUL, bgq_cmplxval1(weyl_xup_v0_c0));
				bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_XUP_KAMUL, bgq_cmplxval2(weyl_xup_v0_c0));

		bgq_su3_accum_weyl_xup(result, weyl_xup);
				bgq_setdesc(BGQREF_XUP_ACCUM, "BGQREF_XUP_ACCUM");
				bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_XUP_ACCUM, bgq_cmplxval1(result_v0_c0));
				bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_XUP_ACCUM, bgq_cmplxval2(result_v0_c0));
	}


	// X- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_xdown);
		if (COMM_X && bgq_direction_isCrossingBorder_collapsed(isOdd_dst, ic_dst, XDOWN)) {
			size_t idxrecv_xdown = g_bgq_collapsed2recvidx[isOdd_dst][ic_dst][XDOWN];
			bgq_weyl_vec *weylptr_xdown = &g_bgq_sec_recv[XDOWN][idxrecv_xdown];
			bgq_su3_weyl_load(weyl_xdown, weylptr_xdown);
					bgq_weylqpx_expect(weyl_xdown, t1_dst, t2_dst, x_dst, y_dst, z_dst, XDOWN, false);
		} else {
			bgq_spinor_vec *spinorptr_xdown = &spinorfield->BGQ_SEC_FULLLAYOUT[bgq_direction_move_collapsed(isOdd_dst, ic_dst, XDOWN)];
			bgq_su3_spinor_decl(spinor_xdown);
			bgq_su3_spinor_load(spinor_xdown, spinorptr_xdown);
					bgq_spinorqpx_direxpect(spinor_xdown, t1_dst, t2_dst, x_dst, y_dst, z_dst, XDOWN);
			bgq_su3_reduce_weyl_xdown(weyl_xdown, spinor_xdown);
		}

		if (gaugemul) {
			bgq_su3_mdecl(gauge_xdown);
			bgq_qvlfduxa(gauge_xdown_c00, gaugesite, 32);
			bgq_qvlfduxa(gauge_xdown_c01, gaugesite, 32);
			bgq_qvlfduxa(gauge_xdown_c02, gaugesite, 32);
			bgq_qvlfduxa(gauge_xdown_c10, gaugesite, 32);
			bgq_qvlfduxa(gauge_xdown_c11, gaugesite, 32);
			bgq_qvlfduxa(gauge_xdown_c12, gaugesite, 32);
			bgq_qvlfduxa(gauge_xdown_c20, gaugesite, 32);
			bgq_qvlfduxa(gauge_xdown_c21, gaugesite, 32);
			bgq_qvlfduxa(gauge_xdown_c22, gaugesite, 32);
					bgq_gaugeqpx_expect(gauge_xdown, t1_dst, t2_dst, x_dst, y_dst, z_dst, XDOWN, false);
			bgq_su3_weyl_mvinvmul(weyl_xdown, gauge_xdown, weyl_xdown);
		}
		if (kamul) {
			bgq_su3_weyl_cmul(weyl_xdown, qka1, weyl_xdown);
		}

		bgq_su3_accum_weyl_xdown(result, weyl_xdown);
				bgq_setdesc(BGQREF_XDOWN_ACCUM, "BGQREF_XDOWN_ACCUM");
				bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_XDOWN_ACCUM, bgq_cmplxval1(result_v0_c0));
				bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_XDOWN_ACCUM, bgq_cmplxval2(result_v0_c0));
	}


	// Y+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_yup);
		if (COMM_Y && bgq_direction_isCrossingBorder_collapsed(isOdd_dst, ic_dst, YUP)) {
			size_t idxrecv_yup = g_bgq_collapsed2recvidx[isOdd_dst][ic_dst][YUP];
			bgq_weyl_vec *weylptr_yup = &g_bgq_sec_recv[YUP][idxrecv_yup];
			bgq_su3_weyl_load(weyl_yup, weylptr_yup);
					bgq_weylqpx_expect(weyl_yup, t1_dst, t2_dst, x_dst, y_dst, z_dst, YUP, false);
		} else {
			bgq_spinor_vec *spinorptr_yup = &spinorfield->BGQ_SEC_FULLLAYOUT[bgq_direction_move_collapsed(isOdd_dst, ic_dst, YUP)];
			bgq_su3_spinor_decl(spinor_yup);
			bgq_su3_spinor_load(spinor_yup, spinorptr_yup);
					bgq_spinorqpx_direxpect(spinor_yup, t1_dst, t2_dst, x_dst, y_dst, z_dst, YUP);
			bgq_su3_reduce_weyl_yup(weyl_yup, spinor_yup);
		}

		if (gaugemul) {
				bgq_su3_mdecl(gauge_yup);
				bgq_qvlfduxa(gauge_yup_c00, gaugesite, 32);
				bgq_qvlfduxa(gauge_yup_c01, gaugesite, 32);
				bgq_qvlfduxa(gauge_yup_c02, gaugesite, 32);
				bgq_qvlfduxa(gauge_yup_c10, gaugesite, 32);
				bgq_qvlfduxa(gauge_yup_c11, gaugesite, 32);
				bgq_qvlfduxa(gauge_yup_c12, gaugesite, 32);
				bgq_qvlfduxa(gauge_yup_c20, gaugesite, 32);
				bgq_qvlfduxa(gauge_yup_c21, gaugesite, 32);
				bgq_qvlfduxa(gauge_yup_c22, gaugesite, 32);
						bgq_gaugeqpx_expect(gauge_yup, t1_dst, t2_dst, x_dst, y_dst, z_dst, YUP, false);
				bgq_su3_weyl_mvmul(weyl_yup, gauge_yup, weyl_yup);
		}
		if (kamul) {
			bgq_su3_weyl_cmul(weyl_yup, qka2, weyl_yup);
		}

		bgq_su3_accum_weyl_yup(result, weyl_yup);
				bgq_setdesc(BGQREF_YUP_ACCUM, "BGQREF_YUP_ACCUM");
				bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_YUP_ACCUM, bgq_cmplxval1(result_v0_c0));
				bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_YUP_ACCUM, bgq_cmplxval2(result_v0_c0));
	}


	// Y- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_ydown);
		if (COMM_Y && bgq_direction_isCrossingBorder_collapsed(isOdd_dst, ic_dst, YDOWN)) {
			size_t idxrecv_ydown = g_bgq_collapsed2recvidx[isOdd_dst][ic_dst][YDOWN];
			bgq_weyl_vec *weylptr_ydown = &g_bgq_sec_recv[YDOWN][idxrecv_ydown];
			bgq_su3_weyl_load(weyl_ydown, weylptr_ydown);
					bgq_weylqpx_expect(weyl_ydown, t1_dst, t2_dst, x_dst, y_dst, z_dst, YDOWN, false);
		} else {
			bgq_spinor_vec *spinorptr_ydown = &spinorfield->BGQ_SEC_FULLLAYOUT[bgq_direction_move_collapsed(isOdd_dst, ic_dst, YDOWN)];
			bgq_su3_spinor_decl(spinor_ydown);
			bgq_su3_spinor_load(spinor_ydown, spinorptr_ydown);
					bgq_spinorqpx_direxpect(spinor_ydown, t1_dst, t2_dst, x_dst, y_dst, z_dst, YDOWN);
			bgq_su3_reduce_weyl_ydown(weyl_ydown, spinor_ydown);
		}

		if (gaugemul) {
				bgq_su3_mdecl(gauge_ydown);
				bgq_qvlfduxa(gauge_ydown_c00, gaugesite, 32);
				bgq_qvlfduxa(gauge_ydown_c01, gaugesite, 32);
				bgq_qvlfduxa(gauge_ydown_c02, gaugesite, 32);
				bgq_qvlfduxa(gauge_ydown_c10, gaugesite, 32);
				bgq_qvlfduxa(gauge_ydown_c11, gaugesite, 32);
				bgq_qvlfduxa(gauge_ydown_c12, gaugesite, 32);
				bgq_qvlfduxa(gauge_ydown_c20, gaugesite, 32);
				bgq_qvlfduxa(gauge_ydown_c21, gaugesite, 32);
				bgq_qvlfduxa(gauge_ydown_c22, gaugesite, 32);
						bgq_gaugeqpx_expect(gauge_ydown, t1_dst, t2_dst, x_dst, y_dst, z_dst, YDOWN, false);
				bgq_su3_weyl_mvinvmul(weyl_ydown, gauge_ydown, weyl_ydown);
		}
		if (kamul) {
			bgq_su3_weyl_cmul(weyl_ydown, qka2, weyl_ydown);
		}

		bgq_su3_accum_weyl_ydown(result, weyl_ydown);
				bgq_setdesc(BGQREF_YDOWN_ACCUM, "BGQREF_YDOWN_ACCUM");
				bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_YDOWN_ACCUM, bgq_cmplxval1(result_v0_c0));
				bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_YDOWN_ACCUM, bgq_cmplxval2(result_v0_c0));
	}


	// Z+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zup);
		if (COMM_Z && bgq_direction_isCrossingBorder_collapsed(isOdd_dst, ic_dst, ZUP)) {
			size_t idxrecv_zup = g_bgq_collapsed2recvidx[isOdd_dst][ic_dst][ZUP];
			bgq_weyl_vec *weylptr_zup = &g_bgq_sec_recv[ZUP][idxrecv_zup];
			bgq_su3_weyl_load(weyl_zup, weylptr_zup);
					bgq_weylqpx_expect(weyl_zup, t1_dst, t2_dst, x_dst, y_dst, z_dst, ZUP, false);
		} else {
			bgq_spinor_vec *spinorptr_zup = &spinorfield->BGQ_SEC_FULLLAYOUT[bgq_direction_move_collapsed(isOdd_dst, ic_dst, ZUP)];
			bgq_su3_spinor_decl(spinor_zup);
			bgq_su3_spinor_load(spinor_zup, spinorptr_zup);
					bgq_spinorqpx_direxpect(spinor_zup, t1_dst, t2_dst, x_dst, y_dst, z_dst, ZUP);
			bgq_su3_reduce_weyl_zup(weyl_zup, spinor_zup);
		}

		if (gaugemul) {
				bgq_su3_mdecl(gauge_zup);
				bgq_qvlfduxa(gauge_zup_c00, gaugesite, 32);
				bgq_qvlfduxa(gauge_zup_c01, gaugesite, 32);
				bgq_qvlfduxa(gauge_zup_c02, gaugesite, 32);
				bgq_qvlfduxa(gauge_zup_c10, gaugesite, 32);
				bgq_qvlfduxa(gauge_zup_c11, gaugesite, 32);
				bgq_qvlfduxa(gauge_zup_c12, gaugesite, 32);
				bgq_qvlfduxa(gauge_zup_c20, gaugesite, 32);
				bgq_qvlfduxa(gauge_zup_c21, gaugesite, 32);
				bgq_qvlfduxa(gauge_zup_c22, gaugesite, 32);
						bgq_gaugeqpx_expect(gauge_zup, t1_dst, t2_dst, x_dst, y_dst, z_dst, ZUP, false);
				bgq_su3_weyl_mvmul(weyl_zup, gauge_zup, weyl_zup);
		}
		if (kamul) {
			bgq_su3_weyl_cmul(weyl_zup, qka1, weyl_zup);
		}

		bgq_su3_accum_weyl_zup(result, weyl_zup);
				bgq_setdesc(BGQREF_ZUP_ACCUM, "BGQREF_ZUP_ACCUM");
				bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_ZUP_ACCUM, bgq_cmplxval1(result_v0_c0));
				bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_ZUP_ACCUM, bgq_cmplxval2(result_v0_c0));
	}


	// Z- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weyl_zdown);
		if (COMM_Z && bgq_direction_isCrossingBorder_collapsed(isOdd_dst, ic_dst, ZDOWN)) {
			size_t idxrecv_zdown = g_bgq_collapsed2recvidx[isOdd_dst][ic_dst][ZDOWN];
			bgq_weyl_vec *weylptr_zdown = &g_bgq_sec_recv[ZDOWN][idxrecv_zdown];
			bgq_su3_weyl_load(weyl_zdown, weylptr_zdown);
					bgq_weylqpx_expect(weyl_zdown, t1_dst, t2_dst, x_dst, y_dst, z_dst, ZDOWN, false);
		} else {
			bgq_spinor_vec *spinorptr_zdown = &spinorfield->BGQ_SEC_FULLLAYOUT[bgq_direction_move_collapsed(isOdd_dst, ic_dst, ZDOWN)];
			bgq_su3_spinor_decl(spinor_zdown);
			bgq_su3_spinor_load(spinor_zdown, spinorptr_zdown);
					bgq_spinorqpx_direxpect(spinor_zdown, t1_dst, t2_dst, x_dst, y_dst, z_dst, ZDOWN);
			bgq_su3_reduce_weyl_zdown(weyl_zdown, spinor_zdown);
		}

		if (gaugemul) {
			bgq_su3_mdecl(gauge_zdown);
			bgq_qvlfduxa(gauge_zdown_c00, gaugesite, 32);
			bgq_qvlfduxa(gauge_zdown_c01, gaugesite, 32);
			bgq_qvlfduxa(gauge_zdown_c02, gaugesite, 32);
			bgq_qvlfduxa(gauge_zdown_c10, gaugesite, 32);
			bgq_qvlfduxa(gauge_zdown_c11, gaugesite, 32);
			bgq_qvlfduxa(gauge_zdown_c12, gaugesite, 32);
			bgq_qvlfduxa(gauge_zdown_c20, gaugesite, 32);
			bgq_qvlfduxa(gauge_zdown_c21, gaugesite, 32);
			bgq_qvlfduxa(gauge_zdown_c22, gaugesite, 32);
					bgq_gaugeqpx_expect(gauge_zdown, t1_dst, t2_dst, x_dst, y_dst, z_dst, ZDOWN, false);
			bgq_su3_weyl_mvinvmul(weyl_zdown, gauge_zdown, weyl_zdown);
		}
		if (kamul) {
			bgq_su3_weyl_cmul(weyl_zdown, qka1, weyl_zdown);
		}

		bgq_su3_accum_weyl_zdown(result, weyl_zdown);
				bgq_setdesc(BGQREF_ACCUM, "BGQREF_ACCUM");
				bgq_setbgqvalue(t1_dst, x_dst, y_dst, z_dst, BGQREF_ACCUM, bgq_cmplxval1(result_v0_c0));
				bgq_setbgqvalue(t2_dst, x_dst, y_dst, z_dst, BGQREF_ACCUM, bgq_cmplxval2(result_v0_c0));
	}

			bgq_spinorqpx_written(result, t1_dst, t2_dst, x_dst, y_dst, z_dst);
	bgq_su3_spinor_store(target, result);
}
#undef BGQ_HOPPINGMATRIXSTENCIL_INC_
