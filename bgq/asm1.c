
#include "bgq_qpx.h"
#include "bgq_HoppingMatrix.h"
#include "bgq_comm.h"




void *somewhere;


bgq_spinor bgq_spinorfield_getspinor(bgq_weylfield_controlblock *field, ucoord t, ucoord x, ucoord y, ucoord z) {
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);

	assert(bgq_local2isOdd(t,x,y,z)==field->isOdd);
	//ucoord ic = bgq_local2collapsed(t,x,y,z);
	ucoord k = bgq_local2k(t,x,y,z);

		bgq_spinorfield_setup(field, field->isOdd, false, false, true, false);
		bgq_su3_spinor_decl(spinor);
		//bgq_HoppingMatrix_loadWeyllayout(spinor,&field->sec_collapsed[ic], bgq_t2t(t,0), bgq_t2t(t,1), x, y, z);

	bgq_weylsite *spinorsite = &field->sec_collapsed[0];

	for (ucoord ic = 0; ic < t; ic+=1) {
		ucoord t1 = bgq_t2t(t,0);
		ucoord t2 = bgq_t2t(t,1);

		 {
			bgq_su3_spinor_decl(result);
			bgq_su3_spinor_zero(spinor);

			asm volatile ("" : : : "memory");

			bgq_su3_weyl_prefetch_double(&spinorsite->d[TDOWN]);
			asm volatile ("" : : : "memory");
			// T+ //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_tup);
				bgq_su3_weyl_load_double(weyl_tup, &spinorsite->d[TUP]);
				asm volatile ("" : : : "memory");
				bgq_weylqpx_expect(weyl_tup, t1, t2, x, y, z, TUP, false);
						bgq_setdesc(BGQREF_TUP_RECV, "BGQREF_TUP_RECV");
						bgq_setbgqvalue(t1, x, y, z, BGQREF_TUP_RECV, bgq_cmplxval1(weyl_tup_v1_c0));
						bgq_setbgqvalue(t2, x, y, z, BGQREF_TUP_RECV, bgq_cmplxval2(weyl_tup_v1_c0));
				bgq_su3_expand_weyl_tup(result, weyl_tup);
			}
			asm volatile ("" : : : "memory");
			bgq_su3_weyl_prefetch_double(&spinorsite->d[XUP]);
			asm volatile ("" : : : "memory");
			// T- //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_tdown);
				bgq_su3_weyl_load_double(weyl_tdown, &spinorsite->d[TDOWN]);
				asm volatile ("" : : : "memory");
				bgq_weylqpx_expect(weyl_tdown, t1, t2, x, y, z, TDOWN, false);
						bgq_setdesc(BGQREF_TDOWN_RECV, "BGQREF_TDOWN_RECV");
						bgq_setbgqvalue(t1, x, y, z, BGQREF_TDOWN_RECV, bgq_cmplxval1(weyl_tdown_v1_c0));
						bgq_setbgqvalue(t2, x, y, z, BGQREF_TDOWN_RECV, bgq_cmplxval2(weyl_tdown_v1_c0));
				bgq_su3_accum_weyl_tdown(result, weyl_tdown);
						bgq_setdesc(BGQREF_TDOWN_ACCUM, "BGQREF_TDOWN_ACCUM");
						bgq_setbgqvalue(t1, x, y, z, BGQREF_TDOWN_ACCUM, bgq_cmplxval1(result_v1_c0));
						bgq_setbgqvalue(t2, x, y, z, BGQREF_TDOWN_ACCUM, bgq_cmplxval2(result_v1_c0));
			}

			asm volatile ("" : : : "memory");
			bgq_su3_weyl_prefetch_double(&spinorsite->d[XDOWN]);
			asm volatile ("" : : : "memory");

			// X+ //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_xup);
				bgq_su3_weyl_load_double(weyl_xup, &spinorsite->d[XUP]);
				asm volatile ("" : : : "memory");
				bgq_weylqpx_expect(weyl_xup, t1, t2, x, y, z, XUP, false);
					bgq_setdesc(BGQREF_XUP_RECV, "BGQREF_XUP_RECV");
					bgq_setbgqvalue(t1, x, y, z, BGQREF_XUP_RECV, bgq_cmplxval1(weyl_xup_v1_c0));
					bgq_setbgqvalue(t2, x, y, z, BGQREF_XUP_RECV, bgq_cmplxval2(weyl_xup_v1_c0));
				bgq_su3_accum_weyl_xup(result, weyl_xup);
					bgq_setdesc(BGQREF_XUP_ACCUM, "BGQREF_XUP_ACCUM");
					bgq_setbgqvalue(t1, x, y, z, BGQREF_XUP_ACCUM, bgq_cmplxval1(result_v1_c0));
					bgq_setbgqvalue(t2, x, y, z, BGQREF_XUP_ACCUM, bgq_cmplxval2(result_v1_c0));
			}
			asm volatile ("" : : : "memory");
			bgq_su3_weyl_prefetch_double(&spinorsite->d[YUP]);
			asm volatile ("" : : : "memory");

#if 0

			// X- //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_xdown);
				bgq_su3_weyl_load_double(weyl_xdown, &spinorsite->d[XDOWN]);
				asm volatile ("" : : : "memory");
				bgq_weylqpx_expect(weyl_xdown, t1, t2, x, y, z, XDOWN, false);
					bgq_setdesc(BGQREF_XDOWN_RECV, "BGQREF_XDOWN_RECV");
					bgq_setbgqvalue(t1, x, y, z, BGQREF_XDOWN_RECV, bgq_cmplxval1(weyl_xdown_v1_c0));
					bgq_setbgqvalue(t2, x, y, z, BGQREF_XDOWN_RECV, bgq_cmplxval2(weyl_xdown_v1_c0));
				bgq_su3_accum_weyl_xdown(result, weyl_xdown);
					bgq_setdesc(BGQREF_XDOWN_ACCUM, "BGQREF_XDOWN_ACCUM");
					bgq_setbgqvalue(t1, x, y, z, BGQREF_XDOWN_ACCUM, bgq_cmplxval1(result_v1_c0));
					bgq_setbgqvalue(t2, x, y, z, BGQREF_XDOWN_ACCUM, bgq_cmplxval2(result_v1_c0));
			}
			asm volatile ("" : : : "memory");
			bgq_su3_weyl_prefetch_double(&spinorsite->d[YDOWN]);
			asm volatile ("" : : : "memory");
			// Y+ //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_yup);
				bgq_su3_weyl_load_double(weyl_yup, &spinorsite->d[YUP]);
				asm volatile ("" : : : "memory");
				bgq_weylqpx_expect(weyl_yup, t1, t2, x, y, z, YUP,false);
					bgq_setdesc(BGQREF_YUP_RECV, "BGQREF_YUP_RECV");
					bgq_setbgqvalue(t1, x, y, z, BGQREF_YUP_RECV, bgq_cmplxval1(weyl_yup_v1_c0));
					bgq_setbgqvalue(t2, x, y, z, BGQREF_YUP_RECV, bgq_cmplxval2(weyl_yup_v1_c0));
				bgq_su3_accum_weyl_yup(result, weyl_yup);
					bgq_setdesc(BGQREF_YUP_ACCUM, "BGQREF_YUP_ACCUM");
					bgq_setbgqvalue(t1, x, y, z, BGQREF_YUP_ACCUM, bgq_cmplxval1(result_v1_c0));
					bgq_setbgqvalue(t2, x, y, z, BGQREF_YUP_ACCUM, bgq_cmplxval2(result_v1_c0));
			}
			asm volatile ("" : : : "memory");
			bgq_su3_weyl_prefetch_double(&spinorsite->d[ZUP]);
			asm volatile ("" : : : "memory");
			// Y- //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_ydown);
				bgq_su3_weyl_load_double(weyl_ydown, &spinorsite->d[YDOWN]);
				asm volatile ("" : : : "memory");
				bgq_weylqpx_expect(weyl_ydown, t1, t2, x, y, z, YDOWN, false);
				bgq_su3_accum_weyl_ydown(result, weyl_ydown);
					bgq_setdesc(BGQREF_YDOWN_ACCUM, "BGQREF_YDOWN_ACCUM");
					bgq_setbgqvalue(t1, x, y, z, BGQREF_YDOWN_ACCUM, bgq_cmplxval1(result_v1_c0));
					bgq_setbgqvalue(t2, x, y, z, BGQREF_YDOWN_ACCUM, bgq_cmplxval2(result_v1_c0));
			}
			asm volatile ("" : : : "memory");
			bgq_su3_weyl_prefetch_double(&spinorsite->d[ZDOWN]);
			asm volatile ("" : : : "memory");
			// Z+ //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_zup);
				bgq_su3_weyl_load_double(weyl_zup, &spinorsite->d[ZUP]);
				asm volatile ("" : : : "memory");
				bgq_weylqpx_expect(weyl_zup, t1, t2, x, y, z, ZUP, false);
					bgq_setdesc(BGQREF_ZUP_RECV, "BGQREF_ZUP_RECV");
					bgq_setbgqvalue(t1, x, y, z, BGQREF_ZUP_RECV, bgq_cmplxval1(weyl_zup_v1_c0));
					bgq_setbgqvalue(t2, x, y, z, BGQREF_ZUP_RECV, bgq_cmplxval2(weyl_zup_v1_c0));
				bgq_su3_accum_weyl_zup(result, weyl_zup);
					bgq_setdesc(BGQREF_ZUP_ACCUM, "BGQREF_ZUP_ACCUM");
					bgq_setbgqvalue(t1, x, y, z, BGQREF_ZUP_ACCUM, bgq_cmplxval1(result_v1_c0));
					bgq_setbgqvalue(t2, x, y, z, BGQREF_ZUP_ACCUM, bgq_cmplxval2(result_v1_c0));
			}
			asm volatile ("" : : : "memory");
			bgq_su3_weyl_prefetch_double(&spinorsite->d[TUP]);
			asm volatile ("" : : : "memory");
			// Z- //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_zdown);
				bgq_su3_weyl_load_double(weyl_zdown, &spinorsite->d[ZDOWN]);
				asm volatile ("" : : : "memory");
				bgq_weylqpx_expect(weyl_zdown, t1, t2, x, y, z, ZDOWN,false);
				bgq_su3_accum_weyl_zdown(result, weyl_zdown);
			}


					bgq_setdesc(BGQREF_ACCUM, "BGQREF_ACCUM");
					bgq_setbgqvalue(t1, x, y, z, BGQREF_ACCUM, bgq_cmplxval1(result_v1_c0));
					bgq_setbgqvalue(t2, x, y, z, BGQREF_ACCUM, bgq_cmplxval2(result_v1_c0));
					asm volatile ("" : : : "memory");
#endif
			bgq_su3_spinor_mov(spinor, result);
			asm volatile ("" : : : "memory");
		}


		//return bgq_spinor_fromqpx(spinor,k);
		 bgq_su3_spinor_store_double(somewhere,spinor);
	}

}
