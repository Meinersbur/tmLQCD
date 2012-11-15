/*
 * IBM xlc refuses to inline functions that big, therefore we have to include the function using the preprocessor
 */

#if 0
					bgq_weyl_vec *targetsite = targetptrs->d[TDOWN];
					asm (
						QVLFDUXA(GAUGE_C00,GAUGESITE,C32)
						QVLFDUXA(GAUGE_C01,GAUGESITE,C32)
						QVLFDUXA(GAUGE_C02,GAUGESITE,C32)
						QVLFDUXA(GAUGE_C10,GAUGESITE,C32)
						QVLFDUXA(GAUGE_C11,GAUGESITE,C32)
						QVLFDUXA(GAUGE_C12,GAUGESITE,C32)
						QVLFDUXA(GAUGE_C20,GAUGESITE,C32)
						QVLFDUXA(GAUGE_C21,GAUGESITE,C32)
						QVLFDUXA(GAUGE_C22,GAUGESITE,C32)

						QVFSUB(WEYL_V0_C0,SPINOR_V0_C0,SPINOR_V2_C0)
						QVFSUB(WEYL_V0_C1,SPINOR_V0_C1,SPINOR_V2_C1)
						QVFSUB(WEYL_V0_C2,SPINOR_V0_C2,SPINOR_V2_C2)

						QVFXMUL(HALF_V0_C0,GAUGE_C00,WEYL_V0_C0)
						QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C0,GAUGE_C00,HALF_V0_C0)
						QVFXMADD(HALF_V0_C0,GAUGE_C10,WEYL_V0_C1,HALF_V0_C0)
						QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C1,GAUGE_C10,HALF_V0_C0)
						QVFXMADD(HALF_V0_C0,GAUGE_C20,WEYL_V0_C2,HALF_V0_C0)
						QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C2,GAUGE_C20,HALF_V0_C0)

						QVFXMUL(HALF_V0_C1,GAUGE_C01,WEYL_V0_C0)
						QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C0,GAUGE_C01,HALF_V0_C1)
						QVFXMADD(HALF_V0_C1,GAUGE_C11,WEYL_V0_C1,HALF_V0_C1)
						QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C1,GAUGE_C11,HALF_V0_C1)
						QVFXMADD(HALF_V0_C1,GAUGE_C21,WEYL_V0_C2,HALF_V0_C1)
						QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C2,GAUGE_C21,HALF_V0_C1)

						QVFXMUL(HALF_V0_C2,GAUGE_C02,WEYL_V0_C0)
						QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C0,GAUGE_C02,HALF_V0_C2) /* end lifetime: WEYL_V0_C0 */
						QVFXMADD(HALF_V0_C2,GAUGE_C12,WEYL_V0_C1,HALF_V0_C2)
						QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C1,GAUGE_C12,HALF_V0_C2) /* end lifetime: WEYL_V0_C1 */
						QVFXMADD(HALF_V0_C2,GAUGE_C22,WEYL_V0_C2,HALF_V0_C2)
						QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C2,GAUGE_C22,HALF_V0_C2) /* end lifetime: WEYL_V0_C2 */

			#if KAMUL
						QVFXMUL    (WEYL_V0_C0, QKA0      , HALF_V0_C0)
						QVFXXNPMADD(HALF_V0_C0, HALF_V0_C0, QKA0      , WEYL_V0_C0)
						QVFXMUL    (WEYL_V0_C1, QKA0      , HALF_V0_C1)
						QVFXXNPMADD(HALF_V0_C1, HALF_V0_C1, QKA0      , WEYL_V0_C1)
						QVFXMUL    (WEYL_V0_C2, QKA0      , HALF_V0_C2)
						QVFXXNPMADD(HALF_V0_C2, HALF_V0_C2, QKA0      , WEYL_V0_C2)
			#endif

						//DCBZ_0(TARGETSITE)
						//DCBZ(C64,TARGETSITE)
						//DCBZ(C128,TARGETSITE)
						QVSTFDXA_0(HALF_V0_C0,TARGETSITE) /* end lifetime: HALF_V0_C0 */
						QVSTFDUXA(HALF_V0_C1,TARGETSITE,C32) /* end lifetime: HALF_V0_C1 */
						QVSTFDUXA(HALF_V0_C2,TARGETSITE,C32) /* end lifetime: HALF_V0_C2 */

						QVFSUB(WEYL_V1_C0,SPINOR_V1_C0,SPINOR_V3_C0)
						QVFSUB(WEYL_V1_C1,SPINOR_V1_C1,SPINOR_V3_C1)
						QVFSUB(WEYL_V1_C2,SPINOR_V1_C2,SPINOR_V3_C2)

						QVFXMUL(HALF_V1_C0,GAUGE_C00,WEYL_V1_C0)
						QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C0,GAUGE_C00,HALF_V1_C0) /* end lifetime: GAUGE_C00 */
						QVFXMADD(HALF_V1_C0,GAUGE_C10,WEYL_V1_C1,HALF_V1_C0)
						QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C1,GAUGE_C10,HALF_V1_C0) /* end lifetime: GAUGE_C10 */
						QVFXMADD(HALF_V1_C0,GAUGE_C20,WEYL_V1_C2,HALF_V1_C0)
						QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C2,GAUGE_C20,HALF_V1_C0) /* end lifetime: GAUGE_C20 */

						QVFXMUL(HALF_V1_C1,GAUGE_C01,WEYL_V1_C0)
						QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C0,GAUGE_C01,HALF_V1_C1) /* end lifetime: GAUGE_C01 */
						QVFXMADD(HALF_V1_C1,GAUGE_C11,WEYL_V1_C1,HALF_V1_C1)
						QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C1,GAUGE_C11,HALF_V1_C1) /* end lifetime: GAUGE_C11 */
						QVFXMADD(HALF_V1_C1,GAUGE_C21,WEYL_V1_C2,HALF_V1_C1)
						QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C2,GAUGE_C21,HALF_V1_C1) /* end lifetime: GAUGE_C21 */

						QVFXMUL(HALF_V1_C2,GAUGE_C02,WEYL_V1_C0)
						QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C0,GAUGE_C02,HALF_V1_C2) /* end lifetime: GAUGE_C02, WEYL_V1_C0 */
						QVFXMADD(HALF_V1_C2,GAUGE_C12,WEYL_V1_C1,HALF_V1_C2)
						QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C1,GAUGE_C12,HALF_V1_C2) /* end lifetime: GAUGE_C12, WEYL_V1_C1 */
						QVFXMADD(HALF_V1_C2,GAUGE_C22,WEYL_V1_C2,HALF_V1_C2)
						QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C2,GAUGE_C22,HALF_V1_C2) /* end lifetime: GAUGE_C22, WEYL_V1_C2 */

			#if KAMUL
						QVFXMUL    (WEYL_V1_C0, QKA0      , HALF_V1_C0)
						QVFXXNPMADD(HALF_V1_C0, HALF_V1_C0, QKA0      , WEYL_V1_C0)
						QVFXMUL    (WEYL_V1_C1, QKA0      , HALF_V1_C1)
						QVFXXNPMADD(HALF_V1_C1, HALF_V1_C1, QKA0      , WEYL_V1_C1)
						QVFXMUL    (WEYL_V1_C2, QKA0      , HALF_V1_C2)
						QVFXXNPMADD(HALF_V1_C2, HALF_V1_C2, QKA0      , WEYL_V1_C2)
			#endif

						QVSTFDUXA(HALF_V1_C0,TARGETSITE,C32) /* end lifetime: HALF_V1_C0 */
						QVSTFDUXA(HALF_V1_C1,TARGETSITE,C32) /* end lifetime: HALF_V1_C1 */
						QVSTFDUXA(HALF_V1_C2,TARGETSITE,C32) /* end lifetime: HALF_V1_C2 */

						: /* out, inout */
							[gaugesite] "+b" (gaugesite),
							[targetsite] "+b" (targetsite)
						: /* input */
							[spinor_v0_c0] "v" (spinor_v0_c0),
							[spinor_v0_c1] "v" (spinor_v0_c1),
							[spinor_v0_c2] "v" (spinor_v0_c2),
							[spinor_v1_c0] "v" (spinor_v1_c0),
							[spinor_v1_c1] "v" (spinor_v1_c1),
							[spinor_v1_c2] "v" (spinor_v1_c2),
							[spinor_v2_c0] "v" (spinor_v2_c0),
							[spinor_v2_c1] "v" (spinor_v2_c1),
							[spinor_v2_c2] "v" (spinor_v2_c2),
							[spinor_v3_c0] "v" (spinor_v3_c0),
							[spinor_v3_c1] "v" (spinor_v3_c1),
							[spinor_v3_c2] "v" (spinor_v3_c2),
			#if KAMUL
							[qka0] "v" (qka0),
			#endif
							[c32]  "r" (32),
							[c64]  "r" (64),
							[c128] "r" (128)
						: /*clobber*/
							"f" GAUGE_C00,
							"f" GAUGE_C01,
							"f" GAUGE_C02,
							"f" GAUGE_C10,
							"f" GAUGE_C11,
							"f" GAUGE_C12,
							"f" GAUGE_C20,
							"f" GAUGE_C21,
							"f" GAUGE_C22,
							"f" WEYL_V0_C0,
							"f" WEYL_V0_C1,
							"f" WEYL_V0_C2,
							"f" HALF_V0_C0,
							"f" HALF_V0_C1,
							"f" HALF_V0_C2
					);
#endif

#define SPINOR_V0_C0 "%[spinor_v0_c0]"
#define SPINOR_V0_C1 "%[spinor_v0_c1]"
#define SPINOR_V0_C2 "%[spinor_v0_c2]"
#define SPINOR_V1_C0 "%[spinor_v1_c0]"
#define SPINOR_V1_C1 "%[spinor_v1_c1]"
#define SPINOR_V1_C2 "%[spinor_v1_c2]"
#define SPINOR_V2_C0 "%[spinor_v2_c0]"
#define SPINOR_V2_C1 "%[spinor_v2_c1]"
#define SPINOR_V2_C2 "%[spinor_v2_c2]"
#define SPINOR_V3_C0 "%[spinor_v3_c0]"
#define SPINOR_V3_C1 "%[spinor_v3_c1]"
#define SPINOR_V3_C2 "%[spinor_v3_c2]"
#define GAUGE_C00 "12"
#define GAUGE_C01 "13"
#define GAUGE_C02 "14"
#define GAUGE_C10 "15"
#define GAUGE_C11 "16"
#define GAUGE_C12 "17"
#define GAUGE_C20 "18"
#define GAUGE_C21 "19"
#define GAUGE_C22 "20"
#define WEYL_V0_C0 "21"
#define WEYL_V0_C1 "22"
#define WEYL_V0_C2 "23"
#define WEYL_V1_C0 "21"
#define WEYL_V1_C1 "22"
#define WEYL_V1_C2 "23"
#define HALF_V0_C0 "24"
#define HALF_V0_C1 "25"
#define HALF_V0_C2 "26"
#define HALF_V1_C0 "24"
#define HALF_V1_C1 "25"
#define HALF_V1_C2 "26"
#define GAUGESITE "%[gaugesite]"
#define TARGETSITE "%[targetsite]"
#define QKA0 "%[qka0]"
#define C32 "%[c32]"
#define C64 "%[c64]"
#define C128 "%[c128]"
#define QVFSUB(dst,arg1,arg2) "qvfsub " dst "," arg1 "," arg2 " \n"
#define QVLFDUXA(dst,addr,offset) "qvlfduxa " dst "," addr "," offset " \n"
#define QVFXMUL(dst,arg1,arg2) "qvfxmul " dst "," arg1 "," arg2 " \n"
#define QVFXXCPNMADD(dst,arg1,arg2,accum) "qvfxxcpnmadd " dst "," arg1 "," arg2 "," accum "\n"
#define QVFXXNPMADD(dst,arg1,arg2,accum) "qvfxxnpmadd " dst "," arg1 "," arg2 "," accum "\n"
#define QVFXMADD(dst,arg1,arg2,accum) "qvfxmadd " dst "," arg1 "," arg2 "," accum "\n"
#define QVSTFDXA_0(src,addr) "qvstfdxa " src "," "0" "," addr " \n"
#define QVSTFDUXA(src,addr,offset) "qvstfduxa " src "," addr "," offset " \n"
#define DCBZ(offset,addr) "dcbz " offset "," addr " \n"
#define DCBZ_0(addr) "dcbz " "0" "," addr " \n"

#ifndef BGQ_COMPUTEWEYL_INC_
#include "bgq_qpx.h"
#include "bgq_spinorfield.h"
#include "bgq_gaugefield.h"

#include <stdbool.h>

#define PRECISION double
#define KAMUL 1

void bgq_HoppingMatrix_compute_storeWeyllayout_raw(bgq_weyl_ptr_t *targetptrs, bgq_gaugesite *gaugesite, bgq_su3_spinor_params(spinor), bgq_params(qka0), bgq_params(qka1), bgq_params(qka2), bgq_params(qka3),ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z, bool kamul)
#endif
{


#if KAMUL
#define IFKAMUL(...) __VA_ARGS__
#else
#define IFKAMUL(...)
#endif


	bgq_su3_matrixnext_prefetch(gaugesite);

// T+ /////////////////////////////////////////////////////////////////////////
				{
#if 1
		bgq_weyl_vec *targetsite = targetptrs->d[TUP];
		asm (
			QVLFDUXA(GAUGE_C00,GAUGESITE,C32)
			QVFSUB(WEYL_V0_C0,SPINOR_V0_C0,SPINOR_V2_C0)
			QVFSUB(WEYL_V0_C1,SPINOR_V0_C1,SPINOR_V2_C1)
			QVFSUB(WEYL_V0_C2,SPINOR_V0_C2,SPINOR_V2_C2)

			QVLFDUXA(GAUGE_C01,GAUGESITE,C32)
			QVFXMUL(HALF_V0_C0,GAUGE_C00,WEYL_V0_C0)



			QVLFDUXA(GAUGE_C02,GAUGESITE,C32)
			QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C0,GAUGE_C00,HALF_V0_C0)
			QVFXMUL(HALF_V0_C1,GAUGE_C01,WEYL_V0_C0)


			QVLFDUXA(GAUGE_C10,GAUGESITE,C32)
			QVFXMUL(HALF_V0_C2,GAUGE_C02,WEYL_V0_C0)
			QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C0,GAUGE_C01,HALF_V0_C1)


			QVLFDUXA(GAUGE_C11,GAUGESITE,C32)
			QVFXMADD(HALF_V0_C0,GAUGE_C10,WEYL_V0_C1,HALF_V0_C0)
			QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C0,GAUGE_C02,HALF_V0_C2) /* end lifetime: WEYL_V0_C0 */


			QVLFDUXA(GAUGE_C12,GAUGESITE,C32)
			QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C1,GAUGE_C10,HALF_V0_C0)
			QVFXMADD(HALF_V0_C1,GAUGE_C11,WEYL_V0_C1,HALF_V0_C1)


			QVLFDUXA(GAUGE_C20,GAUGESITE,C32)
			QVFXMADD(HALF_V0_C2,GAUGE_C12,WEYL_V0_C1,HALF_V0_C2)
			QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C1,GAUGE_C11,HALF_V0_C1)


			QVLFDUXA(GAUGE_C21,GAUGESITE,C32)
			QVFXMADD(HALF_V0_C0,GAUGE_C20,WEYL_V0_C2,HALF_V0_C0)
			QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C1,GAUGE_C12,HALF_V0_C2) /* end lifetime: WEYL_V0_C1 */


			QVLFDUXA(GAUGE_C22,GAUGESITE,C32)
			QVFXMADD(HALF_V0_C1,GAUGE_C21,WEYL_V0_C2,HALF_V0_C1)
			QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C2,GAUGE_C20,HALF_V0_C0)


			QVFXMADD(HALF_V0_C2,GAUGE_C22,WEYL_V0_C2,HALF_V0_C2)
			QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C2,GAUGE_C21,HALF_V0_C1)
			IFKAMUL(QVFXMUL    (WEYL_V0_C0, QKA0      , HALF_V0_C0))


			QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C2,GAUGE_C22,HALF_V0_C2) /* end lifetime: WEYL_V0_C2 */
			IFKAMUL(QVFXXNPMADD(HALF_V0_C0, HALF_V0_C0, QKA0      , WEYL_V0_C0))
			IFKAMUL(QVFXMUL    (WEYL_V0_C1, QKA0      , HALF_V0_C1))
			QVFSUB(WEYL_V1_C0,SPINOR_V1_C0,SPINOR_V3_C0)
			IFKAMUL(QVFXMUL    (WEYL_V0_C2, QKA0      , HALF_V0_C2))
			QVSTFDXA_0(HALF_V0_C0,TARGETSITE) /* end lifetime: HALF_V0_C0 */

			IFKAMUL(QVFXXNPMADD(HALF_V0_C1, HALF_V0_C1, QKA0      , WEYL_V0_C1))
			IFKAMUL(QVFXXNPMADD(HALF_V0_C2, HALF_V0_C2, QKA0      , WEYL_V0_C2))


			QVFSUB(WEYL_V1_C1,SPINOR_V1_C1,SPINOR_V3_C1)
			QVSTFDUXA(HALF_V0_C1,TARGETSITE,C32) /* end lifetime: HALF_V0_C1 */
			QVSTFDUXA(HALF_V0_C2,TARGETSITE,C32) /* end lifetime: HALF_V0_C2 */

			QVFSUB(WEYL_V1_C2,SPINOR_V1_C2,SPINOR_V3_C2)




			QVFXMUL(HALF_V1_C0,GAUGE_C00,WEYL_V1_C0)
			QVFXMUL(HALF_V1_C1,GAUGE_C01,WEYL_V1_C0)
			QVFXMUL(HALF_V1_C2,GAUGE_C02,WEYL_V1_C0)
			QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C0,GAUGE_C00,HALF_V1_C0) /* end lifetime: GAUGE_C00 */
			QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C0,GAUGE_C01,HALF_V1_C1) /* end lifetime: GAUGE_C01 */
			QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C0,GAUGE_C02,HALF_V1_C2) /* end lifetime: GAUGE_C02, WEYL_V1_C0 */
			QVFXMADD(HALF_V1_C0,GAUGE_C10,WEYL_V1_C1,HALF_V1_C0)
			QVFXMADD(HALF_V1_C1,GAUGE_C11,WEYL_V1_C1,HALF_V1_C1)
			QVFXMADD(HALF_V1_C2,GAUGE_C12,WEYL_V1_C1,HALF_V1_C2)
			QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C1,GAUGE_C10,HALF_V1_C0) /* end lifetime: GAUGE_C10 */
			QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C1,GAUGE_C11,HALF_V1_C1) /* end lifetime: GAUGE_C11 */
			QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C1,GAUGE_C12,HALF_V1_C2) /* end lifetime: GAUGE_C12, WEYL_V1_C1 */
			QVFXMADD(HALF_V1_C0,GAUGE_C20,WEYL_V1_C2,HALF_V1_C0)
			QVFXMADD(HALF_V1_C1,GAUGE_C21,WEYL_V1_C2,HALF_V1_C1)
			QVFXMADD(HALF_V1_C2,GAUGE_C22,WEYL_V1_C2,HALF_V1_C2)
			QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C2,GAUGE_C20,HALF_V1_C0) /* end lifetime: GAUGE_C20 */
			QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C2,GAUGE_C21,HALF_V1_C1) /* end lifetime: GAUGE_C21 */
			QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C2,GAUGE_C22,HALF_V1_C2) /* end lifetime: GAUGE_C22, WEYL_V1_C2 */

#if KAMUL
			QVFXMUL    (WEYL_V1_C0, QKA0      , HALF_V1_C0)
			QVFXMUL    (WEYL_V1_C1, QKA0      , HALF_V1_C1)
			QVFXMUL    (WEYL_V1_C2, QKA0      , HALF_V1_C2)
			QVFXXNPMADD(HALF_V1_C0, HALF_V1_C0, QKA0      , WEYL_V1_C0)
			QVFXXNPMADD(HALF_V1_C1, HALF_V1_C1, QKA0      , WEYL_V1_C1)
			QVFXXNPMADD(HALF_V1_C2, HALF_V1_C2, QKA0      , WEYL_V1_C2)
#endif

			QVSTFDUXA(HALF_V1_C0,TARGETSITE,C32) /* end lifetime: HALF_V1_C0 */
			QVSTFDUXA(HALF_V1_C1,TARGETSITE,C32) /* end lifetime: HALF_V1_C1 */
			QVSTFDUXA(HALF_V1_C2,TARGETSITE,C32) /* end lifetime: HALF_V1_C2 */

			: /* out, inout */
				[gaugesite] "+b" (gaugesite),
				[targetsite] "+b" (targetsite)
			: /* input */
				[spinor_v0_c0] "v" (spinor_v0_c0),
				[spinor_v0_c1] "v" (spinor_v0_c1),
				[spinor_v0_c2] "v" (spinor_v0_c2),
				[spinor_v1_c0] "v" (spinor_v1_c0),
				[spinor_v1_c1] "v" (spinor_v1_c1),
				[spinor_v1_c2] "v" (spinor_v1_c2),
				[spinor_v2_c0] "v" (spinor_v2_c0),
				[spinor_v2_c1] "v" (spinor_v2_c1),
				[spinor_v2_c2] "v" (spinor_v2_c2),
				[spinor_v3_c0] "v" (spinor_v3_c0),
				[spinor_v3_c1] "v" (spinor_v3_c1),
				[spinor_v3_c2] "v" (spinor_v3_c2),
#if KAMUL
				[qka0] "v" (qka0),
#endif
				[c32]  "r" (32),
				[c64]  "r" (64),
				[c128] "r" (128)
			: /*clobber*/
				"f" GAUGE_C00,
				"f" GAUGE_C01,
				"f" GAUGE_C02,
				"f" GAUGE_C10,
				"f" GAUGE_C11,
				"f" GAUGE_C12,
				"f" GAUGE_C20,
				"f" GAUGE_C21,
				"f" GAUGE_C22,
				"f" WEYL_V0_C0,
				"f" WEYL_V0_C1,
				"f" WEYL_V0_C2,
				"f" HALF_V0_C0,
				"f" HALF_V0_C1,
				"f" HALF_V0_C2
		);
//#else
					bgq_su3_weyl_decl(weyl_tup);
					bgq_su3_reduce_weyl_tdown(weyl_tup, spinor);

					bgq_su3_mdecl(gauge_tup);
					bgq_qvlfuxa(gauge_tup_c00, gaugesite, 32);
					bgq_qvlfuxa(gauge_tup_c01, gaugesite, 32);
					bgq_qvlfuxa(gauge_tup_c02, gaugesite, 32);
					bgq_qvlfuxa(gauge_tup_c10, gaugesite, 32);
					bgq_qvlfuxa(gauge_tup_c11, gaugesite, 32);
					bgq_qvlfuxa(gauge_tup_c12, gaugesite, 32);
					bgq_qvlfuxa(gauge_tup_c20, gaugesite, 32);
					bgq_qvlfuxa(gauge_tup_c21, gaugesite, 32);
					bgq_qvlfuxa(gauge_tup_c22, gaugesite, 32);
							bgq_gaugeqpx_expect(gauge_tup, t1, t2, x, y, z, TUP, true);
							bgq_setdesc(BGQREF_TDOWN_GAUGE, "BGQREF_TDOWN_GAUGE");
							bgq_setbgqvalue_src(t1, x, y, z, TUP, BGQREF_TDOWN_GAUGE, bgq_cmplxval1(gauge_tup_c00));
							bgq_setbgqvalue_src(t2, x, y, z, TUP, BGQREF_TDOWN_GAUGE, bgq_cmplxval2(gauge_tup_c00));

					bgq_su3_weyl_mvinvmul(weyl_tup, gauge_tup, weyl_tup);
					if (kamul) {
						bgq_su3_weyl_cmul(weyl_tup, qka0, weyl_tup);
					}
							bgq_setdesc(BGQREF_TDOWN_KAMUL,"BGQREF_TDOWN_KAMUL");
							bgq_setbgqvalue_src(t1, x, y, z, TUP, BGQREF_TDOWN_KAMUL, bgq_cmplxval1(weyl_tup_v1_c0));
							bgq_setbgqvalue_src(t2, x, y, z, TUP, BGQREF_TDOWN_KAMUL, bgq_cmplxval2(weyl_tup_v1_c0));

					bgq_su3_weyl_zeroload(targetptrs->d[TUP]);
					bgq_su3_weyl_store(targetptrs->d[TUP], weyl_tup);
					bgq_weylvec_written(targetptrs->d[TUP], t1, t2, x, y, z, TUP, true);
#endif
				}

				bgq_su3_matrixnext_prefetch(gaugesite);

// T- /////////////////////////////////////////////////////////////////////////
					{
#if 1
					bgq_weyl_vec *targetsite = targetptrs->d[TDOWN];
					asm (
						QVLFDUXA(GAUGE_C00,GAUGESITE,C32)
						QVLFDUXA(GAUGE_C01,GAUGESITE,C32)
						QVLFDUXA(GAUGE_C02,GAUGESITE,C32)
						QVLFDUXA(GAUGE_C10,GAUGESITE,C32)
						QVLFDUXA(GAUGE_C11,GAUGESITE,C32)
						QVLFDUXA(GAUGE_C12,GAUGESITE,C32)
						QVLFDUXA(GAUGE_C20,GAUGESITE,C32)
						QVLFDUXA(GAUGE_C21,GAUGESITE,C32)
						QVLFDUXA(GAUGE_C22,GAUGESITE,C32)

						QVFSUB(WEYL_V0_C0,SPINOR_V0_C0,SPINOR_V2_C0)
						QVFSUB(WEYL_V0_C1,SPINOR_V0_C1,SPINOR_V2_C1)
						QVFSUB(WEYL_V0_C2,SPINOR_V0_C2,SPINOR_V2_C2)

						QVFXMUL(HALF_V0_C0,GAUGE_C00,WEYL_V0_C0)
						QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C0,GAUGE_C00,HALF_V0_C0)
						QVFXMADD(HALF_V0_C0,GAUGE_C10,WEYL_V0_C1,HALF_V0_C0)
						QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C1,GAUGE_C10,HALF_V0_C0)
						QVFXMADD(HALF_V0_C0,GAUGE_C20,WEYL_V0_C2,HALF_V0_C0)
						QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C2,GAUGE_C20,HALF_V0_C0)

						QVFXMUL(HALF_V0_C1,GAUGE_C01,WEYL_V0_C0)
						QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C0,GAUGE_C01,HALF_V0_C1)
						QVFXMADD(HALF_V0_C1,GAUGE_C11,WEYL_V0_C1,HALF_V0_C1)
						QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C1,GAUGE_C11,HALF_V0_C1)
						QVFXMADD(HALF_V0_C1,GAUGE_C21,WEYL_V0_C2,HALF_V0_C1)
						QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C2,GAUGE_C21,HALF_V0_C1)

						QVFXMUL(HALF_V0_C2,GAUGE_C02,WEYL_V0_C0)
						QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C0,GAUGE_C02,HALF_V0_C2) /* end lifetime: WEYL_V0_C0 */
						QVFXMADD(HALF_V0_C2,GAUGE_C12,WEYL_V0_C1,HALF_V0_C2)
						QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C1,GAUGE_C12,HALF_V0_C2) /* end lifetime: WEYL_V0_C1 */
						QVFXMADD(HALF_V0_C2,GAUGE_C22,WEYL_V0_C2,HALF_V0_C2)
						QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C2,GAUGE_C22,HALF_V0_C2) /* end lifetime: WEYL_V0_C2 */

			#if KAMUL
						QVFXMUL    (WEYL_V0_C0, QKA0      , HALF_V0_C0)
						QVFXXNPMADD(HALF_V0_C0, HALF_V0_C0, QKA0      , WEYL_V0_C0)
						QVFXMUL    (WEYL_V0_C1, QKA0      , HALF_V0_C1)
						QVFXXNPMADD(HALF_V0_C1, HALF_V0_C1, QKA0      , WEYL_V0_C1)
						QVFXMUL    (WEYL_V0_C2, QKA0      , HALF_V0_C2)
						QVFXXNPMADD(HALF_V0_C2, HALF_V0_C2, QKA0      , WEYL_V0_C2)
			#endif

						//DCBZ_0(TARGETSITE)
						//DCBZ(C64,TARGETSITE)
						//DCBZ(C128,TARGETSITE)
						QVSTFDXA_0(HALF_V0_C0,TARGETSITE) /* end lifetime: HALF_V0_C0 */
						QVSTFDUXA(HALF_V0_C1,TARGETSITE,C32) /* end lifetime: HALF_V0_C1 */
						QVSTFDUXA(HALF_V0_C2,TARGETSITE,C32) /* end lifetime: HALF_V0_C2 */

						QVFSUB(WEYL_V1_C0,SPINOR_V1_C0,SPINOR_V3_C0)
						QVFSUB(WEYL_V1_C1,SPINOR_V1_C1,SPINOR_V3_C1)
						QVFSUB(WEYL_V1_C2,SPINOR_V1_C2,SPINOR_V3_C2)

						QVFXMUL(HALF_V1_C0,GAUGE_C00,WEYL_V1_C0)
						QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C0,GAUGE_C00,HALF_V1_C0) /* end lifetime: GAUGE_C00 */
						QVFXMADD(HALF_V1_C0,GAUGE_C10,WEYL_V1_C1,HALF_V1_C0)
						QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C1,GAUGE_C10,HALF_V1_C0) /* end lifetime: GAUGE_C10 */
						QVFXMADD(HALF_V1_C0,GAUGE_C20,WEYL_V1_C2,HALF_V1_C0)
						QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C2,GAUGE_C20,HALF_V1_C0) /* end lifetime: GAUGE_C20 */

						QVFXMUL(HALF_V1_C1,GAUGE_C01,WEYL_V1_C0)
						QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C0,GAUGE_C01,HALF_V1_C1) /* end lifetime: GAUGE_C01 */
						QVFXMADD(HALF_V1_C1,GAUGE_C11,WEYL_V1_C1,HALF_V1_C1)
						QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C1,GAUGE_C11,HALF_V1_C1) /* end lifetime: GAUGE_C11 */
						QVFXMADD(HALF_V1_C1,GAUGE_C21,WEYL_V1_C2,HALF_V1_C1)
						QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C2,GAUGE_C21,HALF_V1_C1) /* end lifetime: GAUGE_C21 */

						QVFXMUL(HALF_V1_C2,GAUGE_C02,WEYL_V1_C0)
						QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C0,GAUGE_C02,HALF_V1_C2) /* end lifetime: GAUGE_C02, WEYL_V1_C0 */
						QVFXMADD(HALF_V1_C2,GAUGE_C12,WEYL_V1_C1,HALF_V1_C2)
						QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C1,GAUGE_C12,HALF_V1_C2) /* end lifetime: GAUGE_C12, WEYL_V1_C1 */
						QVFXMADD(HALF_V1_C2,GAUGE_C22,WEYL_V1_C2,HALF_V1_C2)
						QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C2,GAUGE_C22,HALF_V1_C2) /* end lifetime: GAUGE_C22, WEYL_V1_C2 */

			#if KAMUL
						QVFXMUL    (WEYL_V1_C0, QKA0      , HALF_V1_C0)
						QVFXXNPMADD(HALF_V1_C0, HALF_V1_C0, QKA0      , WEYL_V1_C0)
						QVFXMUL    (WEYL_V1_C1, QKA0      , HALF_V1_C1)
						QVFXXNPMADD(HALF_V1_C1, HALF_V1_C1, QKA0      , WEYL_V1_C1)
						QVFXMUL    (WEYL_V1_C2, QKA0      , HALF_V1_C2)
						QVFXXNPMADD(HALF_V1_C2, HALF_V1_C2, QKA0      , WEYL_V1_C2)
			#endif

						QVSTFDUXA(HALF_V1_C0,TARGETSITE,C32) /* end lifetime: HALF_V1_C0 */
						QVSTFDUXA(HALF_V1_C1,TARGETSITE,C32) /* end lifetime: HALF_V1_C1 */
						QVSTFDUXA(HALF_V1_C2,TARGETSITE,C32) /* end lifetime: HALF_V1_C2 */

						: /* out, inout */
							[gaugesite] "+b" (gaugesite),
							[targetsite] "+b" (targetsite)
						: /* input */
							[spinor_v0_c0] "v" (spinor_v0_c0),
							[spinor_v0_c1] "v" (spinor_v0_c1),
							[spinor_v0_c2] "v" (spinor_v0_c2),
							[spinor_v1_c0] "v" (spinor_v1_c0),
							[spinor_v1_c1] "v" (spinor_v1_c1),
							[spinor_v1_c2] "v" (spinor_v1_c2),
							[spinor_v2_c0] "v" (spinor_v2_c0),
							[spinor_v2_c1] "v" (spinor_v2_c1),
							[spinor_v2_c2] "v" (spinor_v2_c2),
							[spinor_v3_c0] "v" (spinor_v3_c0),
							[spinor_v3_c1] "v" (spinor_v3_c1),
							[spinor_v3_c2] "v" (spinor_v3_c2),
			#if KAMUL
							[qka0] "v" (qka0),
			#endif
							[c32]  "r" (32),
							[c64]  "r" (64),
							[c128] "r" (128)
						: /*clobber*/
							"f" GAUGE_C00,
							"f" GAUGE_C01,
							"f" GAUGE_C02,
							"f" GAUGE_C10,
							"f" GAUGE_C11,
							"f" GAUGE_C12,
							"f" GAUGE_C20,
							"f" GAUGE_C21,
							"f" GAUGE_C22,
							"f" WEYL_V0_C0,
							"f" WEYL_V0_C1,
							"f" WEYL_V0_C2,
							"f" HALF_V0_C0,
							"f" HALF_V0_C1,
							"f" HALF_V0_C2
					);
#else
								bgq_setdesc(BGQREF_TUP_SOURCE,"BGQREF_TUP_SOURCE");
								bgq_setbgqvalue_src(t1, x, y, z, TDOWN, BGQREF_TUP_SOURCE, bgq_cmplxval1(spinor_v1_c0));
								bgq_setbgqvalue_src(t2, x, y, z, TDOWN, BGQREF_TUP_SOURCE, bgq_cmplxval2(spinor_v1_c0));

						bgq_su3_weyl_decl(weyl_tdown);
						bgq_su3_reduce_weyl_tup(weyl_tdown, spinor);
								bgq_setdesc(BGQREF_TUP, "BGQREF_TUP");
								bgq_setbgqvalue_src(t1, x, y, z, TDOWN, BGQREF_TUP, bgq_cmplxval1(weyl_tdown_v1_c0));
								bgq_setbgqvalue_src(t2, x, y, z, TDOWN, BGQREF_TUP, bgq_cmplxval2(weyl_tdown_v1_c0));

						bgq_su3_mdecl(gauge_tdown);
						bgq_qvlfuxa(gauge_tdown_c00, gaugesite, 32);
						bgq_qvlfuxa(gauge_tdown_c01, gaugesite, 32);
						bgq_qvlfuxa(gauge_tdown_c02, gaugesite, 32);
						bgq_qvlfuxa(gauge_tdown_c10, gaugesite, 32);
						bgq_qvlfuxa(gauge_tdown_c11, gaugesite, 32);
						bgq_qvlfuxa(gauge_tdown_c12, gaugesite, 32);
						bgq_qvlfuxa(gauge_tdown_c20, gaugesite, 32);
						bgq_qvlfuxa(gauge_tdown_c21, gaugesite, 32);
						bgq_qvlfuxa(gauge_tdown_c22, gaugesite, 32);
								bgq_gaugeqpx_expect(gauge_tdown, t1, t2, x, y, z, TDOWN, true);
								bgq_setdesc(BGQREF_TUP_GAUGE, "BGQREF_TUP_GAUGE");
								bgq_setbgqvalue_src(t1, x, y, z, TDOWN, BGQREF_TUP_GAUGE, bgq_cmplxval1(gauge_tdown_c00));
								bgq_setbgqvalue_src(t2, x, y, z, TDOWN, BGQREF_TUP_GAUGE, bgq_cmplxval2(gauge_tdown_c00));

						bgq_su3_weyl_mvmul(weyl_tdown, gauge_tdown, weyl_tdown);
								bgq_setdesc(BGQREF_TUP_WEYL,"BGQREF_TUP_WEYL");
								bgq_setbgqvalue_src(t1, x, y, z, TDOWN, BGQREF_TUP_WEYL, bgq_cmplxval1(weyl_tdown_v1_c0));
								bgq_setbgqvalue_src(t2, x, y, z, TDOWN, BGQREF_TUP_WEYL, bgq_cmplxval2(weyl_tdown_v1_c0));
						if (kamul) {
							bgq_su3_weyl_cmul(weyl_tdown, qka0, weyl_tdown);
						}
								bgq_setdesc(BGQREF_TUP_KAMUL,"BGQREF_TUP_KAMUL");
								bgq_setbgqvalue_src(t1, x, y, z, TDOWN, BGQREF_TUP_KAMUL, bgq_cmplxval1(weyl_tdown_v1_c0));
								bgq_setbgqvalue_src(t2, x, y, z, TDOWN, BGQREF_TUP_KAMUL, bgq_cmplxval2(weyl_tdown_v1_c0));

						bgq_su3_weyl_zeroload(targetptrs->d[TDOWN]);
						bgq_su3_weyl_store(targetptrs->d[TDOWN], weyl_tdown);
						bgq_weylvec_written(targetptrs->d[TDOWN], t1, t2, x, y, z, TDOWN, true);
#endif
					}

					bgq_su3_matrixnext_prefetch(gaugesite);

					// X+ /////////////////////////////////////////////////////////////////////////
					{
#if 1
						bgq_weyl_vec *targetsite = targetptrs->d[XUP];
						asm (
							QVLFDUXA(GAUGE_C00,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C01,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C02,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C10,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C11,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C12,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C20,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C21,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C22,GAUGESITE,C32)

							QVFSUB(WEYL_V0_C0,SPINOR_V0_C0,SPINOR_V2_C0)
							QVFSUB(WEYL_V0_C1,SPINOR_V0_C1,SPINOR_V2_C1)
							QVFSUB(WEYL_V0_C2,SPINOR_V0_C2,SPINOR_V2_C2)

							QVFXMUL(HALF_V0_C0,GAUGE_C00,WEYL_V0_C0)
							QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C0,GAUGE_C00,HALF_V0_C0)
							QVFXMADD(HALF_V0_C0,GAUGE_C10,WEYL_V0_C1,HALF_V0_C0)
							QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C1,GAUGE_C10,HALF_V0_C0)
							QVFXMADD(HALF_V0_C0,GAUGE_C20,WEYL_V0_C2,HALF_V0_C0)
							QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C2,GAUGE_C20,HALF_V0_C0)

							QVFXMUL(HALF_V0_C1,GAUGE_C01,WEYL_V0_C0)
							QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C0,GAUGE_C01,HALF_V0_C1)
							QVFXMADD(HALF_V0_C1,GAUGE_C11,WEYL_V0_C1,HALF_V0_C1)
							QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C1,GAUGE_C11,HALF_V0_C1)
							QVFXMADD(HALF_V0_C1,GAUGE_C21,WEYL_V0_C2,HALF_V0_C1)
							QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C2,GAUGE_C21,HALF_V0_C1)

							QVFXMUL(HALF_V0_C2,GAUGE_C02,WEYL_V0_C0)
							QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C0,GAUGE_C02,HALF_V0_C2) /* end lifetime: WEYL_V0_C0 */
							QVFXMADD(HALF_V0_C2,GAUGE_C12,WEYL_V0_C1,HALF_V0_C2)
							QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C1,GAUGE_C12,HALF_V0_C2) /* end lifetime: WEYL_V0_C1 */
							QVFXMADD(HALF_V0_C2,GAUGE_C22,WEYL_V0_C2,HALF_V0_C2)
							QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C2,GAUGE_C22,HALF_V0_C2) /* end lifetime: WEYL_V0_C2 */

				#if KAMUL
							QVFXMUL    (WEYL_V0_C0, QKA0      , HALF_V0_C0)
							QVFXXNPMADD(HALF_V0_C0, HALF_V0_C0, QKA0      , WEYL_V0_C0)
							QVFXMUL    (WEYL_V0_C1, QKA0      , HALF_V0_C1)
							QVFXXNPMADD(HALF_V0_C1, HALF_V0_C1, QKA0      , WEYL_V0_C1)
							QVFXMUL    (WEYL_V0_C2, QKA0      , HALF_V0_C2)
							QVFXXNPMADD(HALF_V0_C2, HALF_V0_C2, QKA0      , WEYL_V0_C2)
				#endif

							//DCBZ_0(TARGETSITE)
							//DCBZ(C64,TARGETSITE)
							//DCBZ(C128,TARGETSITE)
							QVSTFDXA_0(HALF_V0_C0,TARGETSITE) /* end lifetime: HALF_V0_C0 */
							QVSTFDUXA(HALF_V0_C1,TARGETSITE,C32) /* end lifetime: HALF_V0_C1 */
							QVSTFDUXA(HALF_V0_C2,TARGETSITE,C32) /* end lifetime: HALF_V0_C2 */

							QVFSUB(WEYL_V1_C0,SPINOR_V1_C0,SPINOR_V3_C0)
							QVFSUB(WEYL_V1_C1,SPINOR_V1_C1,SPINOR_V3_C1)
							QVFSUB(WEYL_V1_C2,SPINOR_V1_C2,SPINOR_V3_C2)

							QVFXMUL(HALF_V1_C0,GAUGE_C00,WEYL_V1_C0)
							QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C0,GAUGE_C00,HALF_V1_C0) /* end lifetime: GAUGE_C00 */
							QVFXMADD(HALF_V1_C0,GAUGE_C10,WEYL_V1_C1,HALF_V1_C0)
							QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C1,GAUGE_C10,HALF_V1_C0) /* end lifetime: GAUGE_C10 */
							QVFXMADD(HALF_V1_C0,GAUGE_C20,WEYL_V1_C2,HALF_V1_C0)
							QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C2,GAUGE_C20,HALF_V1_C0) /* end lifetime: GAUGE_C20 */

							QVFXMUL(HALF_V1_C1,GAUGE_C01,WEYL_V1_C0)
							QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C0,GAUGE_C01,HALF_V1_C1) /* end lifetime: GAUGE_C01 */
							QVFXMADD(HALF_V1_C1,GAUGE_C11,WEYL_V1_C1,HALF_V1_C1)
							QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C1,GAUGE_C11,HALF_V1_C1) /* end lifetime: GAUGE_C11 */
							QVFXMADD(HALF_V1_C1,GAUGE_C21,WEYL_V1_C2,HALF_V1_C1)
							QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C2,GAUGE_C21,HALF_V1_C1) /* end lifetime: GAUGE_C21 */

							QVFXMUL(HALF_V1_C2,GAUGE_C02,WEYL_V1_C0)
							QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C0,GAUGE_C02,HALF_V1_C2) /* end lifetime: GAUGE_C02, WEYL_V1_C0 */
							QVFXMADD(HALF_V1_C2,GAUGE_C12,WEYL_V1_C1,HALF_V1_C2)
							QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C1,GAUGE_C12,HALF_V1_C2) /* end lifetime: GAUGE_C12, WEYL_V1_C1 */
							QVFXMADD(HALF_V1_C2,GAUGE_C22,WEYL_V1_C2,HALF_V1_C2)
							QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C2,GAUGE_C22,HALF_V1_C2) /* end lifetime: GAUGE_C22, WEYL_V1_C2 */

				#if KAMUL
							QVFXMUL    (WEYL_V1_C0, QKA0      , HALF_V1_C0)
							QVFXXNPMADD(HALF_V1_C0, HALF_V1_C0, QKA0      , WEYL_V1_C0)
							QVFXMUL    (WEYL_V1_C1, QKA0      , HALF_V1_C1)
							QVFXXNPMADD(HALF_V1_C1, HALF_V1_C1, QKA0      , WEYL_V1_C1)
							QVFXMUL    (WEYL_V1_C2, QKA0      , HALF_V1_C2)
							QVFXXNPMADD(HALF_V1_C2, HALF_V1_C2, QKA0      , WEYL_V1_C2)
				#endif

							QVSTFDUXA(HALF_V1_C0,TARGETSITE,C32) /* end lifetime: HALF_V1_C0 */
							QVSTFDUXA(HALF_V1_C1,TARGETSITE,C32) /* end lifetime: HALF_V1_C1 */
							QVSTFDUXA(HALF_V1_C2,TARGETSITE,C32) /* end lifetime: HALF_V1_C2 */

							: /* out, inout */
								[gaugesite] "+b" (gaugesite),
								[targetsite] "+b" (targetsite)
							: /* input */
								[spinor_v0_c0] "v" (spinor_v0_c0),
								[spinor_v0_c1] "v" (spinor_v0_c1),
								[spinor_v0_c2] "v" (spinor_v0_c2),
								[spinor_v1_c0] "v" (spinor_v1_c0),
								[spinor_v1_c1] "v" (spinor_v1_c1),
								[spinor_v1_c2] "v" (spinor_v1_c2),
								[spinor_v2_c0] "v" (spinor_v2_c0),
								[spinor_v2_c1] "v" (spinor_v2_c1),
								[spinor_v2_c2] "v" (spinor_v2_c2),
								[spinor_v3_c0] "v" (spinor_v3_c0),
								[spinor_v3_c1] "v" (spinor_v3_c1),
								[spinor_v3_c2] "v" (spinor_v3_c2),
				#if KAMUL
								[qka0] "v" (qka0),
				#endif
								[c32]  "r" (32),
								[c64]  "r" (64),
								[c128] "r" (128)
							: /*clobber*/
								"f" GAUGE_C00,
								"f" GAUGE_C01,
								"f" GAUGE_C02,
								"f" GAUGE_C10,
								"f" GAUGE_C11,
								"f" GAUGE_C12,
								"f" GAUGE_C20,
								"f" GAUGE_C21,
								"f" GAUGE_C22,
								"f" WEYL_V0_C0,
								"f" WEYL_V0_C1,
								"f" WEYL_V0_C2,
								"f" HALF_V0_C0,
								"f" HALF_V0_C1,
								"f" HALF_V0_C2
						);
#else
						bgq_su3_weyl_decl(weyl_xup);
						bgq_su3_reduce_weyl_xdown(weyl_xup, spinor);

						bgq_su3_mdecl(gauge_xup);
						bgq_qvlfuxa(gauge_xup_c00, gaugesite, 32);
						bgq_qvlfuxa(gauge_xup_c01, gaugesite, 32);
						bgq_qvlfuxa(gauge_xup_c02, gaugesite, 32);
						bgq_qvlfuxa(gauge_xup_c10, gaugesite, 32);
						bgq_qvlfuxa(gauge_xup_c11, gaugesite, 32);
						bgq_qvlfuxa(gauge_xup_c12, gaugesite, 32);
						bgq_qvlfuxa(gauge_xup_c20, gaugesite, 32);
						bgq_qvlfuxa(gauge_xup_c21, gaugesite, 32);
						bgq_qvlfuxa(gauge_xup_c22, gaugesite, 32);

						bgq_su3_weyl_mvinvmul(weyl_xup, gauge_xup, weyl_xup);
						if (kamul) {
							bgq_su3_weyl_cmul(weyl_xup, qka1, weyl_xup);
						}

						bgq_su3_weyl_zeroload(targetptrs->d[XUP]);
						bgq_su3_weyl_store(targetptrs->d[XUP], weyl_xup);
						bgq_weylvec_written(targetptrs->d[XUP], t1, t2, x,y,z,XUP, true);
#endif
					}

					bgq_su3_matrixnext_prefetch(gaugesite);

					// X- /////////////////////////////////////////////////////////////////////////
					{
#if 1
						bgq_weyl_vec *targetsite = targetptrs->d[XDOWN];
						asm (
							QVLFDUXA(GAUGE_C00,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C01,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C02,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C10,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C11,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C12,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C20,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C21,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C22,GAUGESITE,C32)

							QVFSUB(WEYL_V0_C0,SPINOR_V0_C0,SPINOR_V2_C0)
							QVFSUB(WEYL_V0_C1,SPINOR_V0_C1,SPINOR_V2_C1)
							QVFSUB(WEYL_V0_C2,SPINOR_V0_C2,SPINOR_V2_C2)

							QVFXMUL(HALF_V0_C0,GAUGE_C00,WEYL_V0_C0)
							QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C0,GAUGE_C00,HALF_V0_C0)
							QVFXMADD(HALF_V0_C0,GAUGE_C10,WEYL_V0_C1,HALF_V0_C0)
							QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C1,GAUGE_C10,HALF_V0_C0)
							QVFXMADD(HALF_V0_C0,GAUGE_C20,WEYL_V0_C2,HALF_V0_C0)
							QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C2,GAUGE_C20,HALF_V0_C0)

							QVFXMUL(HALF_V0_C1,GAUGE_C01,WEYL_V0_C0)
							QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C0,GAUGE_C01,HALF_V0_C1)
							QVFXMADD(HALF_V0_C1,GAUGE_C11,WEYL_V0_C1,HALF_V0_C1)
							QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C1,GAUGE_C11,HALF_V0_C1)
							QVFXMADD(HALF_V0_C1,GAUGE_C21,WEYL_V0_C2,HALF_V0_C1)
							QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C2,GAUGE_C21,HALF_V0_C1)

							QVFXMUL(HALF_V0_C2,GAUGE_C02,WEYL_V0_C0)
							QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C0,GAUGE_C02,HALF_V0_C2) /* end lifetime: WEYL_V0_C0 */
							QVFXMADD(HALF_V0_C2,GAUGE_C12,WEYL_V0_C1,HALF_V0_C2)
							QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C1,GAUGE_C12,HALF_V0_C2) /* end lifetime: WEYL_V0_C1 */
							QVFXMADD(HALF_V0_C2,GAUGE_C22,WEYL_V0_C2,HALF_V0_C2)
							QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C2,GAUGE_C22,HALF_V0_C2) /* end lifetime: WEYL_V0_C2 */

				#if KAMUL
							QVFXMUL    (WEYL_V0_C0, QKA0      , HALF_V0_C0)
							QVFXXNPMADD(HALF_V0_C0, HALF_V0_C0, QKA0      , WEYL_V0_C0)
							QVFXMUL    (WEYL_V0_C1, QKA0      , HALF_V0_C1)
							QVFXXNPMADD(HALF_V0_C1, HALF_V0_C1, QKA0      , WEYL_V0_C1)
							QVFXMUL    (WEYL_V0_C2, QKA0      , HALF_V0_C2)
							QVFXXNPMADD(HALF_V0_C2, HALF_V0_C2, QKA0      , WEYL_V0_C2)
				#endif

							//DCBZ_0(TARGETSITE)
							//DCBZ(C64,TARGETSITE)
							//DCBZ(C128,TARGETSITE)
							QVSTFDXA_0(HALF_V0_C0,TARGETSITE) /* end lifetime: HALF_V0_C0 */
							QVSTFDUXA(HALF_V0_C1,TARGETSITE,C32) /* end lifetime: HALF_V0_C1 */
							QVSTFDUXA(HALF_V0_C2,TARGETSITE,C32) /* end lifetime: HALF_V0_C2 */

							QVFSUB(WEYL_V1_C0,SPINOR_V1_C0,SPINOR_V3_C0)
							QVFSUB(WEYL_V1_C1,SPINOR_V1_C1,SPINOR_V3_C1)
							QVFSUB(WEYL_V1_C2,SPINOR_V1_C2,SPINOR_V3_C2)

							QVFXMUL(HALF_V1_C0,GAUGE_C00,WEYL_V1_C0)
							QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C0,GAUGE_C00,HALF_V1_C0) /* end lifetime: GAUGE_C00 */
							QVFXMADD(HALF_V1_C0,GAUGE_C10,WEYL_V1_C1,HALF_V1_C0)
							QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C1,GAUGE_C10,HALF_V1_C0) /* end lifetime: GAUGE_C10 */
							QVFXMADD(HALF_V1_C0,GAUGE_C20,WEYL_V1_C2,HALF_V1_C0)
							QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C2,GAUGE_C20,HALF_V1_C0) /* end lifetime: GAUGE_C20 */

							QVFXMUL(HALF_V1_C1,GAUGE_C01,WEYL_V1_C0)
							QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C0,GAUGE_C01,HALF_V1_C1) /* end lifetime: GAUGE_C01 */
							QVFXMADD(HALF_V1_C1,GAUGE_C11,WEYL_V1_C1,HALF_V1_C1)
							QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C1,GAUGE_C11,HALF_V1_C1) /* end lifetime: GAUGE_C11 */
							QVFXMADD(HALF_V1_C1,GAUGE_C21,WEYL_V1_C2,HALF_V1_C1)
							QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C2,GAUGE_C21,HALF_V1_C1) /* end lifetime: GAUGE_C21 */

							QVFXMUL(HALF_V1_C2,GAUGE_C02,WEYL_V1_C0)
							QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C0,GAUGE_C02,HALF_V1_C2) /* end lifetime: GAUGE_C02, WEYL_V1_C0 */
							QVFXMADD(HALF_V1_C2,GAUGE_C12,WEYL_V1_C1,HALF_V1_C2)
							QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C1,GAUGE_C12,HALF_V1_C2) /* end lifetime: GAUGE_C12, WEYL_V1_C1 */
							QVFXMADD(HALF_V1_C2,GAUGE_C22,WEYL_V1_C2,HALF_V1_C2)
							QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C2,GAUGE_C22,HALF_V1_C2) /* end lifetime: GAUGE_C22, WEYL_V1_C2 */

				#if KAMUL
							QVFXMUL    (WEYL_V1_C0, QKA0      , HALF_V1_C0)
							QVFXXNPMADD(HALF_V1_C0, HALF_V1_C0, QKA0      , WEYL_V1_C0)
							QVFXMUL    (WEYL_V1_C1, QKA0      , HALF_V1_C1)
							QVFXXNPMADD(HALF_V1_C1, HALF_V1_C1, QKA0      , WEYL_V1_C1)
							QVFXMUL    (WEYL_V1_C2, QKA0      , HALF_V1_C2)
							QVFXXNPMADD(HALF_V1_C2, HALF_V1_C2, QKA0      , WEYL_V1_C2)
				#endif

							QVSTFDUXA(HALF_V1_C0,TARGETSITE,C32) /* end lifetime: HALF_V1_C0 */
							QVSTFDUXA(HALF_V1_C1,TARGETSITE,C32) /* end lifetime: HALF_V1_C1 */
							QVSTFDUXA(HALF_V1_C2,TARGETSITE,C32) /* end lifetime: HALF_V1_C2 */

							: /* out, inout */
								[gaugesite] "+b" (gaugesite),
								[targetsite] "+b" (targetsite)
							: /* input */
								[spinor_v0_c0] "v" (spinor_v0_c0),
								[spinor_v0_c1] "v" (spinor_v0_c1),
								[spinor_v0_c2] "v" (spinor_v0_c2),
								[spinor_v1_c0] "v" (spinor_v1_c0),
								[spinor_v1_c1] "v" (spinor_v1_c1),
								[spinor_v1_c2] "v" (spinor_v1_c2),
								[spinor_v2_c0] "v" (spinor_v2_c0),
								[spinor_v2_c1] "v" (spinor_v2_c1),
								[spinor_v2_c2] "v" (spinor_v2_c2),
								[spinor_v3_c0] "v" (spinor_v3_c0),
								[spinor_v3_c1] "v" (spinor_v3_c1),
								[spinor_v3_c2] "v" (spinor_v3_c2),
				#if KAMUL
								[qka0] "v" (qka0),
				#endif
								[c32]  "r" (32),
								[c64]  "r" (64),
								[c128] "r" (128)
							: /*clobber*/
								"f" GAUGE_C00,
								"f" GAUGE_C01,
								"f" GAUGE_C02,
								"f" GAUGE_C10,
								"f" GAUGE_C11,
								"f" GAUGE_C12,
								"f" GAUGE_C20,
								"f" GAUGE_C21,
								"f" GAUGE_C22,
								"f" WEYL_V0_C0,
								"f" WEYL_V0_C1,
								"f" WEYL_V0_C2,
								"f" HALF_V0_C0,
								"f" HALF_V0_C1,
								"f" HALF_V0_C2
						);
#else
						bgq_su3_weyl_decl(weyl_xdown);
						bgq_su3_reduce_weyl_xup(weyl_xdown, spinor);

						bgq_su3_mdecl(gauge_xdown);
						bgq_qvlfuxa(gauge_xdown_c00, gaugesite, 32);
						bgq_qvlfuxa(gauge_xdown_c01, gaugesite, 32);
						bgq_qvlfuxa(gauge_xdown_c02, gaugesite, 32);
						bgq_qvlfuxa(gauge_xdown_c10, gaugesite, 32);
						bgq_qvlfuxa(gauge_xdown_c11, gaugesite, 32);
						bgq_qvlfuxa(gauge_xdown_c12, gaugesite, 32);
						bgq_qvlfuxa(gauge_xdown_c20, gaugesite, 32);
						bgq_qvlfuxa(gauge_xdown_c21, gaugesite, 32);
						bgq_qvlfuxa(gauge_xdown_c22, gaugesite, 32);

						bgq_su3_weyl_mvmul(weyl_xdown, gauge_xdown, weyl_xdown);
						if (kamul) {
							bgq_su3_weyl_cmul(weyl_xdown, qka1, weyl_xdown);
						}

						bgq_su3_weyl_zeroload(targetptrs->d[XDOWN]);
						bgq_su3_weyl_store(targetptrs->d[XDOWN], weyl_xdown);
						bgq_weylvec_written(targetptrs->d[XDOWN], t1, t2, x,y,z,XDOWN, true);
#endif
					}

					bgq_su3_matrixnext_prefetch(gaugesite);

					// Y+ /////////////////////////////////////////////////////////////////////////
					{
#if 1
						bgq_weyl_vec *targetsite = targetptrs->d[YUP];
						asm (
							QVLFDUXA(GAUGE_C00,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C01,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C02,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C10,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C11,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C12,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C20,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C21,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C22,GAUGESITE,C32)

							QVFSUB(WEYL_V0_C0,SPINOR_V0_C0,SPINOR_V2_C0)
							QVFSUB(WEYL_V0_C1,SPINOR_V0_C1,SPINOR_V2_C1)
							QVFSUB(WEYL_V0_C2,SPINOR_V0_C2,SPINOR_V2_C2)

							QVFXMUL(HALF_V0_C0,GAUGE_C00,WEYL_V0_C0)
							QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C0,GAUGE_C00,HALF_V0_C0)
							QVFXMADD(HALF_V0_C0,GAUGE_C10,WEYL_V0_C1,HALF_V0_C0)
							QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C1,GAUGE_C10,HALF_V0_C0)
							QVFXMADD(HALF_V0_C0,GAUGE_C20,WEYL_V0_C2,HALF_V0_C0)
							QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C2,GAUGE_C20,HALF_V0_C0)

							QVFXMUL(HALF_V0_C1,GAUGE_C01,WEYL_V0_C0)
							QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C0,GAUGE_C01,HALF_V0_C1)
							QVFXMADD(HALF_V0_C1,GAUGE_C11,WEYL_V0_C1,HALF_V0_C1)
							QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C1,GAUGE_C11,HALF_V0_C1)
							QVFXMADD(HALF_V0_C1,GAUGE_C21,WEYL_V0_C2,HALF_V0_C1)
							QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C2,GAUGE_C21,HALF_V0_C1)

							QVFXMUL(HALF_V0_C2,GAUGE_C02,WEYL_V0_C0)
							QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C0,GAUGE_C02,HALF_V0_C2) /* end lifetime: WEYL_V0_C0 */
							QVFXMADD(HALF_V0_C2,GAUGE_C12,WEYL_V0_C1,HALF_V0_C2)
							QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C1,GAUGE_C12,HALF_V0_C2) /* end lifetime: WEYL_V0_C1 */
							QVFXMADD(HALF_V0_C2,GAUGE_C22,WEYL_V0_C2,HALF_V0_C2)
							QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C2,GAUGE_C22,HALF_V0_C2) /* end lifetime: WEYL_V0_C2 */

				#if KAMUL
							QVFXMUL    (WEYL_V0_C0, QKA0      , HALF_V0_C0)
							QVFXXNPMADD(HALF_V0_C0, HALF_V0_C0, QKA0      , WEYL_V0_C0)
							QVFXMUL    (WEYL_V0_C1, QKA0      , HALF_V0_C1)
							QVFXXNPMADD(HALF_V0_C1, HALF_V0_C1, QKA0      , WEYL_V0_C1)
							QVFXMUL    (WEYL_V0_C2, QKA0      , HALF_V0_C2)
							QVFXXNPMADD(HALF_V0_C2, HALF_V0_C2, QKA0      , WEYL_V0_C2)
				#endif

							//DCBZ_0(TARGETSITE)
							//DCBZ(C64,TARGETSITE)
							//DCBZ(C128,TARGETSITE)
							QVSTFDXA_0(HALF_V0_C0,TARGETSITE) /* end lifetime: HALF_V0_C0 */
							QVSTFDUXA(HALF_V0_C1,TARGETSITE,C32) /* end lifetime: HALF_V0_C1 */
							QVSTFDUXA(HALF_V0_C2,TARGETSITE,C32) /* end lifetime: HALF_V0_C2 */

							QVFSUB(WEYL_V1_C0,SPINOR_V1_C0,SPINOR_V3_C0)
							QVFSUB(WEYL_V1_C1,SPINOR_V1_C1,SPINOR_V3_C1)
							QVFSUB(WEYL_V1_C2,SPINOR_V1_C2,SPINOR_V3_C2)

							QVFXMUL(HALF_V1_C0,GAUGE_C00,WEYL_V1_C0)
							QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C0,GAUGE_C00,HALF_V1_C0) /* end lifetime: GAUGE_C00 */
							QVFXMADD(HALF_V1_C0,GAUGE_C10,WEYL_V1_C1,HALF_V1_C0)
							QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C1,GAUGE_C10,HALF_V1_C0) /* end lifetime: GAUGE_C10 */
							QVFXMADD(HALF_V1_C0,GAUGE_C20,WEYL_V1_C2,HALF_V1_C0)
							QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C2,GAUGE_C20,HALF_V1_C0) /* end lifetime: GAUGE_C20 */

							QVFXMUL(HALF_V1_C1,GAUGE_C01,WEYL_V1_C0)
							QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C0,GAUGE_C01,HALF_V1_C1) /* end lifetime: GAUGE_C01 */
							QVFXMADD(HALF_V1_C1,GAUGE_C11,WEYL_V1_C1,HALF_V1_C1)
							QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C1,GAUGE_C11,HALF_V1_C1) /* end lifetime: GAUGE_C11 */
							QVFXMADD(HALF_V1_C1,GAUGE_C21,WEYL_V1_C2,HALF_V1_C1)
							QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C2,GAUGE_C21,HALF_V1_C1) /* end lifetime: GAUGE_C21 */

							QVFXMUL(HALF_V1_C2,GAUGE_C02,WEYL_V1_C0)
							QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C0,GAUGE_C02,HALF_V1_C2) /* end lifetime: GAUGE_C02, WEYL_V1_C0 */
							QVFXMADD(HALF_V1_C2,GAUGE_C12,WEYL_V1_C1,HALF_V1_C2)
							QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C1,GAUGE_C12,HALF_V1_C2) /* end lifetime: GAUGE_C12, WEYL_V1_C1 */
							QVFXMADD(HALF_V1_C2,GAUGE_C22,WEYL_V1_C2,HALF_V1_C2)
							QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C2,GAUGE_C22,HALF_V1_C2) /* end lifetime: GAUGE_C22, WEYL_V1_C2 */

				#if KAMUL
							QVFXMUL    (WEYL_V1_C0, QKA0      , HALF_V1_C0)
							QVFXXNPMADD(HALF_V1_C0, HALF_V1_C0, QKA0      , WEYL_V1_C0)
							QVFXMUL    (WEYL_V1_C1, QKA0      , HALF_V1_C1)
							QVFXXNPMADD(HALF_V1_C1, HALF_V1_C1, QKA0      , WEYL_V1_C1)
							QVFXMUL    (WEYL_V1_C2, QKA0      , HALF_V1_C2)
							QVFXXNPMADD(HALF_V1_C2, HALF_V1_C2, QKA0      , WEYL_V1_C2)
				#endif

							QVSTFDUXA(HALF_V1_C0,TARGETSITE,C32) /* end lifetime: HALF_V1_C0 */
							QVSTFDUXA(HALF_V1_C1,TARGETSITE,C32) /* end lifetime: HALF_V1_C1 */
							QVSTFDUXA(HALF_V1_C2,TARGETSITE,C32) /* end lifetime: HALF_V1_C2 */

							: /* out, inout */
								[gaugesite] "+b" (gaugesite),
								[targetsite] "+b" (targetsite)
							: /* input */
								[spinor_v0_c0] "v" (spinor_v0_c0),
								[spinor_v0_c1] "v" (spinor_v0_c1),
								[spinor_v0_c2] "v" (spinor_v0_c2),
								[spinor_v1_c0] "v" (spinor_v1_c0),
								[spinor_v1_c1] "v" (spinor_v1_c1),
								[spinor_v1_c2] "v" (spinor_v1_c2),
								[spinor_v2_c0] "v" (spinor_v2_c0),
								[spinor_v2_c1] "v" (spinor_v2_c1),
								[spinor_v2_c2] "v" (spinor_v2_c2),
								[spinor_v3_c0] "v" (spinor_v3_c0),
								[spinor_v3_c1] "v" (spinor_v3_c1),
								[spinor_v3_c2] "v" (spinor_v3_c2),
				#if KAMUL
								[qka0] "v" (qka0),
				#endif
								[c32]  "r" (32),
								[c64]  "r" (64),
								[c128] "r" (128)
							: /*clobber*/
								"f" GAUGE_C00,
								"f" GAUGE_C01,
								"f" GAUGE_C02,
								"f" GAUGE_C10,
								"f" GAUGE_C11,
								"f" GAUGE_C12,
								"f" GAUGE_C20,
								"f" GAUGE_C21,
								"f" GAUGE_C22,
								"f" WEYL_V0_C0,
								"f" WEYL_V0_C1,
								"f" WEYL_V0_C2,
								"f" HALF_V0_C0,
								"f" HALF_V0_C1,
								"f" HALF_V0_C2
						);
#else
						bgq_su3_weyl_decl(weyl_yup);
						bgq_su3_reduce_weyl_ydown(weyl_yup, spinor);

						bgq_su3_mdecl(gauge_yup);
						bgq_qvlfuxa(gauge_yup_c00, gaugesite, 32);
						bgq_qvlfuxa(gauge_yup_c01, gaugesite, 32);
						bgq_qvlfuxa(gauge_yup_c02, gaugesite, 32);
						bgq_qvlfuxa(gauge_yup_c10, gaugesite, 32);
						bgq_qvlfuxa(gauge_yup_c11, gaugesite, 32);
						bgq_qvlfuxa(gauge_yup_c12, gaugesite, 32);
						bgq_qvlfuxa(gauge_yup_c20, gaugesite, 32);
						bgq_qvlfuxa(gauge_yup_c21, gaugesite, 32);
						bgq_qvlfuxa(gauge_yup_c22, gaugesite, 32);

						bgq_su3_weyl_mvinvmul(weyl_yup, gauge_yup, weyl_yup);
						if (kamul) {
							bgq_su3_weyl_cmul(weyl_yup, qka2, weyl_yup);
						}

						bgq_su3_weyl_zeroload(targetptrs->d[YUP]);
						bgq_su3_weyl_store(targetptrs->d[YUP], weyl_yup);
						bgq_weylvec_written(targetptrs->d[YUP], t1, t2, x,y,z,YUP, true);
#endif
					}

					bgq_su3_matrixnext_prefetch(gaugesite);

					// Y- /////////////////////////////////////////////////////////////////////////
					{
#if 1
						bgq_weyl_vec *targetsite = targetptrs->d[YDOWN];
						asm (
							QVLFDUXA(GAUGE_C00,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C01,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C02,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C10,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C11,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C12,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C20,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C21,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C22,GAUGESITE,C32)

							QVFSUB(WEYL_V0_C0,SPINOR_V0_C0,SPINOR_V2_C0)
							QVFSUB(WEYL_V0_C1,SPINOR_V0_C1,SPINOR_V2_C1)
							QVFSUB(WEYL_V0_C2,SPINOR_V0_C2,SPINOR_V2_C2)

							QVFXMUL(HALF_V0_C0,GAUGE_C00,WEYL_V0_C0)
							QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C0,GAUGE_C00,HALF_V0_C0)
							QVFXMADD(HALF_V0_C0,GAUGE_C10,WEYL_V0_C1,HALF_V0_C0)
							QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C1,GAUGE_C10,HALF_V0_C0)
							QVFXMADD(HALF_V0_C0,GAUGE_C20,WEYL_V0_C2,HALF_V0_C0)
							QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C2,GAUGE_C20,HALF_V0_C0)

							QVFXMUL(HALF_V0_C1,GAUGE_C01,WEYL_V0_C0)
							QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C0,GAUGE_C01,HALF_V0_C1)
							QVFXMADD(HALF_V0_C1,GAUGE_C11,WEYL_V0_C1,HALF_V0_C1)
							QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C1,GAUGE_C11,HALF_V0_C1)
							QVFXMADD(HALF_V0_C1,GAUGE_C21,WEYL_V0_C2,HALF_V0_C1)
							QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C2,GAUGE_C21,HALF_V0_C1)

							QVFXMUL(HALF_V0_C2,GAUGE_C02,WEYL_V0_C0)
							QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C0,GAUGE_C02,HALF_V0_C2) /* end lifetime: WEYL_V0_C0 */
							QVFXMADD(HALF_V0_C2,GAUGE_C12,WEYL_V0_C1,HALF_V0_C2)
							QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C1,GAUGE_C12,HALF_V0_C2) /* end lifetime: WEYL_V0_C1 */
							QVFXMADD(HALF_V0_C2,GAUGE_C22,WEYL_V0_C2,HALF_V0_C2)
							QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C2,GAUGE_C22,HALF_V0_C2) /* end lifetime: WEYL_V0_C2 */

				#if KAMUL
							QVFXMUL    (WEYL_V0_C0, QKA0      , HALF_V0_C0)
							QVFXXNPMADD(HALF_V0_C0, HALF_V0_C0, QKA0      , WEYL_V0_C0)
							QVFXMUL    (WEYL_V0_C1, QKA0      , HALF_V0_C1)
							QVFXXNPMADD(HALF_V0_C1, HALF_V0_C1, QKA0      , WEYL_V0_C1)
							QVFXMUL    (WEYL_V0_C2, QKA0      , HALF_V0_C2)
							QVFXXNPMADD(HALF_V0_C2, HALF_V0_C2, QKA0      , WEYL_V0_C2)
				#endif

							//DCBZ_0(TARGETSITE)
							//DCBZ(C64,TARGETSITE)
							//DCBZ(C128,TARGETSITE)
							QVSTFDXA_0(HALF_V0_C0,TARGETSITE) /* end lifetime: HALF_V0_C0 */
							QVSTFDUXA(HALF_V0_C1,TARGETSITE,C32) /* end lifetime: HALF_V0_C1 */
							QVSTFDUXA(HALF_V0_C2,TARGETSITE,C32) /* end lifetime: HALF_V0_C2 */

							QVFSUB(WEYL_V1_C0,SPINOR_V1_C0,SPINOR_V3_C0)
							QVFSUB(WEYL_V1_C1,SPINOR_V1_C1,SPINOR_V3_C1)
							QVFSUB(WEYL_V1_C2,SPINOR_V1_C2,SPINOR_V3_C2)

							QVFXMUL(HALF_V1_C0,GAUGE_C00,WEYL_V1_C0)
							QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C0,GAUGE_C00,HALF_V1_C0) /* end lifetime: GAUGE_C00 */
							QVFXMADD(HALF_V1_C0,GAUGE_C10,WEYL_V1_C1,HALF_V1_C0)
							QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C1,GAUGE_C10,HALF_V1_C0) /* end lifetime: GAUGE_C10 */
							QVFXMADD(HALF_V1_C0,GAUGE_C20,WEYL_V1_C2,HALF_V1_C0)
							QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C2,GAUGE_C20,HALF_V1_C0) /* end lifetime: GAUGE_C20 */

							QVFXMUL(HALF_V1_C1,GAUGE_C01,WEYL_V1_C0)
							QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C0,GAUGE_C01,HALF_V1_C1) /* end lifetime: GAUGE_C01 */
							QVFXMADD(HALF_V1_C1,GAUGE_C11,WEYL_V1_C1,HALF_V1_C1)
							QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C1,GAUGE_C11,HALF_V1_C1) /* end lifetime: GAUGE_C11 */
							QVFXMADD(HALF_V1_C1,GAUGE_C21,WEYL_V1_C2,HALF_V1_C1)
							QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C2,GAUGE_C21,HALF_V1_C1) /* end lifetime: GAUGE_C21 */

							QVFXMUL(HALF_V1_C2,GAUGE_C02,WEYL_V1_C0)
							QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C0,GAUGE_C02,HALF_V1_C2) /* end lifetime: GAUGE_C02, WEYL_V1_C0 */
							QVFXMADD(HALF_V1_C2,GAUGE_C12,WEYL_V1_C1,HALF_V1_C2)
							QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C1,GAUGE_C12,HALF_V1_C2) /* end lifetime: GAUGE_C12, WEYL_V1_C1 */
							QVFXMADD(HALF_V1_C2,GAUGE_C22,WEYL_V1_C2,HALF_V1_C2)
							QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C2,GAUGE_C22,HALF_V1_C2) /* end lifetime: GAUGE_C22, WEYL_V1_C2 */

				#if KAMUL
							QVFXMUL    (WEYL_V1_C0, QKA0      , HALF_V1_C0)
							QVFXXNPMADD(HALF_V1_C0, HALF_V1_C0, QKA0      , WEYL_V1_C0)
							QVFXMUL    (WEYL_V1_C1, QKA0      , HALF_V1_C1)
							QVFXXNPMADD(HALF_V1_C1, HALF_V1_C1, QKA0      , WEYL_V1_C1)
							QVFXMUL    (WEYL_V1_C2, QKA0      , HALF_V1_C2)
							QVFXXNPMADD(HALF_V1_C2, HALF_V1_C2, QKA0      , WEYL_V1_C2)
				#endif

							QVSTFDUXA(HALF_V1_C0,TARGETSITE,C32) /* end lifetime: HALF_V1_C0 */
							QVSTFDUXA(HALF_V1_C1,TARGETSITE,C32) /* end lifetime: HALF_V1_C1 */
							QVSTFDUXA(HALF_V1_C2,TARGETSITE,C32) /* end lifetime: HALF_V1_C2 */

							: /* out, inout */
								[gaugesite] "+b" (gaugesite),
								[targetsite] "+b" (targetsite)
							: /* input */
								[spinor_v0_c0] "v" (spinor_v0_c0),
								[spinor_v0_c1] "v" (spinor_v0_c1),
								[spinor_v0_c2] "v" (spinor_v0_c2),
								[spinor_v1_c0] "v" (spinor_v1_c0),
								[spinor_v1_c1] "v" (spinor_v1_c1),
								[spinor_v1_c2] "v" (spinor_v1_c2),
								[spinor_v2_c0] "v" (spinor_v2_c0),
								[spinor_v2_c1] "v" (spinor_v2_c1),
								[spinor_v2_c2] "v" (spinor_v2_c2),
								[spinor_v3_c0] "v" (spinor_v3_c0),
								[spinor_v3_c1] "v" (spinor_v3_c1),
								[spinor_v3_c2] "v" (spinor_v3_c2),
				#if KAMUL
								[qka0] "v" (qka0),
				#endif
								[c32]  "r" (32),
								[c64]  "r" (64),
								[c128] "r" (128)
							: /*clobber*/
								"f" GAUGE_C00,
								"f" GAUGE_C01,
								"f" GAUGE_C02,
								"f" GAUGE_C10,
								"f" GAUGE_C11,
								"f" GAUGE_C12,
								"f" GAUGE_C20,
								"f" GAUGE_C21,
								"f" GAUGE_C22,
								"f" WEYL_V0_C0,
								"f" WEYL_V0_C1,
								"f" WEYL_V0_C2,
								"f" HALF_V0_C0,
								"f" HALF_V0_C1,
								"f" HALF_V0_C2
						);
#else
						bgq_su3_weyl_decl(weyl_ydown);
						bgq_su3_reduce_weyl_yup(weyl_ydown, spinor);

						bgq_su3_mdecl(gauge_ydown);
						bgq_qvlfuxa(gauge_ydown_c00, gaugesite, 32);
						bgq_qvlfuxa(gauge_ydown_c01, gaugesite, 32);
						bgq_qvlfuxa(gauge_ydown_c02, gaugesite, 32);
						bgq_qvlfuxa(gauge_ydown_c10, gaugesite, 32);
						bgq_qvlfuxa(gauge_ydown_c11, gaugesite, 32);
						bgq_qvlfuxa(gauge_ydown_c12, gaugesite, 32);
						bgq_qvlfuxa(gauge_ydown_c20, gaugesite, 32);
						bgq_qvlfuxa(gauge_ydown_c21, gaugesite, 32);
						bgq_qvlfuxa(gauge_ydown_c22, gaugesite, 32);

						bgq_su3_weyl_mvmul(weyl_ydown, gauge_ydown, weyl_ydown);
						if (kamul) {
							bgq_su3_weyl_cmul(weyl_ydown, qka2, weyl_ydown);
						}

						bgq_su3_weyl_zeroload(targetptrs->d[YDOWN]);
						bgq_su3_weyl_store(targetptrs->d[YDOWN], weyl_ydown);
						bgq_weylvec_written(targetptrs->d[YDOWN], t1, t2, x,y,z,YDOWN, true);
#endif
					}

					bgq_su3_matrixnext_prefetch(gaugesite);

					// Z+ /////////////////////////////////////////////////////////////////////////
					{
#if 1
						bgq_weyl_vec *targetsite = targetptrs->d[ZUP];
						asm (
							QVLFDUXA(GAUGE_C00,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C01,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C02,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C10,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C11,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C12,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C20,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C21,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C22,GAUGESITE,C32)

							QVFSUB(WEYL_V0_C0,SPINOR_V0_C0,SPINOR_V2_C0)
							QVFSUB(WEYL_V0_C1,SPINOR_V0_C1,SPINOR_V2_C1)
							QVFSUB(WEYL_V0_C2,SPINOR_V0_C2,SPINOR_V2_C2)

							QVFXMUL(HALF_V0_C0,GAUGE_C00,WEYL_V0_C0)
							QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C0,GAUGE_C00,HALF_V0_C0)
							QVFXMADD(HALF_V0_C0,GAUGE_C10,WEYL_V0_C1,HALF_V0_C0)
							QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C1,GAUGE_C10,HALF_V0_C0)
							QVFXMADD(HALF_V0_C0,GAUGE_C20,WEYL_V0_C2,HALF_V0_C0)
							QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C2,GAUGE_C20,HALF_V0_C0)

							QVFXMUL(HALF_V0_C1,GAUGE_C01,WEYL_V0_C0)
							QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C0,GAUGE_C01,HALF_V0_C1)
							QVFXMADD(HALF_V0_C1,GAUGE_C11,WEYL_V0_C1,HALF_V0_C1)
							QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C1,GAUGE_C11,HALF_V0_C1)
							QVFXMADD(HALF_V0_C1,GAUGE_C21,WEYL_V0_C2,HALF_V0_C1)
							QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C2,GAUGE_C21,HALF_V0_C1)

							QVFXMUL(HALF_V0_C2,GAUGE_C02,WEYL_V0_C0)
							QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C0,GAUGE_C02,HALF_V0_C2) /* end lifetime: WEYL_V0_C0 */
							QVFXMADD(HALF_V0_C2,GAUGE_C12,WEYL_V0_C1,HALF_V0_C2)
							QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C1,GAUGE_C12,HALF_V0_C2) /* end lifetime: WEYL_V0_C1 */
							QVFXMADD(HALF_V0_C2,GAUGE_C22,WEYL_V0_C2,HALF_V0_C2)
							QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C2,GAUGE_C22,HALF_V0_C2) /* end lifetime: WEYL_V0_C2 */

				#if KAMUL
							QVFXMUL    (WEYL_V0_C0, QKA0      , HALF_V0_C0)
							QVFXXNPMADD(HALF_V0_C0, HALF_V0_C0, QKA0      , WEYL_V0_C0)
							QVFXMUL    (WEYL_V0_C1, QKA0      , HALF_V0_C1)
							QVFXXNPMADD(HALF_V0_C1, HALF_V0_C1, QKA0      , WEYL_V0_C1)
							QVFXMUL    (WEYL_V0_C2, QKA0      , HALF_V0_C2)
							QVFXXNPMADD(HALF_V0_C2, HALF_V0_C2, QKA0      , WEYL_V0_C2)
				#endif

							//DCBZ_0(TARGETSITE)
							//DCBZ(C64,TARGETSITE)
							//DCBZ(C128,TARGETSITE)
							QVSTFDXA_0(HALF_V0_C0,TARGETSITE) /* end lifetime: HALF_V0_C0 */
							QVSTFDUXA(HALF_V0_C1,TARGETSITE,C32) /* end lifetime: HALF_V0_C1 */
							QVSTFDUXA(HALF_V0_C2,TARGETSITE,C32) /* end lifetime: HALF_V0_C2 */

							QVFSUB(WEYL_V1_C0,SPINOR_V1_C0,SPINOR_V3_C0)
							QVFSUB(WEYL_V1_C1,SPINOR_V1_C1,SPINOR_V3_C1)
							QVFSUB(WEYL_V1_C2,SPINOR_V1_C2,SPINOR_V3_C2)

							QVFXMUL(HALF_V1_C0,GAUGE_C00,WEYL_V1_C0)
							QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C0,GAUGE_C00,HALF_V1_C0) /* end lifetime: GAUGE_C00 */
							QVFXMADD(HALF_V1_C0,GAUGE_C10,WEYL_V1_C1,HALF_V1_C0)
							QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C1,GAUGE_C10,HALF_V1_C0) /* end lifetime: GAUGE_C10 */
							QVFXMADD(HALF_V1_C0,GAUGE_C20,WEYL_V1_C2,HALF_V1_C0)
							QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C2,GAUGE_C20,HALF_V1_C0) /* end lifetime: GAUGE_C20 */

							QVFXMUL(HALF_V1_C1,GAUGE_C01,WEYL_V1_C0)
							QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C0,GAUGE_C01,HALF_V1_C1) /* end lifetime: GAUGE_C01 */
							QVFXMADD(HALF_V1_C1,GAUGE_C11,WEYL_V1_C1,HALF_V1_C1)
							QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C1,GAUGE_C11,HALF_V1_C1) /* end lifetime: GAUGE_C11 */
							QVFXMADD(HALF_V1_C1,GAUGE_C21,WEYL_V1_C2,HALF_V1_C1)
							QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C2,GAUGE_C21,HALF_V1_C1) /* end lifetime: GAUGE_C21 */

							QVFXMUL(HALF_V1_C2,GAUGE_C02,WEYL_V1_C0)
							QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C0,GAUGE_C02,HALF_V1_C2) /* end lifetime: GAUGE_C02, WEYL_V1_C0 */
							QVFXMADD(HALF_V1_C2,GAUGE_C12,WEYL_V1_C1,HALF_V1_C2)
							QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C1,GAUGE_C12,HALF_V1_C2) /* end lifetime: GAUGE_C12, WEYL_V1_C1 */
							QVFXMADD(HALF_V1_C2,GAUGE_C22,WEYL_V1_C2,HALF_V1_C2)
							QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C2,GAUGE_C22,HALF_V1_C2) /* end lifetime: GAUGE_C22, WEYL_V1_C2 */

				#if KAMUL
							QVFXMUL    (WEYL_V1_C0, QKA0      , HALF_V1_C0)
							QVFXXNPMADD(HALF_V1_C0, HALF_V1_C0, QKA0      , WEYL_V1_C0)
							QVFXMUL    (WEYL_V1_C1, QKA0      , HALF_V1_C1)
							QVFXXNPMADD(HALF_V1_C1, HALF_V1_C1, QKA0      , WEYL_V1_C1)
							QVFXMUL    (WEYL_V1_C2, QKA0      , HALF_V1_C2)
							QVFXXNPMADD(HALF_V1_C2, HALF_V1_C2, QKA0      , WEYL_V1_C2)
				#endif

							QVSTFDUXA(HALF_V1_C0,TARGETSITE,C32) /* end lifetime: HALF_V1_C0 */
							QVSTFDUXA(HALF_V1_C1,TARGETSITE,C32) /* end lifetime: HALF_V1_C1 */
							QVSTFDUXA(HALF_V1_C2,TARGETSITE,C32) /* end lifetime: HALF_V1_C2 */

							: /* out, inout */
								[gaugesite] "+b" (gaugesite),
								[targetsite] "+b" (targetsite)
							: /* input */
								[spinor_v0_c0] "v" (spinor_v0_c0),
								[spinor_v0_c1] "v" (spinor_v0_c1),
								[spinor_v0_c2] "v" (spinor_v0_c2),
								[spinor_v1_c0] "v" (spinor_v1_c0),
								[spinor_v1_c1] "v" (spinor_v1_c1),
								[spinor_v1_c2] "v" (spinor_v1_c2),
								[spinor_v2_c0] "v" (spinor_v2_c0),
								[spinor_v2_c1] "v" (spinor_v2_c1),
								[spinor_v2_c2] "v" (spinor_v2_c2),
								[spinor_v3_c0] "v" (spinor_v3_c0),
								[spinor_v3_c1] "v" (spinor_v3_c1),
								[spinor_v3_c2] "v" (spinor_v3_c2),
				#if KAMUL
								[qka0] "v" (qka0),
				#endif
								[c32]  "r" (32),
								[c64]  "r" (64),
								[c128] "r" (128)
							: /*clobber*/
								"f" GAUGE_C00,
								"f" GAUGE_C01,
								"f" GAUGE_C02,
								"f" GAUGE_C10,
								"f" GAUGE_C11,
								"f" GAUGE_C12,
								"f" GAUGE_C20,
								"f" GAUGE_C21,
								"f" GAUGE_C22,
								"f" WEYL_V0_C0,
								"f" WEYL_V0_C1,
								"f" WEYL_V0_C2,
								"f" HALF_V0_C0,
								"f" HALF_V0_C1,
								"f" HALF_V0_C2
						);
#else
						bgq_su3_weyl_decl(weyl_zup);
						bgq_su3_reduce_weyl_zdown(weyl_zup, spinor);

						bgq_su3_mdecl(gauge_zup);
						bgq_qvlfuxa(gauge_zup_c00, gaugesite, 32);
						bgq_qvlfuxa(gauge_zup_c01, gaugesite, 32);
						bgq_qvlfuxa(gauge_zup_c02, gaugesite, 32);
						bgq_qvlfuxa(gauge_zup_c10, gaugesite, 32);
						bgq_qvlfuxa(gauge_zup_c11, gaugesite, 32);
						bgq_qvlfuxa(gauge_zup_c12, gaugesite, 32);
						bgq_qvlfuxa(gauge_zup_c20, gaugesite, 32);
						bgq_qvlfuxa(gauge_zup_c21, gaugesite, 32);
						bgq_qvlfuxa(gauge_zup_c22, gaugesite, 32);

						bgq_su3_weyl_mvinvmul(weyl_zup, gauge_zup, weyl_zup);
						if (kamul) {
							bgq_su3_weyl_cmul(weyl_zup, qka3, weyl_zup);
						}

						bgq_su3_weyl_zeroload(targetptrs->d[ZUP]);
						bgq_su3_weyl_store(targetptrs->d[ZUP], weyl_zup);
						bgq_weylvec_written(targetptrs->d[ZUP], t1, t2, x,y,z,ZUP, true);
#endif
					}

					//bgq_su3_matrixnext_prefetch(gaugesite);

					// Z- /////////////////////////////////////////////////////////////////////////
					{
#if 1
						bgq_weyl_vec *targetsite = targetptrs->d[ZDOWN];
						asm (
							QVLFDUXA(GAUGE_C00,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C01,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C02,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C10,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C11,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C12,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C20,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C21,GAUGESITE,C32)
							QVLFDUXA(GAUGE_C22,GAUGESITE,C32)

							QVFSUB(WEYL_V0_C0,SPINOR_V0_C0,SPINOR_V2_C0)
							QVFSUB(WEYL_V0_C1,SPINOR_V0_C1,SPINOR_V2_C1)
							QVFSUB(WEYL_V0_C2,SPINOR_V0_C2,SPINOR_V2_C2)

							QVFXMUL(HALF_V0_C0,GAUGE_C00,WEYL_V0_C0)
							QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C0,GAUGE_C00,HALF_V0_C0)
							QVFXMADD(HALF_V0_C0,GAUGE_C10,WEYL_V0_C1,HALF_V0_C0)
							QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C1,GAUGE_C10,HALF_V0_C0)
							QVFXMADD(HALF_V0_C0,GAUGE_C20,WEYL_V0_C2,HALF_V0_C0)
							QVFXXCPNMADD(HALF_V0_C0,WEYL_V0_C2,GAUGE_C20,HALF_V0_C0)

							QVFXMUL(HALF_V0_C1,GAUGE_C01,WEYL_V0_C0)
							QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C0,GAUGE_C01,HALF_V0_C1)
							QVFXMADD(HALF_V0_C1,GAUGE_C11,WEYL_V0_C1,HALF_V0_C1)
							QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C1,GAUGE_C11,HALF_V0_C1)
							QVFXMADD(HALF_V0_C1,GAUGE_C21,WEYL_V0_C2,HALF_V0_C1)
							QVFXXCPNMADD(HALF_V0_C1,WEYL_V0_C2,GAUGE_C21,HALF_V0_C1)

							QVFXMUL(HALF_V0_C2,GAUGE_C02,WEYL_V0_C0)
							QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C0,GAUGE_C02,HALF_V0_C2) /* end lifetime: WEYL_V0_C0 */
							QVFXMADD(HALF_V0_C2,GAUGE_C12,WEYL_V0_C1,HALF_V0_C2)
							QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C1,GAUGE_C12,HALF_V0_C2) /* end lifetime: WEYL_V0_C1 */
							QVFXMADD(HALF_V0_C2,GAUGE_C22,WEYL_V0_C2,HALF_V0_C2)
							QVFXXCPNMADD(HALF_V0_C2,WEYL_V0_C2,GAUGE_C22,HALF_V0_C2) /* end lifetime: WEYL_V0_C2 */

				#if KAMUL
							QVFXMUL    (WEYL_V0_C0, QKA0      , HALF_V0_C0)
							QVFXXNPMADD(HALF_V0_C0, HALF_V0_C0, QKA0      , WEYL_V0_C0)
							QVFXMUL    (WEYL_V0_C1, QKA0      , HALF_V0_C1)
							QVFXXNPMADD(HALF_V0_C1, HALF_V0_C1, QKA0      , WEYL_V0_C1)
							QVFXMUL    (WEYL_V0_C2, QKA0      , HALF_V0_C2)
							QVFXXNPMADD(HALF_V0_C2, HALF_V0_C2, QKA0      , WEYL_V0_C2)
				#endif

							//DCBZ_0(TARGETSITE)
							//DCBZ(C64,TARGETSITE)
							//DCBZ(C128,TARGETSITE)
							QVSTFDXA_0(HALF_V0_C0,TARGETSITE) /* end lifetime: HALF_V0_C0 */
							QVSTFDUXA(HALF_V0_C1,TARGETSITE,C32) /* end lifetime: HALF_V0_C1 */
							QVSTFDUXA(HALF_V0_C2,TARGETSITE,C32) /* end lifetime: HALF_V0_C2 */

							QVFSUB(WEYL_V1_C0,SPINOR_V1_C0,SPINOR_V3_C0)
							QVFSUB(WEYL_V1_C1,SPINOR_V1_C1,SPINOR_V3_C1)
							QVFSUB(WEYL_V1_C2,SPINOR_V1_C2,SPINOR_V3_C2)

							QVFXMUL(HALF_V1_C0,GAUGE_C00,WEYL_V1_C0)
							QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C0,GAUGE_C00,HALF_V1_C0) /* end lifetime: GAUGE_C00 */
							QVFXMADD(HALF_V1_C0,GAUGE_C10,WEYL_V1_C1,HALF_V1_C0)
							QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C1,GAUGE_C10,HALF_V1_C0) /* end lifetime: GAUGE_C10 */
							QVFXMADD(HALF_V1_C0,GAUGE_C20,WEYL_V1_C2,HALF_V1_C0)
							QVFXXCPNMADD(HALF_V1_C0,WEYL_V1_C2,GAUGE_C20,HALF_V1_C0) /* end lifetime: GAUGE_C20 */

							QVFXMUL(HALF_V1_C1,GAUGE_C01,WEYL_V1_C0)
							QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C0,GAUGE_C01,HALF_V1_C1) /* end lifetime: GAUGE_C01 */
							QVFXMADD(HALF_V1_C1,GAUGE_C11,WEYL_V1_C1,HALF_V1_C1)
							QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C1,GAUGE_C11,HALF_V1_C1) /* end lifetime: GAUGE_C11 */
							QVFXMADD(HALF_V1_C1,GAUGE_C21,WEYL_V1_C2,HALF_V1_C1)
							QVFXXCPNMADD(HALF_V1_C1,WEYL_V1_C2,GAUGE_C21,HALF_V1_C1) /* end lifetime: GAUGE_C21 */

							QVFXMUL(HALF_V1_C2,GAUGE_C02,WEYL_V1_C0)
							QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C0,GAUGE_C02,HALF_V1_C2) /* end lifetime: GAUGE_C02, WEYL_V1_C0 */
							QVFXMADD(HALF_V1_C2,GAUGE_C12,WEYL_V1_C1,HALF_V1_C2)
							QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C1,GAUGE_C12,HALF_V1_C2) /* end lifetime: GAUGE_C12, WEYL_V1_C1 */
							QVFXMADD(HALF_V1_C2,GAUGE_C22,WEYL_V1_C2,HALF_V1_C2)
							QVFXXCPNMADD(HALF_V1_C2,WEYL_V1_C2,GAUGE_C22,HALF_V1_C2) /* end lifetime: GAUGE_C22, WEYL_V1_C2 */

				#if KAMUL
							QVFXMUL    (WEYL_V1_C0, QKA0      , HALF_V1_C0)
							QVFXXNPMADD(HALF_V1_C0, HALF_V1_C0, QKA0      , WEYL_V1_C0)
							QVFXMUL    (WEYL_V1_C1, QKA0      , HALF_V1_C1)
							QVFXXNPMADD(HALF_V1_C1, HALF_V1_C1, QKA0      , WEYL_V1_C1)
							QVFXMUL    (WEYL_V1_C2, QKA0      , HALF_V1_C2)
							QVFXXNPMADD(HALF_V1_C2, HALF_V1_C2, QKA0      , WEYL_V1_C2)
				#endif

							QVSTFDUXA(HALF_V1_C0,TARGETSITE,C32) /* end lifetime: HALF_V1_C0 */
							QVSTFDUXA(HALF_V1_C1,TARGETSITE,C32) /* end lifetime: HALF_V1_C1 */
							QVSTFDUXA(HALF_V1_C2,TARGETSITE,C32) /* end lifetime: HALF_V1_C2 */

							: /* out, inout */
								[gaugesite] "+b" (gaugesite),
								[targetsite] "+b" (targetsite)
							: /* input */
								[spinor_v0_c0] "v" (spinor_v0_c0),
								[spinor_v0_c1] "v" (spinor_v0_c1),
								[spinor_v0_c2] "v" (spinor_v0_c2),
								[spinor_v1_c0] "v" (spinor_v1_c0),
								[spinor_v1_c1] "v" (spinor_v1_c1),
								[spinor_v1_c2] "v" (spinor_v1_c2),
								[spinor_v2_c0] "v" (spinor_v2_c0),
								[spinor_v2_c1] "v" (spinor_v2_c1),
								[spinor_v2_c2] "v" (spinor_v2_c2),
								[spinor_v3_c0] "v" (spinor_v3_c0),
								[spinor_v3_c1] "v" (spinor_v3_c1),
								[spinor_v3_c2] "v" (spinor_v3_c2),
				#if KAMUL
								[qka0] "v" (qka0),
				#endif
								[c32]  "r" (32),
								[c64]  "r" (64),
								[c128] "r" (128)
							: /*clobber*/
								"f" GAUGE_C00,
								"f" GAUGE_C01,
								"f" GAUGE_C02,
								"f" GAUGE_C10,
								"f" GAUGE_C11,
								"f" GAUGE_C12,
								"f" GAUGE_C20,
								"f" GAUGE_C21,
								"f" GAUGE_C22,
								"f" WEYL_V0_C0,
								"f" WEYL_V0_C1,
								"f" WEYL_V0_C2,
								"f" HALF_V0_C0,
								"f" HALF_V0_C1,
								"f" HALF_V0_C2
						);
#else
						bgq_su3_weyl_decl(weyl_zdown);
						bgq_su3_reduce_weyl_zup(weyl_zdown, spinor);

						bgq_su3_mdecl(gauge_zdown);
						bgq_qvlfuxa(gauge_zdown_c00, gaugesite, 32);
						bgq_qvlfuxa(gauge_zdown_c01, gaugesite, 32);
						bgq_qvlfuxa(gauge_zdown_c02, gaugesite, 32);
						bgq_qvlfuxa(gauge_zdown_c10, gaugesite, 32);
						bgq_qvlfuxa(gauge_zdown_c11, gaugesite, 32);
						bgq_qvlfuxa(gauge_zdown_c12, gaugesite, 32);
						bgq_qvlfuxa(gauge_zdown_c20, gaugesite, 32);
						bgq_qvlfuxa(gauge_zdown_c21, gaugesite, 32);
						bgq_qvlfuxa(gauge_zdown_c22, gaugesite, 32);

						bgq_su3_weyl_mvmul(weyl_zdown, gauge_zdown, weyl_zdown);
						if (kamul) {
							bgq_su3_weyl_cmul(weyl_zdown, qka3, weyl_zdown);
						}

						bgq_su3_weyl_zeroload(targetptrs->d[ZDOWN]);
						bgq_su3_weyl_store(targetptrs->d[ZDOWN], weyl_zdown);
						bgq_weylvec_written(targetptrs->d[ZDOWN], t1, t2, x,y,z,ZDOWN, true);
#endif
					}
}

#undef KAMUL
#undef IFKAMUL
#undef BGQ_COMPUTEWEYL_INC_
