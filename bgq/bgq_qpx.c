/*
 * bgq_qpx.c
 *
 *  Created on: Oct 26, 2012
 *      Author: meinersbur
 */

#define BGQ_QPX_C_
#include "bgq_qpx.h"

#include <stdbool.h>
#include <stdint.h>

#define QPX_EXPECT(q, r0, i0, r1, i1) \
	do { \
		if (((int)bgq_elem0(q) != (int)r0) || ((int)bgq_elem1(q) != (int)i0) || ((int)bgq_elem2(q) != (int)r1) || ((int)bgq_elem3(q) != (int)r1)) {  \
			failed = true; \
			fprintf(stderr, "QPX fail: (%d,%d,%d,%d) found, (%d,%d,%d,%d) expected at %s:%d\n",(int)bgq_elem0(q),(int)bgq_elem1(q),(int)bgq_elem2(q),(int)bgq_elem3(q),(int)(r0),(int)(i0),(int)(r1),(int)(i1), __FILE__,__LINE__); \
		} \
	} while (0)


void bgq_qpx_unittest(void) {
	printf("Running QPX unittest...\n");

	bool failed = false;
	complexdouble data[2*3*3] = {
			1 + 2*_Complex_I,
			3 + 4*_Complex_I,
			5 + 6*_Complex_I,
			7 + 8*_Complex_I,
			9 + 10*_Complex_I,
			11 + 12*_Complex_I,
			13 + 14*_Complex_I,
			15 + 16*_Complex_I,
			17 + 18*_Complex_I,
			19 + 20*_Complex_I,
			21 + 22*_Complex_I,
			23 + 24*_Complex_I,
			25 + 26*_Complex_I,
			27 + 28*_Complex_I,
			29 + 39*_Complex_I,
			31 + 32*_Complex_I,
			33 + 34*_Complex_I,
			35 + 36*_Complex_I
	};

	bgq_vector4double_decl(v0);
	bgq_lda_double(v0, 0, &data);
	QPX_EXPECT(v0,1,2,3,4);

	bgq_vector4double_decl(v1);
	bgq_lda_double(v1, 32, &data);
	QPX_EXPECT(v1,5,6,7,8);

	bgq_vector4double_decl(v3);
	bgq_ld2a_double(v3, 0, &data);
	QPX_EXPECT(v3,1,2,1,2);

	bgq_vector4double_decl(vmerge2);
	bgq_merge2(vmerge2, v0, v1);
	QPX_EXPECT(vmerge2,3,4,5,6);

	bgq_vector4double_decl(vsplat);
	bgq_complxval_splat(vsplat, 1 + 2*_Complex_I);
	QPX_EXPECT(vsplat,1,2,1,2);

	{
		bgq_vector4double_decl(qvlfdxa);
		const uint8_t *addr =  (void*)&data[0];
		bgq_qvlfdxa(qvlfdxa, addr, 32);
		QPX_EXPECT(qvlfdxa,5,6,7,8);
	}

	{
		bgq_vector4double_decl(qvlfduxa);
		uint8_t *addr =  (void*)&data[0];
		bgq_qvlfdxa(qvlfduxa, addr, 32);
		QPX_EXPECT(qvlfduxa,5,6,7,8);
		if ((void*)addr != (void*)&data[2]) {
			failed = true;
			fprintf(stderr, "QPX fail: (%lu) found, (%lu) expected at %s:%d\n", (uintptr_t)addr, (uintptr_t)(&data[2]), __FILE__, __LINE__);
		}
	}

	bgq_su3_weyl_decl(weyl);
	bgq_su3_weyl_load_double(weyl, &data);
	QPX_EXPECT(weyl_v0_c0,1,2,3,4);
	QPX_EXPECT(weyl_v0_c1,5,6,7,8);
	QPX_EXPECT(weyl_v0_c2,9,10,11,12);
	QPX_EXPECT(weyl_v1_c0,13,14,15,16);
	QPX_EXPECT(weyl_v1_c1,17,18,19,20);
	QPX_EXPECT(weyl_v1_c2,21,22,23,24);



	double store[24] = {0};
	bgq_su3_weyl_store_double(&store,weyl);
	for (unsigned i = 0; i < lengthof(store);i+=1) {
		if ((int)store[i] != (int)(i+1)) {
			failed = true;
			fprintf(stderr, "QPX fail i=%d: (%d) found, (%d) expected at %s:%d\n", i, (int)store[i], (int)(i+1), __FILE__, __LINE__);
		}
	}

	bgq_su3_mdecl(matrix);
	bgq_su3_matrix_load_double(matrix, &data);
	QPX_EXPECT(matrix_c00,1,2,3,4);
	QPX_EXPECT(matrix_c01,5,6,7,8);
	QPX_EXPECT(matrix_c02,9,10,11,12);
	QPX_EXPECT(matrix_c10,13,14,15,16);
	QPX_EXPECT(matrix_c11,17,18,19,20);
	QPX_EXPECT(matrix_c12,21,22,23,24);
	QPX_EXPECT(matrix_c20,25,26,27,28);
	QPX_EXPECT(matrix_c21,29,30,31,32);
	QPX_EXPECT(matrix_c22,33,34,35,36);


	if (!failed) {
		printf("QPX unittest passed\n");
	} else {
		//abort();
	}
}
