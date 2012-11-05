
#include "bgq_HoppingMatrix.h"

bool io;
bgq_weylfield_controlblock *tf;
bgq_weylfield_controlblock *sf;
vector4double x;

void callme(int x);




void test(ucoord tid, ucoord threads) {
	const bool isOdd = io;
	const bgq_weylfield_controlblock * restrict targetfield = tf;
	const bgq_weylfield_controlblock * restrict spinorfield = sf;
	bool kamul = false;

#if 1
	bgq_vector4double_decl(qka0);
	bgq_complxval_splat(qka0,ka0);
	bgq_vector4double_decl(qka1);
	bgq_complxval_splat(qka1,ka1);
	bgq_vector4double_decl(qka2);
	bgq_complxval_splat(qka2,ka2);
	bgq_vector4double_decl(qka3);
	bgq_complxval_splat(qka3,ka3);
#endif

	const size_t workload = PHYSICAL_SURFACE;
	const size_t threadload = (workload+threads-1)/threads;
	const size_t begin = tid*threadload;
	const size_t end = min_sizet(workload, begin+threadload);
	for (ucoord is = begin; is<end; is+=1) {//TODO: Check removal
		ucoord ih = bgq_surface2halfvolume(isOdd, is);
		ucoord ic = bgq_surface2collapsed(is);
		ucoord t1 = bgq_halfvolume2t1(isOdd,ih);
		ucoord t2 = bgq_halfvolume2t2(isOdd,ih);
		ucoord x = bgq_halfvolume2x(ih);
		ucoord y = bgq_halfvolume2y(ih);
		ucoord z = bgq_halfvolume2z(ih);


		//TODO: Check strength reduction
		bgq_spinorsite * restrict spinorsite = &spinorfield->sec_fullspinor_surface[is];
		bgq_gaugesite * restrict gaugesite = &g_bgq_gaugefield_fromCollapsed[isOdd][is];
		bgq_weyl_ptr_t * restrict targetptrs = &targetfield->sendptr[ic];

		//TODO: prefetching
		//TODO: Check inlining
		bgq_su3_spinor_decl(spinor);

		bgq_su3_spinor_load_double(spinor, spinorsite);
		bgq_spinorqpx_expect(spinor, t1, t2, x, y, z);
		//bgq_HoppingMatrix_loadFulllayout(spinor, spinorsite, t1, t2, x, y, z);

#define BGQ_COMPUTEWEYL_INC_
//#include "bgq_ComputeWeyl.inc.c"
#undef BGQ_COMPUTEWEYL_INC_
		//bgq_HoppingMatrix_compute_storeWeyllayout(destptrs, gaugesite, spinor, t1, t2, x, y, z, true);

		bgq_HoppingMatrix_compute_storeWeyllayout_alldir(targetptrs, gaugesite, spinor, t1, t2, x, y, z, qka0,qka1,qka2,qka3,kamul);

	}
}


