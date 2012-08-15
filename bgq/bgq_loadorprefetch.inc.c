/*
 * bgq_loadorprefetch.inc.c
 *
 *  Created on: Aug 15, 2012
 *      Author: meinersbur
 */

#undef bgq_su3_spinor_double_load_loadorprefetch
#undef bgq_su3_spinor_double_load_left_loadorprefetch
#undef bgq_su3_spinor_double_load_right_loadorprefetch
#undef bgq_su3_weyl_double_load_loadorprefetch
#undef bgq_su3_weyl_double_load_left_loadorprefetch
#undef bgq_su3_matrix_double_load_loadorprefetch


#if BGQ_LOADORPREFETCH_PREFETCH


#define bgq_su3_spinor_double_load_loadorprefetch(dst,addr) \
	bgq_su3_spinor_prefetch(addr)

#define bgq_su3_spinor_double_load_left_loadorprefetch(dst,addr) \
	bgq_su3_spinor_prefetch(addr)

#define bgq_su3_spinor_double_load_right_loadorprefetch(dst,addr) \
	bgq_su3_spinor_prefetch(addr)

#define bgq_su3_weyl_double_load_loadorprefetch(dst,addr) \
	bgq_su3_spinor_prefetch(addr)

#define bgq_su3_weyl_double_load_left_loadorprefetch(dst,addr) \
	bgq_su3_weyl_prefetch(addr)

#define bgq_su3_matrix_double_load_loadorprefetch(dst,addr) \
	bgq_su3_weyl_prefetch(addr)


#elif BGQ_LOADORPREFETCH_LOAD


#define bgq_su3_spinor_double_load_loadorprefetch(dst,addr) \
	bgq_su3_spinor_double_load(dst,addr);

#define bgq_su3_spinor_double_load_left_loadorprefetch(dst,addr) \
	bgq_su3_spinor_double_load_left(dst,addr)

#define bgq_su3_spinor_double_load_right_loadorprefetch(dst,addr) \
	bgq_su3_spinor_double_load_right(dst,addr)

#define bgq_su3_weyl_double_load_loadorprefetch(dst,addr) \
	bgq_su3_weyl_double_load(dst,addr)

#define bgq_su3_weyl_double_load_left_loadorprefetch(dst,addr) \
	bgq_su3_weyl_double_load_left(dst,addr)

#define bgq_su3_matrix_double_load_loadorprefetch(dst,addr) \
	bgq_su3_matrix_double_load(dst,addr)


#else
// Just undef everything
#endif


#undef BGQ_LOADORPREFETCH_PREFETCH
#undef BGQ_LOADORPREFETCH_LOAD
