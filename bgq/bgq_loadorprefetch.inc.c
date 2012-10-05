/*
 * bgq_loadorprefetch.inc.c
 *
 *  Created on: Aug 15, 2012
 *      Author: meinersbur
 */

#undef bgq_su3_spinor_loadorprefetch
#undef bgq_su3_spinor_loadorprefetch_left
#undef bgq_su3_spinor_loadorprefetch_right
#undef bgq_su3_weyl_loadorprefetch
#undef bgq_su3_weyl_loadorprefetch_left
#undef bgq_su3_weyl_loadorprefetch_right
#undef bgq_su3_matrix_loadorprefetch


#ifndef BGQ_LOADORPREFETCH_PREFETCH
#define BGQ_LOADORPREFETCH_PREFETCH 0
#endif

#ifndef BGQ_LOADORPREFETCH_LOAD
#define BGQ_LOADORPREFETCH_LOAD 0
#endif



#if BGQ_LOADORPREFETCH_PREFETCH

#if BGQ_DCBT_DISABLE

#define bgq_su3_spinor_loadorprefetch(dst,addr)
#define bgq_su3_spinor_loadorprefetch_left(dst,addr)
#define bgq_su3_spinor_loadorprefetch_right(dst,addr)
#define bgq_su3_weyl_loadorprefetch(dst,addr)
#define bgq_su3_weyl_loadorprefetch_left(dst,addr)
#define bgq_su3_weyl_loadorprefetch_right(dst,addr)
#define bgq_su3_matrix_loadorprefetch(dst,addr)

#else

#define bgq_su3_spinor_loadorprefetch(dst,addr) \
	bgq_su3_spinor_prefetch(addr)

#define bgq_su3_spinor_loadorprefetch_left(dst,addr) \
	bgq_su3_spinor_prefetch(addr)

#define bgq_su3_spinor_loadorprefetch_right(dst,addr) \
	bgq_su3_spinor_prefetch(addr)

#define bgq_su3_weyl_loadorprefetch(dst,addr) \
	bgq_su3_weyl_prefetch(addr)

#define bgq_su3_weyl_loadorprefetch_left(dst,addr) \
	bgq_su3_weyl_prefetch(addr)

#define bgq_su3_weyl_loadorprefetch_right(dst,addr) \
	bgq_su3_weyl_prefetch(addr)

#define bgq_su3_matrix_loadorprefetch(dst,addr) \
	bgq_su3_matrix_prefetch(addr)

#endif


#elif BGQ_LOADORPREFETCH_LOAD


#define bgq_su3_spinor_loadorprefetch(dst,addr) \
	bgq_su3_spinor_load(dst,addr);

#define bgq_su3_spinor_loadorprefetch_left(dst,addr) \
	bgq_su3_spinor_load_left(dst,addr)

#define bgq_su3_spinor_loadorprefetch_right(dst,addr) \
	bgq_su3_spinor_load_right(dst,addr)

#define bgq_su3_weyl_loadorprefetch(dst,addr) \
	bgq_su3_weyl_load(dst,addr)

#define bgq_su3_weyl_loadorprefetch_left(dst,addr) \
	bgq_su3_weyl_load_left(dst,addr)

#define bgq_su3_weyl_loadorprefetch_right(dst,addr) \
	bgq_su3_weyl_load_right(dst,addr)

#define bgq_su3_matrix_loadorprefetch(dst,addr) \
	bgq_su3_matrix_load(dst,addr)


#else
// Just undef everything
#endif

#undef BGQ_LOADORPREFETCH_PREFETCH
#undef BGQ_LOADORPREFETCH_LOAD
