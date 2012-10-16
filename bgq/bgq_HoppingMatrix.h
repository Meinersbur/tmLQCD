/*
 * bgq_HoppingMatrix.h
 *
 *  Created on: Oct 13, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_HOPPINGMATRIX_H_
#define BGQ_HOPPINGMATRIX_H_

typedef enum {
	hm_nocom = 1 << 0,
	hm_nooverlap = 1 << 1,
	hm_nokamul = 1 << 2,
	hm_fixedoddness = 1 << 3,

	hm_noprefetchexplicit = 1 << 4,
	hm_noprefetchlist = 1 << 5,
	hm_noprefetchstream = 1 << 6,

	hm_noweylsend = 1 << 7,
	hm_nobody = 1 << 8,
	hm_nosurface = 1 << 9,

	hm_l1pnonstoprecord = 1 << 10,
	hm_experimental = 1 << 11,

	hm_prefetchimplicitdisable = 1 << 12,
	hm_prefetchimplicitoptimistic = 2 << 12,
	hm_prefetchimplicitconfirmed = 3 << 12
} bgq_hmflags;

#endif /* BGQ_HOPPINGMATRIX_H_ */
