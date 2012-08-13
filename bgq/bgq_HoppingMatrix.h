/*
 * bgq_HoppingMatrix.h
 *
 *  Created on: Aug 1, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_HOPPINGMATRIX_H_
#define BGQ_HOPPINGMATRIX_H_

////////////////////////////////////////////////////////////////////////////////
// Weylfields

void bgq_hm_init();
void bgq_hm_free();

bool assert_weylval_t(bgq_weylfield_double weylfield, bool isOdd, int t, int x, int y, int z, int xv, int k, int v, int c, bool isRead, bool isWrite);
bool assert_weylfield_t(bgq_weylfield_double weylfield, bool isOdd, int t, int x, int y, int z, int xv, int k, bool isRead, bool isWrite);
#define BGQ_WEYLSITE_T(weylfield, isOdd, t, xv, y, z, x1, x2, isRead, isWrite)  \
	(assert(assert_weylfield_t(weylfield,isOdd,t,x1,y,z,xv,0, isRead, isWrite)), \
	 assert(assert_weylfield_t(weylfield,isOdd,t,x2,y,z,xv,1, isRead, isWrite)), \
	 &weylfield[((xv)*PHYSICAL_LY + (y))*PHYSICAL_LZ + (z)])

#define BGQ_WEYLSITE_T_LEFT(weylfield, isOdd, t, xv, y, z, x1, x2, isRead, isWrite)  \
	(assert(assert_weylfield_t(weylfield,isOdd,t,x1,y,z,xv,0, isRead, isWrite)), \
	 &weylfield[((xv)*PHYSICAL_LY + (y))*PHYSICAL_LZ + (z)])

#define BGQ_WEYLSITE_T_RIGHT(weylfield, isOdd, t, xv, y, z, x1, x2, isRead, isWrite)  \
	(assert(assert_weylfield_t(weylfield,isOdd,t,x2,y,z,xv,1, isRead, isWrite)), \
	 &weylfield[((xv)*PHYSICAL_LY + (y))*PHYSICAL_LZ + (z)])

#define BGQ_WEYLVAL_T(weylfield, isOdd, t, x, y, z, xv, k, v, c, isRead, isWrite)  \
	(assert(assert_weylfield_t(weylfield,isOdd,t,x,y,z,xv,k, isRead, isWrite)), \
	 &weylfield[((xv)*PHYSICAL_LY + (y))*PHYSICAL_LZ + (z)].s[v][c][k])


bool assert_weylval_x(bgq_weylfield_double weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k, int v, int c, bool isRead, bool isWrite);
bool assert_weylfield_x(bgq_weylfield_double weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k, bool isRead, bool isWrite);
#define BGQ_WEYLSITE_X(weylfield, isOdd, tv, x, y, z, t1, t2, isRead, isWrite)  \
	(assert(assert_weylfield_x(weylfield,isOdd,t1,x,y,z,tv,0, isRead, isWrite)), \
	 assert(assert_weylfield_x(weylfield,isOdd,t2,x,y,z,tv,1, isRead, isWrite)), \
	 &weylfield[((tv)*PHYSICAL_LY + (y))*PHYSICAL_LZ + (z)])

#define BGQ_WEYLVAL_X(weylfield, isOdd, t, x, y, z, tv, k, v, c, isRead, isWrite)  \
	(assert(assert_weylfield_x(weylfield,isOdd,t,x,y,z,tv,k, isRead, isWrite)), \
	 &weylfield[((tv)*PHYSICAL_LY + (y))*PHYSICAL_LZ + (z)].s[v][c][k])


bool assert_weylval_y(bgq_weylfield_double weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k, int v, int c, bool isRead, bool isWrite);
bool assert_weylfield_y(bgq_weylfield_double weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k, bool isRead, bool isWrite);
#define BGQ_WEYLSITE_Y(weylfield, isOdd, tv, x, y, z, t1, t2, isRead, isWrite)  \
	(assert(assert_weylfield_y(weylfield,isOdd,t1,x,y,z,tv,0, isRead, isWrite)), \
	 assert(assert_weylfield_y(weylfield,isOdd,t2,x,y,z,tv,1, isRead, isWrite)), \
	 &weylfield[((tv)*PHYSICAL_LX + (x))*PHYSICAL_LZ + (z)])

#define BGQ_WEYLVAL_Y(weylfield, isOdd, t, x, y, z, tv, k, v, c, isRead, isWrite)  \
	(assert(assert_weylfield_y(weylfield,isOdd,t,x,y,z,tv,k, isRead, isWrite)), \
	 &weylfield[((tv)*PHYSICAL_LX + (x))*PHYSICAL_LZ + (z)].s[v][c][k])


#endif /* BGQ_HOPPINGMATRIX_H_ */
