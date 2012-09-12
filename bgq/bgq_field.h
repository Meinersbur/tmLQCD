#ifndef BGQ_FIELD_H_
#define BGQ_FIELD_H_

#include "bgq_utils.h"
#include "complex_c99.h"
#include "../global.h"
#include "../read_input.h"
#include <stdbool.h>
#include <assert.h>


#ifndef BGQ_FIELD_C_H_
#define EXTERN_INLINE inline
#define EXTERN_FIELD extern
#else
#define EXTERN_INLINE extern inline
#define EXTERN_FIELD
#endif


#define LOCAL_LP 2 /* Even/Odd */
#define LOCAL_LT T
#define LOCAL_LX LX
#define LOCAL_LY LY
#define LOCAL_LZ LZ

#define PHYSICAL_LP 2 /* Even/Odd */
#define PHYSICAL_LTV (LOCAL_LT/(PHYSICAL_LP*PHYSICAL_LK))
#define PHYSICAL_LX LOCAL_LX
#define PHYSICAL_LXV (LOCAL_LX/(PHYSICAL_LP*PHYSICAL_LK))
#define PHYSICAL_LY LOCAL_LY
#define PHYSICAL_LZ LOCAL_LZ
#define PHYSICAL_LK 2 /* Vector unit width (2 complex = 4 reals) */
#define PHYSICAL_LD 5 /* No of directions (for gauge field) = 4 up-directions + 1 for shifted t-line (note that it is one vector-element longer for wraparound values) */



// faces, 2 of each
#define SURFACE_T_SITES ((1)*(PHYSICAL_LX-2)*(PHYSICAL_LY-2)*(PHYSICAL_LZ))
#define SURFACE_X_SITES ((PHYSICAL_LTV-2)*(1)*(PHYSICAL_LY-2)*(PHYSICAL_LZ))
#define SURFACE_Y_SITES ((PHYSICAL_LTV-2)*(PHYSICAL_LX-2)*(1)*(PHYSICAL_LZ))
#define SURFACE_FACE_SITES (2*SURFACE_T_SITES + 2*SURFACE_X_SITES + 2*SURFACE_Y_SITES)

#define SURFACE_T_ZLINES ((1)*(PHYSICAL_LX-2)*(PHYSICAL_LY-2))
#define SURFACE_X_ZLINES ((PHYSICAL_LTV-2)*(1)*(PHYSICAL_LY-2))
#define SURFACE_Y_ZLINES ((PHYSICAL_LTV-2)*(PHYSICAL_LX-2)*(1))
#define SURFACE_FACE_ZLINES (2*SURFACE_T_ZLINES + 2*SURFACE_X_ZLINES + 2*SURFACE_Y_ZLINES)

// edges, 4 of each
#define SURFACE_TX_SITES ((1)*(1)*(PHYSICAL_LY-2)*(PHYSICAL_LZ))
#define SURFACE_TY_SITES ((1)*(PHYSICAL_LX-2)*(1)*(PHYSICAL_LZ))
#define SURFACE_XY_SITES ((PHYSICAL_LTV-2)*(1)*(1)*(PHYSICAL_LZ))
#define SURFACE_EDGE_SITES (4*SURFACE_TX_SITES + 4*SURFACE_TY_SITES + 4*SURFACE_XY_SITES)

#define SURFACE_TX_ZLINES ((1)*(1)*(PHYSICAL_LY-2))
#define SURFACE_TY_ZLINES ((1)*(PHYSICAL_LX-2)*(1))
#define SURFACE_XY_ZLINES ((PHYSICAL_LTV-2)*(1)*(1))
#define SURFACE_EDGE_ZLINES (4*SURFACE_TX_ZLINES + 4*SURFACE_TY_ZLINES + 4*SURFACE_XY_ZLINES)

// vertices, 8
#define SURFACE_TXY_SITES ((1)*(1)*(1)*(PHYSICAL_LZ))
#define SURFACE_VERTICE_SITES (8*SURFACE_TXY_SITES)

#define SURFACE_TXY_ZLINES ((1)*(1)*(1))
#define SURFACE_VERTICE_ZLINES (8*SURFACE_TXY_ZLINES)

// total
#define SURFACE_SITES (SURFACE_FACE_SITES+SURFACE_FACE_SITES+SURFACE_EDGE_SITES)
#define BODY_SITES ((PHYSICAL_LTV-2)*(PHYSICAL_LX-2)*(PHYSICAL_LY-2)*(PHYSICAL_LZ))
#define VOLUME_SITES (VOLUME/(PHYSICAL_LP*PHYSICAL_LK))

#define SURFACE_ZLINES (SURFACE_FACE_ZLINES+SURFACE_EDGE_ZLINES+SURFACE_VERTICE_ZLINES)
#define BODY_ZLINES ((PHYSICAL_LTV-2)*(PHYSICAL_LX-2)*(PHYSICAL_LY-2))


#define GAUGE_VOLUME ((LOCAL_LT+1)*(LOCAL_LX+1)*(LOCAL_LY+1)*(LOCAL_LZ)) /* LOCAL volume */
#define GAUGE_EOVOLUME ((PHYSICAL_LTV+1)*(PHYSICAL_LX+1)*(PHYSICAL_LY+1)*(PHYSICAL_LZ+1/*for wraparound*/)) /* This is not tight/hole-free; we could also define an accessor function per direction*/


typedef enum {
	TUP = 0,
	TDOWN = 1,
	XUP = 2,
	XDOWN = 3,
	YUP = 4,
	YDOWN = 5,
	ZUP = 6,
	ZDOWN = 7,

	TUP_SHIFT = 8,
	TDOWN_SHIFT = 9,

	DIR_UP = 0,
	DIR_DOWN = 1
} direction;

typedef _Complex double complexdouble;
typedef _Complex float complexfloat;



#define BGQ_SPINORSITE_ACCESS(spinorfield, isOdd, tv, x, y, z) \
	(&(spinorfield)[(((tv)*PHYSICAL_LX + (x))*PHYSICAL_LY + (y))*PHYSICAL_LZ + (z)])

#define BGQ_SPINORSITE(spinorfield, isOdd, tv, x, y, z, t1, t2, isRead, isWrite)                         \
		(assert(assert_spinorcoord(spinorfield, isOdd, t1, x, y, z, tv, 0, isRead, isWrite)), \
		 assert(assert_spinorcoord(spinorfield, isOdd, t2, x, y, z, tv, 1, isRead, isWrite)), \
		 BGQ_SPINORSITE_ACCESS(spinorfield, isOdd, tv, x, y, z))

#define BGQ_SPINORSITE_LEFT(spinorfield, isOdd, tv, x, y, z, t1, t2, isRead, isWrite)                         \
		(assert(assert_spinorcoord(spinorfield, isOdd, t1, x, y, z, tv, 0, isRead, isWrite)), \
		 BGQ_SPINORSITE_ACCESS(spinorfield, isOdd, tv, x, y, z))

#define BGQ_SPINORSITE_RIGHT(spinorfield, isOdd, tv, x, y, z, t1, t2, isRead, isWrite)                         \
		(assert(assert_spinorcoord(spinorfield, isOdd, t2, x, y, z, tv, 1, isRead, isWrite)), \
		 BGQ_SPINORSITE_ACCESS(spinorfield, isOdd, tv, x, y, z))

#define BGQ_SPINORVAL(spinorfield,isOdd,t,x,y,z,tv,k, v,c, isRead,isWrite) \
	(assert(assert_spinorfield_coord(spinorfield,isOdd,t,x,y,z,tv,k,v,c,isRead,isWrite)), \
	 &(BGQ_SPINORSITE_ACCESS(spinorfield,isOdd,tv,x,y,z)->s[v][c][k]))


////////////////////////////////////////////////////////////////////////////////
// Gaugefield

#define BGQ_GAUGESITE_ACCESS(gaugefield,isOdd,tv,x,y,z,dir) \
	(assert(-1 <= tv && tv < PHYSICAL_LTV),        \
	 assert(-1 <= x && x < PHYSICAL_LX),        \
	 assert(-1 <= y && y < PHYSICAL_LY),        \
	 assert(-1 <= z && z < PHYSICAL_LZ),        \
	 &(gaugefield)->eodir[(isOdd)][(dir)/2][((((tv)+1)*(PHYSICAL_LX+1) + ((x)+1))*(PHYSICAL_LY+1) + ((y)+1))*(PHYSICAL_LZ+1) + ((z)+1)])

#define BGQ_GAUGESITE(gaugefield,isOdd,tv,x,y,z,dir,t1,t2,isRead,isWrite)         \
	(assert(assert_gaugesite(gaugefield,isOdd,t1,x,y,z,tv,0,dir,isRead,isWrite)),  \
	 assert(assert_gaugesite(gaugefield,isOdd,t2,x,y,z,tv,1,dir,isRead,isWrite)),  \
	 BGQ_GAUGESITE_ACCESS(gaugefield,isOdd,tv,x,y,z,dir))

#define BGQ_GAUGESITE_LEFT(gaugefield,isOdd,tv,x,y,z,dir,t1,t2,isRead,isWrite)         \
	(assert(assert_gaugesite(gaugefield,isOdd,t1,x,y,z,tv,0,dir,isRead,isWrite)),  \
	 BGQ_GAUGESITE_ACCESS(gaugefield,isOdd,tv,x,y,z,dir))

#define BGQ_GAUGESITE_RIGHT(gaugefield,isOdd,tv,x,y,z,dir,t1,t2,isRead,isWrite)      \
	(assert(assert_gaugesite(gaugefield,isOdd,t2,x,y,z,tv,1,dir,isRead,isWrite)),  \
	 BGQ_GAUGESITE_ACCESS(gaugefield,isOdd,tv,x,y,z,dir))

#define BGQ_GAUGEVAL(gaugefield,isOdd,t,x,y,z,tv,k,dir,i,l, isRead,isWrite) \
	(assert(assert_gaugeval(gaugefield,isOdd,t,x,y,z,tv,k,dir,i,l,isRead,isWrite)), \
	 &(BGQ_GAUGESITE_ACCESS(gaugefield,isOdd,tv,x,y,z,dir)->c[i][l][k]))


////////////////////////////////////////////////////////////////////////////////
// Weylfields

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


#define BGQ_WEYLSITE_X(weylfield, isOdd, tv, x, y, z, t1, t2, isRead, isWrite)  \
	(assert(assert_weylfield_x(weylfield,isOdd,t1,x,y,z,tv,0, isRead, isWrite)), \
	 assert(assert_weylfield_x(weylfield,isOdd,t2,x,y,z,tv,1, isRead, isWrite)), \
	 &weylfield[((tv)*PHYSICAL_LY + (y))*PHYSICAL_LZ + (z)])

#define BGQ_WEYLVAL_X(weylfield, isOdd, t, x, y, z, tv, k, v, c, isRead, isWrite)  \
	(assert(assert_weylfield_x(weylfield,isOdd,t,x,y,z,tv,k, isRead, isWrite)), \
	 &weylfield[((tv)*PHYSICAL_LY + (y))*PHYSICAL_LZ + (z)].s[v][c][k])


#define BGQ_WEYLSITE_Y(weylfield, isOdd, tv, x, y, z, t1, t2, isRead, isWrite)  \
	(assert(assert_weylfield_y(weylfield,isOdd,t1,x,y,z,tv,0, isRead, isWrite)), \
	 assert(assert_weylfield_y(weylfield,isOdd,t2,x,y,z,tv,1, isRead, isWrite)), \
	 &weylfield[((tv)*PHYSICAL_LX + (x))*PHYSICAL_LZ + (z)])

#define BGQ_WEYLVAL_Y(weylfield, isOdd, t, x, y, z, tv, k, v, c, isRead, isWrite)  \
	(assert(assert_weylfield_y(weylfield,isOdd,t,x,y,z,tv,k, isRead, isWrite)), \
	 &weylfield[((tv)*PHYSICAL_LX + (x))*PHYSICAL_LZ + (z)].s[v][c][k])


////////////////////////////////////////////////////////////////////////////////
// Common Initialization

void bgq_init_gaugefield_allprec();
void bgq_free_gaugefield_allprec();
void bgq_init_spinorfields_allprec(int count);
void bgq_free_spinorfields_allprec();
void bgq_hm_init_allprec();
void bgq_hm_free_allprec();

void bgq_update_backward_gauge();

enum {
	BGQREF_TUP,
	BGQREF_TDOWN,
	BGQREF_XUP,
	BGQREF_XUP_WEYLSEND,
	BGQREF_XDOWN,
	BGQREF_XDOWN_WEYLREAD,
	BGQREF_YUP,
	BGQREF_YDOWN,
	BGQREF_ZUP,
	BGQREF_ZDOWN,
	BGQREF_STORE,
	BGQREF_count
};

void bgq_initbgqref();
void bgq_setrefvalue(int t, int x, int y, int z, int idx, complexdouble val, char *desc);
void bgq_setbgqvalue(int t, int x, int y, int z, int idx, complexdouble val, char *desc);
void bgq_savebgqref();





#undef EXTERN_INLINE
#undef EXTERN_FIELD

#endif /* BGQ_FIELD_H_ */
