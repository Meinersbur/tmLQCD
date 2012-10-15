/*
 * bgq_field.h
 *
 *  Created on: Oct 13, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_FIELD_H_
#define BGQ_FIELD_H_

#ifndef BGQ_FIELD_C_H_
#define EXTERN_INLINE __inline__
#define EXTERN_FIELD extern
#else
#define EXTERN_INLINE extern __inline__
#define EXTERN_FIELD
#endif



/* global logical coordinates */

#define GLOBAL_LT T_global
#define GLOBAL_LX (LX*N_PROC_X)
#define GLOBAL_LY (LY*N_PROC_Y)
#define GLOBAL_LZ (LZ*N_PROC_Z)
#define GLOBAL_VOLUME (GLOBAL_LT*GLOBAL_LX*GLOBAL_LY*GLOBAL_LZ)

/* local logical coordinates */

#define LOCAL_LT T
#define LOCAL_LX LX
#define LOCAL_LY LY
#define LOCAL_LZ LZ
#define LOCAL_VOLUME VOLUME

#define LOCAL_HALO_T (LOCAL_LX*LOCAL_LY*LOCAL_LZ)
#define LOCAL_HALO_X (LOCAL_LT*LOCAL_LY*LOCAL_LZ)
#define LOCAL_HALO_Y (LOCAL_LT*LOCAL_LX*LOCAL_LZ)
#define LOCAL_HALO_Z (LOCAL_LT*LOCAL_LX*LOCAL_LZ)


/* local physical coordinates */

#define PHYSICAL_LP 2 /* Even/Odd */
#define PHYSICAL_LTV (LOCAL_LT/(PHYSICAL_LP*PHYSICAL_LK))
#define PHYSICAL_LX LOCAL_LX
#define PHYSICAL_LY LOCAL_LY
#define PHYSICAL_LZ LOCAL_LZ
#define PHYSICAL_LK 2 /* Vector unit width (2 complex = 4 reals) */
#define PHYSICAL_LD 8 /* Number of directions */
#define PHYSICAL_VOLUME (PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LY*PHYSICAL_LZ) /* (VOLUME/(PHYSICAL_LP*PHYSICAL_LK)) */
EXTERN_FIELD size_t PHYSICAL_BODY;
#define PHYSICAL_SURFACE (PHYSICAL_VOLUME-PHYSICAL_BODY)

#define PHYSICAL_INDEX_LEXICAL(isOdd, tv, t1, t2, x, y, z) (tv + PHYSICAL_LTV*(x + PHYSICAL_LX*(y + PHYSICAL_LY*(z))))
#define PHYSICAL_LEXICAL2TV(isOdd,ih) ((ih)%PHYSICAL_LTV)
#define PHYSICAL_LEXICAL2X(isOdd,ih) (((ih)/PHYSICAL_LTV)%PHYSICAL_LX)
#define PHYSICAL_LEXICAL2Y(isOdd,ih) (((ih)/(PHYSICAL_LTV*PHYSICAL_LX))%PHYSICAL_LY)
#define PHYSICAL_LEXICAL2Z(isOdd,ih) ((ih)/(PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LY))
#define PHYSICAL_LEXICAL2T1(isOdd,ih) (4*PHYSICAL_LEXICAL2TV(isOdd,ih) + ((PHYSICAL_LEXICAL2X(isOdd,ih)+PHYSICAL_LEXICAL2Y(isOdd,ih)+PHYSICAL_LEXICAL2Z(isOdd,ih)+isOdd)&1))
#define PHYSICAL_LEXICAL2T2(isOdd,ih) (PHYSICAL_LEXICAL2T1(isOdd,ih) + 2)

#define PHYSICAL_HALO_T (COMM_T*PHYSICAL_LX*PHYSICAL_LY*PHYSICAL_LZ/PHYSICAL_LP)
#define PHYSICAL_HALO_X (COMM_X*PHYSICAL_LTV*PHYSICAL_LY*PHYSICAL_LZ)
#define PHYSICAL_HALO_Y (COMM_Y*PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LZ)
#define PHYSICAL_HALO_Z (COMM_Z*PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LY)



#if defined(PARALLELT)
#define COMM_T 1
#define COMM_X 0
#define COMM_Y 0
#define COMM_Z 0
#elif defined(PARALLELX)
#define COMM_T 0
#define COMM_X 1
#define COMM_Y 0
#define COMM_Z 0
#elif defined(PARALLELXT)
#define COMM_T 1
#define COMM_X 1
#define COMM_Y 0
#define COMM_Z 0
#elif defined(PARALLELXY)
#define COMM_T 0
#define COMM_X 1
#define COMM_Y 1
#define COMM_Z 0
#elif defined(PARALLELXYT)
#define COMM_T 1
#define COMM_X 1
#define COMM_Y 1
#define COMM_Z 0
#elif defined(PARALLELXYZ)
#define COMM_T 0
#define COMM_X 1
#define COMM_Y 1
#define COMM_Z 1
#elif defined(PARALLELXYZT)
#define COMM_T 1
#define COMM_X 1
#define COMM_Y 1
#define COMM_Z 1
#else
#error Unkonwn inter-node parallelization
#endif





typedef _Complex double complex_double;
typedef _Complex float complex_float;
#define COMPLEX_PRECISION complex_double

typedef enum {
	TUP = 0,
	TDOWN = 1,
	XUP = 2,
	XDOWN = 3,
	YUP = 4,
	YDOWN = 5,
	ZUP = 6,
	ZDOWN = 7,

	DIR_UP = 0,
	DIR_DOWN = 1
} bgq_direction;

typedef enum {
	P_TUP1,
	P_TUP2,
	P_TDOWN1,
	P_TDOWN2,
	P_XUP,
	P_XDOWN,
	P_YUP,
	P_YDOWN,
	P_ZUP,
	P_ZDOWN
} bgq_physical_direction;
#define P_COUNT 10

typedef enum {
	DIM_T,
	DIM_X,
	DIM_Y,
	DIM_Z
} bgq_dimension;


inline bgq_physical2direction(bgq_physical_direction pd) {
	switch (pd) {
	case P_TUP1:
	case P_TUP2:
		return TUP;
	case P_TDOWN1:
	case P_TDOWN2:
		return TDOWN;
	case P_XUP:
		retrun XUP;
	case P_XDOWN:
		return XDOWN;
	case P_YUP:
		return YUP;
	case P_YDOWN:
		return YDOWN;
	case P_ZUP:
		return ZUP;
	case P_ZDOWN:
		return ZDOWN;
	}
}

inline bgq_dimension bgq_direction2dimension(bgq_direction d) {
	return d/2;
}


typedef struct {
	COMPLEX_PRECISION s[4][3][PHYSICAL_LK]; /* 4*3*2*sizeof(COMPLEX_PRECISION) = 384;192 bytes (6;3 L1 cache lines) */
} bgq_spinorsite;
typedef bgq_spinorsite (*bgq_spinorfield);

typedef struct {
	COMPLEX_PRECISION s[2][3][PHYSICAL_LK]; // 96 byte (3 L1 cache lines)
} bgq_weyl_vec;

typedef struct {
	COMPLEX_PRECISION s[2][3]; // 48 byte
	COMPLEX_PRECISION _padding[2]; // 32 byte
} bgq_weyl_nonvec; // 64 byte (One L1 cache line)

typedef struct {
	bgq_weyl_nonvec tup1;
	bgq_weyl_nonvec tup2;
	bgq_weyl_nonvec tdown1;
	bgq_weyl_nonvec tdown2;
	bgq_weyl_vec xup;
	bgq_weyl_vec xdown;
	bgq_weyl_vec yup;
	bgq_weyl_vec ydown;
	bgq_weyl_vec zup;
	bgq_weyl_vec zdown;
} bgq_weylsite;
size_t bgq_offsetof_weylsite[bgq_physical_direction] = {offsetof(bgq_weylsite,tup1),offsetof(bgq_weylsite,tup2),offsetof(bgq_weylsite,tdown1),offsetof(bgq_weylsite,tdown2),
		offsetof(bgq_weylsite,xup),offsetof(bgq_weylsite,xdown),offsetof(bgq_weylsite,yup),offsetof(bgq_weylsite,ydown),offsetof(bgq_weylsite,zup),offsetof(bgq_weylsite,zdown)};

typedef struct {
	COMPLEX_PRECISION s[2];
} bgq_weyl;

typedef enum {
	sec_send_tup,
	sec_send_tdown,
	sec_send_xup,
	sec_send_xdown,
	sec_send_yup,
	sec_send_ydown,
	sec_send_zup,
	sec_send_zdown,
	sec_recv_tup,
	sec_recv_tdown,
	sec_recv_xup,
	sec_recv_xdown,
	sec_recv_yup,
	sec_recv_ydown,
	sec_recv_zup,
	sec_recv_zdown,
	sec_surface,
	sec_body,
	sec_end
} bgq_weylfield_section;


EXTERN_INLINE size_t bgq_weyl_section_offset(bgq_weylfield_section section) {

	size_t body_volume;
	if (PHYSICAL_LTV<=2 || PHYSICAL_LX<=2 || PHYSICAL_LY<=2 || PHYSICAL_Z<=2) {
		body_volume = 0;
	} else {
		body_volume = (PHYSICAL_LTV-2)*(PHYSICAL_LX-2)*(PHYSICAL_LY-2)*(PHYSICAL_Z-2);
	}


	size_t result = 0;
	if (section > sec_send_tup) {
		result += PHYSCIAL_HALO_T * sizeof(bgq_weyl_nonvec);
		result = (result + 127) & ~127;
	}
	if (section > sec_send_tdown) {
		result = PHYSCIAL_HALO_T * sizeof(bgq_weyl_nonvec);
		result = (result + 127) & ~127;
	}
	if (section > sec_send_xup) {
		result += PHYSICAL_HALO_X * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_send_xdown) {
		result += PHYSICAL_HALO_X * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_send_yup) {
		result += PHYSICAL_HALO_Y * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_send_ydown) {
		result += PHYSICAL_HALO_Y * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_send_zup) {
		result += PHYSICAL_HALO_Z * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_send_zdown) {
		result += PHYSICAL_HALO_Z * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_recv_tup) {
		result += PHYSICAL_HALO_T * sizeof(bgq_weyl_nonvec);
		result = (result + 127) & ~127;
	}
	if (section > sec_recv_tdown) {
		result = PHYSICAL_HALO_T * sizeof(bgq_weyl_nonvec);
		result = (result + 127) & ~127;
	}
	if (section > sec_recv_xup) {
		result += PHYSICAL_HALO_X * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_recv_xdown) {
		result += PHYSICAL_HALO_X * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_recv_yup) {
		result += PHYSICAL_HALO_Y * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_recv_ydown) {
		result += PHYSICAL_HALO_Y * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_recv_zup) {
		result += PHYSICAL_HALO_Z * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_recv_zdown) {
		result += PHYSICAL_HALO_Z * sizeof(bgq_weyl_vec);
		result = (result + 127) & ~127;
	}
	if (section > sec_surface) {
		result += (PHYSICAL_VOLUME-body_volume) * sizeof(bgq_weylsite);
		result = (result + 127) & ~127;
	}
	if (section > sec_body) {
		result += body_volume * sizeof(bgq_weylsite);
	}
	return result;
}



static inline bool bgq_isOdd(size_t t, size_t x, size_t y, size_t z) {
	return (t+x+y+z)&1;
}


static inline bool bgq_isCommSurface(size_t t, size_t x, size_t y, size_t z) {
#ifdef COMM_T
if ((t==0) || (t==LOCAL_LT-1))
	return true;
#endif
#ifdef COMM_X
if ((x==0) || (x==LOCAL_LX-1))
	return true;
#endif
#ifdef COMM_Y
if ((y==0) || (y==LOCAL_LY-1))
	return true;
#endif
#ifdef COMM_Z
if ((z==0) || (z==LOCAL_LZ-1))
	return true;
#endif
return false;
}

static inline bgq_isCommBody(size_t t, size_t x, size_t y, size_t z) {
	return !bgq_isCommSurface(t,x,y,z);
}


#undef EXTERN_INLINE
#undef EXTERN_FIELD

#endif /* BGQ_FIELD_H_ */


