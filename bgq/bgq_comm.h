#ifndef BGQ_COMM_H_
#define BGQ_COMM_H_

#include "bgq_spinorfield.h"
#include "bgq_field.h"
#include "bgq_utils.h"

#include <stdint.h>

#ifndef BGQ_COMM_C_
#define EXTERN_INLINE EXTERN_INLINE_DECLARATION
#define EXTERN_FIELD extern
#define EXTERN_INIT(val)
#else
#define EXTERN_INLINE EXTERN_INLINE_DEFINITION
#define EXTERN_FIELD
#define EXTERN_INIT(val) = (val)
#endif


void bgq_comm_mpi_init(void);
void bgq_comm_spi_init(void);

//TODO: inline?
void bgq_comm_recv(bool nospi, bool sloppy, bgq_weylfield_controlblock *targetfield);
void bgq_comm_send(void);
void bgq_comm_wait(void);





EXTERN_FIELD uint8_t *g_bgq_sec_comm;
EXTERN_FIELD uint8_t *g_bgq_sec_comm_float;
EXTERN_FIELD bgq_weyl_vec_double *g_bgq_sec_recv_double[PHYSICAL_LD];
EXTERN_FIELD bgq_weyl_vec_double *g_bgq_sec_send_double[PHYSICAL_LD];

EXTERN_FIELD bgq_weyl_vec_double *g_bgq_sec_temp_tup_double;
EXTERN_FIELD bgq_weyl_vec_double *g_bgq_sec_temp_tdown_double;


EXTERN_FIELD bgq_weyl_vec_float *g_bgq_sec_recv_float[PHYSICAL_LD];
EXTERN_FIELD bgq_weyl_vec_float *g_bgq_sec_send_float[PHYSICAL_LD];

EXTERN_FIELD bgq_weyl_vec_float *g_bgq_sec_temp_tup_float;
EXTERN_FIELD bgq_weyl_vec_float *g_bgq_sec_temp_tdown_float;


#define g_bgq_sec_recv NAME2(g_bgq_sec_recv,PRECISION)
#define g_bgq_sec_send NAME2(g_bgq_sec_send,PRECISION)

#define g_bgq_sec_temp_tup NAME2(g_bgq_sec_temp_tup,PRECISION)
#define g_bgq_sec_temp_tdown NAME2(g_bgq_sec_temp_tdown,PRECISION)

EXTERN_FIELD bgq_weyl_nonvec_double *g_bgq_sec_recv_unvectorized_double[2];
EXTERN_FIELD bgq_weyl_nonvec_float *g_bgq_sec_recv_unvectorized_float[2];
#define g_bgq_sec_recv_unvectorized NAME2(g_bgq_sec_recv_unvectorized,PRECISION)

EXTERN_FIELD bgq_weyl_nonvec_double *g_bgq_sec_send_unvectorized_double[2];
EXTERN_FIELD bgq_weyl_nonvec_float *g_bgq_sec_send_unvectorized_float[2];
#define g_bgq_sec_send_unvectorized NAME2(g_bgq_sec_send_unvectorized,PRECISION)

typedef struct {
	bgq_weyl_nonvec_double *d_tup; // d_dst !!!
	bgq_weyl_nonvec_double *d_tdown;
	bgq_weyl_vec_double *d_xup;
	bgq_weyl_vec_double *d_xdown;
	bgq_weyl_vec_double *d_yup;
	bgq_weyl_vec_double *d_ydown;
	bgq_weyl_vec_double *d_zup;
	bgq_weyl_vec_double *d_zdown;
} bgq_weylsiteptr_unvectorized_double;

typedef struct {
	bgq_weyl_nonvec_float *d_tup; // d_dst !!!
	bgq_weyl_nonvec_float *d_tdown;
	bgq_weyl_vec_float *d_xup;
	bgq_weyl_vec_float *d_xdown;
	bgq_weyl_vec_float *d_yup;
	bgq_weyl_vec_float *d_ydown;
	bgq_weyl_vec_float *d_zup;
	bgq_weyl_vec_float *d_zdown;
} bgq_weylsiteptr_unvectorized_float;


//EXTERN_FIELD size_t (*g_bgq_comm_collapsed2weylsendidx[PHYSICAL_LP])[PHYSICAL_LD];
//EXTERN_FIELD size_t (*g_bgq_comm_collapsed2weylrecvidx[PHYSICAL_LP])[PHYSICAL_LD];
EXTERN_FIELD bgq_weylsiteptr_unvectorized_double *g_comm_collapsed2sendbufptr_unvectorized_double[PHYSICAL_LP]; //TODO: obsolete
EXTERN_FIELD bgq_weylsiteptr_unvectorized_float *g_comm_collapsed2sendbufptr_unvectorized_float[PHYSICAL_LP];//TODO: obsolete
#define g_comm_collapsed2sendbufptr_unvectorized NAME2(g_comm_collapsed2sendbufptr_unvectorized,PRECISION)

EXTERN_FIELD bgq_weylsiteptr_unvectorized_double *g_comm_collapsed2recvbufptr_unvectorized_double[PHYSICAL_LP];//TODO: obsolete
EXTERN_FIELD bgq_weylsiteptr_unvectorized_float *g_comm_collapsed2recvbufptr_unvectorized_float[PHYSICAL_LP];//TODO: obsolete
#define g_comm_collapsed2recvbufptr_unvectorized NAME2(g_comm_collapsed2recvbufptr_unvectorized,PRECISION)



#if 0
#define BGQ_COMM_INC_H_ 1
#define PRECISION double
#include "bgq_comm.inc.h"
#undef PRECISION

#define BGQ_COMM_INC_H_ 1
#define PRECISION float
#include "bgq_comm.inc.h"
#undef PRECISION
#endif


#undef EXTERN_INLINE
#undef EXTERN_FIELD
#undef EXTERN_INIT

#endif /* BGQ_COMM_H_ */
