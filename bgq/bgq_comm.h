#ifndef BGQ_COMM_H_
#define BGQ_COMM_H_

#include "bgq_utils.h"
#include "bgq_spinorfield.h"

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
void bgq_comm_recv(bool nospi);
void bgq_comm_send(bool nospi);
void bgq_comm_wait(bool nospi);

EXTERN_FIELD uint8_t *g_bgq_sec_comm;
EXTERN_FIELD bgq_weyl_vec *g_bgq_sec_recv[PHYSICAL_LD];
EXTERN_FIELD bgq_weyl_vec *g_bgq_sec_send[PHYSICAL_LD];


#undef EXTERN_INLINE
#undef EXTERN_FIELD
#undef EXTERN_INIT

#endif /* BGQ_COMM_H_ */
