/*
 * bgq_comm.h
 *
 *  Created on: Oct 25, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_COMM_H_
#define BGQ_COMM_H_

#include "bgq_utils.h"

#ifndef BGQ_COMM_C_
#define EXTERN_INLINE EXTERN_INLINE_DECLARATION
#define EXTERN_FIELD extern
#define EXTERN_INIT(val)
#else
#define EXTERN_INLINE EXTERN_INLINE_DEFINITION
#define EXTERN_FIELD
#define EXTERN_INIT(val) = (val)
#endif


#undef EXTERN_INLINE
#undef EXTERN_FIELD
#undef EXTERN_INIT

#endif /* BGQ_COMM_H_ */
