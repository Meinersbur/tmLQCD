/*
 * bgq_spinorfield.h
 *
 *  Created on: Oct 26, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_SPINORFIELD_H_
#define BGQ_SPINORFIELD_H_

#include "bgq_utils.h"

#ifndef BGQ_SPINORFIELD_C_
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


#endif /* BGQ_SPINORFIELD_H_ */
