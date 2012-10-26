/*
 * bgq_gaugefield.h
 *
 *  Created on: Oct 26, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_GAUGEFIELD_H_
#define BGQ_GAUGEFIELD_H_

#include "bgq_utils.h"
#include "bgq_field.h"

#ifndef BGQ_GAUGEFIELD_C_
#define EXTERN_INLINE EXTERN_INLINE_DECLARATION
#define EXTERN_FIELD extern
#define EXTERN_INIT(val)
#else
#define EXTERN_INLINE EXTERN_INLINE_DEFINITION
#define EXTERN_FIELD
#define EXTERN_INIT(val) = (val)
#endif


typedef struct {
	COMPLEX_PRECISION c[3][3][PHYSICAL_LK]; /* 3*3*2*sizeof(COMPLEX_PRECISION) = 288;144 bytes (4.5;2.25 L1 cache lines) */
	//COMPLEX_PRECISION padding[6];
} bgq_gaugesu3;
typedef struct {
	bgq_gaugesu3 su3[8];
} bgq_gaugesite;
typedef struct {
	bgq_gaugesu3 *(eodir[PHYSICAL_LP][PHYSICAL_LD]);
} bgq_gaugeeodir;
typedef bgq_gaugeeodir (*bgq_gaugefield);

EXTERN_FIELD bgq_gaugesite *g_bgq_gaugefield_fromHalfvolume[PHYSICAL_LP]; //TODO: Remove, not needed
EXTERN_FIELD bgq_gaugesite *g_bgq_gaugefield_fromSurface[PHYSICAL_LP];
EXTERN_FIELD bgq_gaugesite *g_bgq_gaugefield_fromBody[PHYSICAL_LP];

void bgq_gaugefield_init();
void bgq_gaugefield_transferfrom(su3 **sourcefield);


#undef EXTERN_INLINE
#undef EXTERN_FIELD
#undef EXTERN_INIT

#endif /* BGQ_GAUGEFIELD_H_ */
