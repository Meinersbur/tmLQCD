#ifndef _HOPPING_UNOPT_H
#define _HOPPING_UNOPT_H

#include "Hopping_Matrix.h"

typedef int bool;
#define false 0
#define true 1

void Hopping_unopt_halfspinor(const int ieo, spinor * const l, spinor * const k);
void Hopping_unopt_fullspinor(const int ieo, spinor * const l, spinor * const k);


#endif /* _HOPPING_UNOPT_H */
