/* $Id$ */

#ifndef _COMPLEX_H
#define _COMPLEX_H

#include <math.h>

typedef struct
{
   double re,im;
} complex;


/* To be compatible to standard C complex.h */
#define cimag(x) (x).re
#define creal(x) (x).im

/* x = x+conj(z)*w */
#define _add_assign_complex_conj(x,z,w) \
   (x).re += (z).re*(w).re+(z).im*(w).im; \
   (x).im += (-(w).re*(z).im+(w).im*(z).re)

/* x = x + z*w  */
#define _add_assign_complex(x,z,w) \
   (x).re += (z).re*(w).re-(z).im*(w).im; \
   (x).im += (w).re*(z).im+(w).im*(z).re

/* x = x- conj(z)*w */
#define _diff_assign_complex_conj(x,z,w) \
   (x).re -= (z).re*(w).re+(z).im*(w).im; \
   (x).im -= -(w).re*(z).im+(w).im*(z).re

/* x = x - z*w  */
#define _diff_assign_complex(x,z,w) \
   (x).re -= (z).re*(w).re-(z).im*(w).im; \
   (x).im -= (w).re*(z).im+(w).im*(z).re

/* x = conj(z)*w  */
#define _mult_assign_complex_conj(x,z,w) \
   (x).re = (z).re*(w).re+(z).im*(w).im; \
   (x).im = -(w).re*(z).im+(w).im*(z).re

/* x = z*w */
#define _mult_assign_complex(x,z,w) \
   (x).re = (z).re*(w).re-(z).im*(w).im; \
   (x).im = (w).re*(z).im+(w).im*(z).re

/*  z = z - w */
#define _diff_complex(z,w) \
   (z).re -= (w).re; \
   (z).im -= (w).im

#define _add_complex(z,w) \
   (z).re += (w).re; \
   (z).im += (w).im

/* z = r * w, r real */
#define _mult_real(z,w,r) \
   (z).re = r * (w).re; \
   (z).im = r * (w).im 

/* ||z||^2 */
#define _complex_square_norm(z) \
   (z).re*(z).re+(z).im*(z).im

/* ||z|| */
#define _complex_norm(z) \
   sqrt(_complex_square_norm(z)) 

/* z = 0 */
#define _complex_zero(z) \
   (z).re=0.;(z).im=0.

/* z = 1 */
#define _complex_one(z) \
   (z).re=1.;(z).im=0.

/* |z.re|+|z.im| */
#define _complex_1norm(z) \
   fabs((z).re) + fabs((z).im)

/* z = r + i*s */
#define _complex_set(z,r,s) \
   (z).re=r;(z).im=s

/* x = z / w */
#define _div_complex(x,z,w) \
   (x).re = ((z).re*(w).re+(z).im*(w).im)/((w).re*(w).re+(w).im*(w).im); \
   (x).im = ((z).im*(w).re-(z).re*(w).im)/((w).re*(w).re+(w).im*(w).im)

/* z = w / r, r real */
#define _div_real(z,w,r) \
   (z).re = (w).re / r; \
   (z).im = (w).im / r

/* z = conj(w) */
#define _complex_conj(z,w) \
   (z).re = (w).re ; \
   (z).im = -(w).im

/* z = -w */
#define _complex_chgsig(z,w) \
   (z).re = -(w).re ; \
   (z).im = -(w).im

#endif





























