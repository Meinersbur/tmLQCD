/* $Id$ */
/*******************************************************************************
 *
 *
 * Subroutines related to the lattice geometry
 *
 * The externally accessible function is
 *
 *   void geometry_eo(void)
 *     Computes the index arrays g_ipt, g_iup, g_idn, g_lexic2eo and g_eo2lexic
 *
 * original Version by
 * Author: Martin Luescher <luscher@mail.desy.ch>
 * Date: 24.10.2000
 *
 * Totally abused by M. Hasenbusch, now used for even-odd
 * decomposition of the lattice
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "su3adj.h"


#ifndef _NEW_GEOMETRY
int Index(const int x0, const int x1, const int x2, const int x3) {
  int y0, y1, y2, y3, ix;

  y0 = (x0 + T ) % T; 
  y1 = (x1 + LX) % LX; 
  y2 = (x2 + LY) % LY; 
  y3 = (x3 + LZ) % LZ;
  ix = ((y0*LX + y1)*LY + y2)*LZ + y3;
  
  y0=x0;
#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
  if(x0 == T) {
    ix = VOLUME + y3 + LZ*y2 + LZ*LY*y1;
  }
  /* the slice at time -1 is put to T+1 */
  else if(x0 == -1) {
    ix = VOLUME + LX*LY*LZ + y3 + LZ*y2 + LZ*LY*y1;
  }
#endif
#if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
  if(x1 == LX){
    ix = VOLUME + 2*LX*LY*LZ + y0*LY*LZ + y2*LZ + y3;
  }
  if(x1 == -1){
    ix = VOLUME + 2*LX*LY*LZ + T*LY*LZ + y0*LY*LZ + y2*LZ + y3;
  }   
  /* The edges */
  /* xt-edge */
  if(x0 == T){
    if(x1 == LX){
      ix = VOLUME+RAND+y2*LZ+y3;
    }
    if(x1 == -1){
      ix = VOLUME+RAND+LY*LZ+y2*LZ+y3;
    }
  }
  if(x0 == -1){
    if(x1 == LX){
      ix = VOLUME+RAND+2*LY*LZ+y2*LZ+y3;
    }
    if(x1 == -1){
      ix = VOLUME+RAND+3*LY*LZ+y2*LZ+y3;
    }
  }
  /* endif of PARALLELXT || PARALLELXYT*/
#endif

#if (defined PARALLELXYT || defined PARALLELXYZT)
  /* y-Rand */
  if(x2 == LY) {
    ix = VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + y0*LX*LZ + y1*LZ + y3;
  }
  if(x2 == -1) {
    ix = VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + T*LX*LZ + y0*LX*LZ + y1*LZ + y3;
  }
  /* the edges */
  /* yx-edge  */
  if(x1 == LX) {
    if(x2 == LY) {
      ix = VOLUME + RAND +  4*LY*LZ + y0*LZ + y3;
    }
    if(x2 == -1) {
      ix = VOLUME + RAND +  4*LY*LZ + T*LZ + y0*LZ + y3;
    }
  }
  if(x1 == -1) {
    if(x2 == LY) {
      ix = VOLUME + RAND +  4*LY*LZ + 2*T*LZ + y0*LZ + y3;
    }
    if(x2 == -1) {
      ix = VOLUME + RAND +  4*LY*LZ + 3*T*LZ + y0*LZ + y3;
    }
  }
  /* ty-edge */
  /* Be carefully here! Here we need y first, then t */
  /* this is because the chain is first t dir, then y direction */
  /* this is oposit to the other edges ! */
  if(x2 == LY) {
    if(x0 == T) {
      ix = VOLUME + RAND +  4*LY*LZ + 4*T*LZ + y1*LZ + y3;
    }
    if(x0 == -1) {
      ix = VOLUME + RAND +  4*LY*LZ + 4*T*LZ + LX*LZ + y1*LZ + y3;
    }
  }
  if(x2 == -1) {
    if(x0 == T) {
      ix = VOLUME + RAND +  4*LY*LZ + 4*T*LZ + 2*LX*LZ + y1*LZ + y3;
    }
    if(x0 == -1) {
      ix = VOLUME + RAND +  4*LY*LZ + 4*T*LZ + 3*LX*LZ + y1*LZ + y3;
    }
  }

  /* endif of PARALLELXYT */
#endif
#if defined PARALLELXYZT
  /* z-Rand */
  if(x3 == LZ) {
    ix = VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ 
      + y0*LX*LY + y1*LY + y2;
  }
  if(x3 == -1) {
    ix = VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY 
      + y0*LX*LY + y1*LY + y2;
  }
  /* the edges */
  /* zx-edge  */
  if(x1 == LX) {
    if(x3 == LZ) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ 
	+ y0*LY + y2;
    }
    if(x3 == -1) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + T*LY 
	+ y0*LY + y2;
    }
  }
  if(x1 == -1) {
    if(x3 == LZ) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 2*T*LY 
	+ y0*LY + y2;
    }
    if(x3 == -1) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 3*T*LY 
	+ y0*LY + y2;
    }
  }
  /* tz-edge */
  if(x3 == LZ) {
    if(x0 == T) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY 
	+ y1*LY + y2;
    }
    if(x0 == -1) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + LX*LY 
	+ y1*LY + y2;
    }
  }
  if(x3 == -1) {
    if(x0 == T) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 2*LX*LY 
	+ y1*LY + y2;
    }
    if(x0 == -1) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 3*LX*LY 
	+ y1*LY + y2;
    }
  }
  /* zy-edge */
  if(x3 == LZ) {
    if(x2 == LY) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY 
	+ y0*LX + y1;
    }
    if(x2 == -1) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY + 2*T*LX 
	+ y0*LX + y1;
    }
  }
  if(x3 == -1) {
    if(x2 == LY) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY + T*LX 
	+ y0*LX + y1;
    }
    if(x2 == -1) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY + 3*T*LX 
	+ y0*LX + y1;
    }
  }


  /* endif of PARALLELXYZT */
#endif

  /* The DBW2 stuff --> second boundary slice */
  /* This we put a the very end.              */
#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
  if(x0 == T+1) { 
    ix = VOLUMEPLUSRAND + y3 + LZ*y2 + LZ*LY*y1;
# if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
    /* t2x */
    if(x1 == LX) {
      ix = VOLUMEPLUSRAND + RAND
	+ y2*LZ + y3;
    }
    else if (x1 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 1*LY*LZ
	+ y2*LZ + y3;
    }
# endif
# if (defined PARALLELXYT || defined PARALLELXYZT)
    /* t2y */
    else if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ
	+ y1*LZ + y3;
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 2*LX*LZ
	+ y1*LZ + y3;
    }    
# endif
# if defined PARALLELXYZT
    /* t2z */
    else if(x3 == LZ) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ
	+ y1*LY + y2;
    }
    else if(x3 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 2*LX*LY
	+ y1*LY + y2;
    }
# endif
  }
  /* the slice at time -2 is put behind the one at time T+1 */
  else if(x0 == -2) {
    ix = VOLUMEPLUSRAND + LX*LY*LZ + y3 + LZ*y2 + LZ*LY*y1;
# if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
    /* t2x */
    if(x1 == LX) {
      ix = VOLUMEPLUSRAND + RAND + 2*LY*LZ
	+ y2*LZ + y3;
    }
    else if (x1 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 3*LY*LZ
	+ y2*LZ + y3;
    }
# endif
# if (defined PARALLELXYT || defined PARALLELXYZT)
    /* t2y */
    else if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + LX*LZ
	+ y1*LZ + y3;
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 3*LX*LZ
	+ y1*LZ + y3;
    }    
# endif
# if defined PARALLELXYZT
    /* t2z */
    else if(x3 == LZ) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + LX*LY
	+ y1*LY + y2;
    }
    else if(x3 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 3*LX*LY
	+ y1*LY + y2;
    }
# endif
  }  
#endif
#if ((defined PARALLELXT) || (defined PARALLELXYT) || defined PARALLELXYZT)
  if(x1 == LX+1) {
    if((x0 < T) && (x0 > -1) && (x2 < LY) && (x2 > -1) && (x3 < LZ) && (x3 > -1)) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + y0*LY*LZ + y2*LZ + y3;
    }
    /* x2t */
    else if(x0 == T) {
      ix = VOLUMEPLUSRAND + RAND + 4*LY*LZ
	+ y2*LZ + y3;
    }
    else if (x0 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 6*LY*LZ
	+ y2*LZ + y3;
    }
# if (defined PARALLELXYT || defined PARALLELXYZT)
    /* x2y */
    else if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ
	+ y0*LZ + y3;
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 1*T*LZ 
	+ y0*LZ + y3;
    }
# endif
# if defined PARALLELXYZT
    /* x2z */
    else if(x3 == LZ) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 4*T*LY
	+ y0*LY + y2;      
    }
    else if(x3 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 5*T*LY
	+ y0*LY + y2;      
    }
# endif
  }
  if(x1 == -2) {
    if((x0 < T) && (x0 > -1) && (x2 < LY) && (x2 > -1) && (x3 < LZ) && (x3 > -1)) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + T*LY*LZ + y0*LY*LZ + y2*LZ + y3;
    }
    /* x2t */
    else if(x0 == T) {
      ix = VOLUMEPLUSRAND + RAND + 5*LY*LZ
	+ y2*LZ + y3;
    }
    else if(x0 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 7*LY*LZ
	+ y2*LZ + y3;
    }
# if (defined PARALLELXYT || defined PARALLELXYZT)
    /* x2y */
    else if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 2*T*LZ
	+ y0*LZ + y3;
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 3*T*LZ
	+ y0*LZ + y3;
    }
# endif
# if defined PARALLELXYZT
    /* x2z */
    else if(x3 == LZ) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 6*T*LY
	+ y0*LY + y2;      
    }
    else if(x3 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 7*T*LY
	+ y0*LY + y2;      
    }
# endif
  }   
#endif
#if (defined PARALLELXYT || defined PARALLELXYZT)
  if(x2 == LY+1) {
    if((x0 < T) && (x0 > -1) && (x1 < LX) && (x1 > -1) && (x3 > -1) && (x3 < LZ)) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + y0*LX*LZ + y1*LZ + y3;
    }
    /* y2x */
    else if(x1 == LX) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 4*T*LZ
	+ y0*LZ + y3;
    }
    else if (x1 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 6*T*LZ
	+ y0*LZ + y3;
    }
    /* y2t */
    else if(x0 == T) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 4*LX*LZ
	+ y1*LZ + y3;
    }
    else if(x0 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 5*LX*LZ
	+ y1*LZ + y3;
    }
#  if defined PARALLELXYZT
    /* y2z */
    else if(x3 == LZ) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 4*T*LX
	+ y0*LX + y1;      
    }
    else if(x3 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 5*T*LX
	+ y0*LX + y1;      
    }
#  endif
  }
  if(x2 == -2) {
    if((x0 < T) && (x0 > -1) && (x1 < LX) && (x1 > -1) && (x3 > -1) && (x3 < LZ)) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + T*LX*LZ + y0*LX*LZ + y1*LZ + y3;      
    }
    /* y2x */
    else if(x1 == LX) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 5*T*LZ
	+ y0*LZ + y3;
    }
    else if (x1 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 7*T*LZ
	+ y0*LZ + y3;
    }
    /* y2t */
    else if(x0 == T) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 6*LX*LZ
	+ y1*LZ + y3;
    }
    else if (x0 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 7*LX*LZ
	+ y1*LZ + y3;
    }
#  if defined PARALLELXYZT
    /* y2z */
    else if(x3 == LZ) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 6*T*LX
	+ y0*LX + y1;      
    }
    else if(x3 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 7*T*LX
	+ y0*LX + y1;      
    }
# endif
  }
#endif
#if defined PARALLELXYZT
  /* z2-Rand */
  if(x3 == LZ+1) {
    if((x0 < T) && (x0 > -1) && (x1 < LX) && (x1 > -1) && (x2 > -1) && (x2 < LY)) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + y0*LX*LY + y1*LY + y2;
    }
    /* z2x */
    else if(x1 == LX) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY 
	+ y0*LY + y2;
    }
    else if (x1 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 2*T*LY
	+ y0*LY + y2;
    }
    /* z2t */
    else if(x0 == T) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 4*LX*LY
	+ y1*LY + y2;
    }
    else if(x0 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 5*LX*LY
	+ y1*LY + y2;
    }
    /* z2y */
    else if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY 
	+ y0*LX + y1;      
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 2*T*LX
	+ y0*LX + y1;      
    }
  }
  if(x3 == -2) {
    if((x0 < T) && (x0 > -1) && (x1 < LX) && (x1 > -1) && (x2 > -1) && (x2 < LY)) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY + y0*LX*LY + y1*LY + y2;
    }
    /* z2x */
    else if(x1 == LX) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + T*LY 
	+ y0*LY + y2;
    }
    else if(x1 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 3*T*LY
	+ y0*LY + y2;
    }
    /* z2t */
    else if(x0 == T) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 6*LX*LY
	+ y1*LY + y2;
    }
    else if(x0 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 7*LX*LY
	+ y1*LY + y2;
    }
    /* z2y */
    else if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 1*T*LX
	+ y0*LX + y1;      
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 3*T*LX
	+ y0*LX + y1;      
    }
  }
#endif
/*   if(ix == 372) { */
/*     printf("## %d %d %d %d ix = %d, %d %d %d %d\n", x0, x1, x2, x3, ix, T, LX, LY, LZ); */
/*   } */
  return(ix);
}

#else
/* now for even/odd geometry in the gauge fields. */
/* this works so far only up to 2dim parallelisation */
int Index(const int x0, const int x1, const int x2, const int x3)
{
  int y0, y1, y2, y3, ix, bndt=0, bndx=0, odd, bndy=0, bndz=0;

#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
  y0 = x0;
  /* the slice at time -1 is put to T+1 */
  if(x0 == -1) y0 = T+1;
  if(x0 == -1 || x0 == T) bndt = 1;
  if(x0 == -2) y0 = 1;
  if(x0 == T+1) y0 = 0;
  if(x0 == -2 || x0 == T+1) bndt = 2;
#else
  y0 = (x0+T) % T;
#endif
  
#if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
  y1 = x1;
  /* the slice at x -1 is put to LX+1 */
  if(x1 == -1) y1=LX+1;
  if(x1 == -1 || x1 == LX) bndx = 1;
  /* the slice at x -2 is put to LX+2 */
  if(x1 == -2) y1 = 1;
  if(x1 == LX+1) y1 = 0;
  if(x1 == -2 || x1 == LX+1) bndx = 2;
#else
  y1 = (x1+LX) % LX;
#endif

#if (defined PARALLELXYT || (defined PARALLELXYZT))
  y2 = x2;
  /* the slice at y -1 is put to LY+1 */
  if(x2 == -1) y2=LY+1;
  if(x2 == -1 || x2 == LY) bndy = 1;
  /* the slice at y -2 is put to LY+2 */
  if(x2 == -2) y2 = 1;
  if(x2 == LY+1) y2 = 0;
  if(x2 == -2 || x2 == LY+1) bndy = 2;
#else
   y2 = (x2+LY) % LY; 
#endif

#if defined PARALLELXYZT
  y3 = x3;
  /* the slice at z -1 is put to LZ+1 */
  if(x3 == -1) y3=LZ+1;
  if(x3 == -1 || x3 == LZ) bndz = 1;
  /* the slice at z -2 is put to LZ+2 */
  if(x3 == -2) y3 = 1;
  if(x3 == LZ+1) y3 = 0;
  if(x3 == -2 || x3 == LY+1) bndz = 2;
#else
   y3 = (x3+LZ) % LZ; 
#endif

   /* even or odd point? */
   if((x0+x1+x2+x3+g_proc_coords[0]*T + g_proc_coords[1]*LX +
       g_proc_coords[2]*LY + g_proc_coords[3]*LZ)%2 == 0) {
     odd = 0;
   } 
   else {
     odd = 1;
   }

   /* The local volume */
   if(bndt == 0 && bndx == 0 && bndy == 0) {
     ix = (y3 + LZ*y2 + LY*LZ*y1 + LX*LY*LZ*y0)/2 + (odd*(VOLUME))/2;
   }
   /* The time boundary */
   else if(bndt == 1 && bndx == 0 && bndy == 0) {
     ix = y0*LX*LY*LZ+(y3 + LZ*y2 + LY*LZ*y1)/2 + (odd*(LX*LY*LZ))/2;
   }
   /* The x boundary */
   else if(bndt == 0 && bndx == 1 && bndy == 0) {
     ix = VOLUME + 2*LX*LY*LZ 
       + (y1-LX)*T*LY*LZ + (y0*LY*LZ + y3 + LZ*y2)/2 + (odd*(LY*LZ*T))/2;
   }
   /* the y boundary */
   else if(bndt == 0 && bndx == 0 && bndy == 1) {
     ix = VOLUME + 2*LZ*(LX*LY+T*LY) 
       + (y2-LY)*T*LX*LZ + (y0*LX*LZ + y3 + LZ*y1)/2 + (odd*(LX*LZ*T))/2;
   }
   /* the xt edges */
   else if(bndt == 1 && bndx == 1 && bndy == 0) {
     ix = VOLUME + RAND 
       + ((y0-T)*2 + (y1-LX))*(LY*LZ) + (y3 + LZ*y2)/2 + (odd*(LY*LZ))/2;
   }
   /* the xy edges */
   else if(bndt == 0 && bndx == 1 && bndy == 1) {
     ix = VOLUME + RAND + 4*LY*LZ
       + ((y1-LX)*2 + (y2-LY))*(T*LZ) + (y3 + LZ*y0)/2 + (odd*(T*LZ))/2;
   }
   /* the ty edges */
   else if(bndt == 1 && bndx == 0 && bndy == 1) {
     ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ
       + ((y0-T)*2 + (y2-LY))*(LX*LZ) + (y3 + LZ*y1)/2 + (odd*(LX*LZ))/2;
   }
   /* now comes rectangular gauge action stuff */
   else if(bndt == 2 && bndx == 0 && bndy == 0) {
     ix = VOLUMEPLUSRAND 
       + y0*LX*LY*LZ + (y3 + LZ*y2 + LY*LZ*y1)/2 + (odd*(LX*LY*LZ))/2;
   }
   else if(bndt == 0 && bndx == 2 && bndy == 0) {
     ix = VOLUMEPLUSRAND + 2*LX*LY*LZ 
       + y1*T*LY*LZ + (y0*LY*LZ + y3 + LZ*y2)/2 + (odd*(LY*LZ*T))/2;
   }
   else if(bndt == 2 && bndx == 1 && bndy == 0) {
     ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ 
       + (2*y0 + (y1-LX))*(LY*LZ) + (y3 + LZ*y2)/2 + (odd*(LY*LZ))/2;
   }
   else if(bndt == 1 && bndx == 2 && bndy == 0) {
     ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 4*LY*LZ 
       + ((y0-T)*2 + y1)*(LY*LZ) + (y3 + LZ*y2)/2 + (odd*(LY*LZ))/2;
   }
   else if(bndt == 2 && bndx == 2) {
     printf("Should not happen in index routine!\n");
     printf("%d %d %d %d\n", x0, x1, x2 ,x3);
     ix = -1;
   }
   else if(bndt == 1 && bndx == 1 && bndy == 1) {
     printf("Should not be on three boundaries in index routine!\n");
     printf("%d %d %d %d\n", x0, x1, x2 ,x3);
     ix = -1;
   }
   else {
     printf("Error in index routine!\n");
     exit(1);
   }
   return( ix );
}

#endif

void geometry(){
  
  int x0,x1,x2,x3,ix;
  int y0, y1, y2, y3;
  int bndcnt=0;
  int i_even,i_odd;
  int startvaluet = 0;
  int startvaluex = 0;
  int startvaluey = 0;
  int startvaluez = 0;
  int * xeven;
  
  xeven = malloc(VOLUMEPLUSRAND*sizeof(int));

#if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  startvaluet = 1;
#endif
#if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  startvaluex = 1;
#endif
#if (defined PARALLELXYT || defined PARALLELXYZT)
  startvaluey = 1;
#endif
#if defined PARALLELXYZT
  startvaluez = 1;
#endif
  /* extended for boundary slices */
  for (x0 = -startvaluet; x0 < (T+startvaluet); x0++){
    for (x1 = -startvaluex; x1 < (LX+startvaluex); x1++){
      for (x2 = -startvaluey; x2 < (LY+startvaluey); x2++){
	for (x3 = -startvaluez; x3 < (LZ+startvaluez); x3++){
	  bndcnt = 0;
	  if(x0 < 0 || x0 > T-1) bndcnt++;
	  if(x1 < 0 || x1 > LX-1) bndcnt++;
	  if(x2 < 0 || x2 > LY-1) bndcnt++;
	  if(x3 < 0 || x3 > LZ-1) bndcnt++;

	  y0=x0; y1=x1; y2=x2; y3=x3;
	  if(x0 == -1) {
	    y0 = T+1;
	  }
	  if(x1 == -1) {
	    y1 = LX+1;
	  }
	  if(x2 == -1) {
	    y2 = LY+1;
	  }
	  if(x3 == -1) {
	    y3 = LZ+1;
	  }
	  if(bndcnt > 2) {
	    /* Should not be needed, set it to -1 */
	    g_ipt[y0][y1][y2][y3] = -1;
	  }
	  else {
	    ix=Index(x0, x1, x2, x3);
	    g_ipt[y0][y1][y2][y3] = ix;
	    /* g_proc_id*T|LX|LY|LZ is added to allow for odd T|LX|LY|LZ when the number of 
	       nodes is even */	
	    if((x0 + x1 + x2 + x3 + 
		g_proc_coords[0]*T + g_proc_coords[1]*LX + 
		g_proc_coords[2]*LY + g_proc_coords[3]*LZ)%2==0) {
	      xeven[ix]=1;
	    } 
	    else {
	      xeven[ix]=0; 
	    }

	    g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	    g_idn[ix][0] = Index(x0-1, x1, x2, x3);
	    
	    g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	    g_idn[ix][1] = Index(x0, x1-1, x2, x3);
	    
	    g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	    g_idn[ix][2] = Index(x0, x1, x2-1, x3);
	    
	    g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	    g_idn[ix][3] = Index(x0, x1, x2, x3-1);

	  }
	}
      }
    }
  }

  i_even=0;
  i_odd=0;
  /*For the spinor fields we need only till VOLUME+RAND */
  for (ix = 0; ix < (VOLUME+RAND); ix++){
    if(xeven[ix]==1){
      g_lexic2eo[ix] = i_even;
      g_lexic2eosub[ix] = i_even;
      g_eo2lexic[i_even] = ix;
      i_even++;
    }
    else{
      g_lexic2eo[ix] = (VOLUME+RAND)/2+i_odd;
      g_lexic2eosub[ix] = i_odd;
      g_eo2lexic[(VOLUME+RAND)/2+i_odd] = ix;
      i_odd++;
    }
  }

#if defined PARALLELXYZT
  ix = 0;
  for(x0 = 0; x0 < T; x0++) {
    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	if((x0 + x1 + x2 +  
	    g_proc_coords[0]*T + g_proc_coords[1]*LX +  
	    g_proc_coords[2]*LY + g_proc_coords[3]*LZ)%2==0) { 
	  g_field_z_ipt_even[ix] = g_lexic2eosub[ g_ipt[x0][x1][x2][0]];
	  ix++;
	}
      }
    }
  }
  for(x0 = 0; x0 < T; x0++) {
    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	if((x0 + x1 + x2 + (LZ-1) + 
	    g_proc_coords[0]*T + g_proc_coords[1]*LX +  
	    g_proc_coords[2]*LY + g_proc_coords[3]*LZ)%2==0) { 
	  g_field_z_ipt_even[ix] = g_lexic2eosub[ g_ipt[x0][x1][x2][LZ-1]];
	  ix++;
	}
      }
    }
  }
  ix = 0;
  for(x0 = 0; x0 < T; x0++) {
    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	if((x0 + x1 + x2 +  
	    g_proc_coords[0]*T + g_proc_coords[1]*LX +  
	    g_proc_coords[2]*LY + g_proc_coords[3]*LZ)%2==1) { 
	  g_field_z_ipt_odd[ix] = g_lexic2eosub[ g_ipt[x0][x1][x2][0]];
	  ix++;
	}
      }
    }
  }
  for(x0 = 0; x0 < T; x0++) {
    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	if((x0 + x1 + x2 + (LZ-1) + 
	    g_proc_coords[0]*T + g_proc_coords[1]*LX +  
	    g_proc_coords[2]*LY + g_proc_coords[3]*LZ)%2==1) { 
	  g_field_z_ipt_odd[ix] = g_lexic2eosub[ g_ipt[x0][x1][x2][LZ-1]];
	  ix++;
	}
      }
    }
  }
#endif
  /* The rectangular gauge action part */
  /* Everything is stored behind VOLUMEPLUSRAND-1 !*/
  if(g_dbw2rand != 0) {
#if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
    if(g_proc_id == 0) {
      printf("# Initialising rectangular gauge action stuff\n");
      fflush(stdout);
    }
    for (x1 = -startvaluex; x1 < (LX+startvaluex); x1++){
      for (x2 = -startvaluey; x2 < (LY+startvaluey); x2++) {
	for (x3 = -startvaluez; x3 < (LZ+startvaluez); x3++) {
	  bndcnt = 0;
	  if(x1 < 0 || x1 > LX-1) bndcnt++;
	  if(x2 < 0 || x2 > LY-1) bndcnt++;
	  if(x3 < 0 || x3 > LZ-1) bndcnt++;	  
	  if(bndcnt < 2) {
	    /* t2 Rand and t2x and t2y */
	    x0 = -2;
	    ix = Index(x0, x1, x2, x3);
	    if(ix < VOLUMEPLUSRAND) {
	      printf("#### -2t %d %d %d %d\n",x0, x1, x2, x3);
	    }

	    g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	    g_idn[ix][0] = -1;

	    if(x1 < LX) g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	    if(x1 > -1) g_idn[ix][1] = Index(x0, x1-1, x2, x3);

	    if(x2 < LY) g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	    if(x2 > -1) g_idn[ix][2] = Index(x0, x1, x2-1, x3);

	    if(x3 < LZ) g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	    if(x3 > -1) g_idn[ix][3] = Index(x0, x1, x2, x3-1);

	    x0 = T+1;
	    ix = Index(x0, x1, x2, x3);
	    if(ix < VOLUMEPLUSRAND) {
	      printf("#### +2t %d %d %d %d\n",x0, x1, x2, x3);
	    }
	    g_iup[ix][0] = -1;
	    g_idn[ix][0] = Index(x0-1, x1, x2, x3);

	    if(x1 < LX) g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	    if(x1 > -1) g_idn[ix][1] = Index(x0, x1-1, x2, x3);

	    if(x2 < LY) g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	    if(x2 > -1) g_idn[ix][2] = Index(x0, x1, x2-1, x3);

	    if(x3 < LZ) g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	    if(x3 > -1) g_idn[ix][3] = Index(x0, x1, x2, x3-1);

	  }
	}
      }
    }    
#endif
#if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
    for (x0 = -startvaluet; x0 < (T+startvaluet); x0++){
      for (x2 = -startvaluey; x2 < (LY+startvaluey); x2++) {
	for (x3 = -startvaluez; x3 < (LZ+startvaluez); x3++) {
	  bndcnt = 0;
	  if(x0 < 0 || x0 > T-1) bndcnt++;
	  if(x2 < 0 || x2 > LY-1) bndcnt++;
	  if(x3 < 0 || x3 > LZ-1) bndcnt++;	  
	  if(bndcnt < 2) {
	    /* x2-Rand and x2t and x2y */
	    x1 = -2;
	    ix = Index(x0, x1, x2, x3);
	    if(ix < VOLUMEPLUSRAND) {
	      printf("#### -2x %d %d %d %d\n",x0, x1, x2, x3);
	    }
	    if(x0 < T) g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	    if(x0 > -1) g_idn[ix][0] = Index(x0-1, x1, x2, x3);

	    g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	    g_idn[ix][1] = -1;

	    if(x2 < LY) g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	    if(x2 > -1) g_idn[ix][2] = Index(x0, x1, x2-1, x3);

	    if(x3 < LZ) g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	    if(x3 > -1) g_idn[ix][3] = Index(x0, x1, x2, x3-1);

	    x1 = LX+1;
	    ix = Index(x0, x1, x2, x3);
	    if(ix < VOLUMEPLUSRAND) {
	      printf("#### +2x %d %d %d %d\n",x0, x1, x2, x3);
	    }
	    if(x0 < T) g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	    if(x0 > -1) g_idn[ix][0] = Index(x0-1, x1, x2, x3);

	    g_iup[ix][1] = -1;
	    g_idn[ix][1] = Index(x0, x1-1, x2, x3);

	    if(x2 < LY) g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	    if(x2 > -1) g_idn[ix][2] = Index(x0, x1, x2-1, x3);

	    if(x3 < LZ) g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	    if(x3 > -1) g_idn[ix][3] = Index(x0, x1, x2, x3-1);
	  }
	}
      }
    }
#endif
#if (defined PARALLELXYT || defined PARALLELXYZT)
    for (x0 = -startvaluet; x0 < (T+startvaluet); x0++){
      for (x1 = -startvaluex; x1 < (LX+startvaluex); x1++) {
	for (x3 = -startvaluez; x3 < (LZ+startvaluez); x3++) {
	  bndcnt = 0;
	  if(x0 < 0 || x0 > T-1) bndcnt++;
	  if(x1 < 0 || x1 > LX-1) bndcnt++;
	  if(x3 < 0 || x3 > LZ-1) bndcnt++;	  
	  if(bndcnt < 2) {
	    /* y2-Rand y2t and y2x */
	    x2 = -2;
	    ix = Index(x0, x1, x2, x3);
	    if(ix < VOLUMEPLUSRAND || ix >= VOLUMEPLUSRAND + g_dbw2rand) {
	      printf("#### -2y %d %d %d %d\n",x0, x1, x2, x3);
	    }
	    if(x0 < T) g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	    if(x0 > -1) g_idn[ix][0] = Index(x0-1, x1, x2, x3);
	    
	    if(x1 < LX) g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	    if(x1 > -1) g_idn[ix][1] = Index(x0, x1-1, x2, x3);
	    
	    g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	    g_idn[ix][2] = -1;
	    
	    if(x3 < LZ) g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	    if(x3 > -1) g_idn[ix][3] = Index(x0, x1, x2, x3-1);
	    
	    x2 = LY+1;
	    ix = Index(x0, x1, x2, x3);
	    if(ix < VOLUMEPLUSRAND || ix >= VOLUMEPLUSRAND + g_dbw2rand) {
	      printf("#### +2y %d %d %d %d\n",x0, x1, x2, x3);
	    }
	    if(x0 < T) g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	    if(x0 > -1) g_idn[ix][0] = Index(x0-1, x1, x2, x3);
	    
	    if(x1 < LX) g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	    if(x1 > -1) g_idn[ix][1] = Index(x0, x1-1, x2, x3);
	    
	    g_iup[ix][2] = -1;
	    g_idn[ix][2] = Index(x0, x1, x2-1, x3);
	    
	    if(x3 < LZ) g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	    if(x3 > -1) g_idn[ix][3] = Index(x0, x1, x2, x3-1);
	  }
	}
      }
    }
#endif
#if (defined PARALLELXYZT)
    for (x0 = -startvaluet; x0 < (T+startvaluet); x0++){
      for (x1 = -startvaluex; x1 < (LX+startvaluex); x1++) {
	for (x2 = -startvaluey; x2 < (LY+startvaluey); x2++) {
	  bndcnt = 0;
	  if(x0 < 0 || x0 > T-1) bndcnt++;
	  if(x1 < 0 || x1 > LX-1) bndcnt++;
	  if(x2 < 0 || x2 > LY-1) bndcnt++;	  
	  if(bndcnt < 2) {
	    /* z2-Rand t2z and z2x z2y*/
	    x3 = -2;
	    ix = Index(x0, x1, x2, x3);
	    if(ix < VOLUMEPLUSRAND || ix >= VOLUMEPLUSRAND + g_dbw2rand) {
	      printf("#### -2z %d %d %d %d %d %d %d %d %d %d %d\n",x0, x1, x2, x3, ix, 
		     VOLUMEPLUSRAND, VOLUMEPLUSRAND + g_dbw2rand, T, LX, LY, LZ);
	    }
	    if(x0 <  T) g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	    if(x0 > -1) g_idn[ix][0] = Index(x0-1, x1, x2, x3);
	    
	    if(x1 < LX) g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	    if(x1 > -1) g_idn[ix][1] = Index(x0, x1-1, x2, x3);
	    
	    if(x2 < LY) g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	    if(x2 > -1) g_idn[ix][2] = Index(x0, x1, x2-1, x3);
	    
	    g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	    g_idn[ix][3] = -1;
	    
	    x3 = LZ+1;
	    ix = Index(x0, x1, x2, x3);
	    if(ix < VOLUMEPLUSRAND || ix >= VOLUMEPLUSRAND + g_dbw2rand) {
	      printf("#### +2z %d %d %d %d\n",x0, x1, x2, x3);
	    }
	    if(x0 <  T) g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	    if(x0 > -1) g_idn[ix][0] = Index(x0-1, x1, x2, x3);
	    
	    if(x1 < LX) g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	    if(x1 > -1) g_idn[ix][1] = Index(x0, x1-1, x2, x3);
	    
	    if(x2 < LY) g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	    if(x2 > -1) g_idn[ix][2] = Index(x0, x1, x2-1, x3);
	    
	    g_iup[ix][3] = -1;
	    g_idn[ix][3] = Index(x0, x1, x2, x3-1);
	  }
	}
      }
    }
#endif
  }
}

static char const rcsid[] = "$Id$";


