/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include "cmalloc.h"

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#if defined BGL && !defined BGP
# include <rts.h>
#endif
#ifdef _USE_SHMEM
# include <mpp/shmem.h>
#endif
#include "global.h"
#include "su3.h"
#include "init_dirac_halfspinor.h"

halfspinor ** NBPointer_;
halfspinor * HalfSpinor_;
halfspinor * HalfSpinor ALIGN;
halfspinor *** NBPointer;

/* The single precision versions */
halfspinor32 ** NBPointer32_;
halfspinor32 * HalfSpinor32_;
halfspinor32 * HalfSpinor32 ALIGN;
halfspinor32 *** NBPointer32;

/* GG */
//halfspinor * GalfSpinor ALIGN;
#define TRANSTST 0

int init_dirac_halfspinor_tile() {
  int ieo=0, i=0, j=0, k;
  int x, y, z, t, mu;
#if (defined BGL && defined _USE_BGLDRAM)
  unsigned int actualSize;
  int rts_return=0;
#endif  

  NBPointer = (halfspinor***) malloc(4*sizeof(halfspinor**));
  NBPointer_ = (halfspinor**) malloc(16*(VOLUME+RAND)*sizeof(halfspinor*));
  NBPointer[0] = NBPointer_;
  NBPointer[1] = NBPointer_ + (8*(VOLUME+RAND)/2);
  NBPointer[2] = NBPointer_ + (16*(VOLUME+RAND)/2);
  NBPointer[3] = NBPointer_ + (24*(VOLUME+RAND)/2);

#ifdef _USE_SHMEM
  HalfSpinor_ = (halfspinor*)shmalloc((8*(VOLUME+RAND)+1)*sizeof(halfspinor));
#elif (defined BGL && defined _USE_BGLDRAM)
  rts_return = rts_get_dram_window(8*((VOLUME+RAND)+1)*sizeof(halfspinor), RTS_STORE_WITHOUT_ALLOCATE, (void**)&HalfSpinor_, &actualSize);
  if(rts_return !=0) {
    return(-1);
  }
#else
  HalfSpinor_ = (halfspinor*)calloc(8*(VOLUME+RAND)+1, sizeof(halfspinor));
#endif
  if(errno == ENOMEM) {
    return(-1);
  }

  /* GG */
#if TRANSTST
  GalfSpinorSrc = (halfspinor*)calloc(1024*1024, sizeof(halfspinor));
  if(errno == ENOMEM) {
    return(-1);
  }
  GalfSpinorDst = (halfspinor*)calloc(1024*1024, sizeof(halfspinor));
  if(errno == ENOMEM) {
    return(-1);
  }
#endif
  

#if ( defined SSE || defined SSE2 || defined SSE3)
  HalfSpinor = (halfspinor*)(((unsigned long int)(HalfSpinor_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  HalfSpinor = HalfSpinor_;
#endif

  for(ieo = 0; ieo < 2; ieo++) {
    for(mu = 0; mu < 4; mu++) {  /* ZL */
      for(i = 0; i < VOLUME/2; i++) {
        j = g_eo2lexic[i + ((ieo+1)%2)*(VOLUME+RAND)/2];
        /* get (t,x,y,z) from j */
        t = j/(LX*LY*LZ);
        x = (j-t*(LX*LY*LZ))/(LY*LZ);
        y = (j-t*(LX*LY*LZ)-x*(LY*LZ))/(LZ);
        z = (j-t*(LX*LY*LZ)-x*(LY*LZ) - y*LZ);
	if( g_lexic2eosub[ g_idn[j][mu] ]  < VOLUME/2 ) {
        NBPointer[ieo][VOLUME*mu+2*i+0]=&HalfSpinor[VOLUME*mu+2*g_lexic2eosub[ g_idn[j][mu] ]+0];  
/*	printf("j,g_lexic2eosub dn:%d, %d\n",j,g_lexic2eosub[ g_idn[j][mu] ]); */
	}
	if( g_lexic2eosub[ g_iup[j][mu] ]  < VOLUME/2 ) {
        NBPointer[ieo][VOLUME*mu+2*i+1]=&HalfSpinor[VOLUME*mu+2*g_lexic2eosub[ g_iup[j][mu] ]+1];  /* ZL */
	}
#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
        if(t == 0 && mu == 0) {
    	  k = 8*VOLUME/2 + (g_lexic2eosub[g_idn[j][0]] - VOLUME/2);
    	  NBPointer[ieo][2*i] = &HalfSpinor[ k ];
/*	  printf("PARAT,j, g_lexic2eosub dn: %d,%d\n",j,g_lexic2eosub[
g_idn[j][0] ]);*/
        } 
        if(t == T-1 && mu == 0) {
  	k = 8*VOLUME/2 + (g_lexic2eosub[g_iup[j][0]] - VOLUME/2);
  	NBPointer[ieo][2*i + 1] = &HalfSpinor[ k ];
        }
#endif
#if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT) || (defined PARALLELXYZ))
        if(x == 0 && mu == 1) {
  	k = 8*VOLUME/2 + (g_lexic2eosub[g_idn[j][1]] - VOLUME/2);
  	NBPointer[ieo][VOLUME+2*i + 0] = &HalfSpinor[ k ];
/*	printf("PARAXT,j, g_lexic2eosub dn: %d,%d\n",j, g_lexic2eosub[
g_idn[j][1] ]); */
        }
        if(x == LX-1 && mu == 1) {
  	k = 8*VOLUME/2 + (g_lexic2eosub[g_iup[j][1]] - VOLUME/2);
  	NBPointer[ieo][VOLUME+2*i + 1] = &HalfSpinor[ k ];
        }
#endif
#if ((defined PARALLELXYT) || (defined PARALLELXYZT) || (defined PARALLELXYZ))
        if(y == 0 && mu == 2) {
  	k = 8*VOLUME/2 + (g_lexic2eosub[g_idn[j][2]] - VOLUME/2);
  	NBPointer[ieo][2*VOLUME+2*i + 0] = &HalfSpinor[ k ];
        }
        if(y == LY-1 && mu == 2) {
  	k = 8*VOLUME/2 + (g_lexic2eosub[g_iup[j][2]] - VOLUME/2);
  	NBPointer[ieo][2*VOLUME+2*i + 1] = &HalfSpinor[ k ];
        }
#endif
#if (defined PARALLELXYZT) || (defined PARALLELXYZ)
        if(z == 0 && mu == 3) {
  	k = 8*VOLUME/2 + (g_lexic2eosub[g_idn[j][3]] - VOLUME/2);
  	NBPointer[ieo][3*VOLUME+2*i + 0] = &HalfSpinor[ k ];
        }
        if(z == LZ-1 && mu == 3) {
  	k = 8*VOLUME/2 + (g_lexic2eosub[g_iup[j][3]] - VOLUME/2);
  	NBPointer[ieo][3*VOLUME+2*i + 1] = &HalfSpinor[ k ];
        }
#endif
/*
printf("mu,t,x,y,z,i,VOLUME*mu+2*i+0,%d,%d,%d,%d,%d,%d,%d\n",mu,t,x,y,z,i,VOLUME*mu+2*i+0);*/
      } /* i */
    }/* mu */
    for(i = VOLUME/2; i < (VOLUME+RAND)/2; i++) {
      for(mu = 0; mu < 8; mu++) {
	NBPointer[ieo][8*i + mu] = NBPointer[ieo][0];
      }
    }
#ifdef MPI
#endif
  }/* ieo */
  
  for(ieo = 2; ieo < 4; ieo++) {
   for(mu = 0; mu < 4; mu++) { 
    for(i = 0; i < VOLUME/2; i++) {
      j = g_eo2lexic[i + ((ieo+0)%2)*(VOLUME+RAND)/2];
      /* get (t,x,y,z) from j */
      t = j/(LX*LY*LZ);
      x = (j-t*(LX*LY*LZ))/(LY*LZ);
      y = (j-t*(LX*LY*LZ)-x*(LY*LZ))/(LZ);
      z = (j-t*(LX*LY*LZ)-x*(LY*LZ) - y*LZ);
  
NBPointer[ieo][VOLUME*mu+2*i] = &HalfSpinor[VOLUME*mu+2*i];
NBPointer[ieo][VOLUME*mu+2*i+1] = &HalfSpinor[VOLUME*mu+2*i+1];     
#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(t == T-1 && mu == 0) {
	NBPointer[ieo][2*i]     = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_iup[j][0] ] - VOLUME/2)];
       }
      if(t == 0 && mu == 0) {
	NBPointer[ieo][2*i + 1] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_idn[j][0] ] - VOLUME/2)];
       }
#endif
#if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT) || (defined PARALLELXYZ))
      if(x == LX-1 && mu == 1) { 
	NBPointer[ieo][VOLUME+2*i + 0] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_iup[j][1] ] - VOLUME/2)];
      }
      if(x == 0 && mu == 1) {
	NBPointer[ieo][VOLUME+2*i + 1] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_idn[j][1] ] - VOLUME/2)];
      }
#endif
#if ((defined PARALLELXYT) || (defined PARALLELXYZT) || (defined PARALLELXYZ))
      if(y == LY-1 && mu == 2) {
	NBPointer[ieo][2*VOLUME+2*i + 0] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_iup[j][2] ] - VOLUME/2)];
      }
      if(y == 0 && mu == 2) {
	NBPointer[ieo][2*VOLUME+2*i + 1] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_idn[j][2] ] - VOLUME/2)];
      }
#endif
#if (defined PARALLELXYZT || (defined PARALLELXYZ))
      if(z == LZ-1 && mu == 3) {
	NBPointer[ieo][3*VOLUME+2*i + 0] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_iup[j][3] ] - VOLUME/2)];
      }
      if(z == 0 && mu == 3) {
	NBPointer[ieo][3*VOLUME+2*i + 1] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_idn[j][3] ] - VOLUME/2)];
      }
#endif
    }/* i*/
      }  /*   for(mu = 0; mu < 4; mu++)  OP    */
     for(i = VOLUME/2; i < (VOLUME+RAND)/2; i++) {
      for(mu = 0; mu < 8; mu++) {
	NBPointer[ieo][8*i + mu] = NBPointer[ieo][0];
      }
    }
#ifdef MPI
/*     NBPointer[ieo][4*VOLUME] = NBPointer[ieo][0];  */
#endif
  } /* ieo */
  return(0);
} /* init_dirac_halfspinor_tile */

/*==============================================*/

int init_dirac_halfspinor() {
  int ieo=0, i=0, j=0, k;
  int x, y, z, t, mu;
#if (defined BGL && defined _USE_BGLDRAM && !defined BGP)
  unsigned int actualSize;
  int rts_return=0;
#endif  

  NBPointer = (halfspinor***) calloc(4,sizeof(halfspinor**));
  CMALLOC_ERROR_EXIT(NBPointer);
  NBPointer_ = (halfspinor**) calloc(16,(VOLUME+RAND)*sizeof(halfspinor*));
  CMALLOC_ERROR_EXIT(NBPointer_);
  NBPointer[0] = NBPointer_;
  NBPointer[1] = NBPointer_ + (8*(VOLUME+RAND)/2);
  NBPointer[2] = NBPointer_ + (16*(VOLUME+RAND)/2);
  NBPointer[3] = NBPointer_ + (24*(VOLUME+RAND)/2);

#ifdef _USE_SHMEM
  if((coid*)(HalfSpinor_ = (halfspinor*)shmalloc((8*(VOLUME+RAND)+1)*sizeof(halfspinor))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }
#elif (defined BGL && defined _USE_BGLDRAM && !defined BGP)
  rts_return = rts_get_dram_window(8*((VOLUME+RAND)+1)*sizeof(halfspinor), RTS_STORE_WITHOUT_ALLOCATE, (void**)&HalfSpinor_, &actualSize);
  if(rts_return !=0) {
    return(-1);
  }
#else
  if((void*)(HalfSpinor_ = (halfspinor*)calloc(8*(VOLUME+RAND)+1, sizeof(halfspinor))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }
#endif
#if ( defined SSE || defined SSE2 || defined SSE3)
  HalfSpinor = (halfspinor*)(((unsigned long int)(HalfSpinor_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  HalfSpinor = HalfSpinor_;
#endif

  for(ieo = 0; ieo < 2; ieo++) {
    for(i = 0; i < VOLUME/2; i++) {
      j = g_eo2lexic[i + ((ieo+1)%2)*(VOLUME+RAND)/2];
      /* get (t,x,y,z) from j */
      t = j/(LX*LY*LZ);
      x = (j-t*(LX*LY*LZ))/(LY*LZ);
      y = (j-t*(LX*LY*LZ)-x*(LY*LZ))/(LZ);
      z = (j-t*(LX*LY*LZ)-x*(LY*LZ) - y*LZ);
      for(mu = 0; mu < 4; mu++) {
	NBPointer[ieo][8*i + 2*mu + 0] = &HalfSpinor[ 8*g_lexic2eosub[ g_idn[j][mu] ] + 2*mu + 0];
	NBPointer[ieo][8*i + 2*mu + 1] = &HalfSpinor[ 8*g_lexic2eosub[ g_iup[j][mu] ] + 2*mu + 1];
      }
#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(t == 0) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_idn[j][0]] - VOLUME/2);
	NBPointer[ieo][8*i] = &HalfSpinor[ k ];
      }
      if(t == T-1) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_iup[j][0]] - VOLUME/2);
	NBPointer[ieo][8*i + 1] = &HalfSpinor[ k ];
      }
#endif
#if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT) || (defined PARALLELXYZ))
      if(x == 0) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_idn[j][1]] - VOLUME/2);
	NBPointer[ieo][8*i + 2] = &HalfSpinor[ k ];
      }
      if(x == LX-1) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_iup[j][1]] - VOLUME/2);
	NBPointer[ieo][8*i + 3] = &HalfSpinor[ k ];
      }
#endif
#if ((defined PARALLELXYT) || (defined PARALLELXYZT) || (defined PARALLELXYZ))
      if(y == 0) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_idn[j][2]] - VOLUME/2);
	NBPointer[ieo][8*i + 4] = &HalfSpinor[ k ];
      }
      if(y == LY-1) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_iup[j][2]] - VOLUME/2);
	NBPointer[ieo][8*i + 5] = &HalfSpinor[ k ];
      }
#endif
#if (defined PARALLELXYZT) || (defined PARALLELXYZ)
      if(z == 0) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_idn[j][3]] - VOLUME/2);
	NBPointer[ieo][8*i + 6] = &HalfSpinor[ k ];
      }
      if(z == LZ-1) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_iup[j][3]] - VOLUME/2);
	NBPointer[ieo][8*i + 7] = &HalfSpinor[ k ];
      }
#endif
    }
    for(i = VOLUME/2; i < (VOLUME+RAND)/2; i++) {
      for(mu = 0; mu < 8; mu++) {
	NBPointer[ieo][8*i + mu] = NBPointer[ieo][0];
      }
    }
#ifdef MPI
#endif
  }
  for(ieo = 2; ieo < 4; ieo++) {
    for(i = 0; i < VOLUME/2; i++) {
      j = g_eo2lexic[i + ((ieo+0)%2)*(VOLUME+RAND)/2];
      /* get (t,x,y,z) from j */
      t = j/(LX*LY*LZ);
      x = (j-t*(LX*LY*LZ))/(LY*LZ);
      y = (j-t*(LX*LY*LZ)-x*(LY*LZ))/(LZ);
      z = (j-t*(LX*LY*LZ)-x*(LY*LZ) - y*LZ);
      for(mu = 0; mu < 8; mu++) {
	NBPointer[ieo][8*i + mu] = &HalfSpinor[8*i + mu];
      }
#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(t == T-1) {
	NBPointer[ieo][8*i]     = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_iup[j][0] ] - VOLUME/2)];
      }
      if(t == 0) {
	NBPointer[ieo][8*i + 1] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_idn[j][0] ] - VOLUME/2)];
      }
#endif
#if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT) || (defined PARALLELXYZ))
      if(x == LX-1) { 
	NBPointer[ieo][8*i + 2] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_iup[j][1] ] - VOLUME/2)];
      }
      if(x == 0) {
	NBPointer[ieo][8*i + 3] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_idn[j][1] ] - VOLUME/2)];
      }
#endif
#if ((defined PARALLELXYT) || (defined PARALLELXYZT) || (defined PARALLELXYZ))
      if(y == LY-1) {
	NBPointer[ieo][8*i + 4] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_iup[j][2] ] - VOLUME/2)];
      }
      if(y == 0) {
	NBPointer[ieo][8*i + 5] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_idn[j][2] ] - VOLUME/2)];
      }
#endif
#if (defined PARALLELXYZT) || (defined PARALLELXYZ)
      if(z == LZ-1) {
	NBPointer[ieo][8*i + 6] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_iup[j][3] ] - VOLUME/2)];
      }
      if(z == 0) {
	NBPointer[ieo][8*i + 7] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_idn[j][3] ] - VOLUME/2)];
      }
#endif
    }
    for(i = VOLUME/2; i < (VOLUME+RAND)/2; i++) {
      for(mu = 0; mu < 8; mu++) {
	NBPointer[ieo][8*i + mu] = NBPointer[ieo][0];
      }
    }
#ifdef MPI
/*     NBPointer[ieo][4*VOLUME] = NBPointer[ieo][0];  */
#endif
  }
  return(0);
}


int init_dirac_halfspinor32() {
  int ieo=0, i=0, j=0, k;
  int x, y, z, t, mu;
#if (defined BGL && defined _USE_BGLDRAM && !defined BGP)
  unsigned int actualSize;
  int rts_return=0;
#endif  
  
  NBPointer32 = (halfspinor32***) calloc(4,sizeof(halfspinor32**));
  CMALLOC_ERROR_EXIT(NBPointer32);
  NBPointer32_ = (halfspinor32**) calloc(16,(VOLUME+RAND)*sizeof(halfspinor32*));
  CMALLOC_ERROR_EXIT(NBPointer32_);
  NBPointer32[0] = NBPointer32_;
  NBPointer32[1] = NBPointer32_ + (8*(VOLUME+RAND)/2);
  NBPointer32[2] = NBPointer32_ + (16*(VOLUME+RAND)/2);
  NBPointer32[3] = NBPointer32_ + (24*(VOLUME+RAND)/2);

#ifdef _USE_SHMEM
  if((void*)(HalfSpinor32_ = (halfspinor32*)shmalloc((8*(VOLUME+RAND)+1)*sizeof(halfspinor32))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }
#elif (defined BGL && defined _USE_BGLDRAM && !defined BGP)
  rts_return = rts_get_dram_window(8*((VOLUME+RAND)+1)*sizeof(halfspinor32), RTS_STORE_WITHOUT_ALLOCATE, (void**)&HalfSpinor32_, &actualSize);
  if(rts_return !=0) {
    return(-1);
  }
#else
  if((void*)(HalfSpinor32_ = (halfspinor32*)calloc(8*(VOLUME+RAND)+1, sizeof(halfspinor32))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(-1);
  }
#endif

#if ( defined SSE || defined SSE2 || defined SSE3)
  HalfSpinor32 = (halfspinor32*)(((unsigned long int)(HalfSpinor32_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  HalfSpinor32 = HalfSpinor32_;
#endif

  for(ieo = 0; ieo < 2; ieo++) {
    for(i = 0; i < VOLUME/2; i++) {
      j = g_eo2lexic[i + ((ieo+1)%2)*(VOLUME+RAND)/2];
      /* get (t,x,y,z) from j */
      t = j/(LX*LY*LZ);
      x = (j-t*(LX*LY*LZ))/(LY*LZ);
      y = (j-t*(LX*LY*LZ)-x*(LY*LZ))/(LZ);
      z = (j-t*(LX*LY*LZ)-x*(LY*LZ) - y*LZ);
      for(mu = 0; mu < 4; mu++) {
	NBPointer32[ieo][8*i + 2*mu + 0] = &HalfSpinor32[ 8*g_lexic2eosub[ g_idn[j][mu] ] + 2*mu + 0];
	NBPointer32[ieo][8*i + 2*mu + 1] = &HalfSpinor32[ 8*g_lexic2eosub[ g_iup[j][mu] ] + 2*mu + 1];
      }
#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(t == 0) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_idn[j][0]] - VOLUME/2);
	NBPointer32[ieo][8*i] = &HalfSpinor32[ k ];
      }
      if(t == T-1) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_iup[j][0]] - VOLUME/2);
	NBPointer32[ieo][8*i + 1] = &HalfSpinor32[ k ];
      }
#endif
#if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT) || (defined PARALLELXYZ))
      if(x == 0) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_idn[j][1]] - VOLUME/2);
	NBPointer32[ieo][8*i + 2] = &HalfSpinor32[ k ];
      }
      if(x == LX-1) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_iup[j][1]] - VOLUME/2);
	NBPointer32[ieo][8*i + 3] = &HalfSpinor32[ k ];
      }
#endif
#if ((defined PARALLELXYT) || (defined PARALLELXYZT) || (defined PARALLELXYZ))
      if(y == 0) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_idn[j][2]] - VOLUME/2);
	NBPointer32[ieo][8*i + 4] = &HalfSpinor32[ k ];
      }
      if(y == LY-1) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_iup[j][2]] - VOLUME/2);
	NBPointer32[ieo][8*i + 5] = &HalfSpinor32[ k ];
      }
#endif
#if (defined PARALLELXYZT) || (defined PARALLELXYZ)
      if(z == 0) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_idn[j][3]] - VOLUME/2);
	NBPointer32[ieo][8*i + 6] = &HalfSpinor32[ k ];
      }
      if(z == LZ-1) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_iup[j][3]] - VOLUME/2);
	NBPointer32[ieo][8*i + 7] = &HalfSpinor32[ k ];
      }
#endif
    }
#ifdef MPI
/*     NBPointer32[ieo][4*VOLUME] = NBPointer32[ieo][0];  */
#endif
  }
  for(ieo = 2; ieo < 4; ieo++) {
    for(i = 0; i < VOLUME/2; i++) {
      j = g_eo2lexic[i + ((ieo+0)%2)*(VOLUME+RAND)/2];
      /* get (t,x,y,z) from j */
      t = j/(LX*LY*LZ);
      x = (j-t*(LX*LY*LZ))/(LY*LZ);
      y = (j-t*(LX*LY*LZ)-x*(LY*LZ))/(LZ);
      z = (j-t*(LX*LY*LZ)-x*(LY*LZ) - y*LZ);
      for(mu = 0; mu < 8; mu++) {
	NBPointer32[ieo][8*i + mu] = &HalfSpinor32[8*i + mu];
      }
#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(t == T-1) {
	NBPointer32[ieo][8*i]     = &HalfSpinor32[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_iup[j][0] ] - VOLUME/2)];
      }
      if(t == 0) {
	NBPointer32[ieo][8*i + 1] = &HalfSpinor32[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_idn[j][0] ] - VOLUME/2)];
      }
#endif
#if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT) || (defined PARALLELXYZ))
      if(x == LX-1) { 
	NBPointer32[ieo][8*i + 2] = &HalfSpinor32[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_iup[j][1] ] - VOLUME/2)];
      }
      if(x == 0) {
	NBPointer32[ieo][8*i + 3] = &HalfSpinor32[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_idn[j][1] ] - VOLUME/2)];
      }
#endif
#if ((defined PARALLELXYT) || (defined PARALLELXYZT) || (defined PARALLELXYZ))
      if(y == LY-1) {
	NBPointer32[ieo][8*i + 4] = &HalfSpinor32[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_iup[j][2] ] - VOLUME/2)];
      }
      if(y == 0) {
	NBPointer32[ieo][8*i + 5] = &HalfSpinor32[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_idn[j][2] ] - VOLUME/2)];
      }
#endif
#if (defined PARALLELXYZT) || (defined PARALLELXYZ)
      if(z == LZ-1) {
	NBPointer32[ieo][8*i + 6] = &HalfSpinor32[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_iup[j][3] ] - VOLUME/2)];
      }
      if(z == 0) {
	NBPointer32[ieo][8*i + 7] = &HalfSpinor32[ 4*VOLUME + RAND/2 + (g_lexic2eosub[ g_idn[j][3] ] - VOLUME/2)];
      }
#endif
    }
#ifdef MPI
/*     NBPointer32[ieo][4*VOLUME] = NBPointer32[ieo][0];  */
#endif
  }
  return(0);
}
