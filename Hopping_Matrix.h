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

#ifndef _HOPPING_MATRIX_H
#define _HOPPING_MATRIX_H

#define EO 0
#define OE 1
#define OO 1
#define EE 0

#include "su3.h"

void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k);
#if BGQ_REPLACE
void ref_Hopping_Matrix(const int ieo, spinor * const l, spinor * const k);
#endif

void printHopVerMsg();

#endif
