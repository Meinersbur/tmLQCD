/***********************************************************************
 * Copyright (C) 2001 Carsten Urbach
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
#ifndef _LINSOLVE_H
#define _LINSOLVE_H

#include <setjmp.h> //MK

extern int callcount; //MK
extern jmp_buf longjmpenv; //MK


int solve_cg(spinor * const k, spinor * const l, double eps_sq, const int rel_prec);
int bicg(spinor * const k, spinor * const l, double eps_sq, const int rel_prec);

#endif
