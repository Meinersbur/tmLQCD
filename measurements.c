/***********************************************************************
 *  
 * Copyright (C) 2008 Carsten Urbach
 *
 * Adapted from monomial.c by Florian Burger 2009/12/16
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

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include "global.h"
#include "default_input_values.h"
#include "read_input.h"
#include "pion_norm.h"
#include "online_measurement.h"
#include "measurements.h"

measurement measurement_list[max_no_measurements];
int no_measurements = 0;



int add_measurement(const int type) {
 
  if(no_measurements == max_no_measurements) {
    fprintf(stderr, "maximal number of measurementss %d exceeded!\n", max_no_measurements);
    exit(-1);
  }
  measurement_list[no_measurements].measurefunc = &dummy_meas;
  measurement_list[no_measurements].initialised = 1;

  no_measurements++;
  return(no_measurements);
}




int init_measurements(){
 int i;
  for(i = 0; i < no_measurements; i++) {
 
    if(measurement_list[i].type == ONLINE) {
      measurement_list[i].measurefunc = &online_measurement;
      measurement_list[i].max_source_slice = g_nproc_t*T;
    }

    if(measurement_list[i].type == PIONNORM) {
      measurement_list[i].measurefunc = &pion_norm;
      measurement_list[i].max_source_slice = g_nproc_z*LZ;
    }  
    
    measurement_list[i].id = i;
 }
return(0);
}



void free_measurements(){

 return;
}



void dummy_meas(const int traj, const int slice) {
  if(g_proc_id == 0) {
    fprintf(stderr, "dummy_meas was called. Was that really intended?\n");
  }
  return;
}




