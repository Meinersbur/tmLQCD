/* $Id$ */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "su3.h"
#include "su3adj.h"
#include "expo.h"
#include "ranlxd.h"
#include "sse.h"
#include "global.h"
#include "linalg_eo.h"
#include "clover_eo.h"
#include "start.h"
#include "sw.h"
#include "linsolve.h"
#include "xchange.h"
#include "deriv_Sb.h"
#include "Hopping_Matrix.h"
#include "tm_operators.h"
#include "hybrid_update.h"
#include "derivative_psf.h"
#include "ext_integrator.h"

extern int count00, count01, count10, count11, count20, count21;

void ext_leap_frog(int * const n_int, const double tau, const int S, const int halfstep) {
  int i,j;
  double eps, eps0;

  if(S == g_nr_of_psf) {
    /* initialize the counter for the inverter */
    count00=0; count01=0; count10=0; count11=0; count20=0; count21=0;
    eps = tau/((double)n_int[S]);
    
    for(i = 1; i < S+1; i++) {
      update_fermion_momenta(0.5 * eps, S-i);
      eps /= ((double)n_int[S-i]);
    }
    gauge_momenta(0.5 * eps);
  }

  eps = tau/((double)n_int[S]);
  if(S == 1) {
    eps0 = eps/((double)n_int[0]);
    for(j = 1; j < n_int[S]; j++) {
      for(i = 0; i < n_int[0]; i++) {
	update_gauge(eps0); 
	gauge_momenta(eps0);
      }
      update_fermion_momenta(eps, S-1);
    }
    for(i = 1; i < n_int[0]; i++) {
      update_gauge(eps0); 
      gauge_momenta(eps0);
    }
    update_gauge(eps0); 
    if(halfstep != 1) {
      gauge_momenta(eps0);
      update_fermion_momenta(eps, S-1);
    }
  }
  else {
    for(i = 1; i < n_int[S]; i++){
      ext_leap_frog(n_int, eps, S-1, 0);
      update_fermion_momenta(eps, S-1);
    }
    if(S == g_nr_of_psf) {
      ext_leap_frog(n_int, eps, S-1, 1);
    }
    else ext_leap_frog(n_int, eps, S-1, halfstep);
    if(halfstep != 1 && S != g_nr_of_psf) {
      update_fermion_momenta(eps, S-1);
    }
  }

  if(S == g_nr_of_psf) {
    for(i = 1; i < S+1; i++) {
      update_fermion_momenta(0.5 * eps, S-i);
      eps /= ((double)n_int[S-i]);
    }
    gauge_momenta(0.5 * eps);
  }
}


void ext_sexton_weingarten2(int * const n_int, const double tau, const int S, const int halfstep) {
  int i,j;
  double eps, eps0;

  if(S == g_nr_of_psf) {
    /* initialize the counter for the inverter */
    count00=0; count01=0; count10=0; count11=0; count20=0; count21=0;
    eps = tau/((double)n_int[S]);
    
    for(i = 1; i < S; i++) {
      update_fermion_momenta(0.5 * eps, S-i);
      eps /= ((double)n_int[S-i]);
    }
    update_fermion_momenta(eps/6., 0);
    eps /= ((double)n_int[0])*2.;
    gauge_momenta(eps/6.);
  }
  
  eps = tau/((double)n_int[S]);
  if(S == 1) {
    eps0 = eps/((double)n_int[0]);
    
    for(i = 1; i < n_int[S]; i++) {
      for(j = 0; j < n_int[0]; j++) {
	update_gauge(eps0/4.);
	gauge_momenta(eps0/3.);
	update_gauge(eps0/4.);
	gauge_momenta(eps0/6.);
      }
      update_fermion_momenta(2.*eps/3., S-1);
      for(j = 0; j < n_int[0]; j++) {
	update_gauge(eps0/4.);
	gauge_momenta(eps0/3.);
	update_gauge(eps0/4.);
	gauge_momenta(eps0/6.);
      }
      update_fermion_momenta(eps/3., S-1);
    }
    for(j = 0; j < n_int[0]; j++) {
      update_gauge(eps0/4.);
      gauge_momenta(eps0/3.);
      update_gauge(eps0/4.);
      gauge_momenta(eps0/6.);
    }
    update_fermion_momenta(2.*eps/3., S-1);
    for(j = 1; j < n_int[0]; j++) {
      update_gauge(eps0/4.);
      gauge_momenta(eps0/3.);
      update_gauge(eps0/4.);
      gauge_momenta(eps0/6.);
    }
    update_gauge(eps0/4.);
    gauge_momenta(eps0/3.);
    update_gauge(eps0/4.);
    if(halfstep != 1) {
      gauge_momenta(eps0/6.);
      update_fermion_momenta(eps/3., S-1);
    }
  }
  else {
    for(i = 1; i < n_int[S]; i++){
      ext_sexton_weingarten(n_int, eps, S-1, 0);
      update_fermion_momenta(eps, S-1);
    }
    if(S == g_nr_of_psf) {
      ext_sexton_weingarten(n_int, eps, S-1, 1);
    }
    else ext_sexton_weingarten(n_int, eps, S-1, halfstep);
    if(halfstep != 1 && S != g_nr_of_psf) {
      update_fermion_momenta(eps, S-1);
    }
  }

  if(S == g_nr_of_psf) {
    for(i = 1; i < S; i++) {
      update_fermion_momenta(0.5 * eps, S-i);
      eps /= ((double)n_int[S-i]);
    }
    update_fermion_momenta(eps/6., 0);
    eps /= ((double)n_int[0])*2.;
    gauge_momenta(eps/6.);
  }
}


void ext_sexton_weingarten(int * const n_int, const double tau, const int S, const int halfstep) {
  int i,j;
  double eps, eps0;

  if(S == g_nr_of_psf) {
    /* initialize the counter for the inverter */
    count00=0; count01=0; count10=0; count11=0; count20=0; count21=0;

    eps = tau/((double)n_int[S]);

    for(i = 1; i < S; i++) {
      update_fermion_momenta(eps/6., S-i);
      eps /= ((double)n_int[S-i])*2;
    }
    update_fermion_momenta(eps/6., 0);
    eps /= ((double)n_int[0])*2;
    gauge_momenta(eps/6.);
    
  }
  
  eps = tau/((double)n_int[S]);
  if(S == 1) {
    eps0 = eps/((double)n_int[0]);
    
    for(i = 1; i < n_int[S]; i++) {
      for(j = 0; j < n_int[0]; j++) {
	update_gauge(eps0/4.);
	gauge_momenta(eps0/3.);
	update_gauge(eps0/4.);
	gauge_momenta(eps0/6.);
      }
      update_fermion_momenta(2.*eps/3., S-1);
      for(j = 0; j < n_int[0]; j++) {
	update_gauge(eps0/4.);
	gauge_momenta(eps0/3.);
	update_gauge(eps0/4.);
	gauge_momenta(eps0/6.);
      }
      update_fermion_momenta(eps/3., S-1);
    }
    for(j = 0; j < n_int[0]; j++) {
      update_gauge(eps0/4.);
      gauge_momenta(eps0/3.);
      update_gauge(eps0/4.);
      gauge_momenta(eps0/6.);
    }
    update_fermion_momenta(2.*eps/3., S-1);
    for(j = 1; j < n_int[0]; j++) {
      update_gauge(eps0/4.);
      gauge_momenta(eps0/3.);
      update_gauge(eps0/4.);
      gauge_momenta(eps0/6.);
    }
    update_gauge(eps0/4.);
    gauge_momenta(eps0/3.);
    update_gauge(eps0/4.);
    if(halfstep != 1) {
      gauge_momenta(eps0/6.);
      update_fermion_momenta(eps/3., S-1);
    }
  }
  else {
    for(i = 1; i < n_int[S]; i++){
      ext_sexton_weingarten(n_int, eps/2., S-1, 0);
      update_fermion_momenta(2.*eps/3., S-1);
      ext_sexton_weingarten(n_int, eps/2., S-1, 0);
      update_fermion_momenta(eps/3., S-1);
    }
    ext_sexton_weingarten(n_int, eps/2., S-1, 0);
    update_fermion_momenta(2.*eps/3., S-1);
    if(S == g_nr_of_psf) {
      ext_sexton_weingarten(n_int, eps/2., S-1, 1);
    }
    else ext_sexton_weingarten(n_int, eps/2., S-1, halfstep);
    if(halfstep != 1 && S != g_nr_of_psf) {
      update_fermion_momenta(eps/3., S-1);
    }
  }

  if(S == g_nr_of_psf) {
    for(i = 1; i < S; i++) {
      update_fermion_momenta(eps/6., S-i);
      eps /= ((double)n_int[S-i])*2;
    }
    update_fermion_momenta(eps/6., 0);
    eps /= ((double)n_int[0])*2;
    gauge_momenta(eps/6.);
  }
}
