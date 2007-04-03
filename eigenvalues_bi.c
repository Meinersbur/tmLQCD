 /*******************************************************************************
 * $Id$
 *
 * Here we compute the nr_of_eigenvalues lowest eigenvalues
 * of (gamma5*D)^2. Therefore we use the arnoldi routines.
 * 
 * The computed eigenvalues are stored in g_eigenvalues
 * and the computed eigenvectors in g_ev
 * 
 * inout:
 *   nr_of_eigenvalues:      input:  Number of eigenvalues to compute
 *                           output: Number of computed eigenvalues
 * input:
 *   crylov_space_dimension: Dimension of crylov space dimension
 *                           to be used in the arnoldi routines
 *   operator_flag:          Choose if we want to use D_Wilson
 *                           or D_Overlap
 *
 * Autor: Thomas Chiarappa
 *        Thomas.Chiarappa@mib.infn.it
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "su3.h"
#include "linalg_eo.h"
#include "start.h"
#include "tm_operators.h"
#include "solver/solver.h"
#include "solver/jdher_bi.h"
#ifdef MPI
#include "solver/pjdher_bi.h"
#endif
#include "eigenvalues_bi.h"
#include "Nondegenerate_Matrix.h"


double eigenvalues_bi(int * nr_of_eigenvalues, const int operator_flag, 
		      const int max_iterations, const double precision,
		      const int maxmin) {


  /*
  static spinor * eigenvectors_ = NULL;
  */
  static bispinor * eigenvectors_bi_ = NULL;
  static int allocated = 0;
  bispinor  *temp_field, *temp_field_ = NULL, *aux, *aux_ = NULL;
  bispinor  *eigenvectors_bi = NULL;
  double * eigenvls_bi = NULL;


  /**********************
   * For Jacobi-Davidson 
   **********************/
  int verbosity = g_debug_level, converged = 0, blocksize = 1, blockwise = 0;
  int solver_it_max = 200, j_max, j_min; 
  /*int it_max = 10000;*/
  complex *eigv_ = NULL, *eigv;
  double decay_min = 1.7, decay_max = 1.5, prec,
    threshold_min = 1.e-3, threshold_max = 5.e-2, 
    startvalue, threshold, decay, returnvalue;
/*   static int v0dim = 0; */
  int v0dim = 0;

  /**********************
   * General variables
   **********************/
  int returncode=0;

  int i, iVol, ix;

  FILE *conf_bifile=NULL;
  char * filename = NULL;
  char conf_bifilename[50];

  
  filename = calloc(200, sizeof(char));
  /*  strcpy(filename,optarg);*/

  if(maxmin == JD_MINIMAL) {
    startvalue = 0.;
    threshold = threshold_min;
    decay = decay_min;
    solver_it_max = 200;
  }
  else {
    startvalue = 50.;
    threshold = threshold_max;
    decay = decay_max;
    solver_it_max = 50;
  }

  if(g_proc_id == g_stdio_proc) {
    printf("Number of %s eigenvalues to compute = %d\n",
	   maxmin ? "maximal" : "minimal",(*nr_of_eigenvalues));
    printf("Using Jacobi-Davidson method! \n");
  }

  if((*nr_of_eigenvalues) < 8){
    j_max = 15;
    j_min = 8;
  }
  else{
    j_max = 2*(*nr_of_eigenvalues);
    j_min = (*nr_of_eigenvalues);
  }
  if(precision < 1.e-14){
    prec = 1.e-14;
  }
  else{
    prec = precision;
  }

  if(allocated == 0) {
    allocated = 1;
#if (defined SSE || defined SSE2 || defined SSE3)
    eigenvectors_bi_ = calloc((VOLUME)/2*(*nr_of_eigenvalues)+1, sizeof(bispinor)); 
    eigenvectors_bi = (bispinor *)(((unsigned long int)(eigenvectors_bi_)+ALIGN_BASE)&~ALIGN_BASE);

    temp_field_ = calloc((VOLUME)/2+1, sizeof(bispinor));
    temp_field = (bispinor *)(((unsigned long int)(temp_field_)+ALIGN_BASE)&~ALIGN_BASE);

    aux_ = calloc((VOLUME)/2+1, sizeof(bispinor));
    aux = (bispinor *)(((unsigned long int)(aux_)+ALIGN_BASE)&~ALIGN_BASE);
#else
    eigenvectors_bi_= calloc((VOLUME)/2*(*nr_of_eigenvalues), sizeof(bispinor));

    temp_field_ = calloc((VOLUME)/2, sizeof(bispinor));
    aux_ = calloc((VOLUME)/2, sizeof(bispinor));

    eigenvectors_bi = eigenvectors_bi_;
    temp_field = temp_field_;
    aux = aux_;
#endif
    eigenvls_bi = (double*)malloc((*nr_of_eigenvalues)*sizeof(double));
  }

  /* compute eigenvalues */

  if((g_proc_id==0) && (g_debug_level > 4)) {
    printf(" Values of   mu = %e     mubar = %e     eps = %e     precision = %e  \n \n", g_mu, g_mubar, g_epsbar, precision);
  }
 

#ifdef MPI
  pjdher_bi((VOLUME)/2*sizeof(bispinor)/sizeof(complex), (VOLUMEPLUSRAND)/2*sizeof(bispinor)/sizeof(complex),
	    startvalue, prec, 
	    (*nr_of_eigenvalues), j_max, j_min, 
	    max_iterations, blocksize, blockwise, v0dim, (complex*) eigenvectors_bi,
	    CG, solver_it_max,
	    threshold, decay, verbosity,
	    &converged, (complex*) eigenvectors_bi, eigenvls_bi,
	    &returncode, maxmin, 1,
	    &Q_Qdagger_ND_BI);
  
  /* IN THE LAST LINE, INSERT:
     Q_Qdagger_ND_BI;   Non-degenerate case - on 1 bispinor 
     Q_Qdagger_ND;      Non-degenerate case - on 2 spinors 
     Qtm_pm_psi;        Degenerate case  -  on 1 spinor 
  */

#else
  jdher_bi((VOLUME)/2*sizeof(bispinor)/sizeof(complex),
	   startvalue, prec, 
	   (*nr_of_eigenvalues), j_max, j_min, 
	   max_iterations, blocksize, blockwise, v0dim, (complex*) eigenvectors_bi,
	   CG, solver_it_max,
	   threshold_min, decay_min, verbosity,
	   &converged, (complex*) eigenvectors_bi, eigenvls_bi,
	   &returncode, maxmin, 1,
	   &Q_Qdagger_ND_BI);
  
  /* IN THE LAST LINE, INSERT:
     Q_Qdagger_ND_BI;   Non-degenerate case - on 1 bispinor 
     Q_Qdagger_ND;      Non-degenerate case - on 2 spinors 
     Qtm_pm_psi;        Degenerate case  -  on 1 spinor 
  */
  
#endif

  (*nr_of_eigenvalues) = converged;
  v0dim = converged;

  returnvalue = eigenvls_bi[0];
  return(returnvalue);
}
