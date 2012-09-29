/***********************************************************************
 * $Id$
 *
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
 *
 * invert for twisted mass QCD
 *
 * Author: Carsten Urbach
 *         urbach@physik.fu-berlin.de
 *
 *******************************************************************************/

#define MAIN_PROGRAM

#include"lime.h"
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif

#include <execinfo.h> /* GG */
#include <unistd.h>   /* GG */
#include <sys/types.h> /* GG */
#include <sys/stat.h> /* GG */
#include <fcntl.h> /* GG */

/* #include "mem.h"  BGGP */

#include <stdlib.h>
#include <stdio.h>

/* GG */
extern FILE *fmemopen (void *__s, size_t __len, __const char *__modes) __THROW;

#include <math.h>
#include <time.h>
#include <string.h>
#include <signal.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "global.h"
#include "getopt.h"
#include "linalg_eo.h"
#include "geometry_eo.h"
#include "start.h"
/*#include "eigenvalues.h"*/
#include "observables.h"
#ifdef MPI
#include "xchange.h"
#endif
#include <io/utils.h>
#include "read_input.h"
#include "mpi_init.h"
#include "sighandler.h"
#include "boundary.h"
#include "solver/solver.h"
#include "init_gauge_field.h"
#include "init_geometry_indices.h"
#include "init_spinor_field.h"
#include "init_moment_field.h"
#include "init_dirac_halfspinor.h"
#include "init_bispinor_field.h"
#include "init_chi_spinor_field.h"
#include "xchange_halffield.h"
#include "stout_smear.h"
#include "invert_eo.h"
#include "monomial.h"
#include "ranlxd.h"
#include "phmc.h"
#include "D_psi.h"
#include "little_D.h"
#include "reweighting_factor.h"
#include "linalg/convert_eo_to_lexic.h"
#include "block.h"
#include "operator.h"
#include "sighandler.h"
#include "solver/dfl_projector.h"
#include "solver/generate_dfl_subspace.h"
#include "prepare_source.h"
#include <io/params.h>
#include <io/gauge.h>
#include <io/spinor.h>
#include <io/utils.h>
#include <assert.h>
#ifdef BGQ
#include "bgq/bgq_field_double.h"
#endif

void usage()
{
  fprintf(stdout, "Inversion for EO preconditioned Wilson twisted mass QCD\n");
  fprintf(stdout, "Version %s \n\n", PACKAGE_VERSION);
  fprintf(stdout, "Please send bug reports to %s\n", PACKAGE_BUGREPORT);
  fprintf(stdout, "Usage:   invert [options]\n");
  fprintf(stdout, "Options: [-f input-filename]\n");
  fprintf(stdout, "         [-o output-filename]\n");
  fprintf(stdout, "         [-v] verbose output (default no)\n");
  fprintf(stdout, "         [-h|-? this help]\n");
  exit(0);
}
/* GG */
/* #define GMALLOC yes   */
#if defined GMALLOC
#define __USE_GNU
#include <dlfcn.h>
#include <stdio.h>

void* malloc(size_t sz) {
    void *(*libc_malloc)(size_t) = dlsym(RTLD_NEXT, "malloc");
    printf("malloc\n");
    return libc_malloc(sz);
}

void free(void *p) {
    void (*libc_free)(void*) = dlsym(RTLD_NEXT, "free");
    printf("free\n");
    libc_free(p);
}
#undef __USE_GNU
#endif

/* GG */
/* #define XMALLOC yes  */
#if defined XMALLOC
void * malloc (size_t elsize)
{
  void * new_mem = malloc (elsize);
  
  if (new_mem == NULL)
    {
      printf ("xmalloc: request for 1 element of size %d failed.\n",
               elsize);
      abort ();
    }
  
  return new_mem;
}

void * calloc (size_t nelem, size_t elsize)
{
  void * new_mem = calloc (nelem, elsize);
  
  if (new_mem == NULL)
    {
      printf ("xcalloc: request for %d elements of size %d failed.\n",
               nelem, elsize);
      abort ();
    }
  
  return new_mem;
}
#endif

/* GG */
/* #define BACKTRACE yes */
#if defined BACKTRACE
     void
     print_trace (void)
     {
       /* GG */
       if ( g_proc_id ) 
	 return;

       void *array[20];
       size_t size;
       char **strings;
       size_t i;
     
       size = backtrace (array, 20);
       strings = backtrace_symbols (array, size);
     
       //printf ("Obtained %zd stack frames.\n", size);
     
       for (i = 0; i < size; i++)
	 if ( ! strstr(strings[i], "../invert [") && ! strstr(strings[i], "../invert(print_t") && ! strstr(strings[i], "/lib64/libc") )
	   printf ("%s==", strings[i]);
/*           printf ("%s\n", strings[i]); */
     
       free (strings);
     }
#else
     void print_trace (void) {}
#endif

/* GG */
/* #define SPINORPRT yes */
#if defined SPINORPRT
void print_r(spinor * const R, const int N){
  int ix, is;
  spinor *r; //,*s;
  double *p;

  for (ix = 0; ix < N; ix++){
    r=(spinor *) R + ix;
    for (is = 0; is < 24; is++){
      p = (double *)r + is;
      printf(" % 23.16e", *p);
    }    
    printf("\n");
  }
}
#else
     void print_r (spinor * const R, const int N) {}
#endif

/* GG */
double gatime, getime;

extern int nstore;
int check_geometry();

int main(int argc, char *argv[])
{
  FILE *parameterfile = NULL;
  int c, j, ix = 0, isample = 0, op_id = 0;
  char * filename = NULL;
  char datafilename[50];
  char parameterfilename[50];
  char conf_filename[50];
  char * input_filename = NULL;
  double plaquette_energy;
  double ratime, retime;

  int memlevel = 0;

  /* GG */
  char *mpiRank;
  char catInput[256];
  pid_t catCmdline;
  int mpi_rank = -1;
  int iterl=0;
  double nrm01, nrm45;
  double nrm0 , nrm4 ;
  /* GG */
  int save_prop_flag = 1;
  int traj_count = 696969;
  /* GG */
  int intrig;
  FILE* yyingg = NULL;
  MPI_Status yyStatus;
  void* yybufgg = NULL;
  int yyCount;
  int yyfd;

#ifdef _GAUGE_COPY
  int kb=0;
#endif
  double nrm1, nrm2;
#ifdef MPI
  double atime=0., etime=0.;
#endif
#ifdef _KOJAK_INST
#pragma pomp inst init
#pragma pomp inst begin(main)
#endif
  

#if (defined SSE || defined SSE2 || defined SSE3)
  signal(SIGILL, &catch_ill_inst);
#endif

  DUM_DERI = 8;
  /* DUM_DERI + 2 is enough (not 7) */
  DUM_SOLVER = DUM_DERI + 5;
  DUM_MATRIX = DUM_SOLVER + 8;
  /* DUM_MATRIX + 2 is enough (not 6) */
  NO_OF_SPINORFIELDS = DUM_MATRIX + 2;

  verbose = 0;
  g_use_clover_flag = 0;
  //g_nr_of_psf = 1;

  /* GG */
  buf_size = 1024;
  ilg_max = 1;
  source_node = 0;
  target_node = 1;

#ifdef MPI
  MPI_Init(&argc, &argv);
  /* GG */
  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);
#else
  g_proc_id = 0;
#endif

  /* GG modified ... */
  while ((c = getopt(argc, argv, "vh?f:o:b:i:s:t:")) != -1) {
    switch (c) {
    case 'f':
      input_filename = calloc(200, sizeof(char));
      strcpy(input_filename,optarg);
      break;
    case 'o':
      filename = calloc(200, sizeof(char));
      strcpy(filename,optarg);
      break;

      /* GG */
    case 'v':
      verbose = 1;
      break;
    case 'b':
      buf_size = atoi(optarg);
      break;
    case 'i':
      ilg_max = atoi(optarg);
      break;
    case 's':
      source_node = atoi(optarg);
      break;
    case 't':
      target_node = atoi(optarg);
      break;

    case 'h':
    case '?':
    default:
      usage();
      break;
    }
  }
  if (input_filename == NULL) {
    input_filename = "hmc.input";
  }
  if (filename == NULL) {
    filename = "output";
  } 

  /* GG */
  //system("printenv");
/*   mpiRank = getenv("OMPI_MCA_ns_nds_vpid"); */
/*   mpi_rank = atoi(mpiRank); */
/*     if(mpi_rank > 0) { */
/*       verbose = 0; */
/*     } */

  /* GG */
  if ( g_proc_id ) 
    verbose = 0;

  /* Read the input file */

  /* GG */
#define MPIO yes
#ifndef MPIO
  read_input(input_filename);
#else
  /* New read style with MPI_Bcast */
  yybufgg = (void *) malloc(8192*sizeof(char));
  yyingg = (FILE*) malloc(sizeof(FILE*));
  if (g_proc_id == 0) {
    yyfd = open(input_filename, O_RDONLY);
    yyCount = read(yyfd, yybufgg, 8192);
    intrig = close(yyfd);
  }
  intrig = MPI_Bcast(yybufgg, 8192, MPI_CHAR, 0, MPI_COMM_WORLD);
  yyingg = fmemopen(yybufgg, strlen(yybufgg), "r");
  intrig = read_input_fh(yyingg);
#endif

  int solver_flag = 0; // MK: Error, there can be multiple solvers
  if(solver_flag == 12 && even_odd_flag == 1) {
    even_odd_flag = 0;
    if(g_proc_id == 0) fprintf(stderr, "CGMMS works only without even/odd! Forcing!\n");
  }

  /* this DBW2 stuff is not needed for the inversion ! */
  if (g_dflgcr_flag == 1) {
    even_odd_flag = 0;
  }
  g_rgi_C1 = 0;
  if (Nsave == 0) {
    Nsave = 1;
  }

  if (g_running_phmc) {
    NO_OF_SPINORFIELDS = DUM_MATRIX + 8;
  }

  tmlqcd_mpi_init(argc, argv);

  /* GG */
  /* Bad on BG
  if (g_proc_id == 0) {
    catCmdline = getpid();
    sprintf(catInput, "( cat /proc/%d/cmdline ; echo ) | tr \"\\000\" \"\\t\" ", catCmdline);
    printf ("\n --- command line for pid: %d >%s< \n", catCmdline, catInput);
    system(catInput);
    printf ("\n --- input_file contents: %s --- \n", input_filename);
    sprintf(catInput, "cat -n %s", input_filename);
    system(catInput);
    printf   (" --- input_file contents: END --- \n\n");
  }
  */

  g_dbw2rand = 0;

  /* starts the single and double precision random number */
  /* generator                                            */
  start_ranlux(rlxd_level, random_seed);

#ifndef MPI
  g_dbw2rand = 0;
#endif

#ifdef _GAUGE_COPY
  j = init_gauge_field(VOLUMEPLUSRAND, 1);
#else
  j = init_gauge_field(VOLUMEPLUSRAND, 0);
#endif
  if (j != 0) {
    fprintf(stderr, "Not enough memory for gauge_fields! Aborting...\n");
    exit(-1);
  }
  j = init_geometry_indices(VOLUMEPLUSRAND);
  if (j != 0) {
    fprintf(stderr, "Not enough memory for geometry indices! Aborting...\n");
    exit(-1);
  }
  if (no_monomials > 0) {
    if (even_odd_flag) {
      j = init_monomials(VOLUMEPLUSRAND / 2, even_odd_flag);
    }
    else {
      j = init_monomials(VOLUMEPLUSRAND, even_odd_flag);
    }
    if (j != 0) {
      fprintf(stderr, "Not enough memory for monomial pseudo fermion  fields! Aborting...\n");
      exit(0);
    }
  }
  if (even_odd_flag) {
    j = init_spinor_field(VOLUMEPLUSRAND / 2, NO_OF_SPINORFIELDS);
  }
  else {
    j = init_spinor_field(VOLUMEPLUSRAND, NO_OF_SPINORFIELDS);
  }
  if (j != 0) {
    fprintf(stderr, "Not enough memory for spinor fields! Aborting...\n");
    exit(-1);
  }

  if (g_running_phmc) {
    j = init_chi_spinor_field(VOLUMEPLUSRAND / 2, 20);
    if (j != 0) {
      fprintf(stderr, "Not enough memory for PHMC Chi fields! Aborting...\n");
      exit(0);
    }
  }

  g_mu = g_mu1;
  if (g_cart_id == 0) {
    /*construct the filenames for the observables and the parameters*/
    strcpy(datafilename, filename);
    strcat(datafilename, ".data");
    strcpy(parameterfilename, filename);
    strcat(parameterfilename, ".para");

    parameterfile = fopen(parameterfilename, "w");
    write_first_messages(parameterfile, 1);
    fclose(parameterfile);
  }

  /* define the geometry */
  geometry();

  /* define the boundary conditions for the fermion fields */
  boundary(g_kappa);

  phmc_invmaxev = 1.;

  init_operators();

  /* this could be maybe moved to init_operators */
#ifdef _USE_HALFSPINOR
  /* GG for tiled version */
#  if ((defined SSE2)||(defined SSE3)||(defined BGL && defined XLC))
  j = init_dirac_halfspinor();
  if (j != 0) {
    fprintf(stderr, "Not enough memory for halffield! Aborting...\n");
    exit(-1);
  }
  if (g_sloppy_precision_flag == 1) {
    j = init_dirac_halfspinor32();
    if (j != 0)
    {
      fprintf(stderr, "Not enough memory for 32-Bit halffield! Aborting...\n");
      exit(-1);
    } else
    if(g_proc_id == 0) {
      printf("init_dirac_halfspinor32: fine \n"); fflush(stdout);
    }
  }
# else
  /* Introduced by OP in hmc_tm ... cloned ... */
j = init_dirac_halfspinor_tile(); 
if ( j!= 0) {
    fprintf(stderr, "Problem with tiled halffield! Aborting...\n");
    exit(-1);
  } else
    if(g_proc_id == 0) {
      printf("init_dirac_halfspinor_tile: fine \n"); fflush(stdout);
    }
# endif

#  if (defined _PERSISTENT)
  if (even_odd_flag)
    init_xchange_halffield();
#  endif
#endif


#ifdef BGQ
  	assert(even_odd_flag);
	bgq_init_gaugefield_double();
	bgq_init_spinorfields_double(NO_OF_SPINORFIELDS, g_running_phmc ? 20 : 0);
	bgq_hm_init_double();
#endif


  for (j = 0; j < Nmeas; j++) {
    sprintf(conf_filename, "%s.%.4d", gauge_input_filename, nstore);
    if (g_cart_id == 0) {
      printf("Reading gauge field from file %s\n", conf_filename);
      fflush(stdout);
    }
    //if( (j = read_gauge_field(conf_filename)) !=0) {
   //   fprintf(stderr, "error %d while reading gauge field from %s\n Aborting...\n", j, conf_filename);
    //  exit(-2);
   // }


    if (g_cart_id == 0) {
      printf("done!\n");
      fflush(stdout);
    }
#ifdef MPI
    xchange_gauge();
#endif
#if BGQ
    bgq_update_backward_gauge();
#endif


    /*compute the energy of the gauge field*/
    plaquette_energy = measure_gauge_action();

    if (g_cart_id == 0) {
      printf("The plaquette value is %e\n", plaquette_energy / (6.*VOLUME*g_nproc));
      fflush(stdout);
    }

    if (use_stout_flag == 1) {
      if (stout_smear_gauge_field(stout_rho , stout_no_iter) != 0)
        exit(1) ;

      plaquette_energy = measure_gauge_action();

      if (g_cart_id == 0) {
        printf("The plaquette value after stouting is %e\n", plaquette_energy / (6.*VOLUME*g_nproc));
        fflush(stdout);
      }
    }

    if (reweighting_flag == 1) {
      reweighting_factor(reweighting_samples, nstore);
    }

    /* Compute minimal eigenvalues, if wanted */
    if (compute_evs != 0) {
      eigenvalues(&no_eigenvalues, 5000, eigenvalue_precision,
                  0, compute_evs, nstore, even_odd_flag);
    }
    if (phmc_compute_evs != 0) {
#ifdef MPI
      MPI_Finalize();
#endif
      return(0);
    }

    /* move to operators as well */
    if (g_dflgcr_flag == 1) {
      /* set up deflation blocks */
      init_blocks(1, 1, 2, 1, g_N_s);

      /* the can stay here for now, but later we probably need */
      /* something like init_dfl_solver called somewhere else  */
      /* create set of approximate lowest eigenvectors ("global deflation subspace") */

      /*       g_mu = 0.; */
      /*       boundary(0.125); */
      generate_dfl_subspace(g_N_s, VOLUME);
      /*       boundary(g_kappa); */
      /*       g_mu = g_mu1; */

      /* Compute little Dirac operators */
/*       alt_block_compute_little_D(); */

/* GG 
      if (g_debug_level > 0){
        check_projectors();
      }
*/
      if (g_debug_level > 1){
        check_little_D_inversion();
      }

    }

#if 0
    /* GG Special test deflation*/
    sprintf(conf_filename2,"%s.%.4d", gauge_input_filename, nstore+1);
    if (g_proc_id == 0){
      printf("Reading Gauge field from file2 %s\n", conf_filename2); fflush(stdout);
    }
    read_lime_gauge_field(conf_filename2);
    xlfmessage = read_message(conf_filename2, "xlf-info");
    gaugelfn = read_message(conf_filename2, "ildg-data-lfn");
    gaugecksum = read_message(conf_filename2, "scidac-checksum");
    if (g_proc_id == 0){
      printf("%s \n", gaugecksum);
      printf("done2!\n"); fflush(stdout);
    }
#ifdef MPI
    xchange_gauge();
#endif
    plaquette_energy = measure_gauge_action();
    if(g_proc_id == 0) {
      printf("The plaquette value2 is %e\n", plaquette_energy/(6.*VOLUME*g_nproc)); fflush(stdout);
    }
    /* GG End special test deflation */
#endif

    if(SourceInfo.type == 1) {
      index_start = 0;
      index_end = 1;
    }
    for(isample = 0; isample < no_samples; isample++) {
      for (ix = index_start; ix < index_end; ix++) {
	for(op_id = 0; op_id < no_operators; op_id++) {
	  /* we use g_spinor_field[0-7] for sources and props for the moment */
	  /* 0-3 in case of 1 flavour  */
	  /* 0-7 in case of 2 flavours */
	  prepare_source(nstore, isample, ix, op_id, 
			 read_source_flag,
			 source_location);

#ifdef BGQ
	  for (int i = 0; i < 7; i += 1) {
		  bgq_transfer_spinorfield_double(i&1, g_spinorfields_double[i], g_spinor_field[i]);
	  }
#endif

	  operator_list[op_id].inverter(op_id, index_start);
	}
      }
    }
    nstore += Nsave;
  }
  
#ifdef MPI
  MPI_Finalize();
#endif

  free_blocks();
  free_dfl_subspace();
  free_gauge_field();
  free_geometry_indices();
  free_spinor_field();
  free_moment_field();
  free_chi_spinor_field();
  return(0);
#ifdef _KOJAK_INST
#pragma pomp inst end(main)
#endif
}
