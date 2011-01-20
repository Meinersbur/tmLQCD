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
 *
 * Hybrid-Monte-Carlo for twisted mass QCD
 *
 * Author: Carsten Urbach
 *         urbach@physik.fu-berlin.de
 *
 * Modified by Jenifer Gonzalez Lopez for the Schroedinger Functional
 *
 *******************************************************************************/
#define MAIN_PROGRAM
#include "lime.h"
#if HAVE_CONFIG_H
#include<config.h>
#endif

#include <execinfo.h> /* GG */
#include <unistd.h>   /* GG */
#include <malloc.h>   /* GG */
#include <sys/resource.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

/* GG */
#include <errno.h>
#define ERRNO_PRINT {if ( errno ) {int gerrsv=errno; printf (" errno value: %d at %s:%d \n", gerrsv, __FILE__, __LINE__); fflush(stdout); errno = 0;}}

#include <stdlib.h>
#include <stdio.h>

/* GG Works for both SGI Altix (icc) and BG/L (xlc) */
extern FILE *fmemopen (void *__s, size_t __len, __const char *__modes) __THROW;

#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <signal.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include <io/params.h>
#include <io/gauge.h>
#include "getopt.h"
#include "ranlxd.h"
#include "geometry_eo.h"
#include "start.h"
#include "observables.h"
#include "measure_rectangles.h"
#ifdef MPI
# include "xchange.h"
#endif
#include "read_input.h"
#include "mpi_init.h"
#include "sighandler.h"
#include "update_tm.h"
#include "init_gauge_field.h"
#include "init_geometry_indices.h"
#include "init_spinor_field.h"
#include "init_moment_field.h"
#include "init_gauge_tmp.h"
#include "init_dirac_halfspinor.h"
#include "init_stout_smear_vars.h"
#include "init_bispinor_field.h"
#include "init_chi_spinor_field.h"
#include "xchange_halffield.h"
#include "test/check_geometry.h"
#include "boundary.h"
#include "phmc.h"
#include "solver/solver.h"
#include "monomial.h"
#include "integrator.h"
#include "sighandler.h"
#include "measurements.h"
#include "sf_calc_action.h"
#include "sf_observables.h"
#include "linsolve.h" //MK

//#define __NON_INSTRUMENT_FUNCTION__    __attribute__((__no_instrument_function__))

void usage(){
  fprintf(stdout, "HMC for Wilson twisted mass QCD\n");
  fprintf(stdout, "Version %s \n\n", PACKAGE_VERSION);
  fprintf(stdout, "Please send bug reports to %s\n", PACKAGE_BUGREPORT);
  fprintf(stdout, "Usage:   hmc_tm [options]\n");
  fprintf(stdout, "Options: [-f input-filename]  default: hmc.input\n");
  fprintf(stdout, "         [-o output-filename] default: output\n");
  fprintf(stdout, "         [-v] more verbosity\n");
  fprintf(stdout, "         [-h|-? this help]\n");
  exit(0);
}

/* GG */
int callback_function(const void *pentry, size_t sz, int useflag, int status,
                      const char *filename, size_t line)
{
  /*
  if (_HEAPOK != status) {
    puts("status is not _HEAPOK.");
    exit(status);
  }
  */
  /*
  if (_USEDENTRY == useflag)
    printf("allocated  %p     %u\n", pentry, sz);
  else
    printf("freed      %p     %u\n", pentry, sz);
  */
  printf("allocated  %p useflag %d  size   %u\n", pentry, useflag, sz);

  return 0;
}

/* GG */
//#define BACKTRACE yes
#if defined BACKTRACE
     void
     //__NON_INSTRUMENT_FUNCTION__
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

       //GG printf ("Obtained %zd stack frames.\n", size);

       for (i = 0; i < size; i++)
	 if ( ! strstr(strings[i], "../phmc_tm [") && ! strstr(strings[i], "../phmc_tm(print_t") && ! strstr(strings[i], "/lib64/tls/libc") )
	   printf ("%s==", strings[i]);
       printf ("\n");
/*           printf ("%s\n", strings[i]); */

       free (strings);
     }

     void
     //__NON_INSTRUMENT_FUNCTION__
     print_tracel (char* line, char* file)
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

       //GG printf ("Obtained %zd stack frames.\n", size);

       for (i = 0; i < size; i++)
	 if ( ! strstr(strings[i], "../phmc_tm [") && ! strstr(strings[i], "../phmc_tm(print_t") && ! strstr(strings[i], "/lib64/tls/libc") )
	   printf ("%s==", strings[i]);
       printf (" line %s file %s \n", line, file);
       /*
       printf ("\n");
       */
/*           printf ("%s\n", strings[i]); */

       free (strings);
     }
#else
     void
     //__NON_INSTRUMENT_FUNCTION__
print_trace (void) {}

     void
     //__NON_INSTRUMENT_FUNCTION__
print_tracel (char* a, char* b) {}
#endif

/* GG */
void tmlqcd_mpi_close (void) {
#ifdef MPI
  MPI_Finalize();
#endif
  abort();
}

/* GG */
double gatime=0.0, getime=0.0;

extern int nstore;

int const rlxdsize = 105;

int main(int argc,char *argv[]) {

  FILE *parameterfile=NULL, *countfile=NULL;
  char *filename = NULL;
  char datafilename[50];
  char parameterfilename[50];
  char gauge_filename[50];
  char nstore_filename[50];
  char tmp_filename[50];
  char *input_filename = NULL;

  int j,ix,mu, trajectory_counter=1;
  struct timeval t1;

  /* GG */
  double  atime=0.0,  etime=0.0;
  int intrig;
  FILE* yyingg = NULL;
  MPI_File myyingg = MPI_FILE_NULL;
  MPI_Info minfo = 0;
  MPI_Status yyStatus;
  void* yybufgg = NULL;
  int yyCount;
  int yyfd;

  /* Energy corresponding to the Gauge part */
  double eneg = 0., plaquette_energy = 0., rectangle_energy = 0.;
  /* Acceptance rate */
  int Rate=0;
  /* Do we want to perform reversibility checks */
  /* See also return_check_flag in read_input.h */
  int return_check = 0;
  /* For getopt */
  int c;

  paramsXlfInfo *xlfInfo;

/* For online measurements */
  measurement * meas;
  int imeas;

#ifdef _KOJAK_INST
#pragma pomp inst init
#pragma pomp inst begin(main)
#endif

  /* GG */
  pid_t gpid;
  char ghost[256];
  char gmall[256];
  char gmalc[256];
  char gmalv[256];
  char gmalt[256];

  //GLOB char pnametrajGlob[256];
  gpid = getpid();
  gethostname(ghost, 256);
  snprintf(gmall, 256, "malloc-ptrace.%s.%d", ghost, gpid);
  setenv("MALLOC_TRACE", gmall, 1);

#if (defined SSE || defined SSE2 || SSE3)
  signal(SIGILL,&catch_ill_inst);
#endif

  /* GG */
  //GOOD mtrace();

  strcpy(gauge_filename,"conf.save");
  strcpy(nstore_filename,".nstore_counter");
  strcpy(tmp_filename, ".conf.tmp");

  verbose = 0;
  g_use_clover_flag = 0;

#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);
#else
  g_proc_id = 0;
#endif

  if (g_proc_id == 0) {
    fprintf(stderr, "MK_Static memusage:\n");
    print_memusage(); // MK
  }

  /* GG */
  /*
  if (g_proc_id==0)
    mtrace();
  */
  /*
  if (g_proc_id!=0)
    muntrace();
  */

  while ((c = getopt(argc, argv, "h?vf:o:")) != -1) {
    switch (c) {
    case 'f':
      input_filename = calloc(200, sizeof(char));
      strcpy(input_filename,optarg);
      break;
    case 'o':
      filename = calloc(200, sizeof(char));
      strcpy(filename,optarg);
      break;
    case 'v':
      verbose = 1;
      break;
    case 'h':
    case '?':
    default:
      usage();
      break;
    }
  }
  if(input_filename == NULL){
    input_filename = "hmc.input";
  }
  if(filename == NULL){
    filename = "output";
  }

  /* Read the input file */

  /* GG */
#define MPIO yes
#ifndef MPIO
  if( (j = read_input(input_filename)) != 0) {
    fprintf(stderr, "Could not find input file: %s\nAborting...\n", input_filename);
    exit(-1);
  }
#else
  /* Proto ...
  int MPI_File_open(MPI_Comm comm, char *filename, int amode,
          MPI_Info info, MPI_File *fh)
  */

  if (g_proc_id == 0) {
    system("df -k");
    system("ls -l /tmp");
  }
#if 1
  /* New read style with MPI_Bcast */
    yybufgg = (void *) malloc(8192*sizeof(char));
    yyingg = (FILE*) malloc(sizeof(FILE*));
  if (g_proc_id == 0) {
    yyfd = open(input_filename, O_RDONLY);
    //yybufgg = (void *) malloc(8192*sizeof(char));
    yyCount = read(yyfd, yybufgg, 8192);
    intrig = close(yyfd);
    /*
    intrig = MPI_Bcast(yybufgg, yyCount, MPI_CHAR, 0, MPI_COMM_WORLD);
    yyingg = fmemopen(yybufgg, strlen(yybufgg), "r");
    intrig = read_input_fh(yyingg);
    */
  }
    intrig = MPI_Bcast(yybufgg, 8192, MPI_CHAR, 0, MPI_COMM_WORLD);
    //???? intrig = MPI_Barrier(MPI_COMM_WORLD);
    //intrig = MPI_Recv(yybufgg, 8192, MPI_CHAR, 0, ??tag, MPI_COMM_WORLD, &yyStatus);
    yyingg = fmemopen(yybufgg, strlen(yybufgg), "r");
    intrig = read_input_fh(yyingg);

#else
  yyingg = (MPI_File*)malloc(sizeof(MPI_File));
  yybufgg = (void *) malloc(8192*sizeof(MPI_CHAR));
  //MPI_Comm_dup(MPI_COMM_WORLD,
  //intrig = MPI_File_open(MPI_COMM_WORLD, input_filename, MPI_MODE_RDONLY, MPI_INFO_NULL, yyingg);
  intrig = MPI_File_open(MPI_COMM_WORLD, input_filename, MPI_MODE_RDONLY, MPI_INFO_NULL, yyingg);
  /* Proto
     int MPI_File_read_all(MPI_File fh, void *buf, int count,
          MPI_Datatype datatype, MPI_Status *status)
  */
  intrig = MPI_File_read_all(*yyingg, yybufgg, 8192, MPI_CHAR, &yyStatus);
  //Good?? intrig = MPI_Get_count(&yyStatus, MPI_CHAR, &yyCount);
  intrig = MPI_File_close(yyingg);
  //Good printf(" Control %s bufsize %d \n", gmall, yyCount);
  //system("ls -l /tmp");

  //yyinfg = (FILE *) yyingg;
  //////intrig = read_input_fh((FILE *) yyingg);
  //?? MPI_MODE_SEQUENTIAL
  intrig = MPI_File_open(MPI_COMM_WORLD, gmall, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, yyingg);
  /*Proto
    int MPI_File_write(MPI_File fh, void *buf, int count,
          MPI_Datatype datatype, MPI_Status *status)
  */
  intrig = MPI_File_write(*yyingg, yybufgg, strlen(yybufgg), MPI_CHAR, &yyStatus);
  intrig = MPI_File_close(yyingg);
  if (g_proc_id == 0)
    system(gmalc);
  intrig = read_input(gmall);
#endif

  /* Special for tests */
  /*
  MPI_Finalize();
  return(0);
  */
#endif

  DUM_DERI = 6;
  DUM_SOLVER = DUM_DERI+8;
  DUM_MATRIX = DUM_SOLVER+6;
  if(g_running_phmc) {
    NO_OF_SPINORFIELDS = DUM_MATRIX+8;
    //NO_OF_SPINORFIELDS = DUM_MATRIX+11; /* GG det0 */
    //NO_OF_SPINORFIELDS = DUM_MATRIX+13; /* GG detr1 */
    //NO_OF_SPINORFIELDS = DUM_MATRIX+15; /* GG detr2 */  //MK: Change suggested by GG
  }
  else {
    NO_OF_SPINORFIELDS = DUM_MATRIX+6;
  }
  DUM_BI_DERI = 6;
  DUM_BI_SOLVER = DUM_BI_DERI+7;

  DUM_BI_MATRIX = DUM_BI_SOLVER+6;
  NO_OF_BISPINORFIELDS = DUM_BI_MATRIX+6;

  tmlqcd_mpi_init(argc, argv);

  init_integrator();

  //atexit(&tmlqcd_mpi_close);

  if(nstore == -1) {
    if (g_proc_id == 0) { /* GG */
    countfile = fopen(nstore_filename, "r");
    if(countfile != NULL) {
      j = fscanf(countfile, "%d %d %s\n", &nstore, &trajectory_counter, gauge_input_filename);
      if(j < 1) nstore = 0;
      if(j < 2) trajectory_counter = 0;
      fclose(countfile);
    }
    else {
      nstore = 0;
      trajectory_counter = 0;
    }

    /* GG */
    }
    intrig = MPI_Bcast(gauge_input_filename, 100, MPI_CHAR, 0, MPI_COMM_WORLD);
    intrig = MPI_Bcast(&nstore, 1, MPI_INT, 0, MPI_COMM_WORLD);
    intrig = MPI_Bcast(&trajectory_counter, 1, MPI_INT, 0, MPI_COMM_WORLD);

  }

  if(g_rgi_C1 == 0.) {
    g_dbw2rand = 0;
  }
#ifndef MPI
  g_dbw2rand = 0;
#endif

  if(g_proc_id == 0) {
    for(j = 0; j < no_monomials; j++) {
      printf("# monomial id %d type = %d timescale %d\n", j, monomial_list[j].type, monomial_list[j].timescale);
    }
  }

  g_mu = g_mu1;

#ifdef _GAUGE_COPY
  j = init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 1);
#else
  j = init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 0);
#endif
  if (j != 0) {
    fprintf(stderr, "Not enough memory for gauge_fields! Aborting...\n");
    exit(0);
  }
  j = init_geometry_indices(VOLUMEPLUSRAND + g_dbw2rand);
  if (j != 0) {
    fprintf(stderr, "Not enough memory for geometry_indices! Aborting...\n");
    exit(0);
  }

  /* GG */
  //if( g_rec_ev == 0 ) {  //MK: Change suggested by GG
  if(even_odd_flag) {
    j = init_monomials(VOLUMEPLUSRAND/2, even_odd_flag);
  }
  else {
    j = init_monomials(VOLUMEPLUSRAND, even_odd_flag);
  }
  if (j != 0) {
    fprintf(stderr, "Not enough memory for monomial pseudo fermion  fields! Aborting...\n");
    exit(0);
  }
  //}  //MK: Change suggested by GG
  /* */

  if(even_odd_flag) {
    j = init_spinor_field(VOLUMEPLUSRAND/2, NO_OF_SPINORFIELDS);
  }
  else {
    j = init_spinor_field(VOLUMEPLUSRAND, NO_OF_SPINORFIELDS);
  }
  if (j != 0) {
    fprintf(stderr, "Not enough memory for spinor fields! Aborting...\n");
    exit(0);
  }

  /* GG */
  //if( g_rec_ev == 0 ) { //MK: Change suggested by GG
  if(even_odd_flag) {
    j = init_csg_field(VOLUMEPLUSRAND/2);
  }
  else {
    j = init_csg_field(VOLUMEPLUSRAND);
  }
  if (j != 0) {
    fprintf(stderr, "Not enough memory for csg fields! Aborting...\n");
    exit(0);
  }

  j = init_moment_field(VOLUME, VOLUMEPLUSRAND);
  if (j != 0) {
    fprintf(stderr, "Not enough memory for moment fields! Aborting...\n");
    exit(0);
  }
  //} //MK: Change suggested by GG
  /* */

  if(g_running_phmc) {
    j = init_bispinor_field(VOLUME/2, NO_OF_BISPINORFIELDS);
    if (j!= 0) {
      fprintf(stderr, "Not enough memory for Bispinor fields! Aborting...\n");
      exit(0);
    }
  }


     /* list and initialize measurements*/
   if(g_proc_id == 0) {
    printf("\n");
    for(j = 0; j < no_measurements; j++) {
      printf("# measurement id %d, type = %d: Frequency %d\n", j, measurement_list[j].type, measurement_list[j].freq);
    }
   }
   init_measurements();

  zero_spinor_field(g_spinor_field[DUM_DERI+4],VOLUME);
  zero_spinor_field(g_spinor_field[DUM_DERI+5],VOLUME);
  zero_spinor_field(g_spinor_field[DUM_DERI+6],VOLUME);

  if(use_stout_flag == 1)
    init_stout_smear_vars(VOLUMEPLUSRAND, stout_no_iter);

  /*construct the filenames for the observables and the parameters*/
  strcpy(datafilename,filename);  strcat(datafilename,".data");
  strcpy(parameterfilename,filename);  strcat(parameterfilename,".para");

  if(g_proc_id == 0){
    parameterfile = fopen(parameterfilename, "a");
    write_first_messages(parameterfile, 0);
  }

  /* define the geometry */
  geometry();

  /* define the boundary conditions for the fermion fields */
  boundary(g_kappa);

  check_geometry();

  /* GG */
  //if( g_rec_ev == 0 ) { //MK: Change suggested by GG
#ifdef _USE_HALFSPINOR
  j = init_dirac_halfspinor();
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for halffield! Aborting...\n");
    exit(0);
  }
  if(g_sloppy_precision_flag == 1) {
    init_dirac_halfspinor32();
  }
#  if (defined _PERSISTENT)
  init_xchange_halffield();
#  endif
#endif
//}  //MK: Change suggested by GG
/* */

  /* Initialise random number generator */
  start_ranlux(rlxd_level, random_seed^nstore);

  /* Set up the gauge field */
  /* continue and restart */
  if(startoption==3 || startoption == 2) {
    if(g_proc_id == 0) {
      printf("# Reading Gauge field from file %s in %d Bit\n",
	     gauge_input_filename, gauge_precision_read_flag);
      fflush(stdout);
    }
    read_gauge_field(gauge_input_filename);

    if (g_proc_id == 0){
      printf("# done!\n"); fflush(stdout);
    }
  }
  else if (startoption == 1) {
    /* hot */
    random_gauge_field(reproduce_randomnumber_flag);
  }
  else if(startoption == 0) {
    /* cold */
    unit_g_gauge_field();
  }

    /* GG special test ecriture-lecture */
#if 1

  /*For parallelization: exchange the gaugefield */
#ifdef MPI
  xchange_gauge();
#endif

  /* GG debug on BG */
      if(g_proc_id==0) {
	//_heap_walk(callback_function);
      }

  if(g_running_phmc) init_phmc();

#endif


  /*********************************************************/
  /* impose SF bc in case it was chosen in the input file */
  /*******************************************************/

  if (bc_flag == 1) { /* if SF */
    dirichlet_boundary_conditions(g_Tbsf);
    dirichlet_boundary_conditions_spatial_links_to_one(g_Tbsf);
    /* sf_boundary_conditions_spatially_constant_abelian_field(g_Tbsf, g_eta); */
    fprintf(parameterfile,"# SF put boundary at time slice: g_Tbsf = %d \n",g_Tbsf);

    /* compute the energy of the gauge field for SF */
    if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) {
      /* NOTE: the factor (1./(2.*3.)) is due to the difference between	our normalisation and Carstens's normalisation
	 when defining the plaquette and rectangle functions */
      plaquette_energy = (1./(2.*3.))*measure_plaquette_sf_iwasaki(g_Tbsf, g_Cs, g_Ct, g_rgi_C0);
      rectangle_energy = (1./(2.*3.))*measure_rectangle_sf_iwasaki(g_Tbsf, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);
      eneg = plaquette_energy + rectangle_energy;
      /* print energy for SF */
      if(g_proc_id==0){
	fprintf(parameterfile,"# First plaquette value for SF: %14.12f \n", plaquette_energy/(6.*VOLUME*g_nproc));
	fprintf(parameterfile,"# First rectangle value for SF: %14.12f \n", rectangle_energy/(12.*VOLUME*g_nproc));
	printf("# First plaquette value for SF: %14.12f \n", plaquette_energy/(6.*VOLUME*g_nproc));
	printf("# First rectangle value for SF: %14.12f \n", rectangle_energy/(12.*VOLUME*g_nproc));
      }
    }
    else {
      /* NOTE: the factor (1./(2.*3.)) is due to the difference between	our normalisation and Carstens's normalisation
	 when defining the plaquette and rectangle functions */
      plaquette_energy = (1./(2.*3.))*measure_plaquette_sf_weights_improvement(g_Tbsf, g_Cs, g_Ct);
      eneg = plaquette_energy;
      /* print plaquette energy for SF */
      if(g_proc_id==0){
	fprintf(parameterfile,"# First plaquette value for SF: %14.12f \n", plaquette_energy/(6.*VOLUME*g_nproc));
	printf("# First plaquette value for SF: %14.12f \n", plaquette_energy/(6.*VOLUME*g_nproc));
      }
    }
  }
  else if (bc_flag == 0) { /*if PBC */
    plaquette_energy=measure_gauge_action();
    if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) {
      rectangle_energy = measure_rectangles();
      if(g_proc_id==0){
	fprintf(parameterfile,"# First rectangle value: %14.12f \n",rectangle_energy/(12.*VOLUME*g_nproc));
      }
    }
    eneg = g_rgi_C0 * plaquette_energy + g_rgi_C1 * rectangle_energy;

    if(g_proc_id == 0) {
      fprintf(parameterfile,"# First plaquette value: %14.12f \n", plaquette_energy/(6.*VOLUME*g_nproc));
      printf("# First plaquette value: %14.12f \n", plaquette_energy/(6.*VOLUME*g_nproc));
      fclose(parameterfile);
    }
  }


  /* set ddummy to zero */
  for(ix = 0; ix < VOLUME+RAND; ix++){
    for(mu=0; mu<4; mu++){
      ddummy[ix][mu].d1=0.;
      ddummy[ix][mu].d2=0.;
      ddummy[ix][mu].d3=0.;
      ddummy[ix][mu].d4=0.;
      ddummy[ix][mu].d5=0.;
      ddummy[ix][mu].d6=0.;
      ddummy[ix][mu].d7=0.;
      ddummy[ix][mu].d8=0.;
    }
  }

  if(g_proc_id == 0) {
    gettimeofday(&t1,NULL);
    countfile = fopen("history_hmc_tm", "a");
    fprintf(countfile, "!!! Timestamp %ld, Nsave = %d, g_mu = %e, g_mu1 = %e, g_mu_2 = %e, g_mu3 = %e, beta = %f, kappa = %f, C1 = %f, ",
	    t1.tv_sec, Nsave, g_mu, g_mu1, g_mu2, g_mu3, g_beta, g_kappa, g_rgi_C1);
    for(j = 0; j < Integrator.no_timescales; j++) {
      fprintf(countfile, "n_int[%d] = %d ", j, Integrator.no_mnls_per_ts[j]);
    }
    fprintf(countfile, "\n");
    fclose(countfile);

    /* GG */
#ifdef MPI
      atime = MPI_Wtime();
#endif
  }

  if (g_proc_id == 0) {
    fprintf(stderr, "MK_Init memusage:\n");
    print_memusage(); // MK
  }

  /* Loop for measurements */
  for(j = 0; j < Nmeas; j++) {
    if(g_proc_id == 0) {
      printf("# Starting trajectory no %d\n", trajectory_counter); fflush(stdout);
    }

    /* GG */
    snprintf(pnametrajGlob, 256, "ProcName %s Trajecount %d", ghost, trajectory_counter);
    /* We limit to first trajec only */
    if ( j )
      debug_detailGlob = 0;

    if(return_check_flag == 1 && trajectory_counter%return_check_interval == 0) return_check = 1;
    else return_check = 0;

    if (setjmp(longjmpenv) == 0) {
    /* GG special test ecriture-lecture */
#if 1
    Rate += update_tm(&plaquette_energy, &rectangle_energy, datafilename, return_check, Ntherm<trajectory_counter);
#endif
      if (g_proc_id == 0)
        fprintf(stderr, "MK_update_tm completed normally; Rate=%d\n", Rate);
    } else {
      if (g_proc_id == 0)
        fprintf(stderr, "MK_update_tm partially skipped; Rate=%d\n", Rate);
    }

    /* GG */
    if(g_proc_id == 0) {
#ifdef MPI
      etime = MPI_Wtime();
#endif
      errno = 0;
      printf("# Ending trajectory no %d loop_index %d Rate: %d in %e sec. (MPI_Wtime)\n", trajectory_counter, j+1, Rate, etime-atime); fflush(stdout);
      ERRNO_PRINT;
    }

#if 1 //MK
    if (g_proc_id == 0) {
        fprintf(stderr, "MK_Plaquette energy: %e\n", plaquette_energy/(6.*VOLUME*g_nproc));
        fprintf(stderr, "MK_Rectangle energy: %e\n", rectangle_energy);
    }
#endif

    /* Save gauge configuration all Nsave times */
    if((Nsave !=0) && (trajectory_counter%Nsave == 0) && (trajectory_counter!=0)) {
      sprintf(gauge_filename,"conf.%.4d", nstore);
      if(g_proc_id == 0) {
	errno = 0;
        countfile = fopen("history_hmc_tm", "a");
	ERRNO_PRINT;
	fprintf(countfile, "%.4d, measurement %d of %d, Nsave = %d, Plaquette = %e, trajectory nr = %d\n",
		nstore, j, Nmeas, Nsave, plaquette_energy/(6.*VOLUME*g_nproc),
		trajectory_counter);
	ERRNO_PRINT;
	fclose(countfile);
	ERRNO_PRINT;
      }
      nstore ++;
    }
    else {
      sprintf(gauge_filename,"conf.save");
    }

    //if(((Nsave !=0) && (trajectory_counter%Nsave == 0) && (trajectory_counter!=0)) || (write_cp_flag == 1) || (j >= (Nmeas - 1))) {
        if (g_proc_id == 0)
	  fprintf(stderr, "MK_Writing gauge configuration to %s\n", gauge_filename);
      /* Write the gauge configuration first to a temporary file */
/*       write_gauge_field_time_p( tmp_filename); */

      xlfInfo = construct_paramsXlfInfo(plaquette_energy/(6.*VOLUME*g_nproc), trajectory_counter);
      write_gauge_field( tmp_filename, gauge_precision_write_flag, xlfInfo);
      free(xlfInfo);


      /* Now move it! */
      if(g_proc_id == 0) {
	errno = 0;
  	rename(tmp_filename, gauge_filename);
	ERRNO_PRINT;
        countfile = fopen(nstore_filename, "w");
	ERRNO_PRINT;
        fprintf(countfile, "%d %d %s\n", nstore, trajectory_counter+1, gauge_filename);
	ERRNO_PRINT;
        fclose(countfile);
	ERRNO_PRINT;
      }
    //}

    /* online measurements */
    for(imeas=0; imeas<no_measurements; imeas++){
      meas = &measurement_list[imeas];
      if(trajectory_counter%meas->freq == 0){
        meas->measurefunc(trajectory_counter, imeas);
      }
    }

    if((g_rec_ev !=0) && (trajectory_counter%g_rec_ev == 0) && (g_running_phmc)) {
      /* GG */
      if (nstore%5 == 0)
	phmc_compute_ev(trajectory_counter, plaquette_energy);
        if (g_proc_id == 0)
	  fprintf(stderr, "MK_phmc_compute_ev %e\n", plaquette_energy);
    }

    if(g_proc_id == 0) {
      verbose = 1;
    }

    /* GG */
#if 0
#ifndef MPIO
  ix = reread_input("hmc.reread");
#else
  /* New reread style with MPI_Bcast */
  if (g_proc_id == 0) {
    if ( (yyfd = open("hmc.reread", O_RDONLY)) == NULL ) {
      ix = 2;
    } else {
      //ix = reread_input("hmc.reread");
      yyCount = read(yyfd, yybufgg, 8192);
      intrig = close(yyfd);
    }
    intrig = MPI_Bcast(yybufgg, 8192, MPI_CHAR, 0, MPI_COMM_WORLD);
    yyingg = fmemopen(yybufgg, strlen(yybufgg), "r");
    intrig = reread_input_fh(yyingg);
    ix = intrig;
  }
#endif
#else
  ix = 2;
#endif

    if(g_proc_id == 0) {
      verbose = 0;
    }

#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if(ix == 0 && g_proc_id == 0) {
      countfile = fopen("history_hmc_tm", "a");
      fprintf(countfile, "# Changed parameter according to hmc.reread: measurement %d of %d\n", j, Nmeas);
      fclose(countfile);
      printf("# Changed parameter according to hmc.reread (see stdout): measurement %d of %d\n", j, Nmeas);
      remove("hmc.reread");
    }

#if 0
    /* GG STOP feature */
  if (g_proc_id == 0) {
    strcpy(gmalv, getenv("PBS_JOBDIR"));
    snprintf(gmall, 256, "%s/%s", gmalv, "STOP");
    //yyfd = open(gmall, O_RDONLY);
#if 0
    if ( (yyingg = fopen(gmall, "r")) != (FILE*) NULL ) {
      intrig = fclose(yyingg);
      printf(" FOUND: %s \n", gmall); fflush(stdout);
      break;
    } else
      printf(" NOTFOUND ! \n");
#else
    intrig = MPI_File_open(MPI_COMM_WORLD, gmall, MPI_MODE_RDONLY, minfo, &myyingg);
    if ( myyingg != MPI_FILE_NULL ) {
      intrig = MPI_File_close(&myyingg);
      printf(" FOUND: %s \n", gmall); fflush(stdout);
      break;
    } else
      printf(" NOTFOUND ! \n");
#endif
  }
#endif

    print_memusage(); // MK

    trajectory_counter++;
  } /* end of loop over trajectories */

  if(g_proc_id==0) {
    printf("Acceptance rate was: %3.2f percent\n", 100.*(double)Rate/(double)Nmeas);
    /* GG */
#ifdef MPI
      etime = MPI_Wtime();
      printf("Simulation achieved for %d trajectories with %d accepted in %e sec. (MPI_Wtime)\n", Nmeas, Rate, etime-atime);
#endif
    fflush(stdout);
    parameterfile = fopen(parameterfilename, "a");
    fprintf(parameterfile, "Acceptance Rate was: %3.2f Percent\n", 100.*(double)Rate/(double)Nmeas);
    fclose(parameterfile);
  }

#if 1
  free_gauge_tmp();
  free_gauge_field();
  free_geometry_indices();
  free_spinor_field();
  free_moment_field();
  free_monomials();
  if(g_running_phmc) {
    free_bispinor_field();
    free_chi_spinor_field();
  }
#endif
  /* End IF PHMC */
/*   if(use_stout_flag == 1) */
/*     free_stout_smear_vars(); */

#ifdef MPI
  MPI_Finalize();
#endif
  return(0);
#ifdef _KOJAK_INST
#pragma pomp inst end(main)
#endif
}

static char const rcsid[] = "$Id$";
