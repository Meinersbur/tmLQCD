/*
 * bgq_HoppingMatrix.c
 *
 *  Created on: Jul 27, 2012
 *      Author: meinersbur
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "bgq_field.h"
#include "bgq.h"
#include "global.h"
#include "boundary.h"
#include <mpi.h>
#include <omp.h>

#define BGQ_HOPPINGMATRIX_C_
#include "bgq_HoppingMatrix.h"



#define SETUP_PERSISTENT_SEND(ISODD, DIM, DIR) \
	MPI_CHECK(MPI_Rsend_init(&sendsurface->dim ## DIM[(ISODD)][(-(DIR)+1)/3], sizeof(sendsurface->dim ## DIM[(ISODD)][(-(DIR)+1)/3]), MPI_BYTE, g_cart_shift[joinDimdir((DIM),(DIR))], g_cart_tag[(ISODD)][joinDimdir((DIM),(DIR))], g_cart, &g_cart_send_requests[(ISODD)][joinDimdir((DIM),(DIR))]))

#define SETUP_PERSISTENT_RECV(ISODD, DIM, DIR) \
	MPI_CHECK(MPI_Recv_init(&recvsurface->dim ## DIM[(ISODD)][(-(DIR)+1)/3], sizeof(recvsurface->dim ## DIM[(ISODD)][(-(DIR)+1)/3]), MPI_BYTE, g_cart_shift[joinDimdir((DIM),(DIR))], g_cart_tag[(ISODD)][joinDimdir((DIM),(DIR))], g_cart, &g_cart_recv_requests[(ISODD)][joinDimdir((DIM),(DIR))]))

bgq_weylfield_double weylxchange_recv_double[6];
bgq_weylfield_double weylxchange_send_double[6];
size_t weylxchange_size_double[3];
int weylexchange_destination[6];

void bgq_hm_init() {
	weylxchange_size_double[T_UP / 2] = PHYSICAL_LX * PHYSICAL_LY * PHYSICAL_LZV * sizeof(bgq_weylsite_double);
	weylxchange_size_double[X_UP / 2] = PHYSICAL_LT * PHYSICAL_LY * PHYSICAL_LZV * sizeof(bgq_weylsite_double);
	weylxchange_size_double[Y_UP / 2] = PHYSICAL_LT * PHYSICAL_LX * PHYSICAL_LZV * sizeof(bgq_weylsite_double);

	for (int d = T_UP; d <= Y_DOWN; d += 1) {
		size_t size = weylxchange_size_double[d / 2];
		weylxchange_recv_double[d] = (bgq_weylfield_double) malloc_aligned(size, 128);
		weylxchange_send_double[d] = (bgq_weylfield_double) malloc_aligned(size, 128);
	}

	weylexchange_destination[T_UP] = g_nb_t_up;
	weylexchange_destination[T_DOWN] = g_nb_t_dn;
	weylexchange_destination[X_UP] = g_nb_x_up;
	weylexchange_destination[X_DOWN] = g_nb_x_dn;
	weylexchange_destination[Y_UP] = g_nb_y_up;
	weylexchange_destination[Y_DOWN] = g_nb_y_dn;
}

void bgq_hm_free() {
	for (int d = T_UP; d <= Y_DOWN; d += 1) {
		free(weylxchange_recv_double[d]);
		weylxchange_recv_double[d] = NULL;
		free(weylxchange_send_double[d]);
		weylxchange_send_double[d] = NULL;
	}
}


#define BGQ_HM_BORDER_NOFUNC 1
#define BGQ_HM_BORDERDIST_NOFUNC 1

#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
//#pragma GCC diagnostic ignored "-Wunused-variable"


#define BGQ_HM_SUFFIX double
#include "bgq_HoppingMatrix.inc.c"

#define BGQ_HM_NOCOM 1
#define BGQ_HM_SUFFIX double_nocom
#include "bgq_HoppingMatrix.inc.c"


// Hopping_Matrix.c compatibility layer
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k) {
	bool isOdd = (ieo!=0);
	bgq_spinorfield_double target = bgq_translate_spinorfield(l);
	bgq_spinorfield_double source = bgq_translate_spinorfield(k);

	bgq_HoppingMatrix_double(isOdd, target, source, g_gaugefield_double);
}


void Hopping_Matrix_nocom(const int ieo, spinor * const l, spinor * const k) {
	bool isOdd = (ieo!=0);
	bgq_spinorfield_double target = bgq_translate_spinorfield(l);
	bgq_spinorfield_double source = bgq_translate_spinorfield(k);

	bgq_HoppingMatrix_double_nocom(isOdd, target, source, g_gaugefield_double);
}

