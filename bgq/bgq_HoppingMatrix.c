/*
 * bgq_HoppingMatrix.c
 *
 *  Created on: Jul 27, 2012
 *      Author: meinersbur
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bgq_field.h"
#include "bgq.h"
#include "global.h"
#include "boundary.h"
#include <mpi.h>
#include <omp.h>

#define BGQ_HOPPINGMATRIX_C_
#include "bgq_HoppingMatrix.h"




bgq_weylfield_double weylxchange_recv_double[6];
bgq_weylfield_double weylxchange_send_double[6];
size_t weylxchange_size_double[3];
int weylexchange_destination[6];

void bgq_hm_init() {
	weylxchange_size_double[TUP / 2] = PHYSICAL_LXV * PHYSICAL_LY * PHYSICAL_LZ * sizeof(bgq_weylsite_double);
	weylxchange_size_double[XUP / 2] = PHYSICAL_LTV * PHYSICAL_LY * PHYSICAL_LZ * sizeof(bgq_weylsite_double);
	weylxchange_size_double[YUP / 2] = PHYSICAL_LTV * PHYSICAL_LX * PHYSICAL_LZ * sizeof(bgq_weylsite_double);

	for (int d = TUP; d <= YDOWN; d += 1) {
		size_t size = weylxchange_size_double[d/2];
		weylxchange_recv_double[d] = (bgq_weylfield_double) malloc_aligned(size, 128);
		weylxchange_send_double[d] = (bgq_weylfield_double) malloc_aligned(size, 128);
	}

	weylexchange_destination[TUP] = g_nb_t_up;
	weylexchange_destination[TDOWN] = g_nb_t_dn;
	weylexchange_destination[XUP] = g_nb_x_up;
	weylexchange_destination[XDOWN] = g_nb_x_dn;
	weylexchange_destination[YUP] = g_nb_y_up;
	weylexchange_destination[YDOWN] = g_nb_y_dn;
}

void bgq_hm_free() {
	for (int d = TUP; d <= YDOWN; d += 1) {
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


#define HoppingMatrix bgq_HoppingMatrix_double
#include "bgq_HoppingMatrix.inc.c"

//#define BGQ_HM_NOCOM 1
//#define BGQ_HM_SUFFIX double_nocom
//#include "bgq_HoppingMatrix.inc.c"


// Hopping_Matrix.c compatibility layer
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k) {
	const bool isOdd = (ieo!=0);
	bgq_spinorfield_double target = bgq_translate_spinorfield(l);
	bgq_spinorfield_double source = bgq_translate_spinorfield(k);

	bgq_HoppingMatrix_double(isOdd, target, source, g_gaugefield_double);
}


void Hopping_Matrix_nocom(const int ieo, spinor * const l, spinor * const k) {
	const bool isOdd = (ieo!=0);
	bgq_spinorfield_double target = bgq_translate_spinorfield(l);
	bgq_spinorfield_double source = bgq_translate_spinorfield(k);

	//bgq_HoppingMatrix_double_nocom(isOdd, target, source, g_gaugefield_double);
}

