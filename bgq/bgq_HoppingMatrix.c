/*
 * bgq_HoppingMatrix.c
 *
 *  Created on: Jul 27, 2012
 *      Author: meinersbur
 */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif

#include "bgq_field.h"
#include "bgq.h"
#include "global.h"
#include "boundary.h"
#include <mpi.h>
#include <omp.h>


#define SETUP_PERSISTENT_SEND(ISODD, DIM, DIR) \
	MPI_CHECK(MPI_Rsend_init(&sendsurface->dim ## DIM[(ISODD)][(-(DIR)+1)/3], sizeof(sendsurface->dim ## DIM[(ISODD)][(-(DIR)+1)/3]), MPI_BYTE, g_cart_shift[joinDimdir((DIM),(DIR))], g_cart_tag[(ISODD)][joinDimdir((DIM),(DIR))], g_cart, &g_cart_send_requests[(ISODD)][joinDimdir((DIM),(DIR))]))

#define SETUP_PERSISTENT_RECV(ISODD, DIM, DIR) \
	MPI_CHECK(MPI_Recv_init(&recvsurface->dim ## DIM[(ISODD)][(-(DIR)+1)/3], sizeof(recvsurface->dim ## DIM[(ISODD)][(-(DIR)+1)/3]), MPI_BYTE, g_cart_shift[joinDimdir((DIM),(DIR))], g_cart_tag[(ISODD)][joinDimdir((DIM),(DIR))], g_cart, &g_cart_recv_requests[(ISODD)][joinDimdir((DIM),(DIR))]))


void *malloc_aligned(size_t size, size_t alignment) {
	void *result=NULL;
	int errcode = posix_memalign(&result,alignment, size) ;
	if (errcode != 0) {
		fprintf(stderr ,"malloc returned %d\n", errcode);
		exit(10);
	}
	return result;
}


bgq_weylfield_double weylxchange_xup_send_double;
bgq_weylfield_double weylxchange_xup_recv_double;
bgq_weylfield_double weylxchange_xdown_send_double;
bgq_weylfield_double weylxchange_xdown_recv_double;
bgq_weylfield_double weylxchange_yup_send_double;
bgq_weylfield_double weylxchange_yup_recv_double;
bgq_weylfield_double weylxchange_ydown_send_double;
bgq_weylfield_double weylxchange_ydown_recv_double;
bgq_weylfield_double weylxchange_zup_send_double;
bgq_weylfield_double weylxchange_zup_recv_double;
bgq_weylfield_double weylxchange_zdown_send_double;
bgq_weylfield_double weylxchange_zdown_recv_double;

void bgq_init() {
	weylxchange_xup_send_double = (bgq_weylfield_double)malloc_aligned(sizeof(*weylxchange_xup_send_double) * PHYSICAL_LY *PHYSICAL_LZ * PHYSICAL_LTV,32);
	weylxchange_xup_recv_double = (bgq_weylfield_double)malloc_aligned(sizeof(*weylxchange_xup_recv_double) * PHYSICAL_LY *PHYSICAL_LZ * PHYSICAL_LTV,32);
	weylxchange_xdown_send_double = (bgq_weylfield_double)malloc_aligned(sizeof(*weylxchange_xdown_send_double) * PHYSICAL_LY *PHYSICAL_LZ * PHYSICAL_LTV,32);
	weylxchange_xdown_recv_double = (bgq_weylfield_double)malloc_aligned(sizeof(*weylxchange_xdown_recv_double) * PHYSICAL_LY *PHYSICAL_LZ * PHYSICAL_LTV,32);

	weylxchange_yup_send_double = (bgq_weylfield_double)malloc_aligned(sizeof(*weylxchange_yup_send_double) * PHYSICAL_LX *PHYSICAL_LZ * PHYSICAL_LTV,32);
	weylxchange_yup_recv_double = (bgq_weylfield_double)malloc_aligned(sizeof(*weylxchange_yup_recv_double) * PHYSICAL_LX *PHYSICAL_LZ * PHYSICAL_LTV,32);
	weylxchange_ydown_send_double = (bgq_weylfield_double)malloc_aligned(sizeof(*weylxchange_ydown_send_double) * PHYSICAL_LX *PHYSICAL_LZ * PHYSICAL_LTV,32);
	weylxchange_ydown_recv_double = (bgq_weylfield_double)malloc_aligned(sizeof(*weylxchange_ydown_recv_double) * PHYSICAL_LX *PHYSICAL_LZ * PHYSICAL_LTV,32);

	weylxchange_zup_send_double = (bgq_weylfield_double)malloc_aligned(sizeof(*weylxchange_zup_send_double) * PHYSICAL_LX *PHYSICAL_LY * PHYSICAL_LTV,32);
	weylxchange_zup_recv_double = (bgq_weylfield_double)malloc_aligned(sizeof(*weylxchange_zup_recv_double) * PHYSICAL_LX *PHYSICAL_LY * PHYSICAL_LTV,32);
	weylxchange_zdown_send_double = (bgq_weylfield_double)malloc_aligned(sizeof(*weylxchange_zdown_send_double) * PHYSICAL_LX *PHYSICAL_LY * PHYSICAL_LTV,32);
	weylxchange_zdown_recv_double = (bgq_weylfield_double)malloc_aligned(sizeof(*weylxchange_zdown_recv_double) * PHYSICAL_LX *PHYSICAL_LY * PHYSICAL_LTV,32);
}

void bgq_destroy() {
	free(weylxchange_xup_send_double);
	free(weylxchange_xup_recv_double);
	free(weylxchange_xdown_send_double);
	free(weylxchange_xdown_recv_double);
	free(weylxchange_yup_send_double);
	free(weylxchange_yup_recv_double);
	free(weylxchange_ydown_send_double);
	free(weylxchange_ydown_recv_double);
	free(weylxchange_zup_send_double);
	free(weylxchange_zup_recv_double);
	free(weylxchange_zdown_send_double);
	free(weylxchange_zdown_recv_double);
}




#define BGQ_HM_NOFUNC 1
#define BGQ_HM_TLINE_NOFUNC 1
#define BGQ_HM_SITE_NOFUNC 1
#define BGQ_HM_DIR_NOFUNC 1

#define BGQ_HM_BORDER_NOFUNC 1
#define BGQ_HM_BORDERDIST_NOFUNC 1

//#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

void bgq_HoppingMatrix_double(bool isOdd, bgq_spinorfield_double spinorfield, bgq_spinorfield_double targetfield, bgq_gaugefield_double gaugefield) {
	MPI_Request request_recv_xup;
	MPI_CHECK(MPI_Irecv(weylxchange_xup_recv_double, sizeof(*weylxchange_xup_recv_double) * PHYSICAL_LY * PHYSICAL_LZ * PHYSICAL_LTV, MPI_BYTE, g_nb_x_up, X_UP, MPI_COMM_WORLD, &request_recv_xup));
	MPI_Request request_recv_xdown;
	MPI_CHECK(MPI_Irecv(weylxchange_xdown_recv_double, sizeof(*weylxchange_xdown_recv_double) * PHYSICAL_LY * PHYSICAL_LZ * PHYSICAL_LTV, MPI_BYTE, g_nb_x_dn,X_DOWN, MPI_COMM_WORLD, &request_recv_xdown));

	MPI_Request request_recv_yup;
	MPI_CHECK(MPI_Irecv(weylxchange_yup_recv_double, sizeof(*weylxchange_yup_recv_double) * PHYSICAL_LX * PHYSICAL_LZ * PHYSICAL_LTV, MPI_BYTE,g_nb_y_up, Y_UP, MPI_COMM_WORLD, &request_recv_yup));
	MPI_Request request_recv_ydown;
	MPI_CHECK(MPI_Irecv(weylxchange_ydown_recv_double, sizeof(*weylxchange_ydown_recv_double) * PHYSICAL_LX * PHYSICAL_LZ * PHYSICAL_LTV, MPI_BYTE, g_nb_y_dn, Y_DOWN,MPI_COMM_WORLD, &request_recv_ydown));

	MPI_Request request_recv_zup;
	MPI_CHECK(MPI_Irecv(weylxchange_zup_recv_double, sizeof(*weylxchange_zup_recv_double) * PHYSICAL_LX * PHYSICAL_LY * PHYSICAL_LTV, MPI_BYTE, g_nb_z_up,Z_UP, MPI_COMM_WORLD, &request_recv_zup));
	MPI_Request request_recv_zdown;
	MPI_CHECK(MPI_Irecv(weylxchange_zdown_recv_double, sizeof(*weylxchange_zdown_recv_double) * PHYSICAL_LX * PHYSICAL_LY * PHYSICAL_LTV, MPI_BYTE, g_nb_z_dn,Z_DOWN, MPI_COMM_WORLD, &request_recv_zdown));

	// Load constants
	bgq_vector4double_decl(qka0);
	bgq_cconst(qka0,creal(ka0),cimag(ka0));
	bgq_vector4double_decl(qka1);
	bgq_cconst(qka1,creal(ka1),cimag(ka1));
	bgq_vector4double_decl(qka2);
	bgq_cconst(qka2,creal(ka2),cimag(ka2));
	bgq_vector4double_decl(qka3);
	bgq_cconst(qka3,creal(ka3),cimag(ka3));

	MPI_CHECK(MPI_Barrier(g_cart_grid)); // To ensure that all ranks started the receive requests (necessary? how expensive is this?)




	// xup border weyl computation
#pragma omp parallel for schedule(static)
	for (int xyz = 0; xyz < 1*(PHYSICAL_LY)*(PHYSICAL_LZ)*(PHYSICAL_LTV); xyz += 1) {
		const int x = PHYSICAL_LX;
		const int y = xyz % (PHYSICAL_LY);
		xyz = xyz / (PHYSICAL_LY);
		const int z = xyz % (PHYSICAL_LZ);
		xyz = xyz / (PHYSICAL_LZ);
		const int tv = xyz;

#define BGQ_HM_XUP_WEYL_SEND 1
#include "bgq_HoppingMatrix_xup.inc.c"
	}

	MPI_Request request_send_xup;
	MPI_CHECK(MPI_Isend(weylxchange_xup_send_double, sizeof(*weylxchange_xup_recv_double) * PHYSICAL_LY * PHYSICAL_LZ * PHYSICAL_LTV, MPI_BYTE, g_nb_x_up, X_UP, MPI_COMM_WORLD, &request_send_xup));



	// xdown border weyl computation
#pragma omp parallel for schedule(static)
	for (int xyz = 0; xyz < 1*(PHYSICAL_LY)*(PHYSICAL_LZ)*(PHYSICAL_LTV); xyz += 1) {
		const int x = 0;
		const int y = xyz % (PHYSICAL_LY);
		xyz = xyz / (PHYSICAL_LY);
		const int z = xyz % (PHYSICAL_LZ);
		xyz = xyz / (PHYSICAL_LZ);
		const int tv = xyz;

#define BGQ_HM_XDOWN_WEYL_SEND 1
#include "bgq_HoppingMatrix_xdown.inc.c"
	}

	MPI_Request request_send_xdown;
	MPI_CHECK(MPI_Isend(weylxchange_xdown_send_double, sizeof(*weylxchange_xdown_recv_double) * PHYSICAL_LY * PHYSICAL_LZ * PHYSICAL_LTV, MPI_BYTE, g_nb_x_dn, X_DOWN, MPI_COMM_WORLD, &request_send_xdown));


	// yup border weyl computation
#pragma omp parallel for schedule(static)
	for (int xyz = 0; xyz < (PHYSICAL_LX)*(1)*(PHYSICAL_LZ)*(PHYSICAL_LTV); xyz += 1) {
		const int x = xyz % (PHYSICAL_LX);
		xyz = xyz / PHYSICAL_LX;
		const int y = PHYSICAL_LY-1;
		const int z = xyz % (PHYSICAL_LZ);
		xyz = xyz / (PHYSICAL_LZ);
		const int tv = xyz;

#define BGQ_HM_YUP_WEYL_SEND 1
#include "bgq_HoppingMatrix_yup.inc.c"
	}

	MPI_Request request_send_yup;
	MPI_CHECK(MPI_Isend(weylxchange_yup_send_double, sizeof(*weylxchange_yup_recv_double) * PHYSICAL_LX * 1 * PHYSICAL_LZ * PHYSICAL_LTV, MPI_BYTE, g_nb_y_up, Y_UP, MPI_COMM_WORLD, &request_send_yup));

	// ydown border weyl computation
#pragma omp parallel for schedule(static)
	for (int xyz = 0; xyz < (PHYSICAL_LX)*(1)*(PHYSICAL_LZ)*(PHYSICAL_LTV); xyz += 1) {
		const int x = xyz % (PHYSICAL_LX);
		xyz = xyz / PHYSICAL_LX;
		const int y = 0;
		const int z = xyz % (PHYSICAL_LZ);
		xyz = xyz / (PHYSICAL_LZ);
		const int tv = xyz;

#define BGQ_HM_YDOWN_WEYL_SEND 1
#include "bgq_HoppingMatrix_ydown.inc.c"
	}
	MPI_Request request_send_ydown;
	MPI_CHECK(MPI_Isend(weylxchange_ydown_send_double, sizeof(*weylxchange_ydown_recv_double) * PHYSICAL_LX * 1 * PHYSICAL_LZ * PHYSICAL_LTV, MPI_BYTE, g_nb_y_dn, Y_DOWN, MPI_COMM_WORLD, &request_send_ydown));


	// zup border weyl computation
#pragma omp parallel for schedule(static)
	for (int xyz = 0; xyz < (PHYSICAL_LX)*(PHYSICAL_LY)*(1)*(PHYSICAL_LTV); xyz += 1) {
		const int x = xyz % (PHYSICAL_LX);
		xyz = xyz / PHYSICAL_LX;
		const int y = xyz % (PHYSICAL_LY);
		xyz = xyz / (PHYSICAL_LY);
		const int z = PHYSICAL_LZ-1;
		const int tv = xyz;

#define BGQ_HM_ZUP_WEYL_SEND 1
#include "bgq_HoppingMatrix_zup.inc.c"
	}
	MPI_Request request_send_zup;
	MPI_CHECK(MPI_Isend(weylxchange_zup_send_double, sizeof(*weylxchange_zup_send_double) * PHYSICAL_LX * PHYSICAL_LY * (1) * PHYSICAL_LTV, MPI_BYTE, g_nb_z_up, Z_UP, MPI_COMM_WORLD, &request_send_zup));


	// zdown border weyl computation
#pragma omp parallel for schedule(static)
	for (int xyz = 0; xyz < (PHYSICAL_LX)*(PHYSICAL_LY)*(1)*(PHYSICAL_LTV); xyz += 1) {
		const int x = xyz % (PHYSICAL_LX);
		xyz = xyz / PHYSICAL_LX;
		const int y = xyz % (PHYSICAL_LY);
		xyz = xyz / (PHYSICAL_LY);
		const int z = 0;
		const int tv = xyz;

#define BGQ_HM_ZDOWN_WEYL_SEND 1
#include "bgq_HoppingMatrix_zdown.inc.c"
	}
	MPI_Request request_send_zdown;
	MPI_CHECK(MPI_Isend(weylxchange_zdown_send_double, sizeof(*weylxchange_zdown_send_double) * PHYSICAL_LX * PHYSICAL_LY * (1) * PHYSICAL_LTV, MPI_BYTE, g_nb_z_dn, Z_DOWN, MPI_COMM_WORLD, &request_send_zdown));


	// Body volume kernel
#pragma omp parallel for schedule(static,1)
	for (int xyz = 0; xyz < (PHYSICAL_LX-2)*(PHYSICAL_LY-2)*(PHYSICAL_LZ-2)/2; xyz += 1) {
		const int z = xyz % (PHYSICAL_LZ-2)+1;
		xyz = xyz/(PHYSICAL_LZ-2);
		const int y = xyz % (PHYSICAL_LY-2)+1;
		xyz = xyz/(PHYSICAL_LY-2);
		const int x = xyz*2 + isOdd;

#define BGQ_HM_TLINE_FLUSHLINE 1
		#include "bgq_HoppingMatrix_tline.inc.c"
	}

#pragma omp parallel for schedule(static,1)
	for (int xyz = 0; xyz < (PHYSICAL_LX-2)*(PHYSICAL_LY-2)*(PHYSICAL_LZ-2)/2; xyz += 1) {
		const int z = xyz % (PHYSICAL_LZ-2)+1;
		xyz = xyz/(PHYSICAL_LZ-2);
		const int y = xyz % (PHYSICAL_LY-2)+1;
		xyz = xyz/(PHYSICAL_LY-2);
		const int x = xyz*2 + !isOdd;

#define BGQ_HM_TLINE_RAGGEDLINE 1
		#include "bgq_HoppingMatrix_tline.inc.c"
	}



	MPI_Status weylxchange_xup_recv_status;
	MPI_CHECK(MPI_Wait(&request_recv_xup, &weylxchange_xup_recv_status));
	assert(get_MPI_count(&weylxchange_xup_recv_status) == sizeof(*weylxchange_xup_recv_double)*(1)*PHYSICAL_LY*PHYSICAL_LZ*PHYSICAL_LTV);

	MPI_Status weylxchange_xdown_recv_status;
	MPI_CHECK(MPI_Wait(&request_recv_xdown, &weylxchange_xdown_recv_status));
	assert(get_MPI_count(&weylxchange_xdown_recv_status) == sizeof(*weylxchange_xdown_recv_double)*(1)*PHYSICAL_LY*PHYSICAL_LZ*PHYSICAL_LTV);

	MPI_Status weylxchange_yup_recv_status;
	MPI_CHECK(MPI_Wait(&request_recv_yup, &weylxchange_yup_recv_status));
	assert(get_MPI_count(&weylxchange_yup_recv_status) == sizeof(*weylxchange_yup_recv_double)*PHYSICAL_LX*(1)*(PHYSICAL_LZ)*(PHYSICAL_LTV));

	MPI_Status weylxchange_ydown_recv_status;
	MPI_CHECK(MPI_Wait(&request_recv_ydown, &weylxchange_ydown_recv_status));
	assert(get_MPI_count(&weylxchange_ydown_recv_status) == sizeof(*weylxchange_ydown_recv_double)*PHYSICAL_LX*(1)*(PHYSICAL_LZ)*(PHYSICAL_LTV));

	MPI_Status weylxchange_zup_recv_status;
	MPI_CHECK(MPI_Wait(&request_recv_zup, &weylxchange_zup_recv_status));
	assert(get_MPI_count(&weylxchange_zup_recv_status) == sizeof(*weylxchange_zup_recv_double)*PHYSICAL_LX*(PHYSICAL_LY)*(1)*(PHYSICAL_LTV));

	MPI_Status weylxchange_zdown_recv_status;
	MPI_CHECK(MPI_Wait(&request_recv_zdown, &weylxchange_zdown_recv_status));
	assert(get_MPI_count(&weylxchange_zdown_recv_status) == sizeof(*weylxchange_zdown_recv_double)*PHYSICAL_LX*(PHYSICAL_LY)*(1)*(PHYSICAL_LTV));





#define BGQ_HM_BORDER_XUP 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_XDOWN 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_YUP 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_YDOWN 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_ZUP 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_ZDOWN 1
#include "bgq_HoppingMatrix_border.inc.c"



#define BGQ_HM_BORDER_XUP 1
#define BGQ_HM_BORDER_YUP 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_XUP 1
#define BGQ_HM_BORDER_YDOWN 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_XUP 1
#define BGQ_HM_BORDER_ZUP 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_XUP 1
#define BGQ_HM_BORDER_ZDOWN 1
#include "bgq_HoppingMatrix_border.inc.c"


#define BGQ_HM_BORDER_XDOWN 1
#define BGQ_HM_BORDER_YUP 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_XDOWN 1
#define BGQ_HM_BORDER_YDOWN 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_XDOWN 1
#define BGQ_HM_BORDER_ZUP 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_XDOWN 1
#define BGQ_HM_BORDER_ZDOWN 1
#include "bgq_HoppingMatrix_border.inc.c"


#define BGQ_HM_BORDER_YUP 1
#define BGQ_HM_BORDER_ZUP 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_YUP 1
#define BGQ_HM_BORDER_ZDOWN 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_YDOWN 1
#define BGQ_HM_BORDER_ZUP 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_YDOWN 1
#define BGQ_HM_BORDER_ZDOWN 1
#include "bgq_HoppingMatrix_border.inc.c"


#define BGQ_HM_BORDER_XUP 1
#define BGQ_HM_BORDER_YUP 1
#define BGQ_HM_BORDER_ZUP 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_XUP 1
#define BGQ_HM_BORDER_YUP 1
#define BGQ_HM_BORDER_ZDOWN 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_XUP 1
#define BGQ_HM_BORDER_YDOWN 1
#define BGQ_HM_BORDER_ZUP 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_XUP 1
#define BGQ_HM_BORDER_YDOWN 1
#define BGQ_HM_BORDER_ZDOWN 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_XDOWN 1
#define BGQ_HM_BORDER_YUP 1
#define BGQ_HM_BORDER_ZUP 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_XDOWN 1
#define BGQ_HM_BORDER_YUP 1
#define BGQ_HM_BORDER_ZDOWN 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_XDOWN 1
#define BGQ_HM_BORDER_YDOWN 1
#define BGQ_HM_BORDER_ZUP 1
#include "bgq_HoppingMatrix_border.inc.c"

#define BGQ_HM_BORDER_XDOWN 1
#define BGQ_HM_BORDER_YDOWN 1
#define BGQ_HM_BORDER_ZDOWN 1
#include "bgq_HoppingMatrix_border.inc.c"

}

#undef BGQ_HM_NOFUNC
#undef BGQ_HM_TLINE_NOFUNC
#undef BGQ_HM_SITE_NOFUNC
#undef BGQ_HM_DIR_NOFUNC

