/*
 * bgq_HoppingMatrix.c
 *
 *  Created on: Jul 27, 2012
 *      Author: meinersbur
 */

#include "bgq_field.h"
#include "bgq.h"
#include "global.h"
#include <mpi.h>




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


void bgq_HoppingMatrix_double(bool isOdd, bgq_spinorfield_double spinorfield, bgq_spinorfield_double targetfield) {
	MPI_Request request_recv_xup;
	MPI_CHECK(MPI_Irecv(weylxchange_xup_recv_double, sizeof(*weylxchange_xup_recv_double) * PHYSICAL_LY * PHYSICAL_LZ * PHYSICAL_LTV, MPI_BYTE, g_nb_x_up, MPI_COMM_WORLD, &request_recv_xup));
	MPI_Request request_recv_xdown;
	MPI_CHECK(MPI_Irecv(weylxchange_xdown_recv_double, sizeof(*weylxchange_xdown_recv_double) * PHYSICAL_LY * PHYSICAL_LZ * PHYSICAL_LTV, MPI_BYTE, g_nb_x_dn, MPI_COMM_WORLD, &request_recv_xdown));

	MPI_Request request_recv_yup;
	MPI_CHECK(MPI_Irecv(weylxchange_yup_recv_double, sizeof(*weylxchange_yup_recv_double) * PHYSICAL_LX * PHYSICAL_LZ * PHYSICAL_LTV, MPI_BYTE,g_nb_y_up,  MPI_COMM_WORLD, &request_recv_yup));
	MPI_Request request_recv_ydown;
	MPI_CHECK(MPI_Irecv(weylxchange_ydown_recv_double, sizeof(*weylxchange_ydown_recv_double) * PHYSICAL_LX * PHYSICAL_LZ * PHYSICAL_LTV, MPI_BYTE, g_nb_y_dn, MPI_COMM_WORLD, &request_recv_ydown));

	MPI_Request request_recv_zup;
	MPI_CHECK(MPI_Irecv(weylxchange_zup_recv_double, sizeof(*weylxchange_zup_recv_double) * PHYSICAL_LX * PHYSICAL_LY * PHYSICAL_LTV, MPI_BYTE, g_nb_z_up, MPI_COMM_WORLD, &request_recv_yup));
	MPI_Request request_recv_zdown;
	MPI_CHECK(MPI_Irecv(weylxchange_zdown_recv_double, sizeof(*weylxchange_zdown_recv_double) * PHYSICAL_LX * PHYSICAL_LY * PHYSICAL_LTV, MPI_BYTE, g_nb_z_dn, MPI_COMM_WORLD, &request_recv_ydown));

	MPI_CHECK(MPI_Barrier(g_cart_grid)); // To ensure that all ranks started the receive requests (necessary? how expensive is this?)

#pragma omp parallel

	// xup border weyl computation
#pragma omp for schedule(static) nowait
	for (int xyz = 0; xyz < 1*(PHYSICAL_LY)*(PHYSICAL_LZ)*(PHYSICAL_LTV); xyz += 1) {
		const int x = PHYSICAL_LX;
		const int y = y % (PHYSICAL_LY);
		xyz = xyz / (PHYSICAL_LY);
		const int z = xyz % (PHYSICAL_LZ);
		xyz = xyz / (PHYSICAL_LZ);
		const int tv = xyz;

#include "bgq_HoppingMatrix_xup.inc.c"
	}
#pragma omp barrier
#pragma omp master
	{
	MPI_Request request_send_xup;
	MPI_CHECK(MPI_Isend(weylxchange_xup_send_double, sizeof(*weylxchange_xup_recv_double) * PHYSICAL_LY * PHYSICAL_LZ * PHYSICAL_LTV, MPI_BYTE, g_nb_x_up,  MPI_COMM_WORLD, &request_send_xup));
	}



	// Body volume kernel
#pragma omp parallel for schedule(static)
	for (int xyz = 0; xyz < (PHYSICAL_LX-2)*(PHYSICAL_LY-2)*(PHYSICAL_LZ-2); xyz += 1) {
		const bool p = xyz % 2;
		xyz = xyz / 2;
		const int z = xyz % (PHYSICAL_LZ-2)+1;
		xyz = xyz/(PHYSICAL_LZ-2);
		const int y = xyz % (PHYSICAL_LY-2)+1;
		xyz = xyz/(PHYSICAL_LY-2)
		const int x = xyz*2 +p+1;


#include "bgq_HoppingMatrix_tline.inc.x"
	}
}



#undef BGQ_HM_NOFUNC
