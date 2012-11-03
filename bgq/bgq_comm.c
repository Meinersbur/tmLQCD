/*
 * bgq_comm.c
 *
 *  Created on: Oct 25, 2012
 *      Author: meinersbur
 */

#define BGQ_COMM_C_
#include "bgq_comm.h"

#include "bgq_field.h"
#ifdef SPI
#include "DirectPut.h"
#endif




static MPI_Request g_bgq_request_recv[COMMDIR_COUNT];
static MPI_Request g_bgq_request_send[COMMDIR_COUNT];



static bgq_dimension bgq_commdim2dimension(ucoord commdim) {
	assert(0 <= commdim && commdim < COMMDIR_COUNT);
	if (COMM_T) {
		if (commdim == 0)
			return DIM_T;
		commdim -= 1;
	}
	if (COMM_X) {
		if (commdim == 0)
			return DIM_X;
		commdim -= 1;
	}
	if (COMM_Y) {
		if (commdim == 0)
			return DIM_Y;
		commdim -= 1;
	}
	assert(COMM_Z);
	assert(commdim==0);
	return DIM_Z;
}


static ucoord bgq_direction2commdir(bgq_direction d) {
	ucoord result = 0;
	if (COMM_T) {
		if (d==TUP)
			return result;
		result+=1;
		if (d==TDOWN)
			return result;
		result+=1;
	}
	if (COMM_X) {
		if(d==XUP)
			return result;
		result+=1;
		if (d==XDOWN)
			return result;
		result+=1;
	}
	if (COMM_Y) {
		if(d==YUP)
			return result;
		result+=1;
		if (d==YDOWN)
			return result;
		result+=1;
	}
	if (COMM_Z) {
		if(d==ZUP)
			return result;
		result+=1;
		if (d==ZDOWN)
			return result;
		result+=1;
	}
	UNREACHABLE
}


static ucoord bgq_commdir2direction(ucoord commdir) {
	if (COMM_T) {
		if (commdir==0)
			return TUP;
		commdir-=1;
		if (commdir==0)
			return TDOWN;
		commdir-=1;
	}
	if (COMM_X) {
		if (commdir==0)
			return XUP;
		commdir-=1;
		if (commdir==0)
			return XDOWN;
		commdir-=1;
	}
	if (COMM_Y) {
		if (commdir==0)
			return YUP;
		commdir-=1;
		if (commdir==0)
			return YDOWN;
		commdir-=1;
	}
	if (COMM_Z) {
		if (commdir==0)
			return ZUP;
		commdir-=1;
		if (commdir==0)
			return ZDOWN;
		commdir-=1;
	}
	UNREACHABLE
}


static int bgq_direction2rank(bgq_direction d) {
	switch (d) {
	case TUP:
	return g_nb_t_up;
	case TDOWN:
	return g_nb_t_dn;
	case XUP:
	return g_nb_x_up;
	case XDOWN:
	return g_nb_x_dn;
	case YUP:
	return g_nb_y_up;
	case YDOWN:
	return g_nb_y_dn;
	case ZUP:
	return g_nb_z_up;
	case ZDOWN:
	return g_nb_z_dn;
	default:
	UNREACHABLE
	}
}


void bgq_comm_init() {
	size_t commbufsize = bgq_weyl_section_offset(sec_comm_end) - bgq_weyl_section_offset(sec_comm);

	uint8_t *buf = (uint8_t*)malloc_aligned(commbufsize, BGQ_ALIGNMENT_L2);
	uint8_t *bufend = buf + commbufsize;
	g_bgq_sec_comm = buf;
	for (bgq_direction d = 0; d < PHYSICAL_LD; d+=1) {
		bgq_weylfield_section sec_send = bgq_direction2section(d, true);
		g_bgq_sec_send[d] = (bgq_weyl_vec*)(buf + bgq_weyl_section_offset(sec_send) - bgq_weyl_section_offset(sec_comm));
		assert((uint8_t*)g_bgq_sec_send[d] < bufend);
		bgq_weylfield_section sec_recv = bgq_direction2section(d, false);
		g_bgq_sec_recv[d] = (bgq_weyl_vec*)(buf + bgq_weyl_section_offset(sec_recv) - bgq_weyl_section_offset(sec_comm));
		assert((uint8_t*)g_bgq_sec_recv[d] < bufend);
	}


	for (ucoord d_src = 0; d_src < PHYSICAL_LD; d_src+=1) {
		bgq_direction d_dst = bgq_direction_revert(d_src);
		bgq_dimension dim = bgq_direction2dimension(d_src);
		if (!bgq_direction_isDistributed(d_src))
			continue;

		ucoord commdir_src = bgq_direction2commdir(d_src);
		bgq_weylfield_section sec_recv = bgq_direction2section(d_src, false);
		size_t secsize = bgq_weyl_section_offset(sec_recv+1) - bgq_weyl_section_offset(sec_recv);
		MPI_CHECK(MPI_Recv_init(g_bgq_sec_recv[d_src], secsize / sizeof(double), MPI_DOUBLE, bgq_direction2rank(d_src), d_dst, g_cart_grid, &g_bgq_request_recv[commdir_src]));

		ucoord commdir_dst = bgq_direction2commdir(d_dst);
		bgq_weylfield_section sec_send = bgq_direction2section(d_dst, true);
		assert(secsize == bgq_weyl_section_offset(sec_send+1) - bgq_weyl_section_offset(sec_send));
		MPI_CHECK(MPI_Send_init(g_bgq_sec_send[d_dst], secsize / sizeof(double), MPI_DOUBLE, bgq_direction2rank(d_dst), d_dst, g_cart_grid, &g_bgq_request_send[commdir_dst]));
	}


#ifdef SPI
	  // here comes the SPI initialization
	  uint64_t messageSizes[COMMDIR_COUNT];
	  uint64_t roffsets[COMMDIR_COUNT];
	  uint64_t soffsets[COMMDIR_COUNT];

	  spi_num_dirs = 2*(COMM_T+COMM_X+COMM_Y+COMM_Z);
	  size_t totalMessageSize = 0;
	for (ucoord d_src = 0; d_src < PHYSICAL_LD; d_src+=1) {
		bgq_direction d_dst = bgq_direction_revert(d_src);
		bgq_dimension dim = bgq_direction2dimension(d_src);
		if (!bgq_direction_isDistributed(d_src))
			continue;

		ucoord commdir_src = bgq_direction2commdir(d_src);
		bgq_weylfield_section sec_recv = bgq_direction2section(d_src, false);
		size_t secsize = bgq_weyl_section_offset(sec_recv+1) - bgq_weyl_section_offset(sec_recv);

		ucoord commdir_dst = bgq_direction2commdir(d_dst);
		bgq_weylfield_section sec_send = bgq_direction2section(d_dst, true);
		assert(secsize == bgq_weyl_section_offset(sec_send+1) - bgq_weyl_section_offset(sec_send));

		messageSizes[commdir_src] = secsize;
		soffsets[commdir_src] = totalMessageSize;
		totalMessageSize += secsize;

		// forward here is backward on the right neighbor
		// and the other way around...
		if(commdir_src%2 == 0)
		  roffsets[commdir_src] = soffsets[commdir_src] + messageSizes[commdir_src];
		else
		  roffsets[commdir_src] = soffsets[commdir_src] - messageSizes[commdir_src-1];
	}




	  Personality_t pers;
	  int rc = 0;
	  // get the CNK personality
	  Kernel_GetPersonality(&pers, sizeof(pers));
	  int mypers[6];
	  mypers[0] = pers.Network_Config.Acoord;
	  mypers[1] = pers.Network_Config.Bcoord;
	  mypers[2] = pers.Network_Config.Ccoord;
	  mypers[3] = pers.Network_Config.Dcoord;
	  mypers[4] = pers.Network_Config.Ecoord;

	  get_destinations(mypers);

	  // adjust the SPI pointers to the send and receive buffers
	  SPIrecvBuffers = (char*)recvBuffer;
	  SPIsendBuffers = (char*)sendBuffer;

	  // Setup the FIFO handles
	  rc = msg_InjFifoInit (&injFifoHandle,
				 0,                      /* startingSubgroupId */
				 0,                      /* startingFifoId     */
				 spi_num_dirs,           /* numFifos   */
				 INJ_MEMORY_FIFO_SIZE+1, /* fifoSize */
				 NULL                    /* Use default attributes */
				 );
	  if(rc != 0) {
	    fprintf(stderr, "msg_InjFifoInit failed with rc=%d\n",rc);
	    exit(1);
	  }

	  // Set up base address table for reception counter and buffer
	  setup_mregions_bats_counters(totalMessageSize);

	  // Create descriptors
	  // Injection Direct Put Descriptor, one for each neighbour
	  SPIDescriptors = (MUHWI_Descriptor_t*)(((uint64_t)SPIDescriptorsMemory+64)&~(64-1));
	  create_descriptors(SPIDescriptors, messageSizes, soffsets, roffsets, spi_num_dirs);

	  // test communication
	  for(unsigned int i = 0; i < RAND/2; i++) {
	    sendBuffer[i].s0.c0 = (double)g_cart_id;
	    sendBuffer[i].s0.c1 = (double)g_cart_id;
	    sendBuffer[i].s0.c2 = (double)g_cart_id;
	    sendBuffer[i].s1.c0 = (double)g_cart_id;
	    sendBuffer[i].s1.c1 = (double)g_cart_id;
	    sendBuffer[i].s1.c2 = (double)g_cart_id;
	  }

	  // Initialize the barrier, resetting the hardware.
	  rc = MUSPI_GIBarrierInit ( &GIBarrier, 0 /*comm world class route */);
	  if(rc) {
	    printf("MUSPI_GIBarrierInit returned rc = %d\n", rc);
	    exit(__LINE__);
	  }
	  // reset the recv counter
	  recvCounter = totalMessageSize;
	  global_barrier(); // make sure everybody is set recv counter

	  //#pragma omp for nowait
	  for (unsigned int j = 0; j < spi_num_dirs; j++) {
	    descCount[ j ] =
	      msg_InjFifoInject ( injFifoHandle,
				  j,
				  &SPIDescriptors[j]);
	  }
	  // wait for receive completion
	  while ( recvCounter > 0 );

	  _bgq_msync();

	  j = 0;
	  for(unsigned int i = 0; i < spi_num_dirs; i++) {
	    if(i == 0) k = g_nb_t_up;
	    if(i == 1) k = g_nb_t_dn;
	    if(i == 2) k = g_nb_x_up;
	    if(i == 3) k = g_nb_x_dn;
	    if(i == 4) k = g_nb_y_up;
	    if(i == 5) k = g_nb_y_dn;
	    if(i == 6) k = g_nb_z_up;
	    if(i == 7) k = g_nb_z_dn;
	    for(int mu = 0; mu < messageSizes[i]/sizeof(halfspinor); mu++) {
	      if(k != (int)creal(recvBuffer[ soffsets[i]/sizeof(halfspinor) + mu ].s0.c0) ||
		 k != (int)creal(recvBuffer[ soffsets[i]/sizeof(halfspinor) + mu ].s0.c1) ||
		 k != (int)creal(recvBuffer[ soffsets[i]/sizeof(halfspinor) + mu ].s0.c2) ||
		 k != (int)creal(recvBuffer[ soffsets[i]/sizeof(halfspinor) + mu ].s1.c0) ||
		 k != (int)creal(recvBuffer[ soffsets[i]/sizeof(halfspinor) + mu ].s1.c1) ||
		 k != (int)creal(recvBuffer[ soffsets[i]/sizeof(halfspinor) + mu ].s1.c2)) {
		if(g_cart_id == 0) {
		  printf("SPI exchange doesn't work for dir %d: %d != %d at point %d\n",
			 i, k ,(int)creal(recvBuffer[ soffsets[i]/sizeof(halfspinor) + mu ].s0.c0), mu);
		}
		j++;
	      }
	    }
	  }
	  if(j > 0) {
	    printf("hmm, SPI exchange failed on proc %d...\n!", g_cart_id);
	  }
	  else {
	    if(g_cart_id == 0) printf("# SPI exchange successfully tested\n");
	  }
#endif
}


//TODO: inline?
void bgq_comm_recv() {
//	return;
	MPI_CHECK(MPI_Startall(COMMDIR_COUNT, g_bgq_request_recv));
}


void bgq_comm_send() {
//	return;
	MPI_CHECK(MPI_Startall(COMMDIR_COUNT, g_bgq_request_send));
}


void bgq_comm_wait() {
//	return;
	MPI_Status recv_status[COMMDIR_COUNT];
	MPI_CHECK(MPI_Waitall(COMMDIR_COUNT, g_bgq_request_recv, recv_status));

	MPI_Status send_status[COMMDIR_COUNT];
	MPI_CHECK(MPI_Waitall(COMMDIR_COUNT, g_bgq_request_send, send_status));


#ifndef NDEBUG
	for (ucoord commdir = 0; commdir < COMMDIR_COUNT; commdir += 1) {
		bgq_direction d = bgq_commdir2direction(commdir);
		bgq_weylfield_section sec = bgq_direction2section(d,false);
		size_t size = bgq_weyl_section_offset(sec+1) - bgq_weyl_section_offset(sec);
		assert(get_MPI_count(&recv_status[commdir]) == size);
	}
#endif
}

