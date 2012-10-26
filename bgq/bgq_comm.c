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




void bgq_comm_init() {
#ifdef SPI
	  // here comes the SPI initialization
	  uint64_t messageSizes[PHYSICAL_LD];
	  uint64_t roffsets[PHYSICAL_LD];
	  uint64_t soffsets[PHYSICAL_LD];

	 spi_num_dirs = 2*(COMM_T+COMM_X+COMM_Y+COMM_Z);
	 size_t i = 0;
	 if (COMM_T) {
		 messageSizes[COMM_ORD_TSEND] = LOCAL_HALO_T * sizeof(bgq_weyl_vec) / PHYSICAL_LP;
		 messageSizes[COMM_ORD_TRECV] = LOCAL_HALO_T * sizeof(bgq_weyl_vec) / PHYSICAL_LP;
	 }
	 if (COMM_X) {
		 messageSizes[COMM_ORD_XSEND] = PHYSICAL_HALO_X * sizeof(bgq_weyl_vec);
		 messageSizes[COMM_ORD_XRECV] = PHYSICAL_HALO_X * sizeof(bgq_weyl_vec);
	 }
	 if (COMM_Y) {
		 messageSizes[COMM_ORD_YSEND] = PHYSICAL_HALO_Y * sizeof(bgq_weyl_vec);
		 messageSizes[COMM_ORD_YRECV] = PHYSICAL_HALO_Y * sizeof(bgq_weyl_vec);
	 	 }
	 if (COMM_Z) {
		 messageSizes[COMM_ORD_ZSEND] = PHYSICAL_HALO_Z * sizeof(bgq_weyl_vec);
		 messageSizes[COMM_ORD_ZRECV] = PHYSICAL_HALO_Z * sizeof(bgq_weyl_vec);
	 	 }

	 totalMessageSize = 0;
	 for(size_t i = 0; i < spi_num_dirs; i++) {
		 soffsets[i] = totalMessageSize;
		 totalMessageSize += messageSizes[i];
	 }

	  for(size_t i = 0; i < spi_num_dirs; i++) {
	    // forward here is backward on the right neighbour
	    // and the other way around...
	    if(i%2 == 0) {
	      roffsets[i] = soffsets[i] + messageSizes[i];
	    }
	    else {
	      roffsets[i] = soffsets[i] - messageSizes[i-1];
	    }
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
	  SPIrecvBuffers = (char *)(recvBuffer);
	  SPIsendBuffers = (char *)(sendBuffer);

	  // Setup the FIFO handles
	  rc = msg_InjFifoInit ( &injFifoHandle,
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
	  SPIDescriptors =
	    ( MUHWI_Descriptor_t *)(((uint64_t)SPIDescriptorsMemory+64)&~(64-1));
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
