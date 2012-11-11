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
#include "../DirectPut.h"
#include <upci/upc_atomic.h>
#endif
#include "../global.h"

#include <mpi.h>




static MPI_Request g_bgq_request_recv[PHYSICAL_LD];
static MPI_Request g_bgq_request_send[PHYSICAL_LD];




static bgq_dimension bgq_commdim2dimension(ucoord commdim) {
	assert(0 <= commdim && commdim < PHYSICAL_LD);
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
	return -1;
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
	return -1;
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
	return -1;
}


#if 0
/////////////////////////////////////////
/// Test library common code
////////////////////////////////////////
/* begin_generated_IBM_copyright_prolog                             */
/*                                                                  */
/* This is an automatically generated copyright prolog.             */
/* After initializing,  DO NOT MODIFY OR MOVE                       */
/* ================================================================ */
/*                                                                  */
/* Licensed Materials - Property of IBM                             */
/*                                                                  */
/* Blue Gene/Q                                                      */
/*                                                                  */
/* (C) Copyright IBM Corp.  2008, 2012                              */
/*                                                                  */
/* US Government Users Restricted Rights -                          */
/* Use, duplication or disclosure restricted                        */
/* by GSA ADP Schedule Contract with IBM Corp.                      */
/*                                                                  */
/* This software is available to you under the                      */
/* Eclipse Public License (EPL).                                    */
/*                                                                  */
/* ================================================================ */
/*                                                                  */
/* end_generated_IBM_copyright_prolog                               */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <stdint.h>
#include <malloc.h>
#include <assert.h>
#include <string.h>


//////////////////////////////////////////
// Basic SPI and HWI includes
//////////////////////////////////////////
#include <hwi/include/bqc/A2_core.h>
#include <hwi/include/bqc/A2_inlines.h>
#include <hwi/include/bqc/MU_PacketCommon.h>
#include <firmware/include/personality.h>
#include <spi/include/mu/Descriptor.h>
#include <spi/include/mu/Descriptor_inlines.h>
#include <spi/include/mu/InjFifo.h>
#include <spi/include/mu/RecFifo.h>
#include <spi/include/mu/Addressing.h>
#include <spi/include/mu/Addressing_inlines.h>
#include <spi/include/mu/GIBarrier.h>
#include <spi/include/kernel/MU.h>
#include <spi/include/kernel/process.h>
#include <spi/include/kernel/location.h>
#include <spi/include/kernel/collective.h>

/**
 * \brief Initialize Buffer with Random Data
 *
 * Initialize a buffer with random data, using the APPSEED env var random number seed,
 * the message's source coordinates, the message's destination coordinates, and a
 * specified unique value.
 *
 * The buffer is initialized backwards so it can be checked backwards.
 *
 * \param [in]  bufPtr     Pointer to the buffer to be initialized.
 * \param [in]  size       Size of the buffer, in bytes.
 * \param [in]  sourcePtr  Pointer to source coordinates (may be NULL)
 * \param [in]  destPtr    Pointer to destination coordinates (may be NULL)
 * \param [in]  uniquePtr  Pointer to unique value (may be NULL)
 * \param [in]  chksumPtr  Pointer to a location where the chksum of the buffer
 *                         is stored (may be NULL).  This result is only
 *                         accurate if size is a multiple of 8.
 *
 * \returns  The buffer is initialized with random data.
 *           The chksum result is returned in the location pointed-to by
 *           chksumPtr.
 */
static void msg_InitBuffer ( void                *bufPtr,
                      size_t               size,
                      MUHWI_Destination_t *sourcePtr,
                      MUHWI_Destination_t *destPtr,
                      uint32_t            *uniquePtr,
                      uint64_t            *chksumPtr );

/**
 * \brief Check Buffer with Random Data
 *
 * Check a buffer with random data, using the APPSEED env var random number seed,
 * the message's source coordinates, the message's destination coordinates,
 * and a specified unique value.
 *
 * The buffer is checked backwards.  This is because the last bytes of the
 * buffer are the last to arrive, and we want to ensure those bytes are there
 * at the time we begin checking the buffer.
 *
 * \param [in]  bufPtr     Pointer to the buffer to be checked.
 * \param [in]  size       Size of the buffer, in bytes.
 * \param [in]  sourcePtr  Pointer to source coordinates (may be NULL)
 * \param [in]  destPtr    Pointer to destination coordinates (may be NULL)
 * \param [in]  uniquePtr  Pointer to unique value (may be NULL)
 * \param [in]  textPtr    Pointer to text shown in error messages
 *
 * \retval  0  The buffer checks-out successfully.
 * \retval -1  The buffer failed to check correctly.  A message was printed
 *             showing where the check failed.
 *
 * \see msg_InitBuffer()
 */
 static int msg_CheckBuffer ( void                *bufPtr,
                      size_t               size,
                      MUHWI_Destination_t *sourcePtr,
                      MUHWI_Destination_t *destPtr,
                      uint32_t            *uniquePtr,
                      char                *textPtr );

/**
 * \brief Poll a Reception Counter Until It Hits Zero
 *
 * This function polls the specified reception counter until it hits
 * zero, with no timeout.
 *
 * \param [in]  receptionCounter  Pointer to the reception counter
 *
 */
 static void msg_CounterPoll( volatile uint64_t *receptionCounter );


/**
 * \brief Injection Fifo Handle
 *
 * This is a "handle" returned from msg_InjFifoInit() and passed into subsequent
 * calls to msg_InjFifoXXXX() functions.  It is used internally within the
 * msg_InjFifoXXXX() functions to anchor resources that have been allocated.
 */
//typedef struct {
//  void* pOpaqueObject;
//} msg_InjFifoHandle_t;


/**
 * \brief Poll a Reception Counter Until It Hits Zero or Exceeds Timeout
 *
 * This function polls the specified reception counter until it hits
 * zero, or until the specified time limit is exceeded.
 *
 * \param [in]  receptionCounter  Pointer to the reception counter
 * \param [in]  counterPollBase   The base value returned from a previous
 *                                call to msg_CounterPollWithTimeoutInit().
 * \param [in]  seconds           Time limit (in seconds)
 *
 * \retval  0  Successful poll.  The counter hit zero.
 * \retval -1  Unsuccessful poll.  The time limit was exceeded.
 *
 * \see msg_CounterPollWithTimeoutInit()
 */
static int msg_CounterPollWithTimeout( volatile uint64_t *receptionCounter,
                                uint64_t           counterPollBase,
                                uint64_t           seconds );

static int msg_InjFifoInit2 ( msg_InjFifoHandle_t *injFifoHandlePtr,
                      uint32_t             startingSubgroupId,
                      uint32_t             startingFifoId,
                      uint32_t             numFifos,
                      size_t               fifoSize,
                      Kernel_InjFifoAttributes_t  *injFifoAttrs );


/**
 * \brief Terminate Injection Fifos
 *
 * Terminate the usage of injection fifos.  This deactivates the fifos and
 * frees all of the storage associated with them (previously allocated during
 * msg_InjFifoInit()).
 *
 * \param [in]  injFifoHandle  The handle returned from msg_InjFifoInit().
 *                             It must be passed into this function untouched
 *                             from when it was returned from msg_InjFifoInit().
 *
 * \note After this function returns, no more InjFifo functions should be called
 *       with this injFifoHandle.
 */
static void msg_InjFifoTerm ( msg_InjFifoHandle_t injFifoHandle );


/**
 * \brief Inject Descriptor into Injection Fifo
 *
 * Inject the specified descriptor into the specified injection fifo.
 *
 * \param [in]  injFifoHandle  The handle returned from msg_InjFifoInit().
 *                             It must be passed into this function untouched
 *                             from when it was returned from msg_InjFifoInit().
 * \param [in]  relativeFifoId  The fifo number, relative to the start of
 *                              the fifos managed by this opaque object.
 *                              For example, if msg_InjFifoInit() was called
 *                              to init fifos in subgroup 2, starting with
 *                              fifo Id 3, the relativeFifoNumber of the
 *                              first fifo is 0, not 3.
 * \param [in]  descPtr         Pointer to the descriptor to be injected.
 *
 * \retval  positiveNumber  The descriptor was successfully injected.  The
 *                          returned value is the sequence number of this
 *                          descriptor.
 * \retval  -1              The descriptor was not injected, most likely because
 *                          there is no room in the fifo.
 */
static uint64_t msg_InjFifoInject2 ( msg_InjFifoHandle_t injFifoHandle,
                             uint32_t            relativeFifoId,
                             MUHWI_Descriptor_t *descPtr );


static inline void msg_BuildPt2PtDirectPutInfo
(MUSPI_Pt2PtDirectPutDescriptorInfo_t   * dinfo,
 uint64_t                                 put_offset,
 uint64_t                                 counter_offset,
 uint64_t                                 buffer,
 uint64_t                                 bytes,
 uint64_t                                 torusInjectionFifoMap,
 int                                      vc,
 uint8_t                                  hintsABCD,
 uint8_t                                  hintsE)
{
  dinfo->Base.Pre_Fetch_Only  = 0;
  dinfo->Base.Payload_Address = buffer;
  dinfo->Base.Message_Length  = bytes;
  dinfo->Base.Torus_FIFO_Map  = torusInjectionFifoMap;
  //  dinfo->Base.Dest            = dest;
  dinfo->Pt2Pt.Hints_ABCD = hintsABCD;

  if ( vc == MUHWI_PACKET_VIRTUAL_CHANNEL_DETERMINISTIC ){
    dinfo->Pt2Pt.Misc1      = hintsE                                 |
                              MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE   |
                              MUHWI_PACKET_USE_DETERMINISTIC_ROUTING |
                              MUHWI_PACKET_DO_NOT_DEPOSIT;
    dinfo->Pt2Pt.Misc2      = MUHWI_PACKET_VIRTUAL_CHANNEL_DETERMINISTIC;
  }
  else if (vc == MUHWI_PACKET_VIRTUAL_CHANNEL_DYNAMIC) {
    dinfo->Pt2Pt.Misc1      = hintsE                               |
                              MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE |
                              MUHWI_PACKET_USE_DYNAMIC_ROUTING     |
                              MUHWI_PACKET_DO_NOT_DEPOSIT;
    dinfo->Pt2Pt.Misc2      = MUHWI_PACKET_VIRTUAL_CHANNEL_DYNAMIC;
  }
  dinfo->Pt2Pt.Skip       = 0;
  dinfo->DirectPut.Rec_Payload_Base_Address_Id = 0;
  dinfo->DirectPut.Rec_Payload_Offset          = put_offset;
  dinfo->DirectPut.Rec_Counter_Base_Address_Id = 0;
  dinfo->DirectPut.Rec_Counter_Offset          = counter_offset;

  dinfo->DirectPut.Pacing = MUHWI_PACKET_DIRECT_PUT_IS_NOT_PACED;
}


/* Init software bytes for the packet header. */
typedef union _swb
{
  struct {
    uint8_t reserved[6];
    uint8_t bytes[14];
  };

  struct {
    uint8_t  reserved[6];
    uint8_t  functionid;
    uint8_t  bytes[5];

    uint32_t message_size_in_bytes;
    uint8_t unused[4];
  } BytesStruct;

} SoftwareBytes_t;


static inline void msg_BuildPt2PtMemoryFifoInfo
(MUSPI_Pt2PtMemoryFIFODescriptorInfo_t  * minfo,
 SoftwareBytes_t                          SoftwareBytes,
 uint                                     rfifoid,
 uint64_t                                 put_offset,
 uint64_t                                 buffer,
 uint64_t                                 bytes,
 uint64_t                                 torusInjectionFifoMap,
 int                                      vc,
 uint8_t                                  hintsABCD,
 uint8_t                                  hintsE)
{
  minfo->Pt2Pt.Hints_ABCD = hintsABCD;
  if ( vc ==  MUHWI_PACKET_VIRTUAL_CHANNEL_DETERMINISTIC ) {
    minfo->Pt2Pt.Misc1      = hintsE         |
      MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE   |
      MUHWI_PACKET_USE_DETERMINISTIC_ROUTING |
      MUHWI_PACKET_DO_NOT_DEPOSIT;
    minfo->Pt2Pt.Misc2 = MUHWI_PACKET_VIRTUAL_CHANNEL_DETERMINISTIC;
  }
  else if (vc == MUHWI_PACKET_VIRTUAL_CHANNEL_DYNAMIC) {
    minfo->Pt2Pt.Misc1      = hintsE       |
      MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE |
      MUHWI_PACKET_USE_DYNAMIC_ROUTING     |
      MUHWI_PACKET_DO_NOT_DEPOSIT;
    minfo->Pt2Pt.Misc2      = MUHWI_PACKET_VIRTUAL_CHANNEL_DYNAMIC;
  }
  else if ( vc ==  MUHWI_PACKET_VIRTUAL_CHANNEL_SYSTEM ) {
    minfo->Pt2Pt.Misc1      = hintsE         |
      MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE   |
      MUHWI_PACKET_USE_DETERMINISTIC_ROUTING |
      MUHWI_PACKET_DO_NOT_DEPOSIT;
    minfo->Pt2Pt.Misc2 = MUHWI_PACKET_VIRTUAL_CHANNEL_SYSTEM;
  }
  else { MUSPI_assert(0); }

  minfo->Pt2Pt.Skip  = 0;
  minfo->MemFIFO.Rec_FIFO_Id    = rfifoid;
  minfo->MemFIFO.Rec_Put_Offset = put_offset;
  minfo->Base.Pre_Fetch_Only= MUHWI_DESCRIPTOR_PRE_FETCH_ONLY_NO;
  minfo->MemFIFO.Interrupt  = MUHWI_PACKET_DO_NOT_INTERRUPT_ON_PACKET_ARRIVAL;
  minfo->MemFIFO.SoftwareBit = 0;
  memcpy( minfo->MemFIFO.SoftwareBytes,
          SoftwareBytes.bytes,
          sizeof( minfo->MemFIFO.SoftwareBytes ) );

  //Warning assume V == P
  minfo->Base.Payload_Address = buffer;
  minfo->Base.Message_Length =  bytes;
  minfo->Base.Torus_FIFO_Map = torusInjectionFifoMap;
}



static inline void msg_BuildCollectiveDirectPutAllreduceInfo
(MUSPI_CollectiveDirectPutDescriptorInfo_t  * dinfo,
 uint64_t                                bytes,
 uint64_t                                buffer,
 uint32_t                                op,
 uint32_t                                sizeoftype,
 MUHWI_Destination_t                     dest,
 uint64_t                                put_offset,
 uint64_t                                counter_offset,
 uint32_t                                classRoute,
 int                                     vc)
{
  uint64_t torusInjectionFifoMap = MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CUSER;  // bit 0 ===> A- Torus FIFO

  if (vc == MUHWI_PACKET_VIRTUAL_CHANNEL_SYSTEM_COLLECTIVE)
    torusInjectionFifoMap = MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CSYSTEM;

  dinfo->Base.Pre_Fetch_Only  = 0;
  dinfo->Base.Payload_Address = (uint64_t)buffer;
  dinfo->Base.Message_Length  = bytes;
  dinfo->Base.Torus_FIFO_Map  = torusInjectionFifoMap;
  dinfo->Base.Dest            = dest;

  dinfo->Collective.Op_Code     = op;
  dinfo->Collective.Word_Length = sizeoftype;
  dinfo->Collective.Class_Route = classRoute;
  dinfo->Collective.Misc        = vc | MUHWI_COLLECTIVE_TYPE_ALLREDUCE;
  dinfo->Collective.Skip        = 0;

  dinfo->DirectPut.Rec_Payload_Base_Address_Id = 0;
  dinfo->DirectPut.Rec_Payload_Offset          = put_offset;
  dinfo->DirectPut.Rec_Counter_Base_Address_Id = 0;
  dinfo->DirectPut.Rec_Counter_Offset          = counter_offset;

  dinfo->DirectPut.Pacing = MUHWI_PACKET_DIRECT_PUT_IS_NOT_PACED;
}


static inline void msg_BuildCollectiveMemoryFifoAllreduceInfo
(MUSPI_CollectiveMemoryFIFODescriptorInfo_t  * minfo,
 SoftwareBytes_t                          SoftwareBytes,
 uint64_t                                bytes,
 uint64_t                                buffer,
 uint32_t                                op,
 uint32_t                                sizeoftype,
 MUHWI_Destination_t                     dest,
 uint64_t                                put_offset,
 uint32_t                                rfifoid,
 uint32_t                                classRoute,
 int                                     vc)
{
  uint64_t torusInjectionFifoMap = MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CUSER;  // bit 0 ===> A- Torus FIFO

  if (vc == MUHWI_PACKET_VIRTUAL_CHANNEL_SYSTEM_COLLECTIVE)
    torusInjectionFifoMap = MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CSYSTEM;

  minfo->Base.Pre_Fetch_Only  = 0;
  minfo->Base.Payload_Address = (uint64_t)buffer;
  minfo->Base.Message_Length  = bytes;
  minfo->Base.Torus_FIFO_Map  = torusInjectionFifoMap;
  minfo->Base.Dest            = dest;

  minfo->Collective.Op_Code     = op;
  minfo->Collective.Word_Length = sizeoftype;
  minfo->Collective.Class_Route = classRoute;
  minfo->Collective.Misc        = vc | MUHWI_COLLECTIVE_TYPE_ALLREDUCE;
  minfo->Collective.Skip        = 0;

  minfo->MemFIFO.Rec_FIFO_Id    = rfifoid;
  minfo->MemFIFO.Rec_Put_Offset = put_offset;

  minfo->MemFIFO.Interrupt  = MUHWI_PACKET_DO_NOT_INTERRUPT_ON_PACKET_ARRIVAL;
  minfo->MemFIFO.SoftwareBit = 0;
  memcpy( minfo->MemFIFO.SoftwareBytes,
          SoftwareBytes.bytes,
          sizeof( minfo->MemFIFO.SoftwareBytes ) );
}


#define SEND_BUFFER_ALIGNMENT   128
#define RECV_BUFFER_ALIGNMENT   128
#define MAX_MESSAGE_SIZE       256

// Sub message size
static int window_size  =  256;   //  size of submessages
static long long messageSizeInBytes = MAX_MESSAGE_SIZE;

// Maximum number of ndoes this test will run on
#define MAX_NUM_NODES          32

//#define NUM_INJ_FIFOS          64
//#define INJ_MEMORY_FIFO_SIZE  ((64*MAX_NUM_NODES) -1)

static int num_iterations  = 100;
static int skip_iterations = 1;
static int check_messages  = 1;
//static int do_dynamic      = 0;

// Allocate static memory for descriptors
static char muDescriptorsMemory[ MAX_NUM_NODES * sizeof(MUHWI_Descriptor_t) + 64 ];

// pointer to descriptor array
MUHWI_Descriptor_t *muDescriptors;

// Allocate static memory for send and receive buffers
static char sendBufMemory[MAX_NUM_NODES * MAX_MESSAGE_SIZE+ SEND_BUFFER_ALIGNMENT];
static char recvBufMemory[MAX_NUM_NODES * MAX_MESSAGE_SIZE+ SEND_BUFFER_ALIGNMENT];


// Random permutation array. The alltoall messages will be sent
// to neighbors in a random permutation
static uint32_t randPerm[PHYSICAL_LD];

// rank2dest cache
#if 0
struct {
  MUHWI_Destination_t dest;
  uint8_t             hintsABCD;
  uint8_t             hintsE;
} rank2dest[MAX_NUM_NODES];
#endif


// pointers to send and receive buffers
static char * recvBuffers;
static char * sendBuffers;

// receive counter
static volatile uint64_t recvCounter;

// base addess table slot for receive buffer and counter
static uint32_t recvBufBatId = 0, recvCntrBatId = 1;

// physical address of send buffers
static uint64_t sendBufPAddr;

#if 0
struct {
  MUHWI_Destination_t dest;
  uint8_t             hintsABCD;
  uint8_t             hintsE;
  torus_t torus;
} nb2dest[PHYSICAL_LD];
#endif


// number of participating processes (1 proc/node)
static unsigned numNodes;

// rank of this process
static unsigned myRank;// obsolete

// Enable different zone routing modes
static uint8_t  zoneRoutingMask = 0;
static unsigned zoneRoutingId   = 0;

// stay on bubble bits
static uint8_t stayOnBubbleMask  = 0;
static unsigned stayOnBubbleFlag = 0;

// Call to create the alltoall descriptors
static void create_descriptors2(void);

// Set up hints, fifomaps and destination address values
static void setup_destinations(Personality_t *p);

// Call to set up the base address table id and memory regions for the
// alltoall
static void setup_mregions_bats_counters2(void);

// Create a random permutation
static void create_random_permutation(void);

MUSPI_GIBarrier_t GIBarrier;

static void alltoall_exit (int rc) {
  exit (rc);
}

//Call GI barrier
static void global_barrier2(void)
{
  int rc = 0;
  uint64_t timeoutCycles = 60UL * 1600000000UL; // about 60 sec at 1.6 ghz
  rc = MUSPI_GIBarrierEnter ( &GIBarrier );
  if (rc)
    {
      printf("MUSPI_GIBarrierEnter failed returned rc = %d\n", rc);
      alltoall_exit(1);
    }

  //=======================================================================

  // Poll for completion of the barrier.
  rc = MUSPI_GIBarrierPollWithTimeout ( &GIBarrier, timeoutCycles);
  if ( rc )
    {
      printf("MUSPI_GIBarrierPollWithTimeout failed returned rc = %d\n", rc);
      DelayTimeBase (200000000000UL);
      alltoall_exit(1);
    }
}



/////////////////////////////////////////////////////////////////////////////
/// Read environment variables to program alltoall
/////////////////////////////////////////////////////////////////////////////
static void getEnvVars()
{
  char *envvar;

  envvar = getenv("MSG_ALLTOALL_MESSAGE_SIZE");
  if (envvar)
    {
      messageSizeInBytes = strtoul( envvar, 0, 10 );

      if (  messageSizeInBytes > MAX_MESSAGE_SIZE )
	messageSizeInBytes = MAX_MESSAGE_SIZE;

      if ( (messageSizeInBytes % SEND_BUFFER_ALIGNMENT) ||
	   (messageSizeInBytes % RECV_BUFFER_ALIGNMENT) )
	{
	  printf("ERROR: MESSAGE_SIZE %lld must be multiple of both SEND %d and RECV %d ALIGNMENT\n",
		 messageSizeInBytes,
		 SEND_BUFFER_ALIGNMENT,
		 RECV_BUFFER_ALIGNMENT);

	  alltoall_exit(1);
	}
    }

  envvar = getenv("MSG_ALLTOALL_CHECK_MESSAGES");
  if (envvar)
    {
      check_messages = strtoul( envvar, 0, 10 );
    }

  envvar = getenv("MSG_ALLTOALL_ZONE_ROUTING_ID");
  if (envvar)
    {
      zoneRoutingId = strtoul( envvar, 0, 10 );
      switch (zoneRoutingId)
	{
	case 0: zoneRoutingMask = MUHWI_PACKET_ZONE_ROUTING_0; break;
	case 1: zoneRoutingMask = MUHWI_PACKET_ZONE_ROUTING_1; break;
	case 2: zoneRoutingMask = MUHWI_PACKET_ZONE_ROUTING_2; break;
	case 3: zoneRoutingMask = MUHWI_PACKET_ZONE_ROUTING_3; break;
	default:
	  printf("ERROR: MSG_ALLTOALL_ZONE_ROUTING_ID is not in range [0,3]\n");
	}
    }

  envvar = getenv("MSG_ALLTOALL_STAY_ON_BUBBLE");
  if (envvar)
    {
      stayOnBubbleFlag = strtoul( envvar, 0, 10 );
      if ( stayOnBubbleFlag )
	stayOnBubbleMask = MUHWI_PACKET_STAY_ON_BUBBLE;
    }


  envvar = getenv("MSG_ALLTOALL_DYNAMIC");
  if (envvar)
    {
      do_dynamic = strtoul( envvar, 0, 10 );
    }
}

static unsigned msg_InjFifoCheckCompletion2(msg_InjFifoHandle_t injFifoHandle,
                                    uint32_t            relativeFifoId,
                                    uint64_t            desc_count);

static int main2(int argc, char **argv)
{
    int rc;
    getEnvVars();

    Personality_t pers;
    Kernel_GetPersonality(&pers, sizeof(pers));

    //////////////////////////////////////////////////////
    ///////////  Set up destinations  ////////////////////
    //////////////////////////////////////////////////////

    setup_destinations(&pers);

    if (myRank == 0)
    {
      printf("main(): Alltoall Performance Test\n");
    }

    ///////////////////////////////////////////////////////
    /////////// Set up Injection Fifos  ///////////////////
    ///////////////////////////////////////////////////////
    msg_InjFifoHandle_t injFifoHandle;
    rc = msg_InjFifoInit ( &injFifoHandle,
			   0,        /* startingSubgroupId */
			   0,        /* startingFifoId     */
			   NUM_INJ_FIFOS,       /* numFifos   */
			   INJ_MEMORY_FIFO_SIZE+1, /* fifoSize */
			   NULL      /* Use default attributes */
			 );
    if (rc != 0)
    {
      printf("msg_InjFifoInit failed with rc=%d\n",rc);
      alltoall_exit(1);
    }

    /////////////////////////////////////////////////////////
    ////////  Set up send and recive buffers   //////////////
    /////////////////////////////////////////////////////////
    //recvBuffers = (char *)(((uint64_t)recvBufMemory+RECV_BUFFER_ALIGNMENT)&~(RECV_BUFFER_ALIGNMENT-1));
    recvBuffers = (char*)(g_bgq_sec_recv[0]);
    assert((uintptr_t)recvBuffers % RECV_BUFFER_ALIGNMENT == 0);
    //sendBuffers = (char *)(((uint64_t)sendBufMemory+SEND_BUFFER_ALIGNMENT)&~(SEND_BUFFER_ALIGNMENT-1));
    sendBuffers = (char*)(g_bgq_sec_send[0]);
    assert((uintptr_t)sendBuffers % SEND_BUFFER_ALIGNMENT == 0);

    if (check_messages)
     {
       uint64_t offset;
       int i;
       for (i=0, offset = 0;
	    i < COMMDIR_COUNT;
	    i++, offset += messageSizeInBytes)
	 {
	   if ( i == myRank ) continue; // no self send
	   msg_InitBuffer ( sendBuffers + offset,
			    messageSizeInBytes,
			    (MUHWI_Destination_t *)&myRank, /* sourcePtr */
			    (MUHWI_Destination_t *)&i,      /* destPtr */
			    NULL ,                          /* uniquePtr */
			    NULL                            /* chksumPtr */ );
	 }
     }

    ////////////////////////////////////////////////////////////////////////
    // Set up base address table for reception counter and buffer
    ////////////////////////////////////////////////////////////////////////
    setup_mregions_bats_counters2();

    ////////////////////////////////////////////////////////////////////////
    // Create descriptors
    ////////////////////////////////////////////////////////////////////////

    // Injection Direct Put Descriptor, one for each neighbor
    muDescriptors =
      ( MUHWI_Descriptor_t *)(((uint64_t)muDescriptorsMemory+64)&~(64-1));

    create_descriptors2();

    // Initialize the barrier, resetting the hardware.
    rc = MUSPI_GIBarrierInit ( &GIBarrier, 0 /*comm world class route */);


    if (rc)
      {
	printf("MUSPI_GIBarrierInit returned rc = %d\n", rc);
	alltoall_exit(__LINE__);
      }

    ///////////////////////////////////////////////////////////////////////////
    //
    // Start the timer.
    // Loop num_iterations times.
    // - Set the reception counter to the message length * numNodes
    // - Loop through each destination, injecting the descriptor.
    // - Wait for receive and send completion
    // Stop the iteration timer.
    //
    //////////////////////////////////////////////////////////////////////////

    uint64_t totalCycles=0;
    uint64_t startTime=0;
    uint64_t totalBytesPerIter     =  messageSizeInBytes * COMMDIR_COUNT;
    master_print("totalBytesPerIter=%llu\n",totalBytesPerIter);
    uint64_t totalBytes            =  totalBytesPerIter * num_iterations;
    uint64_t totalNumberOfMessages =  (COMMDIR_COUNT - 1) * num_iterations;
    uint64_t bytesPer1000Cycles;

    assert (window_size <= messageSizeInBytes);
    if ( myRank == 0 )
      {
	printf("processes:%d, #iterations = %d, Message Size = %llu, Send Buffer Alignment = %u, sendBuffer@ = %p, Recv Buffer Alignment = %u, recvBuffer@=%p,zoneId:%d,stayOnBubble:%d,submessage%d\nPerforming communications...\n",
	       numNodes, num_iterations, messageSizeInBytes, SEND_BUFFER_ALIGNMENT,sendBuffers,
	       RECV_BUFFER_ALIGNMENT,recvBuffers,zoneRoutingId,stayOnBubbleFlag, window_size);
      }

    // create random order for sends
    create_random_permutation();

    uint64_t descCount[NUM_INJ_FIFOS];

    int i, j;
    for ( i=0; i<num_iterations+skip_iterations; i++)
      {
	// reset the recv counter
	recvCounter = totalBytesPerIter;

	if (i >= skip_iterations) startTime = GetTimeBase();

	 printf("Here is %d waiting before the barrier\n", g_proc_id);
	global_barrier(); // make sure everybody is set recv counter
	 printf("Here is %d behind the barrier\n", g_proc_id);

	uint64_t bytes = 0;
	for (bytes = 0; bytes < messageSizeInBytes; bytes += window_size) {
	  uint64_t msize = (bytes <= messageSizeInBytes - window_size) ? window_size : (messageSizeInBytes - bytes);
	  for ( j=0; j<COMMDIR_COUNT; j++)
	    {
	      //if ( myRank == randPerm[j] ) continue; // no self send

		  //master_print("injecting j=%d randPerm[j]=%d, msize=%d messageSizeInBytes=%d sendBufPAddr=%ull\n", j, randPerm[j],msize,messageSizeInBytes,(uintptr_t)sendBufPAddr);
		  assert(0<=randPerm[j]);
		  assert(randPerm[j]<COMMDIR_COUNT);
	      muDescriptors[randPerm[j]].Message_Length = msize;
	      muDescriptors[randPerm[j]].Pa_Payload    =  sendBufPAddr;
	      MUSPI_SetRecPutOffset(&muDescriptors[randPerm[j]], 0);//MUSPI_SetRecPutOffset(&muDescriptors[randPerm[j]], myRank * messageSizeInBytes + bytes);

	      descCount[ j % NUM_INJ_FIFOS ] =
		msg_InjFifoInject ( injFifoHandle,
				    j % NUM_INJ_FIFOS,
				    &muDescriptors[randPerm[j]]);
	      //master_print("end inject commdir=%d\n", j);
	    }
	}

	// wait for receive completion
	 printf("Here is %d waiting for receive\n", g_proc_id);
	while (true) {
		uint64_t left = recvCounter;
		if (left==0)
			break;
		printf("Here is %d received, %llu bytes to go...\n", g_proc_id, left);
	}
	 printf("Here is %d received\n", g_proc_id);

	// wait for send completion
	unsigned sendDone;
	unsigned nexp = ((COMMDIR_COUNT-1) < NUM_INJ_FIFOS) ? (COMMDIR_COUNT-1) : NUM_INJ_FIFOS;
	do
	  {
	    sendDone = 0;
	    for ( j = 0; j < NUM_INJ_FIFOS; j++ )
	      sendDone += msg_InjFifoCheckCompletion( injFifoHandle,
						      j,
						      descCount[j]);
	  }
	while ( sendDone < nexp );

	_bgq_msync(); // Ensure data is available to all cores.

	if (i >= skip_iterations) totalCycles += GetTimeBase() - startTime;
      }

    if (myRank == 0)
      printf("Cycles = %llu, Number of Messages = %llu, Cycles per Message = %llu, Cycles per alltoall = %llu\n",
	     (long long unsigned int) totalBytes,
	     (long long unsigned int) totalCycles,
	     (long long unsigned int) totalNumberOfMessages,
	     (long long unsigned int) totalCycles/num_iterations);

    master_print("Done All2All\n");

    if (check_messages)
      {
	uint64_t offset;
	for (i=0, offset = 0;
	     i < COMMDIR_COUNT;
	     i++, offset += messageSizeInBytes)
	  {
	    if ( i == myRank ) continue; // no self send

	    rc = msg_CheckBuffer ( recvBuffers + offset,
				   messageSizeInBytes,
				   (MUHWI_Destination_t *)&i, /* sourcePtr */
				   (MUHWI_Destination_t *)&myRank, /* destPtr */
				   NULL,   /* uniquePtr */
				   "Checking Receive Buffer" );
	    if ( rc )
	      alltoall_exit(1);
	  }

      }

    msg_InjFifoTerm ( injFifoHandle );

    return 0;
}



static void setup_mregions_bats_counters2(void)
{
  const uint64_t buffersSize =  numNodes * messageSizeInBytes;

  // allocate bat entries for the recive buffer and the receive counter

  uint32_t batIds[2] = { recvBufBatId, recvCntrBatId };
  MUSPI_BaseAddressTableSubGroup_t batSubGrp;

  int rc =  Kernel_AllocateBaseAddressTable( 0/*subgrpId*/,
					     &batSubGrp,
					     2,/*nbatids*/
					     batIds,
					     0 /* "User" use */);

  if ( rc != 0 )
    {
      printf("Kernel_AllocateBaseAddressTable failed with rc=%d\n",rc);
      alltoall_exit(1);
    }

  // Receive buffer bat is set to the PA addr of the receive buffer
  Kernel_MemoryRegion_t memRegion;
  rc = Kernel_CreateMemoryRegion ( &memRegion,
				   recvBuffers,
				   buffersSize);
  if ( rc != 0)
    {
      printf("Kernel_CreateMemoryRegion failed with rc=%d\n",rc);
      alltoall_exit(1);
    }

  uint64_t paAddr =
    (uint64_t)recvBuffers -
    (uint64_t)memRegion.BaseVa +
    (uint64_t)memRegion.BasePa;

  rc = MUSPI_SetBaseAddress ( &batSubGrp,
			      recvBufBatId,
			      paAddr );

  if (rc != 0)
    {
      printf("MUSPI_SetBaseAddress failed with rc=%d\n",rc);
      alltoall_exit(1);
    }

  // Receive counter bat is set to the MU style atomic PA addr of the receive counter
  if ( (uint64_t)(&recvCounter) & 0x7 )
    {
      printf("ERROR: recv counter is not 8 byte aligned\n");
      alltoall_exit(1);
    }

  rc = Kernel_CreateMemoryRegion ( &memRegion,
				   (void *)&recvCounter,
				   sizeof(recvCounter));
  if ( rc != 0)
    {
      printf("Kernel_CreateMemoryRegion failed with rc=%d\n",rc);
      alltoall_exit(1);
    }

  paAddr =
    (uint64_t)&recvCounter -
    (uint64_t)memRegion.BaseVa +
    (uint64_t)memRegion.BasePa;

  uint64_t paAddrAtomic =  MUSPI_GetAtomicAddress(paAddr,MUHWI_ATOMIC_OPCODE_STORE_ADD);

  rc = MUSPI_SetBaseAddress ( &batSubGrp,
			      recvCntrBatId,
			      paAddrAtomic );

  if (rc != 0)
    {
      printf("MUSPI_SetBaseAddress failed with rc=%d\n",rc);
      alltoall_exit(1);
    }

  // Get the send buffers physical address

  rc = Kernel_CreateMemoryRegion ( &memRegion,
				   sendBuffers,
				   buffersSize);
  if ( rc != 0)
    {
      printf("Kernel_CreateMemoryRegion failed with rc=%d\n",rc);
      alltoall_exit(1);
    }

  sendBufPAddr =
    (uint64_t)sendBuffers -
    (uint64_t)memRegion.BaseVa +
    (uint64_t)memRegion.BasePa;

}

static void create_random_permutation(void)
{
  unsigned seedPtr = myRank;
  srandom(myRank);

  int i;
  randPerm[0] = 0;
  for ( i = 1; i < COMMDIR_COUNT; i++ )
    {
      //uint32_t r = random();
      // scale r into the range [0,i]
      //r = r % (i+1); // TODO: deal with modulo bias ?
      //randPerm[i] = randPerm[r];
      //randPerm[r] = i;
      randPerm[i] = i;
    }
}

static void create_descriptors2()
{
  uint64_t anyFifoMap =
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AM |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AP |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BM |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BP |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CM |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CP |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DM |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DP |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EM |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EP;

  uint64_t r_offset = myRank * messageSizeInBytes;
  int i;
  uint64_t s_offset;
  static int  did_print =0;

  for ( i = 0, s_offset = 0;
	i < COMMDIR_COUNT;
	i++, s_offset += messageSizeInBytes )
    {


      // Injection Direct Put Descriptor Information Structure
      MUSPI_Pt2PtDirectPutDescriptorInfo_t dinfo;

      memset( (void*)&dinfo, 0x00, sizeof(dinfo) );

//      dinfo.Base.Payload_Address = sendBufPAddr + roffsets[i];
      dinfo.Base.Message_Length  = messageSizeInBytes;

      dinfo.Base.Torus_FIFO_Map  = anyFifoMap;

      dinfo.Base.Dest = nb2dest[i].dest;

      dinfo.Pt2Pt.Hints_ABCD = nb2dest[i].hintsABCD;
      if ( do_dynamic)
      {
	dinfo.Pt2Pt.Misc1 =
			nb2dest[i].hintsE |
	  MUHWI_PACKET_USE_DYNAMIC_ROUTING |
	  MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE;

	dinfo.Pt2Pt.Misc2 =
	  MUHWI_PACKET_VIRTUAL_CHANNEL_DYNAMIC |
	  zoneRoutingMask |
	  stayOnBubbleMask;
	if ( (myRank ==0) && (did_print ==0)) printf(" dyanmic routing  zoneRoutingMask=%d stayOnBubbleMask=%d\n",
						     zoneRoutingMask, stayOnBubbleMask);
      }
      else
      {
	dinfo.Pt2Pt.Misc1 =
			nb2dest[i].hintsE |
	  MUHWI_PACKET_USE_DETERMINISTIC_ROUTING |
	  MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE;

	dinfo.Pt2Pt.Misc2 =
	  MUHWI_PACKET_VIRTUAL_CHANNEL_DETERMINISTIC |
	  zoneRoutingMask |
	  stayOnBubbleMask;
	if ( (myRank ==0) && (did_print ==0)) printf(" deterministic routing\n");
      }
      did_print++;


      dinfo.Pt2Pt.Skip  = 8; // for checksumming, skip the header
      dinfo.DirectPut.Rec_Payload_Base_Address_Id = recvBufBatId;
      dinfo.DirectPut.Rec_Payload_Offset          = r_offset;
      dinfo.DirectPut.Rec_Counter_Base_Address_Id = recvCntrBatId;
      dinfo.DirectPut.Rec_Counter_Offset          = 0;

      dinfo.DirectPut.Pacing = MUHWI_PACKET_DIRECT_PUT_IS_NOT_PACED;

      int rc = MUSPI_CreatePt2PtDirectPutDescriptor(
						    &muDescriptors[i],
						    &dinfo );
      if (rc != 0)
	{
	  printf("MUSPI_CreatePt2PtDirectPutDescriptor failed with rc=%d\n",rc);
	  alltoall_exit(1);
	}

    } // End: Set up descriptors
}


static int get_destinations(unsigned int *mypers) {

  int tmp[6];
#if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  MPI_Status mstatus;
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_t_up, 0,
	       (void*)tmp, 6, MPI_INT, g_nb_t_dn, 0,
	       g_cart_grid, &mstatus); //TDOWN
  MUSPI_SetUpDestination( &nb2dest[1].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_t_dn, 1,
	       (void*)tmp, 6, MPI_INT, g_nb_t_up, 1,
	       g_cart_grid, &mstatus); //TUP
  MUSPI_SetUpDestination( &nb2dest[0].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );
#endif
#if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_x_up, 2,
	       (void*)tmp, 6, MPI_INT, g_nb_x_dn, 2,
	       g_cart_grid, &mstatus); //XDOWN
  MUSPI_SetUpDestination( &nb2dest[3].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_x_dn, 3,
	       (void*)tmp, 6, MPI_INT, g_nb_x_up, 3,
	       g_cart_grid, &mstatus); //XUP
  MUSPI_SetUpDestination( &nb2dest[2].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );
#endif
#if (defined PARALLELXYT || defined PARALLELXYZT)
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_y_up, 4,
	       (void*)tmp, 6, MPI_INT, g_nb_y_dn, 4,
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[5].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_y_dn, 5,
	       (void*)tmp, 6, MPI_INT, g_nb_y_up, 5,
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[4].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );
#endif
#if (defined PARALLELXYZT)
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_z_up, 6,
	       (void*)tmp, 6, MPI_INT, g_nb_z_dn, 6,
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[7].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_z_dn, 7,
	       (void*)tmp, 6, MPI_INT, g_nb_z_up, 7,
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[6].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );
#endif
  return(0);
}










/* begin_generated_IBM_copyright_prolog                             */
/*                                                                  */
/* This is an automatically generated copyright prolog.             */
/* After initializing,  DO NOT MODIFY OR MOVE                       */
/* ================================================================ */
/*                                                                  */
/* Licensed Materials - Property of IBM                             */
/*                                                                  */
/* Blue Gene/Q                                                      */
/*                                                                  */
/* (C) Copyright IBM Corp.  2008, 2012                              */
/*                                                                  */
/* US Government Users Restricted Rights -                          */
/* Use, duplication or disclosure restricted                        */
/* by GSA ADP Schedule Contract with IBM Corp.                      */
/*                                                                  */
/* This software is available to you under the                      */
/* Eclipse Public License (EPL).                                    */
/*                                                                  */
/* ================================================================ */
/*                                                                  */
/* end_generated_IBM_copyright_prolog                               */


#include <sys/types.h>
#include <stdint.h>
#include <hwi/include/bqc/MU_PacketCommon.h>
#include <spi/include/kernel/MU.h>
#include <spi/include/mu/GIBarrier.h>



/**
 * \brief Initialize Buffer with Random Data
 *
 * Initialize a buffer with random data, using the APPSEED env var random number seed,
 * the message's source coordinates, the message's destination coordinates, and a
 * specified unique value.
 *
 * The buffer is initialized backwards so it can be checked backwards.
 *
 * \param [in]  bufPtr     Pointer to the buffer to be initialized.
 * \param [in]  size       Size of the buffer, in bytes.
 * \param [in]  sourcePtr  Pointer to source coordinates (may be NULL)
 * \param [in]  destPtr    Pointer to destination coordinates (may be NULL)
 * \param [in]  uniquePtr  Pointer to unique value (may be NULL)
 * \param [in]  chksumPtr  Pointer to a location where the chksum of the buffer
 *                         is stored (may be NULL).  This result is only
 *                         accurate if size is a multiple of 8.
 *
 * \returns  The buffer is initialized with random data.
 *           The chksum result is returned in the location pointed-to by
 *           chksumPtr.
 */
static void msg_InitBuffer ( void                *bufPtr,
                      size_t               size,
                      MUHWI_Destination_t *sourcePtr,
                      MUHWI_Destination_t *destPtr,
                      uint32_t            *uniquePtr,
                      uint64_t            *chksumPtr )
{
  /* memset (bufPtr, size, 0); *//* Simple initialization test */
  size_t i = 0;
  char *buf = (char *) bufPtr;
  for (i = 0; i < size; ++i) {
    buf[i] = i & 0xFF;
  }
}

/**
 * \brief Check Buffer with Random Data
 *
 * Check a buffer with random data, using the APPSEED env var random number seed,
 * the message's source coordinates, the message's destination coordinates,
 * and a specified unique value.
 *
 * The buffer is checked backwards.  This is because the last bytes of the
 * buffer are the last to arrive, and we want to ensure those bytes are there
 * at the time we begin checking the buffer.
 *
 * \param [in]  bufPtr     Pointer to the buffer to be checked.
 * \param [in]  size       Size of the buffer, in bytes.
 * \param [in]  sourcePtr  Pointer to source coordinates (may be NULL)
 * \param [in]  destPtr    Pointer to destination coordinates (may be NULL)
 * \param [in]  uniquePtr  Pointer to unique value (may be NULL)
 * \param [in]  textPtr    Pointer to text shown in error messages
 *
 * \retval  0  The buffer checks-out successfully.
 * \retval -1  The buffer failed to check correctly.  A message was printed
 *             showing where the check failed.
 *
 * \see msg_InitBuffer()
 */
static int msg_CheckBuffer ( void                *bufPtr,
                      size_t               size,
                      MUHWI_Destination_t *sourcePtr,
                      MUHWI_Destination_t *destPtr,
                      uint32_t            *uniquePtr,
                      char                *textPtr )
{
  size_t i = 0;
  char *buf = (char *) bufPtr;
  for (i = 0; i < size; ++i) {
    if (buf[i] != (i & 0xFF))
      return 1;
  }

  return 0;
}

/**
 * \brief Poll a Reception Counter Until It Hits Zero
 *
 * This function polls the specified reception counter until it hits
 * zero, with no timeout.
 *
 * \param [in]  receptionCounter  Pointer to the reception counter
 *
 */
static void msg_CounterPoll( volatile uint64_t *receptionCounter )
{
  while (1)
    {
      if (*receptionCounter == 0)
        {
          break;
        }
    }
}

/**
 * \brief Injection Fifo Info Structure
 *
 * This structure is used within the implementation of the msg_InjFifoXXXXX()
 * functions.
 */
typedef struct msg_InjFifoInfo
{
  MUSPI_InjFifoSubGroup_t  subgroup[BGQ_MU_NUM_FIFO_SUBGROUPS_PER_NODE];
  uint32_t                 numFifosInSubgroup[BGQ_MU_NUM_FIFO_SUBGROUPS_PER_NODE];
  void                    *fifoMemoryPtr [BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP *
                                          BGQ_MU_NUM_FIFO_SUBGROUPS_PER_NODE];
  void                    *fifoPtr       [BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP *
                                          BGQ_MU_NUM_FIFO_SUBGROUPS_PER_NODE];
  uint32_t                 startingSubgroupId;
  uint32_t                 startingFifoId;
  uint32_t                 numFifos;
  uint32_t                 numSubgroups;
} msg_InjFifoInfo_t;


static int msg_InjFifoInit2( msg_InjFifoHandle_t *injFifoHandlePtr,
                      uint32_t             startingSubgroupId,
                      uint32_t             startingFifoId,
                      uint32_t             numFifos,
                      size_t               fifoSize,
                      Kernel_InjFifoAttributes_t  *injFifoAttrs )
{
  void                *buffer = NULL;
  uint32_t endingFifoId; // Relative to a subgroup
  uint32_t numFifosInSubgroup;
  int rc;
  uint32_t subgroupId = startingSubgroupId;
  uint32_t fifoIds[BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP];
  Kernel_InjFifoAttributes_t attrs[BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP];
  Kernel_InjFifoAttributes_t defaultAttrs;
  unsigned int i;
  uint64_t lock_cache;

  memset ( &defaultAttrs, 0x00, sizeof(defaultAttrs) );
  if ( injFifoAttrs == NULL ) injFifoAttrs = &defaultAttrs;

  // Malloc space for the info structure
  msg_InjFifoInfo_t *info;
  info = (msg_InjFifoInfo_t *) memalign(32, sizeof(msg_InjFifoInfo_t));
  if ( !info ) return -1;

    // Initialize the info structure
  info->startingSubgroupId = startingSubgroupId;
  info->startingFifoId     = startingFifoId;
  info->numFifos           = numFifos;
  info->numSubgroups       = 0;

  // Malloc space for the injection fifos.  They are 64-byte aligned.
  for (i=0; i<numFifos; i++)
    {
      info->fifoPtr[i] = memalign(64, fifoSize);
      if ( !info->fifoPtr[i] ) return -1;
    }

  // Process one subgroup at a time.
  // - Allocate the fifos.
  // - Init the MU MMIO for the fifos.
  // - Activate the fifos.
  while ( numFifos > 0 )
    {
      info->numSubgroups++;

      // startingFifoId is the starting fifo number relative to the
      // subgroup we are working on.
      // Determine endingFifoId, the ending fifo number relative to
      // the subgroup we are working on.
      endingFifoId = startingFifoId + numFifos-1;
      if ( endingFifoId > (BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP-1) )
        endingFifoId = BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP-1;
      numFifosInSubgroup = endingFifoId - startingFifoId + 1;
      info->numFifosInSubgroup[subgroupId] = numFifosInSubgroup;

      // Init structures for allocating the fifos...
      // - fifo Ids
      // - attributes
      for (i=0; i<numFifosInSubgroup; i++)
        {
          fifoIds[i] = startingFifoId + i;
          memcpy(&attrs[i],injFifoAttrs,sizeof(attrs[i]));
/*        printf("Attrs[%u] = 0x%x\n",i,*((uint32_t*)&attrs[i])); */
/*        printf("InjFifoInit: fifoIds[%u]=%u\n",i,fifoIds[i]); */
        }

      // Allocate the fifos
      rc = Kernel_AllocateInjFifos (subgroupId,
                                    &info->subgroup[subgroupId],
                                    numFifosInSubgroup,
                                    fifoIds,
                                    attrs);
      if ( rc ) {
        printf("msg_InjFifoInit: Kernel_AllocateInjFifos failed with rc=%d\n",rc);
        return rc;
      }

      // Init the MU MMIO for the fifos.
      for (i=0; i<numFifosInSubgroup; i++)
        {
          Kernel_MemoryRegion_t memRegion;
          rc = Kernel_CreateMemoryRegion ( &memRegion,
                                           info->fifoPtr[numFifos-i-1],
                                           fifoSize );
          if ( rc ) {
            printf("msg_InjFifoInit: Kernel_CreateMemoryRegion failed with rc=%d\n",rc);
            return rc;
          }

          rc = Kernel_InjFifoInit (&info->subgroup[subgroupId],
                                   fifoIds[i],
                                   &memRegion,
                                   (uint64_t)info->fifoPtr[numFifos-i-1] -
                                   (uint64_t)memRegion.BaseVa,
                                   fifoSize-1);
          if ( rc ) {
            printf("msg_InjFifoInit: Kernel_InjFifoInit failed with rc=%d\n",rc);
            return rc;
          }

/*        TRACE(("HW freespace=%lx\n", MUSPI_getHwFreeSpace (MUSPI_IdToInjFifo (fifoIds[i],&info->subgroup[subgroupId]))))
; */
        }

      // Activate the fifos.
      rc = Kernel_InjFifoActivate (&info->subgroup[subgroupId],
                                   numFifosInSubgroup,
                                   fifoIds,
                                   KERNEL_INJ_FIFO_ACTIVATE);
      if ( rc ) {
        printf("msg_InjFifoInit: Kernel_InjFifoActivate failed with rc=%d\n",rc);
        return rc;
      }

      startingFifoId = 0; // Next subgroup will start at fifo 0.

      subgroupId++;       // Next subgroup.
      numFifos -= numFifosInSubgroup;
    }

  injFifoHandlePtr->pOpaqueObject = (void *)info;
  return 0;
}



/**
 * \brief Terminate Injection Fifos
 *
 * Terminate the usage of injection fifos.  This deactivates the fifos and
 * frees all of the storage associated with them (previously allocated during
 * msg_InjFifoInit()).
 *
 * \param [in]  injFifoHandle  The handle returned from msg_InjFifoInit().
 *                             It must be passed into this function untouched
 *                             from when it was returned from msg_InjFifoInit().
 *
 * \note After this function returns, no more InjFifo functions should be called
 *       with this injFifoHandle.
 */
static void msg_InjFifoTerm ( msg_InjFifoHandle_t injFifoHandle )
{
  return; /*Simple library do nothing! */
}


/**
 * \brief Inject Descriptor into Injection Fifo
 *
 * Inject the specified descriptor into the specified injection fifo.
 *
 * \param [in]  injFifoHandle  The handle returned from msg_InjFifoInit().
 *                             It must be passed into this function untouched
 *                             from when it was returned from msg_InjFifoInit().
 * \param [in]  relativeFifoId  The fifo number, relative to the start of
 *                              the fifos managed by this opaque object.
 *                              For example, if msg_InjFifoInit() was called
 *                              to init fifos in subgroup 2, starting with
 *                              fifo Id 3, the relativeFifoNumber of the
 *                              first fifo is 0, not 3.
 * \param [in]  descPtr         Pointer to the descriptor to be injected.
 *
 * \retval  positiveNumber  The descriptor was successfully injected.  The
 *                          returned value is the sequence number of this
 *                          descriptor.
 * \retval  -1              The descriptor was not injected, most likely because
 *                          there is no room in the fifo.
 */
static uint64_t msg_InjFifoInject2 ( msg_InjFifoHandle_t injFifoHandle,
                             uint32_t            relativeFifoId,
                             MUHWI_Descriptor_t *descPtr )
{
  msg_InjFifoInfo_t *info = (msg_InjFifoInfo_t*)injFifoHandle.pOpaqueObject;

  uint32_t globalFifoId = (info->startingSubgroupId * BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP) +
    info->startingFifoId + relativeFifoId;

  uint32_t subgroupId   = globalFifoId / BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP;
  uint64_t rc = MUSPI_InjFifoInject (MUSPI_IdToInjFifo( globalFifoId % BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP,
							&info->subgroup[subgroupId] ),
				     descPtr);
  return rc;
}


static unsigned msg_InjFifoCheckCompletion2(msg_InjFifoHandle_t injFifoHandle,
                                    uint32_t            relativeFifoId,
                                    uint64_t            desc_count)
{
  msg_InjFifoInfo_t *info = (msg_InjFifoInfo_t*)injFifoHandle.pOpaqueObject;

  uint32_t globalFifoId = (info->startingSubgroupId * BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP) +
    info->startingFifoId + relativeFifoId;

  uint32_t subgroupId   = globalFifoId / BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP;

  return MUSPI_CheckDescComplete(MUSPI_IdToInjFifo( globalFifoId % BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP,
                                                    &info->subgroup[subgroupId] ),
                                 desc_count);
}




#endif





static MPI_Request g_bgq_request_recv[PHYSICAL_LD];
static MPI_Request g_bgq_request_send[PHYSICAL_LD];


#ifdef SPI
static unsigned mySpirank;
static inline unsigned bgq_abcde2spirank(Personality_t *pers, uint8_t a, uint8_t b,uint8_t c,uint8_t d,uint8_t e) {
	assert(pers);
	torus_t tdims = {
		pers->Network_Config.Anodes,
		pers->Network_Config.Bnodes,
		pers->Network_Config.Cnodes,
		pers->Network_Config.Dnodes,
		pers->Network_Config.Enodes
	};
	torus_t dims = {
	  pers->Network_Config.Anodes,
	  pers->Network_Config.Bnodes,
	  pers->Network_Config.Cnodes,
	  pers->Network_Config.Dnodes,
	  pers->Network_Config.Enodes
	};

	unsigned numNodes = tdims.a * tdims.b * tdims.c * tdims.d * tdims.e;
	unsigned result = ((((a)*dims.b + b)*dims.c + c)*dims.d + d)*dims.e + e;
	assert(result < numNodes);
	return result;
}


static void setup_destinations(Personality_t *pers) {
	torus_t tcoords = {
			pers->Network_Config.Acoord,
			pers->Network_Config.Bcoord,
			pers->Network_Config.Ccoord,
			pers->Network_Config.Dcoord,
			pers->Network_Config.Ecoord
	};

	torus_t tdims = {
			pers->Network_Config.Anodes,
			pers->Network_Config.Bnodes,
			pers->Network_Config.Cnodes,
			pers->Network_Config.Dnodes,
			pers->Network_Config.Enodes
	};

	//numNodes = tdims.a * tdims.b * tdims.c * tdims.d * tdims.e;
	mySpirank = bgq_abcde2spirank(pers, tcoords.a, tcoords.b, tcoords.c, tcoords.e, tcoords.d);
	//myRank = g_proc_id;
	//myRank = mySpirank;

	if (Kernel_ProcessCount() > 1) {
		printf("ERROR: this test only works with ppn==1\n");
		abort();
	}

	for (ucoord cd = 0; cd < COMMDIR_COUNT; cd += 1) {
		bgq_direction d_src = bgq_commdir2direction(cd);
		bgq_direction d_dst = bgq_direction_revert(d_src);

		MPI_Status status;
		int mpirank_to = bgq_direction2rank(d_src);
		int mpirank_from = bgq_direction2rank(d_dst);
		torus_t nb = {-1,-1,-1,-1,-1};
		MPI_CHECK(MPI_Sendrecv(&tcoords, sizeof(tcoords), MPI_BYTE, mpirank_to, 0, &nb, sizeof(nb), MPI_BYTE, mpirank_from, 0, g_cart_grid, &status));
		assert(get_MPI_count(&status) == sizeof(tcoords));

		ucoord nbrank = bgq_abcde2spirank(pers, nb.a, nb.b, nb.c, nb.d, nb.e);
		nb2dest[cd].tcoord = nb; // Just for reference
		MUSPI_SetUpDestination(&nb2dest[cd].dest, nb.a, nb.b, nb.c, nb.d, nb.e);
		nb2dest[cd].hintsABCD = 0;
		nb2dest[cd].hintsE = 0;
		printf("node %d: %d(%d,%d,%d,%d,%d)-%d->%d(%d,%d,%d,%d,%d)\n", g_proc_id, mySpirank, tcoords.a, tcoords.b, tcoords.c, tcoords.e, tcoords.d, cd, nbrank, nb.a, nb.b, nb.c, nb.d, nb.e);
	}
}
#endif



static void bgq_comm_test(bool nospi) {
	  // test communication

	for (ucoord cd = 0; cd < COMMDIR_COUNT; cd+=1) {
		bgq_direction d = bgq_commdir2direction(cd);
		bgq_weylfield_section sec = bgq_direction2section(d, true);
		ucoord index_end = bgq_section_size(sec) / sizeof(bgq_weyl_vec);
		for (size_t i = 0; i < index_end; i += 1) {
			for (ucoord v = 0; v < 2; v += 1){
				for (ucoord c = 0; c < 3; c += 1){
					for (ucoord k = 0; k < 2; k += 1){
						g_bgq_sec_send[d][i].s[v][c][k] = (double) g_cart_id + g_proc_id * _Complex_I;
					}
				}
			}
		}
	}

	bgq_comm_recv(nospi);
	bgq_comm_send(nospi);
	bgq_comm_wait(nospi);

	for (ucoord cd = 0; cd < COMMDIR_COUNT; cd+=1) {
		bgq_direction d = bgq_commdir2direction(cd);
		assert(bgq_direction_isDistributed(d));

		bgq_weylfield_section sec = bgq_direction2section(d, false);
		ucoord index_end = bgq_section_size(sec) / sizeof(bgq_weyl_vec);
		int rank_neighbor = bgq_direction2rank(d);
		for (size_t i = 0; i < index_end; i += 1) {
			for (ucoord v = 0; v < 2; v += 1) {
				for (ucoord c = 0; c < 3; c += 1) {
					for (ucoord k = 0; k < 2; k += 1) {
						complexdouble val = g_bgq_sec_recv[d][i].s[v][c][k];
						if (creal(val) != rank_neighbor) {
							printf("Node %d: Exchange doesn't work for dir %d: %d != %f (proc=%f) at point %zu\n", g_proc_id, d, rank_neighbor, creal(val), cimag(val), i);
							exit(1);
						}
					}
				}
			}
		}
	}

	master_print("Communication nospi=%d tested successfully\n", nospi);
}


static bool g_bgq_comm_common_initialized = false;
static void bgq_comm_common_init(void) {
	if (g_bgq_comm_common_initialized)
		return;
	g_bgq_comm_common_initialized = true;
	bgq_indices_init();

	size_t commbufsize = bgq_weyl_section_offset(sec_comm_end) - bgq_weyl_section_offset(sec_comm);
	uint8_t *buf = (uint8_t*)malloc_aligned(commbufsize, BGQ_ALIGNMENT_L2);
	uint8_t *bufend = buf + commbufsize;
	g_bgq_sec_comm = buf;
	for (bgq_direction d = 0; d < PHYSICAL_LD; d+=1) {
		bgq_dimension dim = bgq_direction2dimension(d);

		bgq_weylfield_section sec_send = bgq_direction2section(d, true);
		bgq_weylfield_section sec_recv = bgq_direction2section(d, !bgq_dimension_isDistributed(dim));

		g_bgq_sec_send[d] = (bgq_weyl_vec*)(buf + bgq_weyl_section_offset(sec_send) - bgq_weyl_section_offset(sec_comm));
		assert((uint8_t*)g_bgq_sec_send[d] <= bufend);

		g_bgq_sec_recv[d] = (bgq_weyl_vec*)(buf + bgq_weyl_section_offset(sec_recv) - bgq_weyl_section_offset(sec_comm));
		assert((uint8_t*)g_bgq_sec_recv[d] <= bufend);

		assert((uintptr_t)g_bgq_sec_send[d] % BGQ_ALIGNMENT_L2 == 0);
		assert((uintptr_t)g_bgq_sec_recv[d] % BGQ_ALIGNMENT_L2 == 0);
	}
}


static bool g_bgq_comm_mpi_initialized = false;
void bgq_comm_mpi_init(void) {
	if (g_bgq_comm_mpi_initialized)
		return;
	g_bgq_comm_mpi_initialized = true;
	bgq_comm_common_init();

	for (ucoord cd = 0; cd < COMMDIR_COUNT; cd+=1) {
		bgq_direction d_src = bgq_commdir2direction(cd);
		assert(bgq_direction_isDistributed(d_src));

		//master_print("d_src=%d\n",d_src);
		bgq_direction d_dst = bgq_direction_revert(d_src);
		//master_print("d_dst=%d\n",d_dst);
		bgq_dimension dim = bgq_direction2dimension(d_src);
		//master_print("dim=%d\n",dim);

		ucoord commdir_src = cd;
		//master_print("commdir_src=%d\n",commdir_src);
		bgq_weylfield_section sec_recv = bgq_direction2section(d_src, false);
		size_t secsize = bgq_weyl_section_offset(sec_recv+1) - bgq_weyl_section_offset(sec_recv);
		assert(secsize > 0);
		MPI_CHECK(MPI_Recv_init(g_bgq_sec_recv[d_src], secsize / sizeof(complexdouble), MPI_DOUBLE_COMPLEX, bgq_direction2rank(d_src), d_dst, g_cart_grid, &g_bgq_request_recv[cd]));
		//master_print("MPI_CHECK(MPI_Recv_init(%zu, %zu, %zu, %zu, %zu, %zu, %zu))\n", (size_t)(g_bgq_sec_recv[d_src]), secsize / sizeof(double), (size_t)(MPI_DOUBLE), (size_t)bgq_direction2rank(d_src), (size_t)d_dst, (size_t)g_cart_grid, (size_t)(&g_bgq_request_recv[commdir_src]));

		ucoord commdir_dst = bgq_direction2commdir(d_dst);
		bgq_weylfield_section sec_send = bgq_direction2section(d_dst, true);
		assert(secsize == bgq_weyl_section_offset(sec_send+1) - bgq_weyl_section_offset(sec_send));
		MPI_CHECK(MPI_Send_init(g_bgq_sec_send[d_dst], secsize / sizeof(complexdouble), MPI_DOUBLE_COMPLEX, bgq_direction2rank(d_dst), d_dst, g_cart_grid, &g_bgq_request_send[cd]));
		//master_print("MPI_CHECK(MPI_Send_init(%zu, %zu, %zu, %zu, %zu, %zu, %zu))\n", g_bgq_sec_send[d_dst], secsize / sizeof(double), MPI_DOUBLE, bgq_direction2rank(d_dst), d_dst, g_cart_grid, &g_bgq_request_send[commdir_dst]);
	}


	bgq_comm_test(true);
}



static bool g_bgq_comm_spi_initialized = false;
void bgq_comm_spi_init(void) {
#ifdef SPI
	if (g_bgq_comm_spi_initialized)
		return;
	g_bgq_comm_spi_initialized = true;
	bgq_comm_common_init();


	size_t messageSizes[PHYSICAL_LD];
	size_t roffsets[PHYSICAL_LD];
	size_t soffsets[PHYSICAL_LD];
	//size_t totalMessageSize = 0;

	// here comes the SPI initialization
	int spi_num_dirs = COMMDIR_COUNT;

	for (ucoord cd = 0; cd < COMMDIR_COUNT; cd+=1) {
		bgq_direction d_src = bgq_commdir2direction(cd);
		bgq_direction d_dst = bgq_direction_revert(d_src);
		bgq_dimension dim = bgq_direction2dimension(d_src);
		assert(bgq_direction_isDistributed(d_src));

		ucoord commdir_src = bgq_direction2commdir(d_src);
		bgq_weylfield_section sec_recv = bgq_direction2section(d_src, false);
		size_t secsize = bgq_section_size(sec_recv);
		ucoord commdir = commdir_src;

		ucoord commdir_dst = bgq_direction2commdir(d_dst);
		bgq_weylfield_section sec_send = bgq_direction2section(d_dst, true);
		assert(secsize == bgq_section_size(sec_send));


		messageSizes[commdir] = secsize;
		soffsets[commdir] = bgq_weyl_section_offset(sec_send) - bgq_weyl_section_offset(sec_send_begin);
		assert(soffsets[commdir] % BGQ_ALIGNMENT_L2 == 0);
		assert(soffsets[commdir] + messageSizes[commdir] <= bgq_weyl_section_offset(sec_send_end) - bgq_weyl_section_offset(sec_send_begin));
		roffsets[cd] = bgq_weyl_section_offset(sec_recv) - bgq_weyl_section_offset(sec_recv_begin);
		assert(roffsets[cd] % BGQ_ALIGNMENT_L2 == 0);
		assert((roffsets[cd] + secsize) <= (bgq_weyl_section_offset(sec_recv_end) - bgq_weyl_section_offset(sec_recv_begin)));
		totalMessageSize += secsize;

		master_print("SPI %llu: d=%llu msize=%zu soffset=%zu d_dst=%llu roffset=%zu\n", cd, d_src, messageSizes[commdir], soffsets[commdir], d_dst, roffsets[cd]);
	}
	assert(totalMessageSize == bgq_weyl_section_offset(sec_recv_end) - bgq_weyl_section_offset(sec_recv_begin));

	do_dynamic = 0; // Use static routing (since only neighbor-to-neighbor communication)


	//do_dynamic = 0; // Use static routing (since only neighbor-to-neighbor communication)
	//char *argv[] = {"exec"};
	//main2(1, argv);
	//master_print("Main2 done\n");

	Personality_t pers;
	int rc = 0;
	// get the CNK personality
	Kernel_GetPersonality(&pers, sizeof(pers));
	setup_destinations(&pers);

	// adjust the SPI pointers to the send and receive buffers
	SPIrecvBuffers = (char*)g_bgq_sec_recv[0];
	assert((uintptr_t)SPIrecvBuffers % BGQ_ALIGNMENT_L2 == 0);
	SPIsendBuffers = (char*)g_bgq_sec_send[0];
	assert((uintptr_t)SPIsendBuffers % BGQ_ALIGNMENT_L2 == 0);

	// Setup the FIFO handles
	rc = msg_InjFifoInit(&injFifoHandle,
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
	uint64_t bufferSize = bgq_weyl_section_offset(sec_recv_end) - bgq_weyl_section_offset(sec_recv_begin);
	assert(bufferSize == totalMessageSize);
	setup_mregions_bats_counters(bufferSize);

	// Create descriptors
	// Injection Direct Put Descriptor, one for each neighbour
	SPIDescriptors = (MUHWI_Descriptor_t*)(((uint64_t)SPIDescriptorsMemory+64)&~(64-1));
	create_descriptors(SPIDescriptors, messageSizes, soffsets, roffsets, spi_num_dirs);

    // Initialize the barrier, resetting the hardware.
    rc = MUSPI_GIBarrierInit(&GIBarrier, 0 /*comm world class route */);
    if (rc) {
    	printf("MUSPI_GIBarrierInit returned rc=%d\n",rc);
    	abort();
    }

#if 0
    for (size_t cd = 0; cd < COMMDIR_COUNT; cd+=1) {
    	bgq_direction d = bgq_commdir2direction(cd);
    	bgq_weylfield_section sec_send = bgq_direction2section(d,true);
    	size_t secsize = bgq_section_size(sec_send);

    	SPIDescriptors[cd].Message_Length = messageSizes[cd];
    	SPIDescriptors[cd].Pa_Payload     = sendBufPAddr + soffsets[cd];
    	MUSPI_SetRecPutOffset(&SPIDescriptors[cd], roffsets[cd]);
    }
#endif


	bgq_comm_test(false);
#endif
}


//TODO: inline?
void bgq_comm_recv(bool nospi) {
	assert(omp_get_thread_num()==0);
	master_print("Comm Receiving...\n");
#ifdef SPI
	if (!nospi) {
	    // reset the recv counter
	    recvCounter = totalMessageSize;
		return;
	}
#endif

	MPI_CHECK(MPI_Startall(COMMDIR_COUNT, g_bgq_request_recv));
}


void bgq_comm_send(bool nospi) {
	assert(omp_get_thread_num()==0);
	master_print("Comm Sending...\n");
#ifdef SPI
	if (!nospi) {
		// make sure everybody has reset recvCounter
		global_barrier(); //TODO: Can we get rid of it?
	    for (size_t cd = 0; cd < COMMDIR_COUNT; cd+=1) {
	      descCount[cd] = msg_InjFifoInject(injFifoHandle, cd, &SPIDescriptors[cd]);
	      if (descCount[cd] == -1) {
	    	 printf("msg_InjFifoInject failed, most likely because there is no room in the fifo\n");
	    	 abort();
	      }
	    }
		return;
	}
#endif

	MPI_CHECK(MPI_Startall(COMMDIR_COUNT, g_bgq_request_send));
}


void bgq_comm_wait(bool nospi) {
	assert(omp_get_thread_num()==0);
	master_print("Comm Waiting...\n");

#if BGQ_QPX
	uint64_t ppc32 = mfspr(SPRN_PPR32);
	ThreadPriority_Low(); // If there is some other work to be done on this node, give it priority
#endif

#ifdef SPI
	if (!nospi) {
		uint64_t startTime = 0;

		// Wait for all data is received
		 printf("node %d: %llu bytes to be received\n", g_proc_id, totalMessageSize);
		 while(recvCounter > 0) {
			 // Check range of pending bytes to receive
			 assert(recvCounter <= totalMessageSize);

			 if (GetTimeBase() - startTime >= 1600) {
				 printf("node %d: %llu bytes left\n", g_proc_id, recvCounter);
				 startTime = GetTimeBase();
			 }
		 }
		 printf("node %d: All data received\n", g_proc_id);

		 // Wait for all data sent
		while (true) {
			size_t sendDone = 0;
			for (unsigned j = 0; j < COMMDIR_COUNT; j+=1) {
				sendDone += msg_InjFifoCheckCompletion(injFifoHandle, j, descCount[j]);
			}
			if (sendDone == COMMDIR_COUNT)
				break;
		}

		_bgq_msync();  // Ensure data is available to all cores.
	} else
#endif
	{
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

#if BGQ_QPX
	mtspr(SPRN_PPR32, ppc32); // Restore original priority
#endif
}


