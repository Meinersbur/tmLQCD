#include "mypapi.h"
#include "config.h"


#if PAPI

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
//#include <common/alignment.h>
#include <papi.h> //Provides definitions for base PAPI function
//#include <spi/bgp_SPI.h> //The Blue Gene/P substrate for PAPI is based on the UPC interfaces documented in this header file.
#include <papiStdEventDefs.h> //Provides a listing of all standard PAPI presets. Please note that not all of these are available on the Blue Gene/P platform
//#include <linux-bgp-native-events.h> //Provides a listing of all available counters native to Blue Gene.
#include <omp.h>

//#include <linux-bgq.h>
//#include <fpapi.h>
#include <upci/events.h>
#include <mpi.h>
#include <assert.h>
#include <kernel/location.h>
#include <bgpm/include/bgpm.h>

//void FPUArith(void); //This method does various calculations which should saturate many of the counters
void List_PAPI_Events(const int pEventSet, int* pEvents, int* xNumEvents);
//void Print_Native_Counters();
//void Print_Native_Counters_via_Buffer(const BGP_UPC_Read_Counters_Struct_t* pBuffer);
//void Print_Native_Counters_for_PAPI_Counters(const int pEventSet);
//void Print_Native_Counters_for_PAPI_Counters_From_List(const int* pEvents, const int pNumEvents);
void Print_PAPI_Counters(const int pEventSet, const long long* pCounters);
//void Print_PAPI_Counters_From_List(const int* pEventList, const int pNumEvents, const long long* pCounters);
void Print_Counters(const int pEventSet);
void Print_PAPI_Events(const int pEventSet);
long long getMyPapiValue(const int eventNum);

#define lengthof(X) (sizeof(X)/sizeof((X)[0]))

int PAPI_Events[256];
long long PAPI_Counters[256];
int xEventSet=PAPI_NULL;

long long xCyc;
long long xNsec;
double xNow;
static double xWtime;
static double xOmpTime;

extern int g_proc_id;

static double now2(){
   struct timeval t; double f_t;
   gettimeofday(&t, NULL);
   f_t = t.tv_usec; f_t = f_t/1000000.0; f_t +=t.tv_sec;
   return f_t;
}

#define PAPI_ERROR(cmd)                                                                                      \
	do {                                                                                                       \
		int RC = (cmd);                                                                                       \
		if (RC != PAPI_OK) {                                                                                   \
			 fprintf(stderr, "MK_PAPI call failed with code %d at line %d on MPI rank %d thread %d: %s\n", RC, __LINE__,  g_proc_id, Kernel_ProcessorID(), TOSTRING(cmd)); \
		}                                                                                                        \
	} while (0)


#define BGPM_ERROR(cmd)                                                                                      \
	do {                                                                                                       \
		int RC = (cmd);                                                                                       \
		if (RC) {                                                                                   \
			 fprintf(stderr, "MK_BGPM call failed with code %d at line %d on MPI rank %d thread %d: %s\n", RC, __LINE__,  g_proc_id, Kernel_ProcessorID(), TOSTRING(cmd)); \
		}                                                                                                        \
	} while (0)


static int PAPI_add_native_event(int EventSet, int EventCode) {
	// For some strange reason (i.e. unknown to me), the numbers from events.h are off by one
	return PAPI_add_event(EventSet, PAPI_NATIVE_MASK | (EventCode-1));
}


mypapi_counters mypapi_merge_counters(mypapi_counters *counters1, mypapi_counters *counters2) {
	if (!counters1->init && !counters2->init) {
		mypapi_counters emptyresult = { 0 };
		return emptyresult;
	}

	mypapi_counters result;
	if (!counters2->init) {
		//if (g_proc_id == 0 && omp_get_thread_num() == 0)
		//		fprintf(stderr, "MK Take counter1: %llu\n", counters1->native[PEVT_CYCLES]);
		return *counters1;
	} else if (!counters1->init) {
		//if (g_proc_id == 0 && omp_get_thread_num() == 0)
		//		fprintf(stderr, "MK Take counter2: %llu\n", counters2->native[PEVT_CYCLES]);
		return *counters2;
	} else {
		if (counters1->set == counters2->set)
			result.set = counters1->set;
		else
			result.set = -1;

		if (counters1->eventset == counters2->eventset)
			result.eventset = counters1->eventset;
		else
			result.eventset = -1;

		if (counters1->threadid == counters2->threadid)
			result.threadid = counters1->threadid;
		else
			result.threadid = -1;

		if (counters1->coreid == counters2->coreid)
			result.coreid = counters1->coreid;
		else
			result.coreid = -1;

		if (counters1->smtid == counters2->smtid)
			result.smtid = counters1->smtid;
		else
			result.smtid = -1;

		if (counters1->ompid == counters2->ompid)
			result.ompid = counters1->ompid;
		else
			result.ompid = -1;

		//for (int i = 0; i < lengthof(counters1->preset); i += 1) {
		//	result.preset[i] = counters1->preset[i] + counters2->preset[i];
		//}
		for (int i = 0; i < lengthof(counters1->native); i += 1) {
			uint64_t c1 = counters1->native[i];
			uint64_t c2 = counters2->native[i];
			uint64_t merged;
			switch (i) {
			case PEVT_CYCLES:
				// Count PEVT_CYCLES just once
				if (counters1->set != 0)
					c1 = 0;
				if (counters2->set != 0)
					c2 = 0;
				// Fallthrough
			default:
				merged = c1 + c2;
				break;
			}
			result.native[i] = merged;
		}

		result.corecycles = ((counters1->set == 0) ?  counters1->corecycles : 0) + ((counters2->set == 0) ?  counters2->corecycles : 0);

		//if (g_proc_id == 0 && omp_get_thread_num() == 0)
		//		fprintf(stderr, "MK Merge result: %llu\n", result.native[PEVT_CYCLES]);
		result.init = true;
	}

	assert(result.init);
	return result;
}


static int mypapi_eventsets[64][MYPAPI_SETS];

static unsigned long int mypapi_getthreadid() {
	return Kernel_ProcessorID();
}


static int PuEventSets[MYPAPI_SETS][64] = {0};
static int L2EventSet[MYPAPI_SETS] = {0};


void mypapi_init() {
	assert(omp_get_thread_num() == 0);
	if (g_proc_id == 0)
		fprintf(stderr, "MK_Init mypapi\n");

#pragma omp parallel
	{
		BGPM_ERROR(Bgpm_Init(BGPM_MODE_SWDISTRIB));

		int tid = Kernel_ProcessorID();
		int cid = Kernel_ProcessorCoreID();
		int sid = Kernel_ProcessorThreadID();
		for (int i = 0; i < MYPAPI_SETS; i += 1) {
			PuEventSets[i][tid] = Bgpm_CreateEventSet();
			assert(PuEventSets[i][tid] >= 0);
		}

		if (tid == 0) {
			for (int i = 0; i < MYPAPI_SETS; i += 1) {
				L2EventSet[i] = Bgpm_CreateEventSet();
				assert(L2EventSet[i] >= 0);
			}
		}

		int j = 0;
		{
			int pues = PuEventSets[j][tid];
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_CYCLES));
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_INST_ALL));
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_LSU_COMMIT_LD_MISSES));
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_LSU_COMMIT_CACHEABLE_LDS));

			if (tid == 0) {
				int l2es = L2EventSet[j];
				BGPM_ERROR(Bgpm_DeleteEventSet(l2es));
				L2EventSet[j] = -1;
			}
		}

		j += 1;
		{
			int pues = PuEventSets[j][tid];
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_CYCLES));
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_L1P_BAS_MISS));
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_L1P_BAS_HIT)); // Hits in prefetch directory
			BGPM_ERROR(Bgpm_AddEvent(pues, PEVT_L1P_BAS_STRM_LINE_ESTB)); // Hits in prefetch directory

			if (tid == 0) {
				int l2es = L2EventSet[j];
				BGPM_ERROR(Bgpm_AddEvent(l2es, PEVT_L2_HITS));
				BGPM_ERROR(Bgpm_AddEvent(l2es, PEVT_L2_MISSES));
				BGPM_ERROR(Bgpm_AddEvent(l2es, PEVT_L2_FETCH_LINE)); // L2 lines loaded from main memory
				BGPM_ERROR(Bgpm_AddEvent(l2es, PEVT_L2_STORE_LINE)); // L2 lines stored to main memory
				BGPM_ERROR(Bgpm_AddEvent(l2es, PEVT_L2_PREFETCH));
				BGPM_ERROR(Bgpm_AddEvent(l2es, PEVT_L2_STORE_PARTIAL_LINE));
			}
		}
	}

	return;


	int retval = PAPI_library_init(PAPI_VER_CURRENT);
	if (retval != PAPI_VER_CURRENT) {
		fprintf(stderr, "Unexpected PAPI version\n");
		exit(1);
	}

	PAPI_ERROR(PAPI_thread_init(mypapi_getthreadid));

	for (int i = 0; i < lengthof(mypapi_eventsets); i+=1) {
		for (int j = 0; j < 64; j += 1) {
			mypapi_eventsets[j][i] = PAPI_NULL;
		}
	}

	for (int i = 0; i < 255; i++)
		PAPI_Counters[i]=0;

	const PAPI_hw_info_t *hwinfo = PAPI_get_hardware_info();
	int counter = PAPI_num_counters();
	if (g_proc_id == 0 && omp_get_thread_num() == 0) {
		printf("CPU frequency: %f Mhz\n", hwinfo->mhz);
		printf("TOT_CYC frequency: %f Mhz\n", hwinfo->clock_mhz);
		printf("CPU: %s %s\n", hwinfo->vendor_string, hwinfo->model_string);

		//BGP_UPC_Initialize_Counter_Config(BGP_UPC_MODE_DEFAULT, BGP_UPC_CFG_EDGE_DEFAULT);

		printf("Available hardware counters: %d\n", counter);
	}
	return;



	//PAPI_start_counters(*events, array_length)
	PAPI_ERROR(PAPI_create_eventset(&xEventSet));

	// Functional units
	// AXU = AXU Execution Unit Events (Quad Floating Point Unit)
	// IU = instruction unit (fetch,decode)
	// XU = Execution unit (branching,Load/Store,integer)
	// QFPU = Quad floating point unit
	// LSU = Load/Store unit
	// MMU = Memory Management unit
	// L1P = L1 (stream/list) prefetch unit
	// INST = Commit unit
	// L2 = L2 cache unit
	// MU = IO/message unit
	// NW = network unit
	// CNK = Interrupt controller

	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_LSU_COMMIT_LD_MISSES));
	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_INST_XU_LD));
	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_LSU_COMMIT_CACHEABLE_LDS));

	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_IU_IL1_MISS));


	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_L1P_LIST_MISMATCH)); // core address does not match a list address
	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_L1P_LIST_STARTED)); // List prefetch process was started
	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_L1P_LIST_OVF_MEM)); // Written pattern exceeded allocated buffer
	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_L1P_LIST_CMP_OVRUN_PREFCH)); // core address advances faster than prefetch lines can be established dropping prefetches
	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_L1P_LIST_SKIP));
	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_L1P_LIST_ABANDON));
	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_L1P_STRM_HIT_LIST));

	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_L2_HITS)); // (node) hits in L2, both load and store. Network Polling store operations from core 17 on BG/Q pollute in this count during normal use.
	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_L2_MISSES)); // (node) core address advances faster than prefetch lines can be established dropping prefetches


	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_L1P_STRM_LINE_ESTB)); // lines established for any reason ot thread
	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_L1P_STRM_EVICT_UNUSED)); // cacheline miss in L2 (both loads and stores)

	//PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_MMU_TLB_HIT_DIRECT_DERAT));
	//PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_MMU_TLB_MISS_DIRECT_DERAT));

	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_L1P_BAS_LD)); // Loads
	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_L1P_BAS_HIT)); // Hits in prefetch directory
	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_L1P_BAS_MISS)); // Misses in L1p by prefetchable loads

	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_XU_BR_COMMIT)); // Number of Branches committed
	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_XU_BR_MISPRED_COMMIT)); // Number of mispredicted Branches committed (does not include target address mispredicted)
	//PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_XU_BR_TARG_ADDR_MISPRED_COMMIT)); // Number of Branch Target addresses mispredicted committed

	//PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_RES_STL)); // Cycles stalled on any resource / Register Dependency Stall
	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_IU_IL1_MISS_CYC)); // IL1 Miss cycles

	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_IU_IS1_STALL_CYC)); // Dependency stalls + PEVT_IU_IS2_STALL_CYC
	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_IU_IS2_STALL_CYC)); // Functional unit busy because of SMT
	//PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_IU_RAW_DEP_HIT_CYC));
	//PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_IU_WAW_DEP_HIT_CYC));


	//PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_LSU_COMMIT_CACHE_INHIB_LD_MISSES));
	//PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_LSU_COMMIT_LD_MISSES));

	//PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_LSU_COMMIT_STS)); // Store instructions / Number of completed store commands.
	//PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_LSU_COMMIT_ST_MISSES));




	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_LSU_COMMIT_DCBT_MISSES));
	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_LSU_COMMIT_DCBT_HITS));




	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_TOT_CYC));
	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_CYCLES));
	//PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_TOT_IIS)); // Instructions issued / all instructions issued per thread
	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_INST_ALL)); // All Instruction Completions

	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_INST_QFPU_FPGRP2));
	//PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_INST_QFPU_FPGRP2_INSTR));

	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_L2_FETCH_LINE));
	PAPI_ERROR(PAPI_add_native_event(xEventSet, PEVT_L2_STORE_LINE));
}


static long long getMyPapiNativeValue(const int eventNum) {
	return getMyPapiValue(PAPI_NATIVE_MASK | (eventNum-1));
}


static long long getMyPapiValue(const int eventNum) {
	int i;
	char xName[256];
	//printf("Print_PAPI_Counters: PAPI_Counters*=%p, pCounters*=%p\n",PAPI_Counters, pCounters);

	int pNumEvents = PAPI_num_events(xEventSet);
	if (pNumEvents) {
		List_PAPI_Events(xEventSet, PAPI_Events, &pNumEvents);
		for (i=0; i<pNumEvents; i++) {
			//if (PAPI_event_code_to_name(PAPI_Events[i], xName)) {
			//	printf("PAPI_event_code_to_name failed on event code %d\n",PAPI_Events[i]);
			//	exit(1);
			//}
 			//printf("%20llu	%3d	0x%8.8x %s\n", PAPI_Counters[i], i, PAPI_Events[i],xName);
			//if ((eventNum  & 0xBFFFFFFF) == (PAPI_Events[i] & 0xBFFFFFFF))
			if (eventNum == PAPI_Events[i])
				return PAPI_Counters[i];
		}
	}

	return -1;
}


static double mypapi_wtime() {
	Personality_t personality;
	BGPM_ERROR(Kernel_GetPersonality(&personality, sizeof(Personality_t)));
	double freq = MEGA * personality.Kernel_Config.FreqMHz;
	long long cycles = GetTimeBase();
	return cycles / freq;
}

static int activeEventSet = -1;

void mypapi_start(int i) {
	assert(omp_get_thread_num() == 0);

	xCyc = GetTimeBase();
	//xCyc = PAPI_get_real_cyc();
	xNsec = mypapi_wtime();
	//xNsec = PAPI_get_real_nsec();
	xNow = now2();
	xWtime = MPI_Wtime();
	xOmpTime = omp_get_wtime();

	activeEventSet = i;

#pragma omp parallel
	{
		int tid = Kernel_ProcessorID();
		int cid = Kernel_ProcessorCoreID();
		int sid = Kernel_ProcessorThreadID();

		int pues = PuEventSets[i][tid];
		BGPM_ERROR(Bgpm_Apply(pues));
		if (tid == 0) {
			int l2es  = L2EventSet[i];
			if (l2es > 0) {
				BGPM_ERROR(Bgpm_Apply(l2es));
				BGPM_ERROR(Bgpm_Start(l2es));
			}
		}
		BGPM_ERROR(Bgpm_Start(pues));
	}
	return;


#pragma omp parallel
	{

	int tid = Kernel_ProcessorID();
	int eventset = mypapi_eventsets[tid][i];

	if (eventset == PAPI_NULL) {
		PAPI_ERROR(PAPI_create_eventset(&eventset));
		mypapi_eventsets[tid][i] = eventset;
		assert(eventset != PAPI_NULL);

		PAPI_ERROR(PAPI_add_native_event(eventset, PEVT_CYCLES));
		switch (i) {
		case 0:
			PAPI_ERROR(PAPI_add_native_event(eventset, PEVT_INST_ALL)); // All Instruction Completions

			PAPI_ERROR(PAPI_add_native_event(eventset, PEVT_LSU_COMMIT_LD_MISSES));
			PAPI_ERROR(PAPI_add_native_event(eventset, PEVT_LSU_COMMIT_CACHEABLE_LDS));
			PAPI_ERROR(PAPI_add_native_event(eventset, PEVT_L1P_BAS_MISS)); // Misses in L1p by prefetchable loads
			PAPI_ERROR(PAPI_add_native_event(eventset, PEVT_L1P_BAS_HIT)); // Hits in prefetch directory
			PAPI_ERROR(PAPI_add_native_event(eventset, PEVT_L1P_BAS_STRM_LINE_ESTB)); // Hits in prefetch directory
			break;
		case 1:
			PAPI_ERROR(PAPI_add_native_event(eventset, PEVT_L2_HITS)); // (node) hits in L2, both load and store. Network Polling store operations from core 17 on BG/Q pollute in this count during normal use.
			PAPI_ERROR(PAPI_add_native_event(eventset, PEVT_L2_MISSES)); // (node) core address advances faster than prefetch lines can be established dropping prefetches
			PAPI_ERROR(PAPI_add_native_event(eventset, PEVT_L2_FETCH_LINE)); // L2 lines loaded from memory
			PAPI_ERROR(PAPI_add_native_event(eventset, PEVT_L2_STORE_LINE)); // L2 lines stored to memory
			PAPI_ERROR(PAPI_add_native_event(eventset, PEVT_L2_PREFETCH)); // L2 lines stored to memory
			PAPI_ERROR(PAPI_add_native_event(eventset, PEVT_L2_STORE_PARTIAL_LINE)); // L2 lines stored to memory
			break;
		default:
			break;
		}
	}

#pragma omp master
	{
	xEventSet = eventset;

	xCyc = PAPI_get_real_cyc();
	xNsec = PAPI_get_real_nsec();
	xNow = now2();
	xWtime = MPI_Wtime();
	xOmpTime = omp_get_wtime();
	}

//#pragma omp parallel
	{
		//PAPI_ERROR(PAPI_reset(eventset));
		PAPI_ERROR(PAPI_start(eventset));
	}
	}
}


static mypapi_counters mypapi_readcounters(int pEventSet, long long *counters) {
	mypapi_counters result = {0};
	int events[256];

	result.threadid = Kernel_ProcessorID(); /* 0..63, Kernel_ProcessorCoreID() << 2 + Kernel_ProcessorThreadID() */
	result.coreid = Kernel_ProcessorCoreID(); /* 0..15 */
	result.smtid = Kernel_ProcessorThreadID(); /* 0..3 */
	result.ompid = omp_get_thread_num();

	int pNumEvents = PAPI_num_events(pEventSet);
	PAPI_ERROR(PAPI_list_events(pEventSet, events, &pNumEvents));

	for (int i = 0; i < pNumEvents; i++) {
		int event = events[i];

			event = (event & PAPI_NATIVE_AND_MASK) + 1;
			assert(event >= 0);
			assert(event < lengthof(result.native));
			long long value = counters[i];

			switch (event) {
			case PEVT_L2_HITS:
			case PEVT_L2_MISSES:
			case PEVT_L2_FETCH_LINE:
			case PEVT_L2_STORE_LINE:
			case PEVT_L2_PREFETCH:
			case PEVT_L2_STORE_PARTIAL_LINE:
				// Per node counters, store just one event
				if (result.threadid != 0)
					value = 0;
				break;
			}
	
			result.native[event] = value;
	}

	// Cycles are the some on the whole node, but shared between smt threads
	// Just one of 2 threads on the same core is issued, so we count cycles one per core
	if (result.smtid == 0)
		result.corecycles = result.native[PEVT_CYCLES];
	else
		result.corecycles = 0;

	result.init = true;
	return result;
}


static mypapi_counters mypapi_bgpm_read(int eventset, int set) {
	mypapi_counters result = { 0 };
	result.set = set;
	result.eventset = eventset;
	result.threadid = Kernel_ProcessorID(); /* 0..63, Kernel_ProcessorCoreID() << 2 + Kernel_ProcessorThreadID() */
	result.coreid = Kernel_ProcessorCoreID(); /* 0..15 */
	result.smtid = Kernel_ProcessorThreadID(); /* 0..3 */
	result.ompid = omp_get_thread_num();


	int numEvts = Bgpm_NumEvents(eventset);
	assert(numEvts > 0);

	uint64_t cnt;
	for (int i = 0; i < numEvts; i += 1) {
		int eventid = Bgpm_GetEventId(eventset, i);

		switch (eventid) {
		case PEVT_CYCLES:
			// Cycles are the some on the whole node, but shared between smt threads
			// Just one of 2 threads on the same core is issued, so we count cycles once per core
			if (result.smtid != 0)
				continue;
			break;
		case PEVT_L2_HITS:
		case PEVT_L2_MISSES:
		case PEVT_L2_FETCH_LINE:
		case PEVT_L2_STORE_LINE:
			// Per node counters, store just one event
			if (result.threadid != 0)
				continue;
			break;
		}

		uint64_t cnt;
		BGPM_ERROR(Bgpm_ReadEvent(eventset, i, &cnt));
		result.native[eventid] = cnt;

		if (g_proc_id == 0 && eventid == PEVT_CYCLES) {
		//	fprintf(stderr, "MK read result: %llu\n", cnt);
		}
	}

	result.init = true;
	if (g_proc_id == 0 && omp_get_thread_num() == 0) {
		//mypapi_print_counters(&result);
	}
	return result;
}


void mypapi_print_counters(mypapi_counters *counters) {
	fprintf(stderr, "*******************************************************************************\n");
	fprintf(stderr, "Set=%d eventset=%d thread=%d core=%d smt=%d omp=%d\n", counters->set, counters->eventset, counters->threadid, counters->coreid, counters->smtid, counters->ompid);
	if (counters->corecycles) {
		fprintf(stderr, "%-25s = %10llu\n", "CORECYCLES", counters->corecycles);
	}
	for (int i = 0; i < lengthof(counters->native); i+=1) {
		uint64_t val = counters->native[i];
		if (val != 0) {
			//const char *label = Bgpm_GetEventIdLabel(i);
			Bgpm_EventInfo_t info;
			BGPM_ERROR(Bgpm_GetEventIdInfo(i, &info));

			fprintf(stderr, "%-25s = %10llu (%s)\n", info.label, val, info.desc);
		}
	}
	fprintf(stderr, "*******************************************************************************\n");
}


mypapi_counters mypapi_stop() {
	assert(omp_get_thread_num() == 0);

	mypapi_counters result;
	result.init = false;
#pragma omp parallel shared(result)
	{
		int tid = Kernel_ProcessorID();
		int cid = Kernel_ProcessorCoreID();
		int sid = Kernel_ProcessorThreadID();
		int i = activeEventSet;
		assert(i >= 0);
		assert(i < MYPAPI_SETS);

		int pues = PuEventSets[i][tid];
		BGPM_ERROR(Bgpm_Stop(pues));
		if (tid == 0) {
			int l2es = L2EventSet[i];
			if (l2es >= 0) {
				BGPM_ERROR(Bgpm_Stop(l2es));
			}
		}

		mypapi_counters local_result = mypapi_bgpm_read(pues, i);
		if (tid == 0) {
			int l2es = L2EventSet[i];
			if (l2es >= 0) {
				mypapi_counters local_result_l2 = mypapi_bgpm_read(l2es, i);
				local_result = mypapi_merge_counters(&local_result, &local_result_l2);
			}
		}
#pragma omp critical
		{
			result = mypapi_merge_counters(&result, &local_result);
		}
	}

	assert(result.init);
	return result;


	long long counters[256];
	result.init = 0;

#pragma omp parallel
{
	long long local_counters[256];
	PAPI_ERROR(PAPI_stop(xEventSet, local_counters));

	mypapi_counters local_result = mypapi_readcounters(xEventSet, counters);
	result = mypapi_merge_counters(&result, &local_result);

#pragma omp master
	{
		if (g_proc_id == 0) {
			Print_PAPI_Counters(xEventSet, local_counters);
		}
	}
}

	long long cyc = PAPI_get_real_cyc() - xCyc;
	long long nsec = PAPI_get_real_nsec() - xNsec;
	double dnow = now2() - xNow;
	double dwtime = MPI_Wtime() - xWtime;
	double domptime = omp_get_wtime() - xOmpTime;
	double sec = (double)nsec * NANO;

	result.secs = dwtime;
	return result;

	long long papiCyc = getMyPapiValue(PAPI_TOT_CYC);
	long long bgpmCyc = getMyPapiNativeValue(PEVT_CYCLES);

		//Print_Counters(xEventSet);


		printf("Cycles: %llu (PAPI_get_real_cyc), %llu (PAPI_TOT_CYC), %llu (PEVT_CYCLES)\n", papiCyc, cyc, bgpmCyc);
		printf("Time: %.5f secs (PAPI_get_real_nsec), %.5f  secs (gettimeofday), %.5f secs (PAPI_TOT_CYC), %.5f secs (PAPI_get_real_cyc), %.5f secs (MPI_WTime), %.5f secs (omp_get_wtime)\n", ((double)nsec)/(1000.0 * 1000.0 * 1000.0), dnow, ((double)papiCyc) / (1600.0 * 1000.0 * 1000.0), ((double)cyc)/(1600.0 * 1000.0 * 1000.0), dwtime, domptime);
		double duration = ((double)nsec)/(1000.0 * 1000.0 * 1000.0);

		int nThreads = omp_get_max_threads();
		int nThreadsOnCore = (nThreads > 16) ? (nThreads / 16) : 1;

		//int xNumEvents = PAPI_num_events(xEventSet);
		//List_PAPI_Events(xEventSet, PAPI_Events, &xNumEvents);
		//long long l1_dch = BGP_UPC_Read_Counter_Value(BGP_PU0_DCACHE_HIT, BGP_UPC_READ_EXCLUSIVE);
		//long long l1_dcm = BGP_UPC_Read_Counter_Value(BGP_PU0_DCACHE_MISS, BGP_UPC_READ_EXCLUSIVE);
		//long long l1_dca = l1_dch + l1_dcm;
		//if (l1_dca > 0)
		//printf("L1 hit rate: %.2f%% (%llu Hits, %llu Misses)\n", 100.0 * l1_dch / l1_dca,  l1_dch,l1_dcm);

		//long long flipAdd = BGP_UPC_Read_Counter_Value(BGP_PU0_FPU_ADD_SUB_1, BGP_UPC_READ_EXCLUSIVE);
		//long long flipMul = BGP_UPC_Read_Counter_Value(BGP_PU0_FPU_MULT_1, BGP_UPC_READ_EXCLUSIVE);
		//long long flipFma = BGP_UPC_Read_Counter_Value(BGP_PU0_FPU_FMA_2, BGP_UPC_READ_EXCLUSIVE);
		//long long flipAddDh = BGP_UPC_Read_Counter_Value(BGP_PU0_FPU_ADD_SUB_2, BGP_UPC_READ_EXCLUSIVE);
		//long long flipMulDh = BGP_UPC_Read_Counter_Value(BGP_PU0_FPU_MULT_2, BGP_UPC_READ_EXCLUSIVE);
		//long long flipFmaDh = BGP_UPC_Read_Counter_Value(BGP_PU0_FPU_FMA_4, BGP_UPC_READ_EXCLUSIVE);
		//long long flop = flipAdd + flipMul + 2* flipFma + 2*(flipAddDh + flipMulDh + 2 * flipFmaDh);
		//if (duration > 0)
		//printf("PAPI %.1f MFlop/s\n", ((double) flop) / (duration * 1000000.0));
		//printf("\n");

			//float rtime;
			//float ptime;
			//long long flops;
			//float mflops;
		//PAPI_ERROR(PAPI_flops(&rtime, &ptime, &flops, &mflops));
		//printf("PAPI repost %f mflops/s", mflops);


		long long nFpOps = getMyPapiNativeValue(PEVT_INST_QFPU_FPGRP2);
		printf("PAPI mflop/s (thread): %f\n", (double)nFpOps / (MEGA * sec));
		printf("PAPI mflop/s (thread,adapted): %f\n", (double)nFpOps / (1.3 * MEGA * sec));
		printf("PAPI mflop/s (node): %f\n", (double)nFpOps * omp_get_max_threads() / (MEGA * sec));
		printf("PAPI mflop/s (node,adapted): %f\n", (double)nFpOps * omp_get_max_threads() / (1.3 * MEGA * sec));

		//long long nFpIns = getMyPapiNativeValue(PEVT_INST_QFPU_FPGRP2_INSTR);
		//printf("PAPI OPI: %f\n", (double)nFpOps / (double)nFpIns);

		//long long nStalls = getMyPapiValue(PAPI_STL_ICY);

		long long nInstrCommit1 = getMyPapiNativeValue(PEVT_INST_ALL);
		//long long nInstrCommit2 = getMyPapiNativeValue(PEVT_XU_PPC_COMMIT);
		printf("PAPI CPI: %f\n", (double)papiCyc/(double)nInstrCommit1);
		//printf("PAPI CPI: %f\n", (double)papiCyc/(double)nInstrCommit2);
		long long nInstrFetchMiss = getMyPapiNativeValue(PEVT_IU_IL1_MISS);
		printf("PAPI Instruction load miss: %f\n", (double)nInstrFetchMiss/(double)nInstrCommit1);

		//long long nResStalls = getMyPapiValue(PAPI_RES_STL);
		//printf("PAPI resource stalls: %f%%\n", 100.0*(double)nResStalls/(double)papiCyc);

		//long long nFpuStalls = getMyPapiValue(PAPI_FP_STAL);
		//printf("PAPI stalled for any AXU/FXU dependency (excludes IS2 stall): %f%%\n", 100.0*(double)nFpuStalls/(double)papiCyc);

		long long nL1IStalls = getMyPapiNativeValue(PEVT_IU_IL1_MISS_CYC);
		printf("PAPI L1I stalls: %f%%\n", 100.0*(double)nL1IStalls/(double)nInstrFetchMiss);
		printf("PAPI Instruction Miss latency: %f cycles\n", (double)nL1IStalls/(double)papiCyc);

		long long nDcbtMisses = getMyPapiNativeValue(PEVT_LSU_COMMIT_DCBT_MISSES);
		long long nDcbtHits = getMyPapiNativeValue(PEVT_LSU_COMMIT_DCBT_HITS);
		printf("PAPI LSU DCBT Hit rate: %f%%\n", 100.0*(double)(nDcbtHits)/(double)(nDcbtMisses+nDcbtHits));

		//long long nCachableLds = getMyPapiNativeValue(PEVT_LSU_COMMIT_CACHEABLE_LDS);
		//long long nNonchachableLds = getMyPapiNativeValue(PEVT_LSU_COMMIT_CACHE_INHIB_LD_MISSES);
		//printf("PAPI LSU %% cachable loads: %f%%\n", 100.0*(double)(nCachableLds)/(double)(nCachableLds+nNonchachableLds));
		long long nLoads = getMyPapiNativeValue(PEVT_LSU_COMMIT_CACHEABLE_LDS);

		//long long nTlbHits = getMyPapiNativeValue(PEVT_MMU_TLB_HIT_DIRECT_DERAT);
		//long long nTlbMiss = getMyPapiNativeValue(PEVT_MMU_TLB_MISS_DIRECT_DERAT);
		//printf("PAPI MMU %% TLB hits: %f%%\n", 100.0*(double)(nTlbHits)/(double)(nTlbHits+nTlbMiss));

		//long long nStores = getMyPapiNativeValue(PEVT_LSU_COMMIT_STS);
		//long long nStoreMisses = getMyPapiNativeValue(PEVT_LSU_COMMIT_ST_MISSES);
		//printf("PAPI LSU %% store hits: %f%%\n", 100.0*(double)(nStores-nStoreMisses)/(double)(nStores));


		long long nFUDepStalls = getMyPapiNativeValue(PEVT_IU_IS1_STALL_CYC);
		long long nFUStalls = getMyPapiNativeValue(PEVT_IU_IS2_STALL_CYC);
		//long long nRAWStalls= getMyPapiNativeValue(PEVT_IU_RAW_DEP_HIT_CYC);
		//long long nWAWStalls = getMyPapiNativeValue(PEVT_IU_WAW_DEP_HIT_CYC);
		printf("PAPI IU Dep Stalls: %f%%\n", 100.0*(double)(nFUDepStalls-nFUStalls)/(double)(papiCyc));
		printf("PAPI IU FU Stalls: %f%%\n", 100.0*(double)(nFUStalls)/(double)(papiCyc));
		//printf("PAPI IU Read-after-Write Stalls: %f%%\n", 100.0*(double)(nRAWStalls)/(double)(papiCyc));
		//printf("PAPI IU Write-after-Write Stalls: %f%%\n", 100.0*(double)(nWAWStalls)/(double)(papiCyc));

		long long nBrInstr = getMyPapiNativeValue(PEVT_XU_BR_COMMIT);
		long long nBrCondMispredict = getMyPapiNativeValue(PEVT_XU_BR_MISPRED_COMMIT);
		//long long nBrTargetMispredict = getMyPapiNativeValue(PEVT_XU_BR_TARG_ADDR_MISPRED_COMMIT);
		printf("PAPI XU Conditional branch mispredict: %f%%\n", 100.0*(double)(nBrCondMispredict)/(double)(nBrInstr));
		//printf("PAPI XU Branch target mispredict: %f%%\n", 100.0*(double)(nBrTargetMispredict)/(double)(nBrInstr));


		long long nPrefetchLoads = getMyPapiNativeValue(PEVT_L1P_BAS_LD);
		long long nPrefetchHits = getMyPapiNativeValue(PEVT_L1P_BAS_HIT);
		long long nPrefetchMisses = getMyPapiNativeValue(PEVT_L1P_BAS_MISS);
		printf("PAPI L1P Prefetches to Loads: %f%%\n", 100.0*(double)(nPrefetchLoads)/(double)(nLoads));
		printf("PAPI L1P hit: %f%%\n", 100.0*(double)(nPrefetchHits)/(double)(nPrefetchLoads));
		printf("PAPI L1P miss: %f%%\n", 100.0*(double)(nPrefetchMisses)/(double)(nPrefetchLoads));
		// For some reason hit+miss>100%

		long long nPrefetchStreamEstb = getMyPapiNativeValue(PEVT_L1P_STRM_LINE_ESTB);
		long long nPrefetchStreamUnused = getMyPapiNativeValue(PEVT_L1P_STRM_EVICT_UNUSED);
		printf("PAPI L1P Stream uselessly prefetched lines: %f%%\n", 100.0*(double)(nPrefetchStreamUnused)/(double)(nPrefetchStreamEstb));

		long long nL2Hits = getMyPapiNativeValue(PEVT_L2_HITS); // per node
		long long nL2Misses = getMyPapiNativeValue(PEVT_L2_MISSES); // per node
		printf("PAPI L2 Hit rate: %f%%\n", 100.0*(double)(nL2Hits)/(double)(nL2Hits+nL2Misses));

		long long nL2LinesFetched = getMyPapiNativeValue(PEVT_L2_FETCH_LINE);
		long long nL2LinesStores = getMyPapiNativeValue(PEVT_L2_STORE_LINE);
		printf("PAPI Main memory read bandwidth: %f mbyte/s\n", 128.0 * nL2LinesFetched / (1.0 * sec * 1024 * 1024));
		printf("PAPI Main memory write bandwidth: %f mbyte/s\n", 128.0 * nL2LinesStores / (1.0 * sec * 1024 * 1024));

		long long nLoadL1Misses = getMyPapiNativeValue(PEVT_LSU_COMMIT_LD_MISSES);
		long long nLoadInstrs = getMyPapiNativeValue(PEVT_INST_XU_LD);
		printf("PAPI LSU L1 Hit ratio v1: %f%%\n", 100.0 * (nLoads - nLoadL1Misses) / (double)nLoads);
		printf("PAPI LSU L1 Hit ratio v2: %f%%\n", 100.0 * (nLoadInstrs - nLoadL1Misses) / (double)nLoadInstrs);
}





/* Print_Counters */
void Print_Counters(const int pEventSet) {
	//printf("\n***** Start Print Counter Values *****\n");
	// Print_Native_Counters_via_Buffer((BGP_UPC_Read_Counters_Struct_t*)Native_Buffer);
	//Print_Native_Counters();
	//Print_Native_Counters_for_PAPI_Counters(pEventSet);
	Print_PAPI_Counters(pEventSet, PAPI_Counters);
	//printf("\n***** End Print Counter Values *****\n");
	return;
}


/* Print_PAPI_Counters */
static void Print_PAPI_Counters(const int pEventSet, const long long* pCounters) {
	int i;
	char xName[256];
	printf("***** Start Print of PAPI Counter Values *****\n");
	//printf("Print_PAPI_Counters: PAPI_Counters*=%p, pCounters*=%p\n",PAPI_Counters, pCounters);

	int pNumEvents = PAPI_num_events(pEventSet);
	printf("Number of Counters = %d\n", pNumEvents);
	if (pNumEvents) {
		printf("Calculated Value      Location  Number   Event Name\n");
		printf("--------------------- --------- ---------- ----------------------------------------------\n");
		List_PAPI_Events(pEventSet, PAPI_Events, &pNumEvents);
		for (i=0; i<pNumEvents; i++) {
			if (PAPI_event_code_to_name(PAPI_Events[i], xName)) {
				printf("PAPI_event_code_to_name failed on event code %d\n",PAPI_Events[i]);
				exit(1);
			}
 			printf("%20llu	%3d	0x%8.8x %s\n", pCounters[i], i, PAPI_Events[i],xName);
		}
	}
	printf("***** End Print of PAPI Counter Values *****\n");
	return;
}


/* Print_PAPI_Counters_From_List */
static void Print_PAPI_Counters_From_List(const int* pEventList, const int pNumEvents, const long long* pCounters) {
	int i;
	char xName[256];

	printf("\n***** Start Print of PAPI Counter Values *****\n");
	printf("Number of Counters = %d\n", pNumEvents);

	if (pNumEvents) {
		printf("Calculated Value      Number    Location   Event Name\n");
		printf("--------------------- --------- ---------- ----------------------------------------------\n");
		for (i=0; i<pNumEvents; i++) {
			if (PAPI_event_code_to_name(pEventList[i], xName)) {
				printf("PAPI_event_code_to_name failed on event code %d\n",pEventList[i]);
				exit(1);
			}
 			printf("%20ll	u%3d0x	%8.8x %s\n", pCounters[i], i, pEventList[i], xName);
		}
	}
	printf("***** End Print of PAPI Counter Values *****\n");
	return;
}


/* Print_PAPI_Events */
static void Print_PAPI_Events(const int pEventSet) {
	int i;
	char xName[256];
	int pNumEvents = PAPI_num_events(pEventSet);
	List_PAPI_Events(pEventSet, PAPI_Events, &pNumEvents);

	for (i=0; i<pNumEvents; i++) {
		if (!PAPI_event_code_to_name(PAPI_Events[i], xName))
			printf("PAPI Counter Location %3.3d: 0x%8.8x %s\n", i, PAPI_Events[i],xName);
		else
			printf("PAPI Counter Location %3.3d: Not mapped\n", i);
	}
	return;
}


/* List_PAPI_Events */
static void List_PAPI_Events(const int pEventSet, int* pEvents, int* pNumEvents) {
	int xRC = PAPI_list_events(pEventSet, pEvents, pNumEvents);
	if (xRC != PAPI_OK) {
		printf("FAILURE: PAPI_list_events failed, returned xRC=%d...\n", xRC);
		exit(1);
	}
	return;
}

#else

void mypapi_init(){}
void mypapi_start(int i){}
mypapi_counters mypapi_stop(){}

#endif

