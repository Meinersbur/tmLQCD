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

#define PAPI_ERROR(cmd) if ((RC = (cmd))!=PAPI_OK) { fprintf(stderr, "MK_PAPI call failed with code %d at line %d: %s\n", RC, __LINE__, TOSTRING(cmd));  }

#if 0
static void printNativeEvents() {
	int retval, EventSet = PAPI_NULL;
	unsigned int native = 0x0;
	PAPI_event_info_t info;

	/* Initialize the library */
	//retval = PAPI_library_init(PAPI_VER_CURRENT);
	//if (retval != PAPI_VER_CURRENT) {
	//	printf("PAPI library init error!\n");
	////	exit(1);
	}

	/* Create an EventSet */
	int retval = PAPI_create_eventset(&EventSet);
	if (retval != PAPI_OK) handle_error(retval);

	/* Find the first available native event */
	native = PAPI_NATIVE_MASK | 0;
	retval = PAPI_enum_event(&native, PAPI_ENUM_FIRST);
	if (retval != PAPI_OK) handle_error(retval);

	/* Add it to the eventset */
	retval = PAPI_add_event(EventSet, native);
	if (retval != PAPI_OK) handle_error(retval);

	/* Exit successfully */
	printf("Success!\n");
}
#endif

static int PAPI_add_native_event(int EventSet, int EventCode) {
	// For some strange reason (i.e. unknown to me), the numbers from events.h are off by one
	return PAPI_add_event(EventSet, PAPI_NATIVE_MASK | (EventCode-1));
}

void mypapi_init() {
	if (g_proc_id == 0)
		fprintf(stderr, "MK_Init mypapi\n");

	int i;
	int RC;
	int retval = PAPI_library_init(PAPI_VER_CURRENT);
	if (retval != PAPI_VER_CURRENT) {
		fprintf(stderr, "Unexpected PAPI version\n");
		exit(1);
	}

	if (g_proc_id != 0)
		return;

	for (i = 0; i < 255; i++)
		PAPI_Counters[i]=0;

 		const PAPI_hw_info_t *hwinfo = PAPI_get_hardware_info();
		printf("CPU frequency: %f Mhz\n", hwinfo->mhz);
		printf("TOT_CYC frequency: %f Mhz\n", hwinfo->clock_mhz);
		printf("CPU: %s %s\n", hwinfo->vendor_string, hwinfo->model_string);

	//BGP_UPC_Initialize_Counter_Config(BGP_UPC_MODE_DEFAULT, BGP_UPC_CFG_EDGE_DEFAULT);
	int counter = PAPI_num_counters();
	printf("Available hardware counters: %d\n", counter);
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

void mypapi_start() {
	int RC;
	if (g_proc_id != 0)
		return;
	if (omp_get_thread_num() != 0)
		return;

	for (int i = 0; i < 255; i++)
		PAPI_Counters[i]=0;
	PAPI_ERROR(PAPI_reset(xEventSet));
	//fprintf(stderr, "MK_MyPAPI start\n");
	xCyc = PAPI_get_real_cyc();
	xNsec = PAPI_get_real_nsec();
	xNow = now2();
	xWtime = MPI_Wtime();
	xOmpTime =  omp_get_wtime();
	PAPI_ERROR(PAPI_start(xEventSet));
}


void mypapi_stop() {
	if (g_proc_id != 0)
		return;
	if (omp_get_thread_num() != 0)
		return;

	int i;
	for (i = 0; i < 255; i++)
		PAPI_Counters[i]=0;

	int RC;
	RC = PAPI_stop(xEventSet, PAPI_Counters);
	RC = PAPI_read(xEventSet, PAPI_Counters);

	if (g_proc_id == 0) {
		long long cyc = PAPI_get_real_cyc() - xCyc;
		long long nsec = PAPI_get_real_nsec() - xNsec;
		double dnow = now2() - xNow;
		double dwtime = MPI_Wtime() - xWtime;
		double domptime = omp_get_wtime() - xOmpTime;
		double sec = (double)nsec * NANO;

		Print_Counters(xEventSet);

		long long papiCyc = getMyPapiValue(PAPI_TOT_CYC);
		long long bgpmCyc = getMyPapiNativeValue(PEVT_CYCLES);
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

/* Print_Native_Counters */
//void Print_Native_Counters() {
//	printf("\n***** Start Print of Native Counter Values *****\n");
//	BGP_UPC_Print_Counter_Values(BGP_UPC_READ_EXCLUSIVE);
//	printf("***** End Print of Native Counter Values *****\n");
//	return;
//}

/* Print_Native_Counters_for_PAPI_Counters */
#if 0
void Print_Native_Counters_for_PAPI_Counters(const int pEventSet) {
	printf("\n***** Start Print of Native Counter Values for PAPI Counters*****\n");
	int xNumEvents = PAPI_num_events(pEventSet);
	if (xNumEvents) {
		List_PAPI_Events(pEventSet, PAPI_Events, &xNumEvents);
		Print_Native_Counters_for_PAPI_Counters_From_List(PAPI_Events, xNumEvents);
	} else {
		printf("No events are present in the event set.\n");
	}
	printf("***** End Print of Native Counter Values for PAPI Counters *****\n");
	return;
}

/* Print_Native_Counters_for_PAPI_Counters_From_List */
void Print_Native_Counters_for_PAPI_Counters_From_List(const int* pEvents, const int pNumEvents) {
	int i, j, xRC;
	char xName[256];
	BGP_UPC_Event_Id_t xNativeEventId;
	PAPI_event_info_t xEventInfo;
	//BGP_UPC_Print_Counter_Values(); // DLH
	for (i=0; i<pNumEvents; i++) {
		xRC = PAPI_event_code_to_name(PAPI_Events[i], xName);
		if (!xRC) {
			xRC = PAPI_get_event_info(PAPI_Events[i], &xEventInfo);
			if (xRC) {
				printf("FAILURE: PAPI_get_event_info failed for %s, xRC=%d\n", xName,xRC);
				exit(1);
			}
			printf("\n*** PAPI Counter Location %3.3d: 0x%8.8x %s\n", i,PAPI_Events[i], xName);
			if (PAPI_Events[i] & 0x80000000) {
				// Preset event
				for (j=0; j<xEventInfo.count; j++) {
					xNativeEventId = (BGP_UPC_Event_Id_t)(xEventInfo.code[j]&0xBFFFFFFF);
					//printf("Preset: j=%d, xEventInfo.code[j]=0x%8.8x,xNativeEventId=0x%8.8x\n", j, xEventInfo.code[j], xNativeEventId);
					BGP_UPC_Print_Counter_Value(xNativeEventId, BGP_UPC_READ_EXCLUSIVE);
				}
			} else {
				// Native event
				xNativeEventId = (BGP_UPC_Event_Id_t)(PAPI_Events[i]&0xBFFFFFFF);
				//printf("Native: i=%d, PAPI_Events[i]=0x%8.8x,xNativeEventId=0x%8.8x\n", i, PAPI_Events[i], xNativeEventId);
				BGP_UPC_Print_Counter_Value(xNativeEventId, BGP_UPC_READ_EXCLUSIVE);
			}
		} else {
			printf("\n*** PAPI Counter Location %3.3d: Not mapped\n", i);
		}
	}
}
#endif

/* Print_PAPI_Counters */
void Print_PAPI_Counters(const int pEventSet, const long long* pCounters) {
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
void Print_PAPI_Counters_From_List(const int* pEventList, const int pNumEvents, const long long* pCounters) {
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
void Print_PAPI_Events(const int pEventSet) {
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
void List_PAPI_Events(const int pEventSet, int* pEvents, int* pNumEvents) {
	int xRC = PAPI_list_events(pEventSet, pEvents, pNumEvents);
	if (xRC != PAPI_OK) {
		printf("FAILURE: PAPI_list_events failed, returned xRC=%d...\n", xRC);
		exit(1);
	}
	return;
}

#else

void mypapi_init(){}
void mypapi_start(){}
void mypapi_stop(){}

#endif

