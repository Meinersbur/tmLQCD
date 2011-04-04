#include "mypapi.h"


#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <common/alignment.h>
#include "papi.h" //Provides definitions for base PAPI function
#include <spi/bgp_SPI.h> //The Blue Gene/P substrate for PAPI is based on the UPC interfaces documented in this header file.
#include "papiStdEventDefs.h" //Provides a listing of all standard PAPI presets. Please note that not all of these are available on the Blue Gene/P platform
#include <linux-bgp-native-events.h> //Provides a listing of all available counters native to Blue Gene.


void FPUArith(void); //This method does various calculations which should saturate many of the counters
void List_PAPI_Events(const int pEventSet, int* pEvents, int* xNumEvents);
void Print_Native_Counters();
void Print_Native_Counters_via_Buffer(const BGP_UPC_Read_Counters_Struct_t* pBuffer);
void Print_Native_Counters_for_PAPI_Counters(const int pEventSet);
void Print_Native_Counters_for_PAPI_Counters_From_List(const int* pEvents, const int pNumEvents);
void Print_PAPI_Counters(const int pEventSet, const long long* pCounters);
void Print_PAPI_Counters_From_List(const int* pEventList, const int pNumEvents, const long long* pCounters);
void Print_Counters(const int pEventSet);
void Print_PAPI_Events(const int pEventSet);
long long getMyPapiValue(const int eventNum);

int PAPI_Events[256];
long long PAPI_Counters[256];
int xEventSet=PAPI_NULL;

long long xCyc;
long long xNsec;
double xNow;

extern int g_proc_id;

static double now2(){
   struct timeval t; double f_t;
   gettimeofday(&t, NULL);
   f_t = t.tv_usec; f_t = f_t/1000000.0; f_t +=t.tv_sec;
   return f_t;
}

#define PAPI_ERROR(cmd) if (RC = (cmd)) { fprintf(stderr, "MK_PAPI call failed with code %d at line %d: %s\n", RC, __LINE__, TOSTRING(cmd));  }

void mypapi_init() {
	if (g_proc_id == 0)
		fprintf(stderr, "MK_Init mypapi\n");

	int i;
	int RC;
	PAPI_library_init(PAPI_VER_CURRENT);
	if (g_proc_id != 0)
		return;

	for (i = 0; i < 255; i++)
		PAPI_Counters[i]=0;

 		const PAPI_hw_info_t *hwinfo = PAPI_get_hardware_info();
		printf("CPU frequency: %f Mhz\n", hwinfo->mhz);
		printf("TOT_CYC frequency: %f Mhz\n", hwinfo->clock_mhz);
		printf("CPU: %s %s\n", hwinfo->vendor_string, hwinfo->model_string);

	BGP_UPC_Initialize_Counter_Config(BGP_UPC_MODE_DEFAULT, BGP_UPC_CFG_EDGE_DEFAULT);
	RC = PAPI_create_eventset(&xEventSet);

	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_TOT_CYC));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_FP_INS));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_FP_OPS));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_FMA_INS));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_LD_INS));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_LST_INS));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_SR_INS));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_L1_DCR));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_L1_DCM));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_L1_ICM));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_L1_TCM));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_L1_DCR));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_L1_ICA));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_L1_ICR));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_L1_ICW));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_L1_TCH));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_L1_TCA));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_L1_TCW));
	//PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_L2_DCR));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_L2_DCW));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_L2_DCM));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_L2_DCH));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_L2_DCA));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_L3_LDM));
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_PRF_DM)); // Prefetch misses
	PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_CA_SNP)); // # of cycles with idle FPU units
	//PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_BGP_TS_DPKT));
	//PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_BGP_TS_32B));

	PAPI_ERROR(PAPI_add_event(xEventSet, PNE_BGP_MISC_ELAPSED_TIME));
	PAPI_ERROR(PAPI_add_event(xEventSet, PNE_BGP_PU0_DCACHE_LINEFILLINPROG));
	//PAPI_ERROR(PAPI_add_event(xEventSet, PNE_BGP_PU1_DCACHE_LINEFILLINPROG));
	PAPI_ERROR(PAPI_add_event(xEventSet, PNE_BGP_PU0_L2_VALID_PREFETCH_REQUESTS));
	PAPI_ERROR(PAPI_add_event(xEventSet, PNE_BGP_PU0_L2_PREFETCH_HITS_IN_FILTER));
	PAPI_ERROR(PAPI_add_event(xEventSet, PNE_BGP_PU0_L2_PREFETCH_HITS_IN_STREAM));
	PAPI_ERROR(PAPI_add_event(xEventSet, PNE_BGP_PU0_L2_CYCLES_PREFETCH_PENDING));
	PAPI_ERROR(PAPI_add_event(xEventSet, PNE_BGP_PU0_L2_PREFETCHABLE_REQUESTS));
	PAPI_ERROR(PAPI_add_event(xEventSet, PNE_BGP_PU0_DCACHE_HIT));
	//PAPI_ERROR(PAPI_add_event(xEventSet, PNE_BGP_PU1_DCACHE_HIT));
}

void mypapi_start() {
	int RC;
	if (g_proc_id != 0)
		return;

	for (int i = 0; i < 255; i++)
		PAPI_Counters[i]=0;
	PAPI_ERROR(PAPI_reset(xEventSet));
	fprintf(stderr, "MK_MyPAPI start\n");
	xCyc = PAPI_get_real_cyc();
	xNsec = PAPI_get_real_nsec();
	xNow = now2();
	PAPI_ERROR(PAPI_start(xEventSet));
}


/* Print_Counters */
void Print_Counters(const int pEventSet) {
	printf("\n***** Start Print Counter Values *****\n");
	// Print_Native_Counters_via_Buffer((BGP_UPC_Read_Counters_Struct_t*)Native_Buffer);
	Print_Native_Counters();
	Print_Native_Counters_for_PAPI_Counters(pEventSet);
	Print_PAPI_Counters(pEventSet, PAPI_Counters);
	printf("\n***** End Print Counter Values *****\n");
	return;
}


void mypapi_stop() {
	if (g_proc_id != 0)
		return;

	int i;
	for (i = 0; i < 255; i++)
		PAPI_Counters[i]=0;

	int RC;
	RC = PAPI_stop(xEventSet, PAPI_Counters);
	RC = PAPI_read(xEventSet, PAPI_Counters);
	long long cyc = PAPI_get_real_cyc() - xCyc;
	long long nsec = PAPI_get_real_nsec() - xNsec;
	double dnow = now2() - xNow;

	if (g_proc_id == 0) {
		Print_Counters(xEventSet);

		long long papiCyc = getMyPapiValue(PAPI_TOT_CYC);
		printf("Cycles: %llu (PAPI_get_real_cyc), %llu (PAPI_TOT_CYC)\n", papiCyc, cyc);
		printf("Time: %.5f secs (PAPI_get_real_nsec), %.5f  secs (gettimeofday), %.5f secs (PAPI_TOT_CYC), %.5f secs (PAPI_get_real_cyc)\n", ((double)nsec)/(1000.0 * 1000.0 * 1000.0), dnow, ((double)papiCyc) / (850.0 * 1000.0 * 1000.0), ((double)cyc)/(850.0 * 1000.0 * 1000.0));
		double duration = ((double)nsec)/(1000.0 * 1000.0 * 1000.0);

		//int xNumEvents = PAPI_num_events(xEventSet);
		//List_PAPI_Events(xEventSet, PAPI_Events, &xNumEvents);
		long long l1_dch = BGP_UPC_Read_Counter_Value(BGP_PU0_DCACHE_HIT, BGP_UPC_READ_EXCLUSIVE);
		long long l1_dcm = BGP_UPC_Read_Counter_Value(BGP_PU0_DCACHE_MISS, BGP_UPC_READ_EXCLUSIVE);
		long long l1_dca = l1_dch + l1_dcm;
		if (l1_dca > 0)
		printf("L1 hit rate: %.2f%% (%llu Hits, %llu Misses)\n", 100.0 * l1_dch / l1_dca,  l1_dch,l1_dcm);

		long long flipAdd = BGP_UPC_Read_Counter_Value(BGP_PU0_FPU_ADD_SUB_1, BGP_UPC_READ_EXCLUSIVE);
		long long flipMul = BGP_UPC_Read_Counter_Value(BGP_PU0_FPU_MULT_1, BGP_UPC_READ_EXCLUSIVE);
		long long flipFma = BGP_UPC_Read_Counter_Value(BGP_PU0_FPU_FMA_2, BGP_UPC_READ_EXCLUSIVE);
		long long flipAddDh = BGP_UPC_Read_Counter_Value(BGP_PU0_FPU_ADD_SUB_2, BGP_UPC_READ_EXCLUSIVE);
		long long flipMulDh = BGP_UPC_Read_Counter_Value(BGP_PU0_FPU_MULT_2, BGP_UPC_READ_EXCLUSIVE);
		long long flipFmaDh = BGP_UPC_Read_Counter_Value(BGP_PU0_FPU_FMA_4, BGP_UPC_READ_EXCLUSIVE);
		long long flop = flipAdd + flipMul + 2* flipFma + 2*(flipAddDh + flipMulDh + 2 * flipFmaDh);
		if (duration > 0)
		printf("PAPI %.1f MFlop/s\n", ((double) flop) / (duration * 1000000.0));
		printf("\n");
	}
}

long long getMyPapiValue(const int eventNum) {
	int i;
	char xName[256];
	//printf("Print_PAPI_Counters: PAPI_Counters*=%p, pCounters*=%p\n",PAPI_Counters, pCounters);

	int pNumEvents = PAPI_num_events(xEventSet);
	if (pNumEvents) {
		List_PAPI_Events(xEventSet, PAPI_Events, &pNumEvents);
		for (i=0; i<pNumEvents; i++) {
			if (PAPI_event_code_to_name(PAPI_Events[i], xName)) {
				printf("PAPI_event_code_to_name failed on event code %d\n",PAPI_Events[i]);
				exit(1);
			}
 			//printf("%20llu	%3d	0x%8.8x %s\n", PAPI_Counters[i], i, PAPI_Events[i],xName);
			if ((eventNum  & 0xBFFFFFFF) == (PAPI_Events[i] & 0xBFFFFFFF))
				return PAPI_Counters[i];
		}
	}

	return 0;
}



/* Print_Counters */
void Print_Counters(const int pEventSet) {
	printf("\n***** Start Print Counter Values *****\n");
	// Print_Native_Counters_via_Buffer((BGP_UPC_Read_Counters_Struct_t*)Native_Buffer);
	//Print_Native_Counters();
	Print_Native_Counters_for_PAPI_Counters(pEventSet);
	Print_PAPI_Counters(pEventSet, PAPI_Counters);
	printf("\n***** End Print Counter Values *****\n");
	return;
}

/* Print_Native_Counters */
void Print_Native_Counters() {
	printf("\n***** Start Print of Native Counter Values *****\n");
	BGP_UPC_Print_Counter_Values(BGP_UPC_READ_EXCLUSIVE);
	printf("***** End Print of Native Counter Values *****\n");
	return;
}

/* Print_Native_Counters_for_PAPI_Counters */
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

/* Print_PAPI_Counters */
void Print_PAPI_Counters(const int pEventSet, const long long* pCounters) {
	int i;
	char xName[256];
	printf("\n***** Start Print of PAPI Counter Values *****\n");
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



