
#include <mpi.h>

int g_proc_id;
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)


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

int PAPI_Events[256];
long long PAPI_Counters[256];
int xEventSet=PAPI_NULL;

#define PAPI_ERROR(cmd) if (RC = (cmd)) { fprintf(stderr, "MK_PAPI call failed with code %d at line %d: %s\n", RC, __LINE__, TOSTRING(cmd));  }

void mypapi_init() {
        PAPI_library_init(PAPI_VER_CURRENT);
        if (g_proc_id != 0)
                return;
        fprintf(stderr, "MK_Init mypapi\n");

        int i;
        int RC;

        for (i = 0; i < 255; i++)
                PAPI_Counters[i]=0;

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
        PAPI_ERROR(PAPI_add_event(xEventSet, PAPI_L2_DCR));
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

        fprintf(stderr, "MK_Init mypapi END\n");
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

        if (g_proc_id == 0)
                Print_Counters(xEventSet);
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
                printf("Calculated Value Location Event Number Event Name\n");
                printf("-------------------- -------- --------------------------------------------------------\n");
                List_PAPI_Events(pEventSet, PAPI_Events, &pNumEvents);
                for (i=0; i<pNumEvents; i++) {
                        if (PAPI_event_code_to_name(PAPI_Events[i], xName)) {
                                printf("PAPI_event_code_to_name failed on event code %d\n",PAPI_Events[i]);
                                exit(1);
                        }
                        printf("%20llu  %3d     0x%8.8x %s\n", pCounters[i], i, PAPI_Events[i],xName);
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
                printf("Calculated Value Location Event Number Event Name\n");
                printf("-------------------- -------- --------------------------------------------------------\n");
                for (i=0; i<pNumEvents; i++) {
                        if (PAPI_event_code_to_name(pEventList[i], xName)) {
                                printf("PAPI_event_code_to_name failed on event code %d\n",pEventList[i]);
                                exit(1);
                        }
                        printf("%20ll   u%3d0x  %8.8x %s\n", pCounters[i], i, pEventList[i], xName);
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




void main(int argc,char *argv[]) {
	  MPI_Init(&argc, &argv);
  /* GG */
  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);
  mypapi_init();

}


