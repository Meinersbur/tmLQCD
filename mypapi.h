#ifndef MYPAPI
#define MYPAPI

#include <stdbool.h>
#include <papi.h>
#include <upci/events.h>

#ifndef STRINGIFY
#define STRINGIFY(x) #x
#endif

#ifndef TOSTRING
#define TOSTRING(x) STRINGIFY(x)
#endif 

#define __Pragma(x) _Pragma(#x)

typedef struct {
	bool init;
	int set;
	int eventset;
	int threadid;
	int coreid;
	int smtid;
	int ompid;
	//long long preset[PAPI_END_idx];
	uint64_t native[UPCI_NUM_EVENTS];
	uint64_t corecycles;
	uint64_t nodecycles;
	double secs;
} mypapi_counters;

typedef enum {
	pi_cpi,
	pi_corecpi,
	pi_hitinl1,

	pi_l1phitrate,
	pi_l2hitrate,
	pi_dcbthitrate,
	__pi_COUNT,

	pi_hitinl1p
} mypapi_interpretations;

#define MYPAPI_SETS 2
void mypapi_init();
void mypapi_start(int i);
mypapi_counters mypapi_stop();

mypapi_counters mypapi_merge_counters(mypapi_counters *counters1, mypapi_counters *counters2);
void mypapi_print_counters(mypapi_counters *counters);

#define NANO  (1e-9)
#define MICRO (1e-6)
#define MILLI (1e-3)

#define KILO (1e3)
#define MEGA (1e6)
#define GIGA (1e9)
#define TERA (1e12)


#endif
