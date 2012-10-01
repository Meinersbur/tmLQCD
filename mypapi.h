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
	int threadid;
	int coreid;
	int smtid;
	long long preset[PAPI_END_idx];
	long long native[UPCI_NUM_EVENTS];
	double secs;
} mypapi_counters;

typedef enum {
	pi_hitinl1,
	pi_hitinl1p,
	pi_hitinl2,
	pi_hitinmain,
	__pi_COUNT
} mypapi_interpretations;

#define MYPAPI_SETS 1
void mypapi_init();
void mypapi_start(int i);
mypapi_counters mypapi_stop();

mypapi_counters mypapi_merge_counters(mypapi_counters *counters1, mypapi_counters *counters2);

#define NANO  (1e-9)
#define MICRO (1e-6)
#define MILLI (1e-3)

#define KILO (1e3)
#define MEGA (1e6)
#define GIGA (1e9)
#define TERA (1e12)


#endif
