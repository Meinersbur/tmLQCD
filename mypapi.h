#ifndef MYPAPI
#define MYPAPI

#ifndef STRINGIFY
#define STRINGIFY(x) #x
#endif

#ifndef TOSTRING
#define TOSTRING(x) STRINGIFY(x)
#endif 

#define __Pragma(x) _Pragma(#x)

void mypapi_init();
void mypapi_start();
void mypapi_stop();

#define NANO  (1e-9)
#define MICRO (1e-6)
#define MILLI (1e-3)

#define KILO (1e3)
#define MEGA (1e6)
#define GIGA (1e9)
#define TERA (1e12)



#endif
