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

#endif
