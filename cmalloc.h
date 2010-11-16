/* GG */
#include <errno.h>
#define CMALLOC_ERROR_EXIT(p) {if ( p==NULL ) {printf ("calloc/malloc errno: %d at %s:%d \n", errno, __FILE__, __LINE__); fflush(stdout); errno = 0; exit(69);}}
