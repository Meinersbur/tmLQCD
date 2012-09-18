#include "global.h"
#include <stdio.h>
#include <string.h>
#include <mcheck.h>

#undef malloc
#undef posix_memalign

static long long total = 0;
void *nalloc(size_t size, const char *term, const char *file, const char* func, int line) {
    total += size;
    if (g_proc_id == 0)
        fprintf(stderr, "malloc: %s = %lu bytes, %.2f MB total (%s at %s:%d)\n", term, size, 1.0 * total / (1024*1024), func, file, line);
    void* result = malloc(size);
    if (result) {
    	//memset(result, 0, size);
    } else if (g_proc_id == 0)
    	fprintf(stderr, "malloc: returned NULL\n");
    return result;
}


int nposix_memalign(void **memptr, size_t alignment, size_t size, const char *term, const char *file, const char* func, int line) {
    total += size;
    if (g_proc_id == 0)
        fprintf(stderr, "posix_memalign: %s = %lu bytes, %.2f MB total (%s at %s:%d)\n", term, size, 1.0 * total / (1024*1024), func, file, line);
    int result = posix_memalign(memptr, alignment, size);
    if (result && (g_proc_id == 0))
    	fprintf(stderr, "posix_memalign: status code %d\n", result);
    if (*memptr) {
    	//memset(*memptr, 0, size);
    } else if (g_proc_id == 0)
    	fprintf(stderr, "posix_memalign: returned NULL\n");
    return result;
}
