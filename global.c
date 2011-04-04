#include "global.h"
#include <stdio.h>
#include <string.h>
#include <mcheck.h>

#undef malloc

static int total = 0;
void *nalloc(size_t size, const char *term, const char *file, const char* func, int line) {
    total += size;
    if (g_proc_id == 0)
        fprintf(stderr, "malloc: %s = %d bytes, %.2f MB total (%s at %s:%d)\n", term, size, 1.0 * total / (1024*1024), func, file, line);
    void* result = malloc(size);
    memset(result, 0, size);
    return result;
}

