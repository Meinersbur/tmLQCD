#include "global.h"
#include <stdio.h>
#include <string.h>
#include <mcheck.h>

#undef malloc
void *nalloc(size_t size, const char *term, const char *file, const char* func, int line) {
    if (g_proc_id == 0)
    fprintf(stderr, "malloc: %s = %d bytes (%s at %s:%d)\n", term, size, func, file, line);
    void* result = malloc(size);
    memset(result, 0, size);
    return result;
}

