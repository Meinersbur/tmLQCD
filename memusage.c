#include "memusage.h"

#include <sys/resource.h>
#include <common/bgp_personality.h>
#include <common/bgp_personality_inlines.h>
#include <spi/kernel_interface.h>
#include <spi/bgp_SPI.h>

#include "global.h"

/* returns memory per core in MBytes */
unsigned bg_coreMB()
{
    unsigned procMB, coreMB;
    static _BGP_Personality_t mybgp;
    Kernel_GetPersonality(&mybgp, sizeof(_BGP_Personality_t));
    procMB = BGP_Personality_DDRSizeMB(&mybgp);
    coreMB = procMB/Kernel_ProcessCount();
    return coreMB;
}

/* return maximum memory usage of process in kBytes */
unsigned bg_usedKB()
{
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) != 0)
        return 0;
    return usage.ru_maxrss;
}

void print_memusage()
{
    if (g_proc_id != 0)
        return;

    {
        unsigned memPerCore = bg_coreMB();
        fprintf(stderr, "MK_Memory available per Core:  %10ld MB\n", memPerCore);
    }

    struct rusage usage;
    memset(&usage, 0, sizeof(usage));
    int returncode = 0;
    if (returncode = getrusage(RUSAGE_SELF, &usage))
    {
        fprintf(stderr, "MK_getrusage error; return code %d\n", returncode);
    }
    else
    {
        fprintf(stderr, "MK_User time used:             %d,%d secs\n", usage.ru_utime.tv_sec, usage.ru_utime.tv_usec);
        fprintf(stderr, "MK_System time used:           %d,%d secs\n", usage.ru_stime.tv_sec, usage.ru_stime.tv_usec);
        fprintf(stderr, "MK_Maximum resident set size:  %10ld kB\n", usage.ru_maxrss);
        fprintf(stderr, "MK_Integral shared memory size:%10ld kB\n", usage.ru_ixrss);
        fprintf(stderr, "MK_Integral unshared data size:%10ld kB\n", usage.ru_idrss);
        fprintf(stderr, "MK_Integral unshared stack s.: %10ld kB\n", usage.ru_isrss);
        fprintf(stderr, "MK_Page reclaims:              %10ld\n", usage.ru_minflt);
        fprintf(stderr, "MK_Page faults:                %10ld\n", usage.ru_majflt);
        fprintf(stderr, "MK_Swaps:                      %10ld\n", usage.ru_nswap);
        fprintf(stderr, "MK_Block input operations:     %10ld\n", usage.ru_inblock);
        fprintf(stderr, "MK_Block output operations:    %10ld\n", usage.ru_oublock);
        fprintf(stderr, "MK_Messages sent:              %10ld\n", usage.ru_msgsnd);
        fprintf(stderr, "MK_Messages received:          %10ld\n", usage.ru_msgrcv);
        fprintf(stderr, "MK_Signals received:           %10ld\n", usage.ru_nsignals);
        fprintf(stderr, "MK_Voluntary context switches: %10ld\n", usage.ru_nvcsw);
        fprintf(stderr, "MK_Involuntary context sw.:    %10ld\n", usage.ru_nivcsw);
    }

    {
        unsigned int memory_size = 0;
        Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &memory_size);
        fprintf(stderr, "MK_Memory size HEAP:           %10ld MB\n", (long)memory_size);

        Kernel_GetMemorySize(KERNEL_MEMSIZE_STACK, &memory_size);
        fprintf(stderr, "MK_Memory size STACK:          %10ld MB\n", (long)memory_size);

        Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &memory_size);
        fprintf(stderr, "MK_Memory available HEAP:      %10ld MB\n", (long)memory_size);

        Kernel_GetMemorySize(KERNEL_MEMSIZE_STACKAVAIL, &memory_size);
        fprintf(stderr, "MK_Memory available STACK:     %10ld MB\n", (long)memory_size);

        Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPMAX, &memory_size);
        fprintf(stderr, "MK_Maximum memory HEAP:        %10ld MB\n", (long)memory_size);
    }
}


