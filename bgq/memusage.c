#include "memusage.h"

#include "config.h"
#include <sys/resource.h>

#ifdef BGP
#include <common/bgp_personality.h>
#include <common/bgp_personality_inlines.h>
#include <spi/kernel_interface.h>
#include <spi/bgp_SPI.h>
#endif

#ifdef BGQ
#include <kernel/process.h>
#include <kernel/location.h>
#include <kernel/memory.h>
//#include <kernel/personality.h>
#endif

#include "global.h"
#include <string.h>

#include <stdio.h>
#include <malloc.h>


#if defined(BGQ) || defined(BGP)
/* returns memory per core in MBytes */
unsigned bg_coreMB()
{
    unsigned procMB, coreMB;
    static Personality_t mybgp;
    Kernel_GetPersonality(&mybgp, sizeof(Personality_t));
    //procMB = BGP_Personality_DDRSizeMB(&mybgp);
    coreMB = procMB/Kernel_ProcessCount();
    return coreMB;
}
#endif

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

#if defined(BGQ) || defined(BGP)
    {
        unsigned memPerCore = bg_coreMB();
        fprintf(stderr, "MK_Memory available per Core:  %10ld MB\n", memPerCore);
    }
#endif

    struct rusage usage;
    memset(&usage, 0, sizeof(usage));
    int returncode = 0;
    if ((returncode = getrusage(RUSAGE_SELF, &usage)))
    {
        fprintf(stderr, "MK_getrusage error; return code %d\n", returncode);
    }
    else
    {
        fprintf(stderr, "MK_User time used:             %10ld,%6ld secs\n", (long)(usage.ru_utime.tv_sec), (long)(usage.ru_utime.tv_usec));
        fprintf(stderr, "MK_System time used:           %10ld,%6ld secs\n", (long)(usage.ru_stime.tv_sec), (long)(usage.ru_stime.tv_usec));
        fprintf(stderr, "MK_Maximum resident set size:  %10ld KB\n", usage.ru_maxrss);
        fprintf(stderr, "MK_Signals                     %ld\n", usage.ru_nsignals);
#if 0
        fprintf(stderr, "MK_Integral shared memory size:%10ld MB\n", (long)(usage.ru_ixrss) / 1024);
        fprintf(stderr, "MK_Integral unshared data size:%10ld MB\n", (long)(usage.ru_idrss) / 1024);
        fprintf(stderr, "MK_Integral unshared stack s.: %10ld MB\n", (long)(usage.ru_isrss) / 1024);
        fprintf(stderr, "MK_Page reclaims:              %10ld\n", (long)(usage.ru_minflt));
        fprintf(stderr, "MK_Page faults:                %10ld\n", (long)(usage.ru_majflt));
        fprintf(stderr, "MK_Swaps:                      %10ld\n", (long)(usage.ru_nswap));
        fprintf(stderr, "MK_Block input operations:     %10ld\n", (long)(usage.ru_inblock));
        fprintf(stderr, "MK_Block output operations:    %10ld\n", (long)(usage.ru_oublock));
        fprintf(stderr, "MK_Messages sent:              %10ld\n", (long)(usage.ru_msgsnd));
        fprintf(stderr, "MK_Messages received:          %10ld\n", (long)(usage.ru_msgrcv));
        fprintf(stderr, "MK_Signals received:           %10ld\n", (long)(usage.ru_nsignals));
        fprintf(stderr, "MK_Voluntary context switches: %10ld\n", (long)(usage.ru_nvcsw));
        fprintf(stderr, "MK_Involuntary context sw.:    %10ld\n", (long)(usage.ru_nivcsw));
#endif
    }

#if defined(BGP) || defined(BGQ)
    {
        unsigned long memory_size = 0;
        //Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &memory_size);
        //fprintf(stderr, "MK_Memory size HEAP:           %10ld MB\n", (long)memory_size / (1024 * 1024));

        //Kernel_GetMemorySize(KERNEL_MEMSIZE_STACK, &memory_size);
        //fprintf(stderr, "MK_Memory size STACK:          %10ld MB\n", (long)memory_size / (1024 * 1024));

        //Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &memory_size);
        //fprintf(stderr, "MK_Memory available HEAP:      %10ld MB\n", (long)memory_size / (1024 * 1024));

        //Kernel_GetMemorySize(KERNEL_MEMSIZE_STACKAVAIL, &memory_size);
        //fprintf(stderr, "MK_Memory available STACK:     %10ld MB\n", (long)memory_size / (1024 * 1024));

        Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPMAX, &memory_size);
        fprintf(stderr, "MK_Maximum memory HEAP:        %10ld MB\n", (long)memory_size / (1024 * 1024));

        //Kernel_GetMemorySize(KERNEL_MEMSIZE_GUARD, &memory_size);
        //fprintf(stderr, "MK_Heap guardpage:             %10ld MB\n", (long)memory_size / (1024 * 1024));

        //Kernel_GetMemorySize(KERNEL_MEMSIZE_SHARED, &memory_size);
        //fprintf(stderr, "MK_Shared memory:              %10ld MB\n", (long)memory_size / (1024 * 1024));

#if 0
        Kernel_GetMemorySize(KERNEL_MEMSIZE_PERSIST, &memory_size);
        fprintf(stderr, "MK_Persistent memory:          %10ld MB\n", (long)memory_size / (1024 * 1024));
#endif


#ifdef BGQ
        uint64_t shared, persist, heapavail, stackavail, stack, heap, guard, mmap;
        Kernel_GetMemorySize(KERNEL_MEMSIZE_SHARED, &shared);
        Kernel_GetMemorySize(KERNEL_MEMSIZE_PERSIST, &persist);
        Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &heapavail);
        Kernel_GetMemorySize(KERNEL_MEMSIZE_STACKAVAIL, &stackavail);
        Kernel_GetMemorySize(KERNEL_MEMSIZE_STACK, &stack);
        Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
        Kernel_GetMemorySize(KERNEL_MEMSIZE_GUARD, &guard);
        Kernel_GetMemorySize(KERNEL_MEMSIZE_MMAP, &mmap);

        fprintf(stderr, "MK_Allocated heap: %.2f MB, avail. heap: %.2f MB\n", (double)heap/(1024*1024), (double)heapavail/(1024*1024));
        fprintf(stderr, "MK_Allocated stack: %.2f MB, avail. stack: %.2f MB\n", (double)stack/(1024*1024), (double)stackavail/(1024*1024));
        fprintf(stderr, "MK_Memory: shared: %.2f MB, persist: %.2f MB, guard: %.2f MB, mmap: %.2f MB\n", (double)shared/(1024*1024), (double)persist/(1024*1024), (double)guard/(1024*1024), (double)mmap/(1024*1024));
#endif

        struct mallinfo m;
        m = mallinfo();

        unsigned int arena = m.arena;          /* size to sbrk */
        fprintf(stderr, "MK_arena = %d \n", arena);

        unsigned int uordblks = m.uordblks;     /* chunks in use, in bytes */
        fprintf(stderr, "MK_uordblks = %d \n", uordblks);

        unsigned int hblkhd = m.hblkhd;         /* mmap memory in bytes */
        fprintf(stderr, "MK_hblkhd = %d \n", hblkhd);

        unsigned int total_heap = uordblks + hblkhd;
        fprintf(stderr, "MK_total_heap = %d \n", total_heap);

        fflush(stderr);
    }
#endif

}


