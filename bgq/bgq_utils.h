/*
 * bgq_utils.h
 *
 *  Created on: Aug 4, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_UTILS_H_
#define BGQ_UTILS_H_

#include <mpi.h>
#include "complex_c99.h"
//#include <complex.h>
#include <assert.h>

#ifndef BGQ_UTILS_C_
#define EXTERN_INLINE inline
#define EXTERN_FIELD extern
#define EXTERN_INIT(val)
#else
#define EXTERN_INLINE extern inline
#define EXTERN_FIELD
#define EXTERN_INIT(val) = val
#endif



#define cs2c99(cs) \
	((cs).re + _Complex_I*(cs).im)

EXTERN_INLINE complex c992cs(double _Complex c99) {
	complex result = { creal(c99), cimag(c99) };
	//result.re = __creal(c99);
	//result.im = __cimag(c99);
	return result;
}


#define _CONCAT(X,Y) X##Y
#define CONCAT(X,Y) _CONCAT(X,Y)

#define CONCAT2(s1,s2) CONCAT(s1,s2)
#define CONCAT3(s1,s2,s3) CONCAT(CONCAT2(s1,s2),s3)
#define CONCAT4(s1,s2,s3,s4) CONCAT(CONCAT3(s1,s2,s3),s4)
#define CONCAT5(s1,s2,s3,s4,s5) CONCAT(s1,CONCAT4(s2,s3,s4,s5))

#define NAME2(s1,s2) CONCAT3(s1,_,s2)
#define NAME3(s1,s2,s3) NAME2(CONCAT3(s1,_,s2),s3)

#define MAKENAME CONCAT2(tmp,__LINE__) /* TODO: Use __COUNTER__ if available */
#define MAKENAME1(s1) CONCAT3(MAKENAME, _, s1)
#define MAKENAME2(s1,s2) CONCAT3(MAKENAME1(s1), _, s2)
#define MAKENAME3(s1,s2,s3) CONCAT3(MAKENAME2(s1,s2), _, s3)
#define MAKENAME4(s1,s2,s3,s4) CONCAT3(MAKENAME3(s1,s2,s3), _, s4)
#define MAKENAME5(s1,s2,s3,s4,s5) CONCAT3(MAKENAME4(s1,s2,s3,s4), _, s5)
#define MAKENAME6(s1,s2,s3,s4,s5,s6) CONCAT3(MAKENAME5(s1,s2,s3,s4,s5), _, s6)

#ifndef STRINGIFY
#define STRINGIFY(V) #V
#endif
#ifndef TOSTRING
#define TOSTRING(V) STRINGIFY(V)
#endif




#define MPI_CHECK(RTNCODE)                                                                 \
	do {                                                                                   \
    	int mpi_rtncode = (RTNCODE);                                                       \
    	if (mpi_rtncode != MPI_SUCCESS) {                                                  \
			fprintf(stderr, "MPI call %s at %s:%d failed: errorcode %d\n", #RTNCODE, __FILE__, __LINE__, mpi_rtncode);  \
			assert(!"MPI call " #RTNCODE " failed");                                       \
			abort();                                                                       \
		}                                                                                  \
	} while (0)

#define L1P_CHECK(RTNCODE)                                                                 \
	do {                                                                                   \
    	int mpi_rtncode = (RTNCODE);                                                       \
    	if (mpi_rtncode != MPI_SUCCESS) {                                                  \
			fprintf(stderr, "L1P call %s at %s:%d failed: errorcode %d\n", #RTNCODE, __FILE__, __LINE__, mpi_rtncode);  \
			assert(!"L1P call " #RTNCODE " failed");                                       \
			abort();                                                                       \
		}                                                                                  \
	} while (0)

EXTERN_INLINE int get_MPI_count(MPI_Status *status) {
	int count;
	MPI_Get_count(status, MPI_BYTE, &count);
	return count;
}

#define WORKLOAD_DECL(COUNTER, TOTAL) \
	int xyz_counter = (COUNTER);         \
	int xyz_orig;            \
	int xyz_torig;           \
	int xyz_total = (TOTAL); \
	int xyz_param

// true:  xyz = 0            .. TRUE_COUNT -> 0 .. TRUE_COUNT
// false: xyz = TRUE_COUNT+1 .. xyz_total  -> 0 .. xyz_total-TRUE_COUNT
//TODO: Can also do a variant without conditional
#define WORKLOAD_SPLIT(TRUE_COUNT) \
	((xyz_counter < TRUE_COUNT) ? (                    \
		xyz_total = (TRUE_COUNT),              \
		true                                   \
	) : (                                      \
		xyz_total = xyz_total - (TRUE_COUNT),  \
		xyz_counter = xyz_counter - (TRUE_COUNT),              \
		assert(xyz_total >= 0 && "Expected more to split"),   \
		false                                  \
	))
/*
 xyz_orig		xyz				return
 0				0				1
 1				1				1
 2				2				1
 ...
 TRUE_COUNT-1	TRUE_COUNT-1	1
 TRUE_COUNT		0				0
 1				0
 2				0
 ...
 xyz_torig-1						0
 */

#define WORKLOAD_PARAM(LENGTH)            \
	(assert((mod(xyz_total, (LENGTH))==0) && "Loop bounds must be a multiple of this parameter"), \
	 xyz_param = (LENGTH),                \
	 xyz_orig = xyz_counter,                      \
	 xyz_torig = xyz_total,               \
	 xyz_total = xyz_torig / xyz_param,   \
	 xyz_counter = xyz_orig / xyz_param,          \
	 mod(xyz_orig, xyz_param))
/*
 xyz_orig	xyz			return (=PARAM)
 0			0			0
 1			0			1
 2			0			2
 ...
 LENGTH-1	0			LENGTH-1
 LENGTH		1			0
 1			1
 1			2
 ...
 xyz_torig	xyz_total	LENGTH

 xyz_orig =	xyz*LENGTH + return
 */

#define WORKLOAD_SECTION(SECTIONS)        \
	(xyz_param = (SECTIONS),              \
	 xyz_orig  = xyz_counter,                     \
	 xyz_torig = xyz_total,               \
	 xyz_total = xyz_param,  			  \
	 xyz_counter = xyz_orig / (xyz_torig/xyz_param),          \
	 mod(xyz_orig, (xyz_torig/xyz_param)))
/*
 xyz_orig	xyz			return
 0			0			0
 1			0			1
 2			0			2
 ...
 0
 1			0
 1			1
 1			2
 ...
 xyz_torig-1	SECTIONS-1
 xyz_torig	SECTIONS	xyz_torig/SECTIONS
 (=xyz_total)

 xyz_orig =	xyz*(xyz_torig/SECTIONS) + return
 */

#define WORKLOAD_TILE(TILES)   					   \
		(xyz_param = (TILES),                      \
		 xyz_orig = xyz_counter,                           \
		 xyz_torig = xyz_total,                    \
		 xyz_total = TILES,     				   \
		 xyz_counter = mod(xyz_orig, (xyz_torig / xyz_param)), \
		 xyz_orig / (xyz_torig / xyz_param))
/*
 xyz_orig	xyz	(=TILE)	return
 0			0			0
 1			1			0
 2			2			0
 ...
 TILES-1		TILES-1		0
 TILES		0			1
 1			1
 2			1
 ...
 xyz_torig	TILES		xyz_torig/TILES
 (=xyz_total)

 xyz_orig =	xyz +		return*TILES
 */

#define WORKLOAD_CHUNK(LENGTH)                 \
	(assert((mod(xyz_total, (LENGTH))==0) && "Loop bounds must be a multiple of this parameter"), \
	 xyz_param = (LENGTH),                     \
	 xyz_orig = xyz_counter,                           \
	 xyz_torig = xyz_total,                    \
	 xyz_total = xyz_torig / xyz_param,        \
	 xyz_counter = mod(xyz_orig, xyz_total),               \
	 xyz_orig / xyz_total)
/*
 xyz_orig	xyz			return (=CHUNCK)
 0			0			0
 1			1			0
 2			2			0
 ...
 0
 0			1
 1			1
 2			1
 ...
 xyz_torig-1	xyz_total-1	LENGTH-1
 xyz_torig	xyz_total	LENGTH
 (=xyz_torig/LENGTH)

 xyz_orig =	xyz + 		return*xyz_total
 */

#define WORKLOAD_CHECK \
	assert(xyz_counter==0); \
	assert(xyz_total==1);


#ifndef BGQMOD
#define BGQMOD 0
#endif

// dividend >= -1
// divisor > 0
// 0 <= result < divisor
EXTERN_INLINE int moddown(const int dividend, const int divisor) {
#if BGQMOD==0
	// Compilers can therefore optimize it to bit-operations specific to the target machine
	// This is not possible with the %-operator on signed operands because it required the result to have the sign of the dividend (hence its the remainder, not a modulus)
	unsigned int udividend = dividend + divisor;
	unsigned int udivisor = divisor;
	return udividend % udivisor; // Use unsigned operation here to enable the optimizer to use bit-tricks (like dividend&1 if divisor==2)
#elif BGQMOD==1
	int result = dividend % divisor;
	if (result < 0)
		result += divisor;
	return result;
#elif BGQMOD==2
	if (dividend >= 0)
		return dividend % divisor;
	else
		return dividend % divisor + divisor;
#elif BGQMOD==3
	return ((dividend % divisor) + divisor) % divisor;
#endif
}

 // dividend >= 0
 // divisor > 0
 // 0 <= result < divisor
EXTERN_INLINE int mod(const int dividend, const int divisor) {
  	unsigned int udividend = dividend;
 	unsigned int udivisor = divisor;
 	return udividend % udivisor; // Use unsigned operation here to enable the optimizer to use bit-tricks (like dividend&1 if divisor==2)
}


#ifndef BGQDIV
#define BGQDIV 0
#endif

// Always round towards negative infinity
// dividend >= -1
// divisor > 0
// result >= -1
EXTERN_INLINE int divdown(const int dividend, const int divisor) {
#if BGQDIV==0
	unsigned int udividend = dividend + divisor; // Dividend can be -1
	unsigned int udivisor = divisor;
	return (int)(udividend / udivisor) - 1;
#elif BGQDIV==1
	return (dividend - mod(dividend,divisor)) / divisor;
#elif BGQDIV==2
	if (dividend >= 0)
		return dividend / divisor;
	else
		return (dividend - divisor + 1) / divisor;
#endif
}



#define BGQ_ENTER_FUNC                                     \
	if (g_proc_id==0) {                                    \
		fprintf(stderr, "MK ENTER_FUNC %s\n",  __func__);  \
	}

#define BGQ_BEACON \
	{  \
	static bool CONCAT(beacon, __LINE__); \
	if (!CONCAT(beacon, __LINE__)) \
		master_print("BEACON: File: %s, line %d, func: %s\n", __FILE__, __LINE__, __func__); \
	CONCAT(beacon, __LINE__) = true; \
	}

EXTERN_FIELD double bgq_g_wtick EXTERN_INIT(0);

EXTERN_INLINE void bgq_init_wtime() {
	bgq_g_wtick = MPI_Wtick();
}

// return wall clock time in seconds
EXTERN_INLINE double bgq_wtime() {
	assert(bgq_g_wtick != 0);
	return MPI_Wtime();
}


#define master_print(...)              \
	if (g_proc_id == 0)                 \
		if (omp_get_thread_num() == 0)  \
			fprintf(stderr, __VA_ARGS__)

#define master_error(errcode, ...) \
	do {                            \
		master_print(__VA_ARGS__);  \
		exit(errcode);              \
	} while (0)


void *malloc_aligned(size_t size, size_t alignment);


EXTERN_INLINE double max(double const lhs, double const rhs) {
	if (lhs > rhs)
		return lhs;
	return rhs;
}

EXTERN_INLINE double min(double const lhs, double const rhs) {
	if (lhs > rhs)
		return rhs;
	return lhs;
}


#undef EXTERN_INLINE
#undef EXTERN_FIELD

#endif /* BGQ_UTILS_H_ */
