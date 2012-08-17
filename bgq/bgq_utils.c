/*
 * bgq_utils.h
 *
 *  Created on: Aug 4, 2012
 *      Author: meinersbur
 */


#define BGQ_UTILS_C_
#include "bgq_utils.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


void *malloc_aligned(size_t size, size_t alignment) {
	void *result = NULL;
	int errcode = posix_memalign(&result, alignment, size);
	if (errcode != 0) {
		fprintf(stderr, "malloc returned %d\n", errcode);
		exit(10);
	}
	memset(result, 0, size);
	return result;
}
