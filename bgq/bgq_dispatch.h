/*
 * bgq_dispatch.h
 *
 *  Created on: Oct 15, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_DISPATCH_H_
#define BGQ_DISPATCH_H_

#include <string.h>

typedef void (*bgq_worker_func)(void *arg, size_t tid, size_t threads);
typedef int (*bgq_master_func)(void *arg);

int bgq_parallel(bgq_master_func master_func, void *master_arg);

void bgq_worker();
void bgq_master_call(bgq_worker_func worker_func, void *arg);
void bgq_master_sync();

#endif /* BGQ_DISPATCH_H_ */
