/*
 * bgq_dispatch.h
 *
 *  Created on: Oct 15, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_DISPATCH_H_
#define BGQ_DISPATCH_H_

typedef void (bgq_dispatch_func*)(void *arg, size_t tid, size_t threads);

void bgq_worker();
void bgq_master_call(bgq_dispatch_func func, void *arg);
void bgq_master_nomorework();

#endif /* BGQ_DISPATCH_H_ */
