/*
 * bgq_field.inc.h
 *
 *  Created on: Aug 15, 2012
 *      Author: meinersbur
 */


#include "bgq_field.h"
#include <omp.h>

#if BGQ_PREFETCH_LIST
#include <l1p/pprefetch.h>
#endif

#ifndef BGQ_PRECISION
#warning Specify a precision before including this file
#define BGQ_PRECISION 64
#include "bgq_precisionselect.inc.c"
#endif






// Type names
#define bgq_spinorsite NAME2(bgq_spinorsite,PRECISION)
#define bgq_spinorfield NAME2(bgq_spinorfield,PRECISION)

#define bgq_gaugesite NAME2(bgq_gaugesite,PRECISION)
#define bgq_gaugeeodir NAME2(bgq_gaugeeodir,PRECISION)
#define bgq_gaugefield NAME2(bgq_gaugefield,PRECISION)

#define bgq_weylsite NAME2(bgq_weylsite,PRECISION)
#define bgq_weylfield NAME2(bgq_weylfield,PRECISION)

// Field names
#define g_spinorfields NAME2(g_spinorfields,PRECISION)
#define g_spinorfields_data NAME2(g_spinorfields_data,PRECISION)
#define g_spinorfields_data_coords NAME2(g_spinorfields_data_coords,PRECISION)

#define g_num_total_spinorfields NAME2(g_num_total_spinorfields,PRECISION)
int g_num_total_spinorfields;

typedef struct {
	COMPLEX_PRECISION s[4][3][PHYSICAL_LK]; /* 4*3*2*sizeof(COMPLEX_PRECISION) = 384;192 bytes (6;3 L1 cache lines) */
	//COMPLEX_PRECISION padding[4][1][PHYSICAL_LK];
} bgq_spinorsite;

typedef bgq_spinorsite (*bgq_spinorfield);


typedef struct {
	COMPLEX_PRECISION c[3][3][PHYSICAL_LK]; /* 3*3*2*sizeof(COMPLEX_PRECISION) = 288;144 bytes (4.5;2.25 L1 cache lines) */
	//COMPLEX_PRECISION padding[6];
} bgq_gaugesite;
typedef struct {
	bgq_gaugesite *(eodir[PHYSICAL_LP][PHYSICAL_LD]);
} bgq_gaugeeodir;
typedef bgq_gaugeeodir (*bgq_gaugefield);


typedef struct {
	COMPLEX_PRECISION s[2][3][PHYSICAL_LK]; /* 2*3*2*sizeof(COMPLEX_PRECISION) = 192;96 bytes (3;1.5 L1 cache lines) */
	//COMPLEX_PRECISION padding[2][1][PHYSICAL_LK]; /* 4*16=64 byte To fill complete L2 cache line */
} bgq_weylsite;
typedef bgq_weylsite (*bgq_weylfield);



////////////////////////////////////////////////////////////////////////////////
// Spinorfields



extern bgq_spinorfield *g_spinorfields;

#define bgq_init_spinorfields NAME2(bgq_init_spinorfields,PRECISION)
void bgq_init_spinorfields(int count, int chi_count);
#define bgq_free_spinofields NAME2(bgq_free_spinofields,PRECISION)
void bgq_free_spinofields();

//#define bgq_num_spinorfields NAME2(bgq_num_spinorfields,PRECISION)
//extern int bgq_num_spinorfields;

#define bgq_spinorfield_find_index NAME2(bgq_spinorfield_find_index,PRECISION)
int bgq_spinorfield_find_index(bgq_spinorfield spinorfield);

#define bgq_translate_spinorfield NAME2(bgq_translate_spinorfield,PRECISION)
bgq_spinorfield bgq_translate_spinorfield(spinor * const field);
#define bgq_transfer_spinorfield NAME2(bgq_transfer_spinorfield,PRECISION)
void bgq_transfer_spinorfield(bool isOdd, bgq_spinorfield targetfield, spinor *sourcefield);
#define bgq_spinorfield_resetcoord NAME2(bgq_spinorfield_resetcoord,PRECISION)
void bgq_spinorfield_resetcoord(bgq_spinorfield spinorfield, bool isOdd, int expected_reads_min, int expected_reads_max, int expected_writes_min, int expected_writes_max);

#define bgq_spinorfield_transfer_back NAME2(bgq_spinorfield_transfer_back,PRECISION)
void bgq_spinorfield_transfer_back(const bool isOdd, spinor * const targetfield, bgq_spinorfield sourcefield);

#define assert_spinorfield_coord NAME2(assert_spinorfield_coord,PRECISION)
bool assert_spinorfield_coord(bgq_spinorfield spinorfield, bool isOdd, int t, int x, int y, int z, int tv, int k, int v, int c, bool isRead, bool isWrite);
#define assert_spinorcoord NAME2(assert_spinorcoord,PRECISION)
bool assert_spinorcoord(bgq_spinorfield spinorfield, bool isOdd, int t, int x, int y, int z, int tv, int k, bool isRead, bool isWrite);

#define bgq_spinorfield_get NAME2(bgq_spinorfield_get,PRECISION)
COMPLEX_PRECISION bgq_spinorfield_get(bgq_spinorfield spinorfield, bool isOdd, int t, int x, int y, int z, int v, int c);
#define bgq_spinorfield_set NAME2(bgq_spinorfield_set,PRECISION)
void bgq_spinorfield_set(bgq_spinorfield spinorfield, bool isOdd, int t, int x, int y, int z, int v, int c, COMPLEX_PRECISION value);

#define bgq_spinorfield_compare NAME2(bgq_spinorfield_compare,PRECISION)
double bgq_spinorfield_compare(const bool isOdd, bgq_spinorfield const bgqfield, spinor * const reffield, bool silent);

#define bgq_spinorfield_isOdd NAME2(bgq_spinorfield_isOdd,PRECISION)
bool bgq_spinorfield_isOdd(bgq_spinorfield spinorfield);

#define bgq_spinorfield_setOdd NAME2(bgq_spinorfield_setOdd,PRECISION)
void bgq_spinorfield_setOdd(bgq_spinorfield spinorfield, bool isOdd, bool overwrite);

////////////////////////////////////////////////////////////////////////////////
// Gaugefield

#define g_gaugefield NAME2(g_gaugefield,PRECISION)
extern bgq_gaugefield g_gaugefield;

#define bgq_init_gaugefield NAME2(bgq_init_gaugefield,PRECISION)
void bgq_init_gaugefield();
#define bgq_free_gaugefield NAME2(bgq_free_gaugefield,PRECISION)
void bgq_free_gaugefield();

#define bgq_transfer_gaugefield NAME2(bgq_transfer_gaugefield,PRECISION)
void bgq_transfer_gaugefield(bgq_gaugefield targetfield, su3 **sourcefield);

#define bgq_gaugefield_resetcoord NAME2(bgq_gaugefield_resetcoord,PRECISION)
void bgq_gaugefield_resetcoord(bgq_gaugefield gaugefield, int expected_reads_min, int expected_reads_max, int expected_writes_min, int expected_writes_max);

#define assert_gaugeval NAME2(assert_gaugeval,PRECISION)
bool assert_gaugeval(bgq_gaugefield gaugefield, bool isOdd, int t, int x, int y, int z, int tv, int k, direction dir, int i, int l, bool isRead, bool isWrite);

#define assert_gaugesite NAME2(assert_gaugesite,PRECISION)
bool assert_gaugesite(bgq_gaugefield gaugefield, bool isOdd, int t, int x, int y, int z, int tv, int k, direction dir, bool isRead, bool isWrite);

#define bgq_gaugefield_get NAME2(bgq_gaugefield_get,PRECISION)
COMPLEX_PRECISION bgq_gaugefield_get(bgq_gaugefield gaugefield, bool isOdd, int t, int x, int y, int z, direction d, int i, int l);

#define bgq_gaugefield_set NAME2(bgq_gaugefield_set,PRECISION)
void bgq_gaugefield_set(bgq_gaugefield gaugefield, bool isOdd, int t, int x, int y, int z, direction d, int i, int l, double _Complex value);


////////////////////////////////////////////////////////////////////////////////
// Weylfields

#define weylxchange_recv NAME2(weylxchange_recv,PRECISION)
extern bgq_weylfield weylxchange_recv[6];
#define weylxchange_send NAME2(weylxchange_send,PRECISION)
extern bgq_weylfield weylxchange_send[6];
#define weylxchange_size NAME2(weylxchange_size,PRECISION)
extern size_t weylxchange_size[3];
#define weylexchange_destination NAME2(weylexchange_destination,PRECISION)
extern int weylexchange_destination[6];
#define weylexchange_request_recv NAME2(weylexchange_request_recv,PRECISION)
extern MPI_Request weylexchange_request_recv[6];
#define weylexchange_request_send NAME2(weylexchange_request_send,PRECISION)
extern MPI_Request weylexchange_request_send[6];

#define recvbuf NAME2(recvbuf,PRECISION)
//extern double recvbuf;
#define recvrequest NAME2(recvrequest,PRECISION)
extern MPI_Request recvrequest;
#define sendbuf NAME2(sendbuf,PRECISION)
//extern double sendbuf;
#define sendrequest NAME2(sendrequest,PRECISION)
extern MPI_Request sendrequest;

#define bgq_hm_init NAME2(bgq_hm_init,PRECISION)
void bgq_hm_init();
#define bgq_hm_free NAME2(bgq_hm_free,PRECISION)
void bgq_hm_free();

#define assert_weylval_t NAME2(assert_weylval_t,PRECISION)
bool assert_weylval_t(bgq_weylfield weylfield, bool isOdd, int t, int x, int y, int z, int xv, int k, int v, int c, bool isRead, bool isWrite);
#define assert_weylfield_t NAME2(assert_weylfield_t,PRECISION)
bool assert_weylfield_t(bgq_weylfield weylfield, bool isOdd, int t, int x, int y, int z, int xv, int k, bool isRead, bool isWrite);
#define bgq_weylfield_t_resetcoord NAME2(bgq_weylfield_t_resetcoord,PRECISION)
void bgq_weylfield_t_resetcoord(bgq_weylfield weylfield, int t, bool isOdd, int expected_reads_min, int expected_reads_max, int expected_writes_min, int expected_writes_max);

#define assert_weylval_x NAME2(assert_weylval_x,PRECISION)
bool assert_weylval_x(bgq_weylfield weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k, int v, int c, bool isRead, bool isWrite);
#define assert_weylfield_x NAME2(assert_weylfield_x,PRECISION)
bool assert_weylfield_x(bgq_weylfield weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k, bool isRead, bool isWrite);
#define bgq_weylfield_x_resetcoord NAME2(bgq_weylfield_x_resetcoord,PRECISION)
void bgq_weylfield_x_resetcoord(bgq_weylfield weylfield, int x, bool isOdd, int expected_reads_min, int expected_reads_max, int expected_writes_min, int expected_writes_max);

#define assert_weylval_y NAME2(assert_weylval_y,PRECISION)
bool assert_weylval_y(bgq_weylfield weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k, int v, int c, bool isRead, bool isWrite);
#define assert_weylfield_y NAME2(assert_weylfield_y,PRECISION)
bool assert_weylfield_y(bgq_weylfield weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k, bool isRead, bool isWrite);
#define bgq_weylfield_y_resetcoord NAME2(bgq_weylfield_y_resetcoord,PRECISION)
void bgq_weylfield_y_resetcoord(bgq_weylfield weylfield, int y, bool isOdd, int expected_reads_min, int expected_reads_max, int expected_writes_min, int expected_writes_max);


#define bgq_weylsite_callback NAME2(bgq_weylsite_callback,PRECISION)
typedef void (*bgq_weylsite_callback)(bgq_weylfield weylfield, direction dir, bool isSend, bool isOdd, int t, int x, int y, int z, int v, int c, COMPLEX_PRECISION *val, int tag);

#define bgq_weylfield_foreach NAME2(bgq_weylfield_foreach,PRECISION)
void bgq_weylfield_foreach(bgq_weylfield weylfield, direction dir, bool isSend, bool isOdd, bgq_weylsite_callback callback, int tag);

#if 0
#define bgq_setbgqval NAME2(bgq_setbgqval,PRECISION)
static void bgq_setbgqval(bgq_weylfield weylfield, direction dir, bool isSend, bool isOdd, int t, int x, int y, int z, int v, int c, COMPLEX_PRECISION *val, int tag) {
	if ( (v==1) && (c==2) ) {
		char buf[20];
		snprintf(buf, sizeof(buf), "bgqval%d", tag);
		bgq_setbgqvalue(t,x,y,z,tag,*val,buf);
	}
}
#endif


#if 0
#define bgq_expected NAME2(bgq_expected,PRECISION)
static COMPLEX_PRECISION bgq_expected(direction dir, bool isOdd, int t, int x, int y, int z, int v, int c, bool isSend, int tag) {
	COMPLEX_PRECISION result = t + (tag + 0.1*dir)*_Complex_I;

	if (!isSend) {
		result = conj(result);
	}

	return result;
}
#endif


#if 0
#define bgq_weylfield_setcoordfield NAME2(bgq_weylfield_setcoordfield,PRECISION)
static void bgq_weylfield_setcoordfield(bgq_weylfield weylfield, direction dir, bool isOdd, bool isSend, int tag) {
	if (isSend) {
		switch (dir) {
		case XDOWN:
		case XUP: {
			//for (direction d = XUP; d <= XDOWN; d+=1)
				for (int tv = 0; tv < PHYSICAL_LTV; tv+=1)
					for (int y = 0; y < LOCAL_LY; y+=1)
						for (int z = 0; z < LOCAL_LY; z+=1) {
							int x = (dir==XUP) ? -1 : LOCAL_LX;
							const int t1 = ((isOdd+x+y+z)&1)+tv*PHYSICAL_LP*PHYSICAL_LK;
							const int t2 = t1 + 2;
							bgq_weylsite *weylsite_send = BGQ_WEYLSITE_X(weylfield, !isOdd, tv, x + 2*(dir==XUP)-1, y, z, t1, t2, false, false);
							for (int v = 0; v < 2; v+=1)
								for (int c = 0; c < 3; c+=1) {
									weylsite_send->s[v][c][0] = bgq_expected(dir, isOdd, t1, x, y, z, v, c, isSend,tag);
									weylsite_send->s[v][c][1] = bgq_expected(dir, isOdd, t2, x, y, z, v, c, isSend,tag);
								}
						}
		} break;
		default:
			break;
		}
	} else {
		switch (dir) {
		case XDOWN:
		case XUP: {
			//for (direction d = XUP; d <= XDOWN; d+=1)
						for (int tv = 0; tv < PHYSICAL_LTV; tv+=1)
							for (int y = 0; y < LOCAL_LY; y+=1)
								for (int z = 0; z < LOCAL_LY; z+=1) {
									int x = (dir==XUP) ?  LOCAL_LX-1 : 0;
									const int t1 = ((isOdd+x+y+z)&1)+tv*PHYSICAL_LP*PHYSICAL_LK;
									const int t2 = t1 + 2;
									bgq_weylsite *weylsite_recv = BGQ_WEYLSITE_X(weylfield, !isOdd, tv, x + 2*(dir==XUP)-1, y, z, t1, t2, false, false);
									for (int v = 0; v < 2; v+=1)
										for (int c = 0; c < 3; c+=1) {
											weylsite_recv->s[v][c][0] = bgq_expected(dir, isOdd, t1, x, y, z, v, c, isSend,tag);
											weylsite_recv->s[v][c][1] = bgq_expected(dir, isOdd, t2, x, y, z, v, c, isSend,tag);
										}
								}
		} break;
		default:
			break;
		}
	}
}

#define bgq_weylfield_cmpcoordfield NAME2(bgq_weylfield_cmpcoordfield,PRECISION)
static void bgq_weylfield_cmpcoordfield(bgq_weylfield weylfield, direction dir, bool isOdd, bool isSend, int tag) {
	if (isSend) {
		switch (dir) {
		case XDOWN:
		case XUP: {
				for (int tv = 0; tv < PHYSICAL_LTV; tv+=1)
					for (int y = 0; y < LOCAL_LY; y+=1)
						for (int z = 0; z < LOCAL_LY; z+=1) {
							int x = (dir==XUP) ? -1 : LOCAL_LX;
							const int t1 = ((isOdd+x+y+z)&1)+tv*PHYSICAL_LP*PHYSICAL_LK;
							const int t2 = t1 + 2;
							bgq_weylsite *weylsite_send = BGQ_WEYLSITE_X(weylfield, !isOdd, tv, x + 2*(dir==XUP)-1, y, z, t1, t2, false, false);
							for (int v = 0; v < 2; v+=1)
								for (int c = 0; c < 3; c+=1) {
									COMPLEX_PRECISION expect1 = bgq_expected(dir, isOdd, t1, x, y, z, v, c, isSend,tag);
									COMPLEX_PRECISION real1 = weylsite_send->s[v][c][0];
									if (real1 != expect1) {
										master_print("Mismatch at (%d,%d,%d,%d,dir=%d,isSend=%d) expected=%f + %fi != real=%f + %fi\n", t1, x, y, z, dir, isSend, creal(expect1), cimag(expect1), creal(real1), cimag(real1));
									}

									COMPLEX_PRECISION expect2 = bgq_expected(dir, isOdd, t2, x, y, z, v, c, isSend,tag);
									COMPLEX_PRECISION real2 = weylsite_send->s[v][c][1];
									if (real2 != expect2) {
										master_print("Mismatch at (%d,%d,%d,%d,dir=%d,isSend=%d) expected=%f + %fi != real=%f + %fi\n", t2, x, y, z, dir, isSend, creal(expect2), cimag(expect2), creal(real2), cimag(real2));
									}
								}
						}
		} break;
		default:
			break;
		}
	} else {
		switch (dir) {
		case XDOWN:
		case XUP: {
						for (int tv = 0; tv < PHYSICAL_LTV; tv+=1)
							for (int y = 0; y < LOCAL_LY; y+=1)
								for (int z = 0; z < LOCAL_LY; z+=1) {
									int x = (dir==XUP) ?  LOCAL_LX-1 : 0;
									const int t1 = ((isOdd+x+y+z)&1)+tv*PHYSICAL_LP*PHYSICAL_LK;
									const int t2 = t1 + 2;
									bgq_weylsite *weylsite_recv = BGQ_WEYLSITE_X(weylfield, !isOdd, tv, x + 2*(dir==XUP)-1, y, z, t1, t2, false, false);
									for (int v = 0; v < 2; v+=1)
										for (int c = 0; c < 3; c+=1) {
											COMPLEX_PRECISION expect1 = bgq_expected(dir, isOdd, t1, x, y, z, v, c, isSend,tag);
												COMPLEX_PRECISION real1 = weylsite_recv->s[v][c][0];
												if (real1 != expect1) {
													master_print("Mismatch at (%d,%d,%d,%d,dir=%d,isSend=%d) expected=%f + %fi != real=%f + %fi\n", t1, x, y, z, dir, isSend, creal(expect1), cimag(expect1), creal(real1), cimag(real1));
												}

												COMPLEX_PRECISION expect2 = bgq_expected(dir, isOdd, t2, x, y, z, v, c, isSend,tag);
												COMPLEX_PRECISION real2 = weylsite_recv->s[v][c][1];
												if (real2 != expect2) {
													master_print("Mismatch at (%d,%d,%d,%d,dir=%d,isSend=%d) expected=%f + %fi != real=%f + %fi\n", t2, x, y, z, dir, isSend, creal(expect2), cimag(expect2), creal(real2), cimag(real2));
												}
										}
								}
		} break;
		default:
			break;
		}
	}
}
#endif

#if BGQ_PREFETCH_LIST
#define bgq_listprefetch_handle NAME2(bgq_listprefetch_handle,PRECISION)
#ifndef BGQ_FIELD_INC_C_
extern
#endif
L1P_Pattern_t *(*bgq_listprefetch_handle)[2/*even/odd*/][2/*kamul/nokamul*/][64/*For each thread*/][2/*nobody*/][2/*nosurface*/][6/*total threads*/]
#ifdef BGQ_FIELD_INC_C_
= NULL
#endif
;
#endif

//#ifndef BGQ_FIELD_INC_C_
//#include "bgq_precisionselect.inc.c"
//#endif
