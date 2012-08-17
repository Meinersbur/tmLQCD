/*
 * bgq_field.inc.h
 *
 *  Created on: Aug 15, 2012
 *      Author: meinersbur
 */

#include "bgq_field.h"

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

// Function names
#define bgq_init_spinorfields NAME2(bgq_init_spinorfields,PRECISION)
#define bgq_free_spinofields NAME2(bgq_free_spinofields,PRECISION)
#define bgq_translate_spinorfield NAME2(bgq_translate_spinorfield,PRECISION)
#define bgq_transfer_spinorfield NAME2(bgq_transfer_spinorfield,PRECISION)
#define bgq_spinorfield_resetcoord NAME2(bgq_spinorfield_resetcoord,PRECISION)
#define assert_spinorfield_coord NAME2(assert_spinorfield_coord,PRECISION)
#define assert_spinorcoord NAME2(assert_spinorcoord,PRECISION)



typedef struct {
	COMPLEX_PRECISION s[4][3][PHYSICAL_LK];
} bgq_spinorsite;

typedef bgq_spinorsite (*bgq_spinorfield);


typedef struct {
	COMPLEX_PRECISION c[3][3][PHYSICAL_LK];
} bgq_gaugesite;
typedef struct {
	bgq_gaugesite *(eodir[PHYSICAL_LP][PHYSICAL_LD]);
} bgq_gaugeeodir;
typedef bgq_gaugeeodir (*bgq_gaugefield);


typedef struct {
	COMPLEX_PRECISION s[2][3][PHYSICAL_LK];
} bgq_weylsite;
typedef bgq_weylsite (*bgq_weylfield);



////////////////////////////////////////////////////////////////////////////////
// Spinorfields



extern bgq_spinorfield *g_spinorfields;


void bgq_init_spinorfields(int count);
void bgq_free_spinofields();

bgq_spinorfield bgq_translate_spinorfield(spinor * const field);
void bgq_transfer_spinorfield(bool isOdd, bgq_spinorfield targetfield, spinor *sourcefield);
void bgq_spinorfield_resetcoord(bgq_spinorfield spinorfield, bool isOdd, int expected_reads_min, int expected_reads_max, int expected_writes_min, int expected_writes_max);

bool assert_spinorfield_coord(bgq_spinorfield spinorfield, bool isOdd, int t, int x, int y, int z, int tv, int k, int v, int c, bool isRead, bool isWrite);
bool assert_spinorcoord(bgq_spinorfield spinorfield, bool isOdd, int t, int x, int y, int z, int tv, int k, bool isRead, bool isWrite);

#define bgq_spinorfield_get NAME2(bgq_spinorfield_get,PRECISION)
COMPLEX_PRECISION bgq_spinorfield_get(bgq_spinorfield spinorfield, bool isOdd, int t, int x, int y, int z, int v, int c);
#define bgq_spinorfield_set NAME2(bgq_spinorfield_set,PRECISION)
void bgq_spinorfield_set(bgq_spinorfield spinorfield, bool isOdd, int t, int x, int y, int z, int v, int c, COMPLEX_PRECISION value);



////////////////////////////////////////////////////////////////////////////////
// Gaugefield

#define g_gaugefield NAME2(g_gaugefield,PRECISION)
bgq_gaugefield g_gaugefield;

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



//#ifndef BGQ_FIELD_INC_C_
//#include "bgq_precisionselect.inc.c"
//#endif

