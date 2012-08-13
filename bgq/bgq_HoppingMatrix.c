/*
 * bgq_HoppingMatrix.c
 *
 *  Created on: Jul 27, 2012
 *      Author: meinersbur
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bgq_field.h"
#include "bgq.h"
#include "global.h"
#include "boundary.h"
#include <mpi.h>
#include <omp.h>

#define BGQ_HOPPINGMATRIX_C_
#include "bgq_HoppingMatrix.h"



////////////////////////////////////////////////////////////////////////////////
// Weylfields

typedef union {
	_Complex double val;
	struct {
		unsigned t : 8;
		unsigned x : 8;
		unsigned y : 8;
		unsigned z : 8;
		unsigned writes : 8;
		unsigned reads : 8;
		bool isOdd : 1;
		unsigned v : 1; // 0..1
		unsigned c : 2; // 0..2
		bool init : 1;
		// 52 bits = 7 bytes
	} coord;
} bgq_weylcoord;

bgq_weylfield_double weylxchange_recv_double[6];
bgq_weylfield_double weylxchange_send_double[6];
size_t weylxchange_size_double[3];
int weylexchange_destination[6];

#ifndef NDEBUG
bgq_weylcoord *weylxchange_recv_double_debug[PHYSICAL_LP][6];
bgq_weylcoord *weylxchange_send_double_debug[PHYSICAL_LP][6];
#endif



void bgq_hm_init() {
	weylxchange_size_double[TUP / 2] = PHYSICAL_LXV * PHYSICAL_LY * PHYSICAL_LZ * sizeof(bgq_weylsite_double);
	weylxchange_size_double[XUP / 2] = PHYSICAL_LTV * PHYSICAL_LY * PHYSICAL_LZ * sizeof(bgq_weylsite_double);
	weylxchange_size_double[YUP / 2] = PHYSICAL_LTV * PHYSICAL_LX * PHYSICAL_LZ * sizeof(bgq_weylsite_double);

	for (direction d = TUP; d <= YDOWN; d += 1) {
		size_t size = weylxchange_size_double[d/2];
		weylxchange_recv_double[d] = (bgq_weylfield_double) malloc_aligned(size, 128);
		weylxchange_send_double[d] = (bgq_weylfield_double) malloc_aligned(size, 128);
#ifndef NDEBUG
		for (int isOdd = false; isOdd <= true; isOdd += 1) {
			weylxchange_recv_double_debug[isOdd][d] = (bgq_weylcoord*) malloc_aligned(size, 128);
			weylxchange_send_double_debug[isOdd][d] = (bgq_weylcoord*) malloc_aligned(size, 128);
		}
#endif
	}

	weylexchange_destination[TUP] = g_nb_t_up;
	weylexchange_destination[TDOWN] = g_nb_t_dn;
	weylexchange_destination[XUP] = g_nb_x_up;
	weylexchange_destination[XDOWN] = g_nb_x_dn;
	weylexchange_destination[YUP] = g_nb_y_up;
	weylexchange_destination[YDOWN] = g_nb_y_dn;
}

void bgq_hm_free() {
	for (direction d = TUP; d <= YDOWN; d += 1) {
		free(weylxchange_recv_double[d]);
		weylxchange_recv_double[d] = NULL;
		free(weylxchange_send_double[d]);
		weylxchange_send_double[d] = NULL;
#ifndef NDEBUG
		for (int isOdd = false; isOdd <= true; isOdd += 1) {
			free(weylxchange_recv_double_debug[isOdd][d]);
			weylxchange_recv_double_debug[isOdd][d] = NULL;
			free(weylxchange_send_double_debug[isOdd][d]);
			weylxchange_send_double_debug[isOdd][d] = NULL;
		}
#endif
	}
}


#ifndef NDEBUG
static bgq_weylcoord *bgq_weylfield_coordref(bool isOdd, _Complex double *val) {
	assert( sizeof(bgq_weylcoord) == sizeof(_Complex double) );
	// Logically, odd and even sites are different locations
	// But in HoppingMatrix, only one of them is used in a time, odd and even computation/communication do not interleave
	// Therefore we coalesce even and odd fields to the same physical memory location
	// But logically these sites are still distinct, so we need different debug infos for them

	for (direction d = TUP; d <= YDOWN; d += 1) {
		if ( ((char*)weylxchange_recv_double[d] <= (char*)val) && ((char*)val < (char*)weylxchange_recv_double[d] + weylxchange_size_double[d/2]) ) {
			size_t offset = (char*)val - (char*)(weylxchange_recv_double[d]);
			return (bgq_weylcoord*)((char*)weylxchange_recv_double_debug[isOdd][d] + offset);
		}

		if ( ((char*)weylxchange_send_double[d] <= (char*)val) && ((char*)val < (char *)weylxchange_send_double[d] + weylxchange_size_double[d/2]) ) {
			size_t offset = (char*)val - (char*)(weylxchange_send_double[d]);
			return (bgq_weylcoord*)((char*)weylxchange_send_double_debug[isOdd][d] + offset);
		}
	}

	master_error(1, "Unknown Weylfield\n");
	return NULL;
}

static void bgq_weylfield_checkcoord(bgq_weylfield_double weylfield, bool isOdd, int t, int x, int y, int z, int v, int c, _Complex double *value) {
	bgq_weylcoord *coord = bgq_weylfield_coordref(isOdd, value);
	if (coord->coord.init) {
		assert(coord->coord.isOdd == isOdd);
		assert(coord->coord.t == t+1);
		assert(coord->coord.x == x+1);
		assert(coord->coord.y == y+1);
		assert(coord->coord.z == z);
		assert(coord->coord.v == v);
		assert(coord->coord.c == c);
	} else {
		coord->coord.init = true;
		coord->coord.isOdd = isOdd;
		coord->coord.t = t+1;
		coord->coord.x = x+1;
		coord->coord.y = y+1;
		coord->coord.z = z;
		coord->coord.v = v;
		coord->coord.c = c;
		coord->coord.writes = 0;
		coord->coord.reads = 0;
	}
}



void bgq_weylfield_t_resetcoord(bgq_weylfield_double weylfield, int t, bool isOdd, int expected_reads_min, int expected_reads_max, int expected_writes_min, int expected_writes_max) {
	assert(weylfield);

//#pragma omp parallel for schedule(static)
	for (int xyz = 0; xyz < LOCAL_LX*LOCAL_LY*LOCAL_LZ; xyz += 1) {
		WORKLOAD_DECL(xyz, LOCAL_LX*LOCAL_LY*LOCAL_LZ);
		const int x = WORKLOAD_PARAM(LOCAL_LX);
		const int y = WORKLOAD_PARAM(LOCAL_LY);
		const int z = WORKLOAD_PARAM(LOCAL_LZ);
		WORKLOAD_CHECK

		if ( ((t+x+y+z)&1) != isOdd ) // Very BAD for branch predictor
			continue;

		const int xeo = x / PHYSICAL_LP;
		const int xv = xeo / PHYSICAL_LK;
		const int k = mod(xeo, PHYSICAL_LK);

		for (int v = 0; v < 2; v += 1) {
			for (int c = 0; c < 3; c += 1) {
				_Complex double *value = BGQ_WEYLVAL_T(weylfield,isOdd,t,x,y,z,xv,k,v,c,false,false);
				bgq_weylcoord *coord = bgq_weylfield_coordref(isOdd, value);

				assert(coord->coord.init);
				int reads = coord->coord.reads;
				int writes = coord->coord.writes;

				if (expected_reads_min >= 0)
					assert(reads >= expected_reads_min);
				if (expected_reads_max >= 0)
					assert(reads <= expected_reads_max);
				if (expected_writes_min >= 0)
					assert(writes >= expected_writes_min);
				if (expected_writes_max >= 0)
					assert(writes <= expected_writes_max);

				coord->coord.writes = 0;
				coord->coord.reads = 0;
			}
		}
	}
}



bool assert_weylval_t(bgq_weylfield_double weylfield, bool isOdd, int t, int x, int y, int z, int xv, int k, int v, int c, bool isRead, bool isWrite) {
	assert(weylfield);
	assert(false <= isOdd && isOdd <= true); // bogus
	assert( (t == -1) || (t == 0) || (t == LOCAL_LT-1) || (t == LOCAL_LT) ); /* We are one out of the volume, either up or down */
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	assert(0 <= xv && xv < PHYSICAL_LXV);
	assert(0 <= k && k < PHYSICAL_LK);

	// Only 4 fields of this type actually exist
	assert( (weylfield == weylxchange_recv_double[TUP]) || (weylfield == weylxchange_recv_double[TDOWN])
			|| (weylfield == weylxchange_send_double[TUP]) || (weylfield == weylxchange_send_double[TDOWN]) );

	// Check that the coordinate is really an odd/even coordinate
	assert(((t+x+y+z)&1) == isOdd);

	const int xeo = x/PHYSICAL_LP;

	// Check that zv and k match the coordinate
	assert(xeo/PHYSICAL_LK == xv);
	assert(mod(xeo,PHYSICAL_LK) == k);

	// Get the memory address it points to
	const int idx = ((xv)*PHYSICAL_LY + (y))*PHYSICAL_LZ + (z);
	assert(0 <= idx && idx < PHYSICAL_LXV*PHYSICAL_LY*PHYSICAL_LZ);

	// Validate the memory address
	bgq_weylsite_double *site = &weylfield[idx];
	assert(&weylfield[0] <= site && site < &weylfield[PHYSICAL_LXV*PHYSICAL_LY*PHYSICAL_LZ]);
	assert(mod((size_t)&site->s[v][c][0],32)==0);

	_Complex double *val = &site->s[v][c][k];
	assert(&weylfield[0].s[0][0][0] <= val && val < &weylfield[PHYSICAL_LXV*PHYSICAL_LY*PHYSICAL_LZ].s[0][0][0]);

	bgq_weylcoord *coord = bgq_weylfield_coordref(isOdd, val);
	bgq_weylfield_checkcoord(weylfield, isOdd, t, x, y, z, v, c, val);

	assert(coord->coord.init);
	if (isRead)
		coord->coord.reads += 1;
	if (isWrite)
		coord->coord.writes += 1;

	return true;
}


bool assert_weylfield_t(bgq_weylfield_double weylfield, bool isOdd, int t, int x, int y, int z, int xv, int k, bool isRead, bool isWrite) {
	for (int v = 0; v < 2; v += 1) {
		for (int c = 0; c < 3; c += 1) {
			assert_weylval_t(weylfield,isOdd,t,x,y,z,xv,k,v,c,isRead,isWrite);
		}
	}
	return true;
}


void bgq_weylfield_x_resetcoord(bgq_weylfield_double weylfield, int x, bool isOdd, int expected_reads_min, int expected_reads_max, int expected_writes_min, int expected_writes_max) {
	assert(weylfield);

//#pragma omp parallel for schedule(static)
	for (int tyz = 0; tyz < LOCAL_LT*LOCAL_LY*LOCAL_LZ; tyz += 1) {
		WORKLOAD_DECL(tyz, LOCAL_LT*LOCAL_LY*LOCAL_LZ);
		const int t = WORKLOAD_PARAM(LOCAL_LT);
		const int y = WORKLOAD_PARAM(LOCAL_LY);
		const int z = WORKLOAD_PARAM(LOCAL_LZ);
		WORKLOAD_CHECK

		if ( ((t+x+y+z)&1) != isOdd ) // Very BAD for branch predictor
			continue;

		const int teo = t / PHYSICAL_LP;
		const int tv = teo / PHYSICAL_LK;
		const int k = mod(teo, PHYSICAL_LK);

		for (int v = 0; v < 2; v += 1) {
			for (int c = 0; c < 3; c += 1) {
				_Complex double *value = BGQ_WEYLVAL_X(weylfield,isOdd,t,x,y,z,tv,k,v,c,false,false);
				bgq_weylcoord *coord = bgq_weylfield_coordref(isOdd, value);

				assert(coord->coord.init);
				int reads = coord->coord.reads;
				int writes = coord->coord.writes;

				if (expected_reads_min >= 0)
					assert(reads >= expected_reads_min);
				if (expected_reads_max >= 0)
					assert(reads <= expected_reads_max);
				if (expected_writes_min >= 0)
					assert(writes >= expected_writes_min);
				if (expected_writes_max >= 0)
					assert(writes <= expected_writes_max);

				coord->coord.writes = 0;
				coord->coord.reads = 0;
			}
		}
	}
}



bool assert_weylval_x(bgq_weylfield_double weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k, int v, int c, bool isRead, bool isWrite) {
	assert(weylfield);
	assert(false <= isOdd && isOdd <= true);
	assert(0 <= t && t < LOCAL_LT);
	assert( (x == -1) || (x == 0) || (x == LOCAL_LX-1) || (x == LOCAL_LX) );
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= k && k < PHYSICAL_LK);

	// Only 4 fields of this type actually exist
	assert((weylfield == weylxchange_recv_double[XUP]) || (weylfield == weylxchange_recv_double[XDOWN])
			|| (weylfield == weylxchange_send_double[XUP]) || (weylfield == weylxchange_send_double[XDOWN]) );

	// Check that the coordinate is really an odd/even coordinate
	assert(((t+x+y+z)&1) == isOdd);

	const int teo = t/PHYSICAL_LP;

	// Check that zv and k match the coordinate
	assert(teo/PHYSICAL_LK == tv);
	assert(mod(teo,PHYSICAL_LK) == k);

	// Get the memory address it points to
	const int idx = (tv*PHYSICAL_LY + y)*PHYSICAL_LZ + z;
	assert(0 <= idx && idx < PHYSICAL_LTV*PHYSICAL_LY*PHYSICAL_LZ);

	// Validate the memory address
	bgq_weylsite_double *site = &weylfield[idx];
	assert(&weylfield[0] <= site && site < &weylfield[PHYSICAL_LTV*PHYSICAL_LY*PHYSICAL_LZ]);
	assert(mod((size_t)&site->s[v][c][0],32)==0);

	_Complex double *val = &site->s[v][c][k];
	assert(&weylfield[0].s[0][0][0] <= val && val < &weylfield[PHYSICAL_LTV*PHYSICAL_LY*PHYSICAL_LZ].s[0][0][0]);

	bgq_weylcoord *coord = bgq_weylfield_coordref(isOdd, val);
	bgq_weylfield_checkcoord(weylfield, isOdd, t, x, y, z, v, c, val);

	assert(coord->coord.init);
	if (isRead)
		coord->coord.reads += 1;
	if (isWrite)
		coord->coord.writes += 1;

	return true;
}


bool assert_weylfield_x(bgq_weylfield_double weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k, bool isRead, bool isWrite) {
	for (int v = 0; v < 2; v += 1) {
		for (int c = 0; c < 3; c += 1) {
			assert_weylval_x(weylfield,isOdd,t,x,y,z,tv,k,v,c,isRead,isWrite);
		}
	}
	return true;
}



void bgq_weylfield_y_resetcoord(bgq_weylfield_double weylfield, int y, bool isOdd, int expected_reads_min, int expected_reads_max, int expected_writes_min, int expected_writes_max) {
	assert(weylfield);

//#pragma omp parallel for schedule(static)
	for (int txz = 0; txz < LOCAL_LT*LOCAL_LX*LOCAL_LZ; txz += 1) {
		WORKLOAD_DECL(txz, LOCAL_LT*LOCAL_LX*LOCAL_LZ);
		const int t = WORKLOAD_PARAM(LOCAL_LT);
		const int x = WORKLOAD_PARAM(LOCAL_LX);
		const int z = WORKLOAD_PARAM(LOCAL_LZ);
		WORKLOAD_CHECK

		if ( ((t+x+y+z)&1) != isOdd ) // Very BAD for branch predictor
			continue;

		const int teo = t / PHYSICAL_LP;
		const int tv = teo / PHYSICAL_LK;
		const int k = mod(teo, PHYSICAL_LK);

		for (int v = 0; v < 2; v += 1) {
			for (int c = 0; c < 3; c += 1) {
				_Complex double *value = BGQ_WEYLVAL_Y(weylfield,isOdd,t,x,y,z,tv,k,v,c,false,false);
				bgq_weylcoord *coord = bgq_weylfield_coordref(isOdd, value);

				assert(coord->coord.init);
				int reads = coord->coord.reads;
				int writes = coord->coord.writes;

				if (expected_reads_min >= 0)
					assert(reads >= expected_reads_min);
				if (expected_reads_max >= 0)
					assert(reads <= expected_reads_max);
				if (expected_writes_min >= 0)
					assert(writes >= expected_writes_min);
				if (expected_writes_max >= 0)
					assert(writes <= expected_writes_max);

				coord->coord.writes = 0;
				coord->coord.reads = 0;
			}
		}
	}
}


bool assert_weylval_y(bgq_weylfield_double weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k, int v, int c, bool isRead, bool isWrite) {
	assert(weylfield);
	assert(false <= isOdd && isOdd <= true);
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert( (-1 == y) || (y == 0) || (y == LOCAL_LY-1) || (y == LOCAL_LY) );
	assert(0 <= z && z < LOCAL_LZ);
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= k && k < PHYSICAL_LK);

	// Only 4 fields of this type actually exist
	assert((weylfield == weylxchange_recv_double[YUP]) || (weylfield == weylxchange_recv_double[YDOWN])
			|| (weylfield == weylxchange_send_double[YUP]) || (weylfield == weylxchange_send_double[YDOWN]));

	// Check that the coordinate is really an odd/even coordinate
	assert(((t+x+y+z)&1) == isOdd);
	const int teo = t/PHYSICAL_LP;

	// Check that zv and k match the coordinate
	assert(teo/PHYSICAL_LK == tv);
	assert(mod(teo,PHYSICAL_LK) == k);

	// Get the memory address it points to
	const int idx = (tv*PHYSICAL_LX + x)*PHYSICAL_LZ + z;
	assert(0 <= idx && idx < PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LZ);

	// Validate the memory address
	bgq_weylsite_double *site = &weylfield[idx];
	assert(&weylfield[0] <= site && site < &weylfield[PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LZ]);
	assert(mod((size_t)&site->s[v][c][0],32)==0);

	_Complex double *val = &site->s[v][c][k];
	assert(&weylfield[0].s[0][0][0] <= val && val < &weylfield[PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LZ].s[0][0][0]);

	bgq_weylcoord *coord = bgq_weylfield_coordref(isOdd, val);
	bgq_weylfield_checkcoord(weylfield, isOdd, t, x, y, z, v, c, val);

	assert(coord->coord.init);
	if (isRead)
		coord->coord.reads += 1;
	if (isWrite)
		coord->coord.writes += 1;

	return true;
}

bool assert_weylfield_y(bgq_weylfield_double weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k, bool isRead, bool isWrite) {
	for (int v = 0; v < 2; v += 1) {
		for (int c = 0; c < 3; c += 1) {
			assert_weylval_y(weylfield,isOdd,t,x,y,z,tv,k,v,c,isRead,isWrite);
		}
	}
	return true;
}
#endif

////////////////////////////////////////////////////////////////////////////////







#define BGQ_HM_BORDER_NOFUNC 1
#define BGQ_HM_BORDERDIST_NOFUNC 1

//#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
//#pragma GCC diagnostic ignored "-Wunused-variable"


#define HoppingMatrix bgq_HoppingMatrix_double
#include "bgq_HoppingMatrix.inc.c"
#undef HoppingMatrix


#define BGQ_HM_NOCOM 1
#define HoppingMatrix bgq_HoppingMatrix_double_nocom
#include "bgq_HoppingMatrix.inc.c"
#undef HoppingMatrix

// Hopping_Matrix.c compatibility layer
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k) {
	const bool isOdd = (ieo!=0);
	bgq_spinorfield_double target = bgq_translate_spinorfield(l);
	bgq_spinorfield_double source = bgq_translate_spinorfield(k);

	bgq_HoppingMatrix_double(isOdd, target, source, g_gaugefield_double);
}


void Hopping_Matrix_nocom(const int ieo, spinor * const l, spinor * const k) {
	const bool isOdd = (ieo!=0);
	bgq_spinorfield_double target = bgq_translate_spinorfield(l);
	bgq_spinorfield_double source = bgq_translate_spinorfield(k);

	bgq_HoppingMatrix_double_nocom(isOdd, target, source, g_gaugefield_double);
}


