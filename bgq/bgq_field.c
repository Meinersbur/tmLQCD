
#define BGQ_FIELD_C_H_
#include "bgq_field.h"

#include "bgq.h"
#include "../geometry_eo.h"
#include "../global.h"

#include <string.h>
#include <stdbool.h>

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


////////////////////////////////////////////////////////////////////////////////
// Spinorfield

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
		unsigned v : 2; // 0..3
		unsigned c : 2; // 0..2
		bool init : 1;
		// 52 bits = 7 bytes
	} coord;
} bgq_spinorcoord;
bgq_spinorcoord bgq_spinorcoord_encode(bool isOdd, int t, int x, int y, int z, int v, int c) {
	bgq_spinorcoord result;
	assert(sizeof(result) == sizeof(_Complex double));
	result.val = 0; // Write zeros

	result.coord.isOdd = isOdd;
	result.coord.t = t;
	result.coord.x = x;
	result.coord.y = y;
	result.coord.z = z;
	result.coord.v = v;
	result.coord.c = c;
	result.coord.writes = 0;
	result.coord.reads = 0;

	return result;
}


static bgq_spinorcoord *bgq_spinorfield_coordref(_Complex double *site) {
	size_t offset = (char*)site - (char*)g_spinorfields_doubledata;
	bgq_spinorcoord *result = (bgq_spinorcoord*)((char*)g_spinorfields_doubledata_coords + offset);
	return result;
}


static void bgq_spinorfield_checkcoord(bgq_spinorfield_double spinorfield, bool isOdd, int t, int x, int y, int z, int tv, int k, int v, int c, _Complex double *value) {
	bgq_spinorcoord *coord = bgq_spinorfield_coordref(value);
	if (coord->coord.init) {
		assert(coord->coord.isOdd == isOdd);
		assert(coord->coord.t == t);
		assert(coord->coord.x == x);
		assert(coord->coord.y == y);
		assert(coord->coord.z == z);
		assert(coord->coord.v == v);
		assert(coord->coord.c == c);
	} else {
		coord->coord.init = true;
		coord->coord.isOdd = isOdd;
		coord->coord.t = t;
		coord->coord.x = x;
		coord->coord.y = y;
		coord->coord.z = z;
		coord->coord.v = v;
		coord->coord.c = c;
		coord->coord.writes = 0;
		coord->coord.reads = 0;
	}
}


void bgq_spinorfield_resetcoord(bgq_spinorfield_double spinorfield, bool isOdd, int expected_reads_min, int expected_reads_max, int expected_writes_min, int expected_writes_max) {
	assert(spinorfield);

#pragma omp parallel for schedule(static)
	for (int txyz = 0; txyz < VOLUME; txyz += 1) {
		WORKLOAD_DECL(txyz, VOLUME);
		const int t = WORKLOAD_PARAM(LOCAL_LT);
		const int x = WORKLOAD_PARAM(LOCAL_LX);
		const int y = WORKLOAD_PARAM(LOCAL_LY);
		const int z = WORKLOAD_PARAM(LOCAL_LZ);
		WORKLOAD_CHECK

		if ( ((t+x+y+z)&1) != isOdd ) // Very BAD for branch predictor
			continue;

		const int teo = t / PHYSICAL_LP;
		const int tv = teo / PHYSICAL_LK;
		const int k = mod(teo, PHYSICAL_LK);

		for (int v = 0; v < 4; v += 1) {
			for (int c = 0; c < 3; c += 1) {
				_Complex double *value = BGQ_SPINORVAL(spinorfield,isOdd,t,x,y,z,tv,k,v,c,false,false);
				bgq_spinorcoord *coord = bgq_spinorfield_coordref(value);

				assert (coord->coord.init);
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

				//bgq_spinorfield_checkcoord(spinorfield,isOdd,t,x,y,z,tv,k,v,c,value);
				coord->coord.writes = 0;
				coord->coord.reads = 0;
			}
		}
	}
}



void bgq_init_spinorfields(int count) {
	g_num_spinorfields = count;
	int datasize = count * sizeof(*g_spinorfields_doubledata) * VOLUME_SITES;
	g_spinorfields_doubledata = malloc_aligned(datasize, 128);
#ifndef NDEBUG
	g_spinorfields_doubledata_coords = malloc_aligned(datasize, 128);
#endif

	g_spinorfields_double = malloc(count * sizeof(*g_spinorfields_double));
#ifndef NDEBUG
	g_spinorfield_isOdd = malloc(count * sizeof(*g_spinorfield_isOdd));
#endif
	for (int i = 0; i < count; i += 1) {
		g_spinorfields_double[i] = g_spinorfields_doubledata + i * VOLUME_SITES;
#ifndef NDEBUG
		g_spinorfield_isOdd[i] = -1; // Unknown yet
#endif
	}
}


static void bgq_reset_spinorfield(bgq_spinorfield_double spinorfield, bool isOdd) {
	// Fill the fields with the coordinates
	for (int txyz = 0; txyz < VOLUME; txyz += 1) {
		WORKLOAD_DECL(txyz, GAUGE_VOLUME);
		const int t = WORKLOAD_PARAM(LOCAL_LT);
		const int x = WORKLOAD_PARAM(LOCAL_LX);
		const int y = WORKLOAD_PARAM(LOCAL_LY);
		const int z = WORKLOAD_PARAM(LOCAL_LZ);
		WORKLOAD_CHECK

		if (((t+x+y+z)&1) != isOdd)
			continue;

		for (int v = 0; v < 4; v+=1)
			for (int c = 0; c < 3; c += 1) {
				_Complex double *val = bgq_spinorfield_double_ref(spinorfield, isOdd, t,x,y,z,v,c,false,false);
				bgq_spinorcoord *coord = bgq_spinorfield_coordref(val);
				*coord = bgq_spinorcoord_encode(isOdd,t,x,y,z,v,c);
			}
	}
}


static void bgq_check_spinorfield(bgq_spinorfield_double spinorfield, int expected_reads, int expected_writes) {
	// Find out whether it was an even or odd field part
	// Get the index of the used spinorfield
	const int fieldsize = sizeof(bgq_spinorsite_double) * VOLUME_SITES; // alignment???
	long long offset = (char*)spinorfield - (char*)g_spinorfields_doubledata;
	assert(offset >= 0);
	assert(mod(offset, fieldsize) == 0);
	int index = offset / fieldsize;
	assert(index < g_num_spinorfields);

	assert(g_spinorfield_isOdd[index] != -1);
	const bool isOdd = (g_spinorfield_isOdd[index] != 0);

	for (int txyz = 0; txyz < VOLUME; txyz += 1) {
		WORKLOAD_DECL(txyz, GAUGE_VOLUME);
		const int t = WORKLOAD_PARAM(LOCAL_LT);
		const int x = WORKLOAD_PARAM(LOCAL_LX);
		const int y = WORKLOAD_PARAM(LOCAL_LY);
		const int z = WORKLOAD_PARAM(LOCAL_LZ);
		WORKLOAD_CHECK

		if (((t+x+y+z)&1) != isOdd)
			continue;

		const int teo = t / PHYSICAL_LP;
		const int tv = teo/PHYSICAL_LK;
		const int k = mod(teo,PHYSICAL_LK);

		for (int v = 0; v < 4; v+=1)
			for (int c = 0; c < 3; c += 1) {
				assert(assert_spinorfield_coord(spinorfield, isOdd, t,x,y,z,tv,k,v,c,false,false));
				_Complex double *site = bgq_spinorfield_double_ref(spinorfield, isOdd, t,x,y,z,v,c,false,false);
				bgq_spinorcoord *coord = bgq_spinorfield_coordref(site);

				if (expected_reads >= 0)
					assert(coord->coord.reads == expected_reads);
				if (expected_writes >= 0)
					assert(coord->coord.writes == expected_writes);
			}
	}
}




void bgq_free_spinofields() {
	free(g_spinorfields_double);
	g_spinorfields_double = NULL;
	free(g_spinorfields_doubledata);
	g_spinorfields_doubledata = NULL;
	g_num_spinorfields = 0;
	free(g_spinorfield_isOdd);
	g_spinorfield_isOdd = NULL;
}


typedef struct {
	double _Complex c[3];
} su3_vector_array64;

typedef struct {
	su3_vector_array64 v[4];
} spinor_array64;



void bgq_transfer_spinorfield(const bool isOdd, bgq_spinorfield_double const targetfield, spinor * const sourcefield) {
	assert(sourcefield);
	assert(targetfield);

#ifndef NDEBUG
	bgq_spinorfield_resetcoord(targetfield, isOdd, -1, -1, -1, -1);
#endif

//#pragma omp parallel for schedule(static)
	for (int txyz = 0; txyz < VOLUME; txyz+=1) {
		WORKLOAD_DECL(txyz, VOLUME);
		const int t = WORKLOAD_CHUNK(LOCAL_LT);
		const int x = WORKLOAD_CHUNK(LOCAL_LX);
		const int y = WORKLOAD_CHUNK(LOCAL_LY);
		const int z = WORKLOAD_CHUNK(LOCAL_LZ);
		WORKLOAD_CHECK

		if (((t+x+y+z)&1)!=isOdd)
			continue;

		const int ix = g_ipt[t][x][y][z]; /* lexic coordinate */
		assert(ix == Index(t,x,y,z));

		int iy = g_lexic2eo[ix]; /* even/odd coordinate (even and odd sites in two different fields of size VOLUME/2, first even field followed by odd) */
		assert(0 <= iy && iy < (VOLUME+RAND));
		int icx = g_lexic2eosub[ix]; /*  even/odd coordinate relative to field base */
		assert(0 <= icx && icx < VOLUME/2);
		assert(icx == iy - (isOdd ? (VOLUME+RAND)/2 : 0));

		spinor_array64 *sp = (spinor_array64*)&sourcefield[icx];
		for (int v = 0; v < 4; v+=1) {
			for (int c = 0; c < 3; c+=1) {
				bgq_spinorfield_double_set(targetfield,isOdd,t,x,y,z,v,c,sp->v[v].c[c]);
			}
		}
	}

#ifndef NDEBUG
	bgq_spinorfield_resetcoord(targetfield, isOdd, 0, 0, 1, 1);
#endif
}



bool assert_spinorfield_coord(bgq_spinorfield_double spinorfield, bool isOdd, int t, int x, int y, int z, int tv, int k, int v, int c, bool isRead, bool isWrite) {
	assert(spinorfield);
	assert(false <= isOdd && isOdd <= true);
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= k && k < PHYSICAL_LK);

	// Get the index of the used spinorfield
	const int fieldsize = sizeof(bgq_spinorsite_double) * VOLUME_SITES; // alignment???
	long long offset = (char*)spinorfield - (char*)g_spinorfields_doubledata;
	assert(offset >= 0);
	assert(offset % fieldsize == 0);
	int index = offset / fieldsize;
	assert(index < g_num_spinorfields);

	// Check that the coordinate is really an odd/even coordinate
	assert(((t+x+y+z)&1) == isOdd);

	// Check that field is used as an odd/even field
	if (g_spinorfield_isOdd[index] == -1) {
		// not yet defined, just ensure that at all following uses are the same
		g_spinorfield_isOdd[index] = isOdd;
	} else {
		assert(g_spinorfield_isOdd[index] == isOdd);
	}

	const int teo = t/PHYSICAL_LP;

	// Check that zv and k match the coordinate
	assert(teo/PHYSICAL_LK == tv);
	assert(mod(teo,PHYSICAL_LK) == k);

	// Get the memory address it points to
	const int idx = ((tv*PHYSICAL_LX + x)*PHYSICAL_LY + y)*PHYSICAL_LZ + z;
	assert(0 <= idx && idx < PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LY*PHYSICAL_LZ);

	// Validate the address
	bgq_spinorsite_double *site = &spinorfield[idx];
	assert(g_spinorfields_double[index] <= site && site < g_spinorfields_double[index+1]);
	assert(mod((size_t)site,32)==0);

	_Complex double *val = &site->s[v][c][k];
	assert( (&g_spinorfields_double[index]->s[0][0][0] <= val) && (val < &g_spinorfields_double[index+1]->s[0][0][0]) );
	assert(mod((size_t)site,16)==0);
	assert(mod((size_t)&site->s[v][c][0],32)==0);

	// Get the equivalent address in debug address space, and check it
	bgq_spinorfield_checkcoord(spinorfield,isOdd,t,x,y,z,tv,k,v,c,val);
	bgq_spinorcoord *coord = bgq_spinorfield_coordref(val);

	// Mark uses of this value
	if (isRead)
		coord->coord.reads += 1;
	if (isWrite)
		coord->coord.writes += 1;

	return true; // All checks passed
}

bool assert_spinorcoord(bgq_spinorfield_double spinorfield, bool isOdd, int t, int x, int y, int z, int tv, int k, bool isRead, bool isWrite) {
	for (int v = 0; v < 4; v += 1) {
			for (int c = 0; c < 3; c += 1) {
				assert_spinorfield_coord(spinorfield, isOdd, t, x, y, z, tv, k, v, c, isRead, isWrite);
			}
	}

	return true;
}



////////////////////////////////////////////////////////////////////////////////
// Gaugefield

typedef
	union {
		_Complex double val;
		struct {
			unsigned t : 8;
			unsigned x : 8;
			unsigned y : 8;
			unsigned z : 8;
			unsigned writes : 8;
			unsigned reads : 8;
			unsigned dir : 4; // T_UP..ZUP
			unsigned i : 2; // 0..3
			unsigned l : 2; // 0..3
			bool isOdd : 1;
			bool init : 1;
			// 58 bits = 8 bytes
		} coord;
} bgq_gaugecoord;


static bgq_gaugecoord *bgq_gaugefield_coordref(_Complex double *site) {
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (direction dir = TUP; dir <= TUP_SHIFT; dir += 2) {
			bgq_gaugesite_double *eodata = g_gaugefield_doubledata.eodir[isOdd][dir/2];
			if ((void*)eodata <= (void*)site && (void*)site < (void*)(eodata + GAUGE_EOVOLUME)) {
				// found!, get equivialent debug site
				size_t offset = (char*)site - (char*)eodata;
				char *debugdata = (char*)g_gaugefield_doubledata_debug.eodir[isOdd][dir/2];
				bgq_gaugecoord *result = (bgq_gaugecoord*)(debugdata + offset);
				return result;
			}
		}
	}

	assert(!"Pointer does not point to gaugefield");
	return NULL;
}


static void bgq_gaugefield_checkcoord(bgq_gaugefield_double gaugefield, bool isOdd, int t, int x, int y, int z, int tv, int k, direction dir, int i, int l, _Complex double *value) {
	bgq_gaugecoord *coord = bgq_gaugefield_coordref(value);
	if (coord->coord.init) {
		assert(coord->coord.isOdd == isOdd);
		assert(coord->coord.t == t+1);
		assert(coord->coord.x == x+1);
		assert(coord->coord.y == y+1);
		assert(coord->coord.z == mod(z,LOCAL_LZ));
		assert(coord->coord.dir == dir);
		assert(coord->coord.i == i);
		assert(coord->coord.l == l);
	} else {
		coord->coord.init = true;
		coord->coord.isOdd = isOdd;
		coord->coord.t = t+1;
		coord->coord.x = x+1;
		coord->coord.y = y+1;
		coord->coord.z = mod(z,LOCAL_LZ);
		coord->coord.dir = dir;
		coord->coord.i = i;
		coord->coord.l = l;
		coord->coord.writes = 0;
		coord->coord.reads = 0;
	}
}

static void bgq_gaugefield_resetcoord_checkval(bgq_gaugefield_double gaugefield, bool isOdd, int t, int x, int y, int z, int tv, int k, direction dir, int i, int l,
		_Complex double *val,
		int expected_reads, int expected_writes) {

	bgq_gaugecoord *coord = bgq_gaugefield_coordref(val);
	int reads = 0;
	int writes = 0;
	if (coord->coord.init) {
		reads = coord->coord.reads;
		writes = coord->coord.writes;
	}
	if (expected_reads >= 0)
		assert(reads == reads);
	if (expected_writes >= 0)
		assert(writes == writes);

	bgq_gaugefield_checkcoord(gaugefield,isOdd,t,x,y,z,tv,k,dir,i,l,val);
	coord->coord.writes = 0;
	coord->coord.reads = 0;
}

void bgq_gaugefield_resetcoord(bgq_gaugefield_double gaugefield, int expected_reads, int expected_writes) {
	assert(gaugefield);

#pragma omp parallel for schedule(static)
	for (int txyz = 0; txyz < GAUGE_VOLUME; txyz += 1) {
		WORKLOAD_DECL(txyz, GAUGE_VOLUME);
		const int t = WORKLOAD_PARAM(LOCAL_LT+1) - 1;
		const int x = WORKLOAD_PARAM(LOCAL_LX+1) - 1;
		const int y = WORKLOAD_PARAM(LOCAL_LY+1) - 1;
		const int z = WORKLOAD_PARAM(LOCAL_LZ);
		WORKLOAD_CHECK

		const bool isOdd = (t+x+y+z)&1;
		const int teo = (t+1) / PHYSICAL_LP;
		const int tv = teo / PHYSICAL_LK;
		const int k = mod(teo, PHYSICAL_LK);

		for (direction dir = TUP; dir <= ZUP; dir += 2) {
			// Overflow is only needed for TDOWN, XDOWN, YDOWN into their dimension
			if ((t == -1) && (dir != TUP))
				continue;
			if ((x == -1) && (dir != XUP))
				continue;
			if ((y == -1) && (dir != YUP))
				continue;

			for (int i = 0; i < 3; i += 1) {
				for (int l = 0; l < 3; l += 1) {
					{
						_Complex double *value = BGQ_GAUGEVAL(gaugefield,isOdd,t,x,y,z,tv,k, dir,i,l,false,false);
						bgq_gaugefield_resetcoord_checkval(gaugefield, isOdd, t,x,y,z, tv,k, dir, i, l, value, expected_reads, expected_writes);
					}

					if (dir == ZUP && z == LOCAL_LZ-1) {
						_Complex double *wrapvalue = BGQ_GAUGEVAL(gaugefield,isOdd,t,x,y,-1,tv,k, dir,i,l,false,false);
						bgq_gaugefield_resetcoord_checkval(gaugefield, isOdd, t,x,y,-1, tv,k, dir, i, l, wrapvalue, expected_reads, expected_writes);
					}

					if (dir == TUP) {
						const int teo_shift = mod(teo+1, 1+LOCAL_LT/PHYSICAL_LP);
						const int tv_shift = teo_shift / PHYSICAL_LK;
						const int k_shift = mod(teo_shift, PHYSICAL_LK);

						_Complex double *shiftvalue = BGQ_GAUGEVAL(gaugefield, isOdd, t, x, y, z, tv_shift, k_shift, TUP_SHIFT, i, l, false, false);
						bgq_gaugefield_resetcoord_checkval(gaugefield, isOdd, t, x, y, z, tv_shift, k_shift, TUP_SHIFT, i, l, shiftvalue, expected_reads, expected_writes);
					}
				}
			}
		}
	}
}


void bgq_init_gaugefield() {
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (direction dir = TUP; dir <= TUP_SHIFT; dir += 2) {
			g_gaugefield_doubledata.eodir[isOdd][dir/2] = malloc_aligned(GAUGE_EOVOLUME * sizeof(bgq_gaugesite_double), 128 /*L2 cache line size*/);
		}
	}

	g_gaugefield_double = &g_gaugefield_doubledata;

#ifndef NDEBUG
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (direction dir = TUP; dir <= TUP_SHIFT; dir += 2) {
			g_gaugefield_doubledata_debug.eodir[isOdd][dir/2] = malloc_aligned(GAUGE_EOVOLUME * sizeof(bgq_gaugesite_double), 128 /*L2 cache line size*/);
		}
	}

	bgq_gaugefield_resetcoord(g_gaugefield_double, -1, -1);
#endif
}




void bgq_free_gaugefield() {
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (int dir = TUP; dir <= TUP_SHIFT; dir += 2) {
			free(g_gaugefield_doubledata.eodir[isOdd][dir/2]);
			g_gaugefield_doubledata.eodir[isOdd][dir/2] = NULL;
		}
	}
	g_gaugefield_double = NULL;

#ifndef NDEBUG
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (int dir = TUP; dir <= TUP_SHIFT; dir += 2) {
			free(g_gaugefield_doubledata_debug.eodir[isOdd][dir/2]);
			g_gaugefield_doubledata_debug.eodir[isOdd][dir/2] = NULL;
		}
	}
#endif
}
typedef struct {
	double _Complex c[3][3];
} su3_array64;

void bgq_transfer_gaugefield(bgq_gaugefield_double const targetfield, su3 ** const sourcefield) {
	assert(targetfield);
	assert(sourcefield);

#ifndef NDEBUG
	bgq_gaugefield_resetcoord(targetfield, -1, -1);
#endif

#pragma omp parallel for schedule(static)
	for (int txyz = 0; txyz < GAUGE_VOLUME; txyz += 1) {
		WORKLOAD_DECL(txyz, GAUGE_VOLUME);
		const int t = WORKLOAD_PARAM(LOCAL_LT+1) - 1;
		const int x = WORKLOAD_PARAM(LOCAL_LX+1) - 1;
		const int y = WORKLOAD_PARAM(LOCAL_LY+1) - 1;
		const int z = WORKLOAD_PARAM(LOCAL_LZ);
		WORKLOAD_CHECK

		const bool isOdd = (t+x+y+z)&1;
		const int ix = Index(t, x, y, z); /* lexic coordinate; g_ipt[t][x][y][z] is not defined for -1 coordinates */

		for (direction dir = TUP; dir <= ZUP; dir += 2) {
			// Overflow is only needed for TDOWN, XDOWN, YDOWN into their dimension
			if ((t == -1) && (dir != TUP))
				continue;
			if ((x == -1) && (dir != XUP))
				continue;
			if ((y == -1) && (dir != YUP))
				continue;

			su3_array64 *m = (su3_array64*) &g_gauge_field[ix][dir / 2];
			for (int i = 0; i < 3; i += 1) {
				for (int l = 0; l < 3; l += 1) {
					bgq_gaugefield_double_set(targetfield, isOdd, t, x, y, z, dir, i, l, m->c[i][l]);
				}
			}
		}
	}

#ifndef NDEBUG
	bgq_gaugefield_resetcoord(targetfield, 0, 1);
#endif
}


bool assert_gaugesite(bgq_gaugefield_double gaugefield, bool isOdd, int t, int x, int y, int z, int tv, int k, direction dir, bool isRead, bool isWrite) {
	for (int i = 0; i < 3; i += 1) {
		for (int l = 0; l < 3; l += 1) {
			assert_gaugeval(gaugefield, isOdd, t, x, y, z, tv, k, dir, i, l, isRead, isWrite);
		}
	}

	return true;
}


bool assert_gaugeval(bgq_gaugefield_double gaugefield, bool isOdd, int t, int x, int y, int z, int tv, int k, direction dir, int i, int l, bool isRead, bool isWrite) {
	assert(gaugefield);
	assert(false <= isOdd && isOdd <= true);
	assert(-1 <= t && t < LOCAL_LT);
	assert(-1 <= x && x < LOCAL_LX);
	assert(-1 <= y && y < LOCAL_LY);
	assert(-1 <= z && z < PHYSICAL_LZ);
	assert(0 <= tv && tv <= PHYSICAL_LTV);
	assert(0 <= k && k < PHYSICAL_LK);
	assert(TUP <= dir && dir <= TDOWN_SHIFT);
	assert(0 <= i && i < 3);
	assert(0 <= l && l < 3);

	// There is really just one gaugefield
	assert(gaugefield == g_gaugefield_double);

	// We really only store the up-fields
	assert((dir&1)==0);

	// Check that the coordinate is really an odd/even coordinate
	assert(((t+x+y+z)&1) == isOdd);
	int teo = (t+1)/PHYSICAL_LP; // Because t=-1 is a valid index, shift everything right
	if (dir == TUP_SHIFT) {
		teo = mod(teo+1, 1+LOCAL_LT/PHYSICAL_LP);
	}

	// Check that zv and k match the coordinate
	assert(teo/PHYSICAL_LK == tv);
	assert(mod(teo,PHYSICAL_LK) == k);

	// Get the index
	const int idx = (((tv)*(PHYSICAL_LX+1) + ((x)+1))*(PHYSICAL_LY+1) + ((y)+1))*(PHYSICAL_LZ+1) + ((z)+1);
	assert((0 <= idx) && (idx < GAUGE_EOVOLUME));

	// Get the address
	bgq_gaugesite_double *eofield = gaugefield->eodir[isOdd][dir/2];
	bgq_gaugesite_double *site = &eofield[idx];
	_Complex double *address = &site->c[i][l][k];
	assert(mod((size_t)&site->c[i][l][0],32)==0);


	// get the debug data
	bgq_gaugefield_checkcoord(gaugefield,isOdd,t,x,y,mod(z,LOCAL_LZ),tv,k,dir,i,l,address);
	bgq_gaugecoord *debugdata = bgq_gaugefield_coordref(address);
	if (isRead)
		debugdata->coord.reads += 1;
	if (isWrite)
		debugdata->coord.writes += 1;


	// All checks passed
	return true;
}



////////////////////////////////////////////////////////////////////////////////
// Weylfields


extern bgq_weylfield_double weylxchange_recv_double[6];
extern bgq_weylfield_double weylxchange_send_double[6];

bool assert_weylfield_t(bgq_weylfield_double weylfield, bool isOdd, int t, int x, int y, int z, int xv, int k) {
	assert(weylfield);
	assert(false <= isOdd && isOdd <= true);
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	assert(0 <= xv && xv < PHYSICAL_LXV);
	assert(0 <= k && k < PHYSICAL_LK);

	// Only 4 fields of this type actually exist
	assert((weylfield == weylxchange_recv_double[TUP]) || (weylfield == weylxchange_recv_double[TDOWN])
			|| (weylfield == weylxchange_send_double[TUP]) || (weylfield == weylxchange_send_double[TDOWN]) );

	// Check that the coordinate is really an odd/even coordinate
	assert(((t+x+y+z)&1) == isOdd);

	const int xeo = x/PHYSICAL_LP;

	// Check that zv and k match the coordinate
	assert(xeo/PHYSICAL_LK == xv);
	assert(mod(xeo,PHYSICAL_LK) == k);

	// Get the memory address it points to
	const int idx = (xv*PHYSICAL_LY + y)*PHYSICAL_LZ + z;
	assert(0 <= idx && idx < PHYSICAL_LXV*PHYSICAL_LY*PHYSICAL_LZ);

	// Validate the memory address
	bgq_weylsite_double *address = &weylfield[idx];
	assert(&weylfield[0] <= address && address < &weylfield[PHYSICAL_LXV*PHYSICAL_LY*PHYSICAL_LZ]);
	assert(mod((size_t)address,32)==0);

	return true;
}

bool assert_weylfield_x(bgq_weylfield_double weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k){
	assert(weylfield);
	assert(false <= isOdd && isOdd <= true);
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
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
	bgq_weylsite_double *address = &weylfield[idx];
	assert(&weylfield[0] <= address && address < &weylfield[PHYSICAL_LTV*PHYSICAL_LY*PHYSICAL_LZ]);
	assert(mod((size_t)address,32)==0);

	return true;
}

bool assert_weylfield_y(bgq_weylfield_double weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k) {
	assert(weylfield);
	assert(false <= isOdd && isOdd <= true);
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
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
	bgq_weylsite_double *address = &weylfield[idx];
	assert(&weylfield[0] <= address && address < &weylfield[PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LZ]);
	assert(mod((size_t)address,32)==0);

	return true;
}

