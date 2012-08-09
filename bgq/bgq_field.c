
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
		// 51 bits = 7 bytes
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






void bgq_init_spinorfields(int count) {
	g_num_spinorfields = count;
	int datasize = count * sizeof(*g_spinorfields_doubledata) * VOLUME_SITES;
	g_spinorfields_doubledata = malloc_aligned(datasize, 128);
#ifndef NDEBUG
	g_spinorfields_doubledata_coords = malloc_aligned(datasize, 128);
#endif

	g_spinorfields_double = malloc(count * sizeof(*g_spinorfields_double));
	g_spinorfield_isOdd = malloc(count * sizeof(*g_spinorfield_isOdd));
	for (int i = 0; i < count; i += 1) {
		g_spinorfields_double[i] = g_spinorfields_doubledata + i * VOLUME_SITES;
		g_spinorfield_isOdd[i] = -1; // Unknown yet
	}
}

static bgq_spinorcoord *bgq_spinorfield_coordref(_Complex double *site) {
	size_t offset = (char*)site - (char*)g_spinorfields_doubledata;
	bgq_spinorcoord *result = (bgq_spinorcoord*)((char*)g_spinorfields_doubledata_coords + offset);
	return result;
}

void bgq_reset_spinorfield(bgq_spinorfield_double spinorfield, bool isOdd) {
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
				_Complex double *val = bgq_spinorfield_double_ref(spinorfield, isOdd, t,x,y,z,v,c,false,false,false);
				bgq_spinorcoord *coord = bgq_spinorfield_coordref(val);
				*coord = bgq_spinorcoord_encode(isOdd,t,x,y,z,v,c);
			}
	}
}

void bgq_check_spinorfield(bgq_spinorfield_double spinorfield, int expected_reads, int expected_writes) {
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

		const int teo = t / PHYSICAL_LP;
		const int tv = teo/PHYSICAL_LK;
		const int k = mod(teo,PHYSICAL_LK);

		if (((t+x+y+z)&1) != isOdd)
			continue;

		for (int v = 0; v < 4; v+=1)
			for (int c = 0; c < 3; c += 1) {
				assert(assert_spinorfield_coord(spinorfield, isOdd, t,x,y,z,tv,k,v,c,false,false));
				_Complex double *site = bgq_spinorfield_double_ref(spinorfield, isOdd, t,x,y,z,v,c,false,false,true);
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

#pragma omp parallel for schedule(static)
	for (int txy = 0; txy < VOLUME/2; txy+=1) {
		WORKLOAD_DECL(txy, VOLUME/2);
		const int t = WORKLOAD_CHUNK(LOCAL_LT);
		const int x = WORKLOAD_CHUNK(LOCAL_LX);
		const int y = WORKLOAD_CHUNK(LOCAL_LY);
		const int z = WORKLOAD_CHUNK(LOCAL_LZ/2)*2 + (t+x+y+isOdd)%2;
		WORKLOAD_CHECK

		assert((t+x+y+z)%2==isOdd);

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
#ifndef NDEBUG
				_Complex double *val = bgq_spinorfield_double_ref(targetfield,isOdd,t,x,y,z,v,c,false,false,false);
				bgq_spinorcoord *coord = bgq_spinorfield_coordref(val);

				*coord = bgq_spinorcoord_encode(isOdd,t,x,y,z,v,c);
#endif
			}
		}
	}
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
	bgq_spinorsite_double *address = &spinorfield[idx];
	assert(g_spinorfields_double[index] <= address);
	assert(address <= g_spinorfields_double[VOLUME/2]);
	assert(mod((size_t)address,32)==0);

	// Get the equivalent address in debug address space
	bgq_spinorsite_double *debugsite = (bgq_spinorsite_double*)(((char*)g_spinorfields_doubledata_coords) + offset);
    bgq_spinorcoord *coord = (bgq_spinorcoord*)&debugsite->s[v][c][k];

	// Verify that we are really accessing the correct data
			assert(coord->coord.isOdd == isOdd);
			assert(coord->coord.t == t);
			assert(coord->coord.x == x);
			assert(coord->coord.y == y);
			assert(coord->coord.z == z);
			assert(coord->coord.v == v);
			assert(coord->coord.c == c);

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
			unsigned dir : 4; // T_UP..Z_UP
			unsigned i : 2; // 0..3
			unsigned l : 2; // 0..3
			bool isOdd : 1;
			// 57 bits = 8 bytes
		} coord;
} bgq_gaugecoord;
bgq_gaugecoord bgq_gaugecoord_encode(bool isOdd, int t, int x, int y, int z, direction dir, int i, int l) {
	bgq_gaugecoord result;
	assert(sizeof(result) == sizeof(_Complex double));
	result.val = 0; // Write zeros

	result.coord.isOdd = isOdd;
	result.coord.t = t;
	result.coord.x = x;
	result.coord.y = y;
	result.coord.z = z;
	result.coord.dir = dir;
	result.coord.i = i;
	result.coord.l = l;
	result.coord.writes = 0;
	result.coord.reads = 0;

	return result;
}


void bgq_init_gaugefield() {
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (direction dir = T_UP; dir <= T_UP_SHIFT; dir += 2) {
			g_gaugefield_doubledata.eodir[isOdd][dir/2] = malloc_aligned(GAUGE_EOVOLUME * sizeof(bgq_gaugesite_double), 128 /*L2 cache line size*/);
		}
	}

	g_gaugefield_double = &g_gaugefield_doubledata;

#ifndef NDEBUG
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (direction dir = T_UP; dir <= T_UP_SHIFT; dir += 2) {
			g_gaugefield_doubledata_debug.eodir[isOdd][dir/2] = malloc_aligned(GAUGE_EOVOLUME * sizeof(bgq_gaugesite_double), 128 /*L2 cache line size*/);
		}
	}
#endif
}

static bgq_gaugecoord *bgq_gaugefield_coordref(_Complex double *site) {
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (direction dir = T_UP; dir <= T_UP_SHIFT; dir += 2) {
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


void bgq_free_gaugefield() {
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (int dir = T_UP; dir <= T_UP_SHIFT; dir += 2) {
			free(g_gaugefield_doubledata.eodir[isOdd][dir/2]);
			g_gaugefield_doubledata.eodir[isOdd][dir/2] = NULL;
		}
	}
	g_gaugefield_double = NULL;

#ifndef NDEBUG
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (int dir = T_UP; dir <= T_UP_SHIFT; dir += 2) {
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
	BGQ_ENTER_FUNC
	assert(targetfield);
	assert(sourcefield);

#pragma omp parallel for schedule(static)
	for (int tzy = 0; tzy < GAUGE_VOLUME; tzy += 1) {
		WORKLOAD_DECL(tzy, GAUGE_VOLUME);
		const int z = WORKLOAD_PARAM(LOCAL_LZ);
		const int y = WORKLOAD_PARAM(LOCAL_LY+1) - 1;
		const int x = WORKLOAD_PARAM(LOCAL_LX+1) - 1;
		const int t = WORKLOAD_PARAM(LOCAL_LT+1) - 1;
		WORKLOAD_CHECK

		const bool isOdd = (8 + t + x + y + z) % 2;
		const int ix = Index(t, x, y, z); /* lexic coordinate */

		for (direction d = T_UP; d <= Z_UP; d += 2) {
			su3_array64 *m = (su3_array64*) &g_gauge_field[ix][d / 2];
			for (int i = 0; i < 3; i += 1) {
				for (int l = 0; l < 3; l += 1) {
					//complex *c = ((complex*)m)+i*3+l;
					bgq_gaugefield_double_set(targetfield, isOdd, t, x, y, z, d, i, l, m->c[i][l]);
#ifndef NDEBUG
					_Complex double *val = bgq_gaugefield_double_ref(targetfield, isOdd, t,x,y,z,d,i,l,false,false,false);
					bgq_gaugecoord *debugval = bgq_gaugefield_coordref(val );
					*debugval = bgq_gaugecoord_encode(isOdd,t,x,y,z,d,i,l);
			#endif
				}
			}
		}
	}
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
	assert(-1 <= z && z < LOCAL_LZ);
	assert(-1 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= k && k < PHYSICAL_LK);
	assert(T_UP <= dir && dir <= T_DOWN_SHIFT);
	assert(0 <= i && i < 3);
	assert(0 <= l && l < 3);

	// There is really just one gaugefield
	assert(gaugefield == g_gaugefield_double);

	// Check that the coordinate is really an odd/even coordinate
	assert(((t+x+y+z)&1) == isOdd);
	const int teo = t/PHYSICAL_LP;

	// Check that zv and k match the coordinate
	assert(teo/PHYSICAL_LK == tv);
	assert(mod(teo,PHYSICAL_LK) == k);

	// We really only store the up-fields
	assert((dir&1)==0);

	// Get the index
	const int idx = (((tv)*PHYSICAL_LX + ((x)+1))*PHYSICAL_LY + ((y)+1))*PHYSICAL_LZ + (z);
	assert(0 <= idx && idx < (PHYSICAL_LTV+1)*(PHYSICAL_LX+1)*(PHYSICAL_LY+1)*(PHYSICAL_LZ));

	// Get the address
	bgq_gaugesite_double *eofield = gaugefield->eodir[isOdd][dir/2];
	bgq_gaugesite_double *site = &eofield[idx];
	_Complex double *address = site->c[i][l];
	assert(mod((size_t)address,32)==0);

	// get the debug data
	bgq_gaugecoord *debugdata = bgq_gaugefield_coordref(address);
	assert(debugdata->coord.isOdd == isOdd);
	assert(debugdata->coord.t == t);
	assert(debugdata->coord.x == x);
	assert(debugdata->coord.y == y);
	assert(debugdata->coord.z == z);
	assert(debugdata->coord.dir == dir);
	assert(debugdata->coord.i == i);
	assert(debugdata->coord.l == l);

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
	assert((weylfield == weylxchange_recv_double[T_UP]) || (weylfield == weylxchange_recv_double[T_DOWN])
			|| (weylfield == weylxchange_send_double[T_UP]) || (weylfield == weylxchange_send_double[T_DOWN]) );

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
	assert((weylfield == weylxchange_recv_double[X_UP]) || (weylfield == weylxchange_recv_double[X_DOWN])
			|| (weylfield == weylxchange_send_double[X_UP]) || (weylfield == weylxchange_send_double[X_DOWN]) );

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
	assert((weylfield == weylxchange_recv_double[Y_UP]) || (weylfield == weylxchange_recv_double[Y_DOWN])
			|| (weylfield == weylxchange_send_double[Y_UP]) || (weylfield == weylxchange_send_double[Y_DOWN]));

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

