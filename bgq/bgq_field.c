


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


void bgq_init_spinorfields(int count) {
	BGQ_ENTER_FUNC
	g_num_spinorfields = count;
	int datasize = count * sizeof(*g_spinorfields_doubledata) * VOLUME/2;
	g_spinorfields_doubledata = malloc_aligned(datasize, 128);


	g_spinorfields_double = malloc(count * sizeof(*g_spinorfields_double));
	g_spinorfield_isOdd = malloc(count * sizeof(*g_spinorfield_isOdd));
	for (int i = 0; i < count; i += 1) {
		g_spinorfields_double[i] = g_spinorfields_doubledata + i * VOLUME/2;
		g_spinorfield_isOdd[i] = -1; // Unknown yet
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


void bgq_init_gaugefield() {
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (direction dir = T_UP; dir <= T_UP_SHIFT; dir += 2) {
			g_gaugefield_doubledata.eodir[isOdd][dir/2] = malloc_aligned((PHYSICAL_LTV+1) * (PHYSICAL_LX+1) * (PHYSICAL_LY+1) * PHYSICAL_LZ * sizeof(bgq_gaugesite_double), 128 /*L2 cache line size*/);
		}
	}

	g_gaugefield_double = &g_gaugefield_doubledata;
}

void bgq_free_gaugefield() {
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (int dir = T_UP; dir <= T_UP_SHIFT; dir += 2) {
			free(g_gaugefield_doubledata.eodir[isOdd][dir/2]);
			g_gaugefield_doubledata.eodir[isOdd][dir/2] = NULL;
		}
	}
	g_gaugefield_double = NULL;
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
				}
			}
		}
	}
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
			}
		}
	}
}



bool assert_spinorcoord(bgq_spinorfield_double spinorfield, bool isOdd, int t, int x, int y, int z, int tv, int k) {
	assert(spinorfield);
	assert(false <= isOdd && isOdd <= true);
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= k && k < PHYSICAL_LK);

	// Get the index of the used spinorfield
	const int fieldsize = sizeof(bgq_spinorsite_double) * VOLUME/2; // alignment???
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

	return true; // All checks passed
}



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


bool assert_gaugesite(bgq_gaugefield_double gaugefield, bool isOdd, int t, int x, int y, int z, int tv, int k, direction dir) {
	assert(gaugefield);
	assert(false <= isOdd && isOdd <= true);
	assert(-1 <= t && t < LOCAL_LT);
	assert(-1 <= x && x < LOCAL_LX);
	assert(-1 <= y && y < LOCAL_LY);
	assert(-1 <= z && z < LOCAL_LZ);
	assert(-1 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= k && k < PHYSICAL_LK);
	assert(T_UP <= dir && dir <= T_DOWN_SHIFT);

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
	bgq_gaugesite_double *address = &eofield[idx];
	//assert(&g_gaugefield_doubledata[0] <= address && address < &g_gaugefield_doubledata[(PHYSICAL_LTV+1)*(PHYSICAL_LX+1)*(PHYSICAL_LY+1)*(PHYSICAL_LZ)]);
	assert(mod((size_t)address,32)==0);

	// All checks passed
	return true;
}

