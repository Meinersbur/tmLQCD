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

static int g_num_spinorfields = 0;
void bgq_init_spinorfields(int count) {
	g_num_spinorfields = count;
	int datasize = count * sizeof(*g_spinorfields_doubledata) * VOLUME;
	g_spinorfields_doubledata = malloc_aligned(datasize, 128);
	memset(g_spinorfields_doubledata, 0, datasize);
	g_spinorfields_double = malloc(count * sizeof(*g_spinorfields_double));

	for (int i = 0; i < count; i += 1) {
		g_spinorfields_double[i] = g_spinorfields_doubledata + i * VOLUME;
	}
}

void bgq_free_spinofields() {
	free(g_spinorfields_double);
	g_spinorfields_double = NULL;
	free(g_spinorfields_doubledata);
	g_spinorfields_doubledata = NULL;
	g_num_spinorfields = 0;
}

void bgq_init_gaugefield() {
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (direction dir = T_UP; dir <= Z_UP_SHIFT; dir += 2) {
			int zlinelength = PHYSICAL_LZV;
			if (dir == Z_UP_SHIFT)
				zlinelength += 1;

			g_gaugefield_doubledata.eodir[isOdd][dir / 2] = malloc_aligned((PHYSICAL_LT + 1) * (PHYSICAL_LX + 1) * (PHYSICAL_LY + 1) * zlinelength * sizeof(bgq_gaugesite_double), 128 /*L2 cache line size*/);
		}
	}

	g_gaugefield_double = &g_gaugefield_doubledata;
}

void bgq_free_gaugefield() {
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (int dir = T_UP; dir <= Z_UP_SHIFT; dir += 2) {
			free(g_gaugefield_doubledata.eodir[isOdd][dir / 2]);
			g_gaugefield_doubledata.eodir[isOdd][dir / 2] = NULL;
		}
	}
	g_gaugefield_double = NULL;
}
typedef struct {
	double _Complex c[3][3];
} su3_array64;

void bgq_transfer_gaugefield(bgq_gaugefield_double const targetfield, su3 ** const sourcefield) {
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
		int icx = g_lexic2eosub[ix]; /*  even/odd coordinate relative to field base */
		assert(icx == iy - (isOdd ? (VOLUME+RAND)/2 : 0));

		spinor_array64 *sp = (spinor_array64*)&sourcefield[icx];
		for (int v = 0; v < 4; v+=1) {
			for (int c = 0; c < 3; c+=1) {
				double _Complex *sdata = bgq_spinorfield_double_local_to_physical(targetfield,isOdd,t,x,y,z,v,c);
				*sdata = sp->v[v].c[c];
			}
		}
	}
}






