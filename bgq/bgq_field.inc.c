
#define BGQ_FIELD_INC_C_




#ifndef PRECISION
//#define BGQ_PRECISION 64
//#include "bgq_field.inc.h"
#include "bgq_field_double.h"
#endif

#include "bgq_field.inc.h"
#include "bgq.h"
#include "../geometry_eo.h"
#include "../global.h"

#include <string.h>
#include <stdbool.h>
#include <omp.h>
#include <stddef.h>

#ifdef XLC
#include <l1p/sprefetch.h>
#endif
#if BGQ_PREFETCH_LIST
#include <l1p/pprefetch.h>
#endif

#ifndef BGQ_FIELD_COORDCHECK
#define BGQ_FIELD_COORDCHECK 0
#endif



////////////////////////////////////////////////////////////////////////////////
// Spinorfield

static bgq_spinorsite *g_spinorfields_data = NULL;
bgq_spinorfield *g_spinorfields = NULL;
static int g_num_spinorfields = -1;
static int g_num_chi_spinorfields = -1;
static int g_chi_up_spinorfield_first = -1;
static int g_chi_dn_spinorfield_first = -1;
int g_num_total_spinorfields = -1;


#if BGQ_FIELD_COORDCHECK
static bgq_spinorsite *g_spinorfields_data_coords = NULL;
#endif
static int *g_spinorfield_isOdd = NULL;



bgq_spinorfield bgq_translate_spinorfield(spinor * const field) {
	int index = -1;

	if (g_chi_up_spinor_field) {
		size_t chi_fieldsize = (char*) g_chi_up_spinor_field[1] - (char*) g_chi_up_spinor_field[0];
		assert(chi_fieldsize > 0);
		if ((char*) g_chi_up_spinor_field[0] <= (char*) field && (char*) field < ((char*) g_chi_up_spinor_field[g_num_chi_spinorfields - 1] + chi_fieldsize)) {
			// This is a spinor from g_chi_up_spinor_field
			size_t offset = (char*) field - (char*) g_chi_up_spinor_field[0];
			assert(offset >= 0);
			assert(offset % chi_fieldsize == 0);
			index = g_chi_up_spinorfield_first + offset / chi_fieldsize;
			assert(index >= g_chi_up_spinorfield_first);
		} else if ((char*) g_chi_dn_spinor_field[0] <= (char*) field && (char*) field < ((char*) g_chi_dn_spinor_field[g_num_chi_spinorfields - 1] + chi_fieldsize)) {
			// This is a spinor from g_chi_dn_spinor_field
			size_t offset = (char*) field - (char*) g_chi_dn_spinor_field[0];
			assert(offset >= 0);
			assert(offset % chi_fieldsize == 0);
			index = g_chi_dn_spinorfield_first + offset / chi_fieldsize;
			assert(index >= g_chi_dn_spinorfield_first);
		}
	}

	if (index == -1) {
		int V = even_odd_flag ? VOLUMEPLUSRAND / 2 : VOLUMEPLUSRAND;
		size_t fieldsize = V * sizeof(*field);

		// This computes the original index address of the passed field; be aware that its correctness depends on the implementation of init_spinorfield
		size_t offset = (char*) field - (char*) g_spinor_field[0];
		assert(offset >= 0);
		assert(offset % fieldsize == 0);
		index = offset / fieldsize;
	}

	assert(index >= 0);
	assert(index < g_num_total_spinorfields);
	bgq_spinorfield result = g_spinorfields[index];
	assert(result >= g_spinorfields_data);
	return result;
}


static COMPLEX_PRECISION *bgq_spinorfield_ref(bgq_spinorfield spinorfield, bool isOdd, int t, int x, int y, int z, int v, int c, bool isRead, bool isWrite) {
	const int teo = t / PHYSICAL_LP;
	const int tv = teo/PHYSICAL_LK;
	const int k = mod(teo,PHYSICAL_LK);

	COMPLEX_PRECISION *val = BGQ_SPINORVAL(spinorfield, isOdd, t, x, y, z, tv, k, v, c, isRead, isWrite);
	return val;
}


COMPLEX_PRECISION bgq_spinorfield_get(bgq_spinorfield spinorfield, bool isOdd, int t, int x, int y, int z, int v, int c) {
	COMPLEX_PRECISION *ptr = bgq_spinorfield_ref(spinorfield,isOdd,t,x,y,z,v,c,true,false);
	return *ptr;
}


void bgq_spinorfield_set(bgq_spinorfield spinorfield, bool isOdd, int t, int x, int y, int z, int v, int c, COMPLEX_PRECISION value) {
	COMPLEX_PRECISION *ptr = bgq_spinorfield_ref(spinorfield,isOdd,t,x,y,z,v,c,false,true);
	*ptr = value;
}



#if BGQ_FIELD_COORDCHECK
typedef union {
	COMPLEX_PRECISION val;
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


static bgq_spinorcoord *bgq_spinorfield_coordref(COMPLEX_PRECISION *site) {
	size_t offset = (char*)site - (char*)g_spinorfields_data;
	bgq_spinorcoord *result = (bgq_spinorcoord*)((char*)g_spinorfields_data_coords + offset);
	return result;
}


static void bgq_spinorfield_checkcoord(bgq_spinorfield spinorfield, bool isOdd, int t, int x, int y, int z, int tv, int k, int v, int c, COMPLEX_PRECISION *value) {
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


void bgq_spinorfield_resetcoord(bgq_spinorfield spinorfield, bool isOdd, int expected_reads_min, int expected_reads_max, int expected_writes_min, int expected_writes_max) {
	assert(spinorfield);

	bool expected_reads_min_reached = false;
	bool expected_reads_max_reached = false;
	bool expected_writes_min_reached = false;
	bool expected_writes_max_reached = false;

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
				COMPLEX_PRECISION *value = BGQ_SPINORVAL(spinorfield,isOdd,t,x,y,z,tv,k,v,c,false,false);
				bgq_spinorcoord *coord = bgq_spinorfield_coordref(value);

				assert (coord->coord.init);
				int reads = coord->coord.reads;
				int writes = coord->coord.writes;

				if (expected_reads_min >= 0) {
					assert(reads >= expected_reads_min);
					if (reads == expected_reads_min)
						expected_reads_min_reached = true;
				}
				if (expected_reads_max >= 0) {
					assert(reads <= expected_reads_max);
					if (reads == expected_reads_max)
						expected_reads_max_reached = true;
				}
				if (expected_writes_min >= 0) {
					assert(writes >= expected_writes_min);
					if (writes == expected_writes_min)
						expected_writes_min_reached = true;
				}
				if (expected_writes_max >= 0) {
					assert(writes <= expected_writes_max);
					if (writes == expected_writes_max)
						expected_writes_max_reached = true;
				}


				coord->coord.writes = 0;
				coord->coord.reads = 0;
			}
		}
	}

	assert(expected_reads_min < 0 || expected_reads_min_reached);
	assert(expected_reads_max < 0 || expected_reads_max_reached);
	assert(expected_writes_min < 0 || expected_writes_min_reached);
	assert(expected_writes_max < 0 || expected_writes_max_reached);
}
#endif


void bgq_init_spinorfields(int count, int chi_count) {
	// preconditions
	if (!even_odd_flag)
		master_error(1, "ERROR: even_odd_flag must be set\n");

	if (LOCAL_LT < 8) /* even/odd, 2-complex vectors, left and right border cannot be the same */
		master_error(1, "ERROR: Local T-dimension (=%d) must be at least 8 such that left and right (vectorized) surface do not overlap\n", LOCAL_LT);

	if (mod(LOCAL_LT,4)!=0) /* even/odd, 2-complex vectors, */
		master_error(1, "ERROR: Local T-dimension (=%d) must be a multiple of 4 (vector width * 2 for even/odd)\n", LOCAL_LT);

	if (LOCAL_LX < 2) /* border up- and down- surfaces cannot collapse */
		master_error(1, "ERROR: Local X-dimension must be larger than 2, such that surfaces do not overlap\n");

	if (mod(LOCAL_LX,2)!=0) /* even/odd */
		master_error(1, "ERROR: Local X-dimension must be a multiple of 2 for even/odd\n");

	if (mod(LOCAL_LX,4)!=0)
		master_error(1, "ERROR: Local X-dimension must be a multiple of 4 (vector width * 2 for even/odd) because of vectorized transmission of halfspinors in T-direction\n");

	if (LOCAL_LY < 2) /* border up- and down- surfaces cannot collapse */
		master_error(1, "ERROR: Local Y-dimension must be larger than 2, such that surfaces do not overlap\n");

	if (mod(LOCAL_LY,2)!=0) /* even/odd */
		master_error(1, "ERROR: Local Y-dimension must be a multiple of 2 for even/odd\n");

	if (LOCAL_LZ < 2) /* border up- and down- surfaces cannot collapse */
		master_error(1, "ERROR: Local X-dimension must be larger than 2, such that surfaces do not overlap\n");

	if (mod(LOCAL_LZ,2)!=0) /* even/odd */
		master_error(1, "ERROR: Local X-dimension must be a multiple of 2 for even/odd\n");

	if (BODY_SITES < 0)
		master_error(1, "ERROR: negative body volume\n");

	if (BODY_SITES == 0)
		master_print("WARNING: Local lattice consists of surface only\n");

	if ( pow(PHYSICAL_LTV*LOCAL_LX*LOCAL_LZ,1.0/3.0) >=  LOCAL_LZ)
		master_print("WARNING: Make the local Z-dimension the longest, data reuse has been optimized for it\n");

	master_print("INFO: Body-to-volume site ratio: %f (the higher the better)\n", (double)BODY_SITES / (double)VOLUME_SITES);

	int nThreads = omp_get_max_threads();
	#pragma omp parallel
	{
		#pragma omp master
		{
			assert(nThreads == omp_get_num_threads());
			master_print("INFO: OMP_NUM_THREADS=%d\n", omp_get_num_threads());
		}
	}

	if (PHYSICAL_LTV*LOCAL_LX*LOCAL_LZ < omp_get_num_threads()) {
		master_print("WARNING: Less z-lines than threads to process them\n");
	} else if ((double)mod(PHYSICAL_LTV*LOCAL_LX*LOCAL_LZ, omp_get_num_threads())/(double)omp_get_num_threads() > 0.5) {
		master_print("WARNING: z-lines to threads imbalance\n");
	}


	g_num_spinorfields = count;
	g_num_chi_spinorfields = chi_count;
	g_num_total_spinorfields = count + 2*chi_count;
	g_chi_up_spinorfield_first = count;
	g_chi_dn_spinorfield_first = count + chi_count;

#if 0
	int datasize = g_num_total_spinorfields * READTOTLENGTH * sizeof(PRECISION);
#else
	int datasize = g_num_total_spinorfields * sizeof(*g_spinorfields_data) * VOLUME_SITES;
#endif
	g_spinorfields_data = malloc_aligned(datasize, 128);
	#if BGQ_FIELD_COORDCHECK
		g_spinorfields_data_coords = malloc_aligned(datasize, 128);
	#endif

	g_spinorfields = malloc(g_num_total_spinorfields * sizeof(*g_spinorfields));
	//#if BGQ_FIELD_COORDCHECK
		g_spinorfield_isOdd = malloc(g_num_total_spinorfields * sizeof(*g_spinorfield_isOdd));
	//#endif
	for (int i = 0; i < g_num_total_spinorfields; i += 1) {
		g_spinorfields[i] = g_spinorfields_data + i * VOLUME_SITES;
		//#if BGQ_FIELD_COORDCHECK
			g_spinorfield_isOdd[i] = -1; // Unknown yet
		//#endif
	}

#if BGQ_PREFETCH_LIST
	bgq_listprefetch_handle = malloc(g_num_total_spinorfields * g_num_total_spinorfields * sizeof(*bgq_listprefetch_handle));
	memset(bgq_listprefetch_handle, 0, g_num_total_spinorfields * g_num_total_spinorfields * sizeof(*bgq_listprefetch_handle));
#endif
}


void bgq_free_spinofields() {
	free(g_spinorfields);
	g_spinorfields = NULL;
	free(g_spinorfields_data);
	g_spinorfields_data = NULL;
	g_num_total_spinorfields = 0;
	g_num_spinorfields = 0;
	g_num_chi_spinorfields = 0;

	free(g_spinorfield_isOdd);
	g_spinorfield_isOdd = NULL;
#if BGQ_FIELD_COORDCHECK
	free(g_spinorfields_data_coords);
	g_spinorfields_data_coords = NULL;
#endif
	//TODO: free(bgq_listprefetch_handle), L1P_DeallocatePattern(bgq_listprefetch_handle[...])
}


typedef struct {
	double _Complex c[3];
} su3_vector_array64;

typedef struct {
	su3_vector_array64 v[4];
} spinor_array64;


void bgq_transfer_spinorfield(const bool isOdd, bgq_spinorfield targetfield, spinor * const sourcefield) {
	assert(sourcefield);
	assert(targetfield);

#if BGQ_FIELD_COORDCHECK
	bgq_spinorfield_resetcoord(targetfield, isOdd, -1, -1, -1, -1);
#endif

#pragma omp parallel for schedule(static)
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
				bgq_spinorfield_set(targetfield,isOdd,t,x,y,z,v,c,sp->v[v].c[c]);
			}
		}
	}

#if BGQ_FIELD_COORDCHECK
	bgq_spinorfield_resetcoord(targetfield, isOdd, 0, 0, 1, 1);
#endif

	bgq_spinorfield_setOdd(targetfield, isOdd, true);
}


void bgq_spinorfield_transfer_back(const bool isOdd, spinor * const targetfield, bgq_spinorfield const sourcefield) {
	assert(sourcefield);
	assert(targetfield);
	bgq_spinorfield_setOdd(sourcefield, isOdd, false);

#pragma omp parallel for schedule(static)
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

		spinor_array64 *sp = (spinor_array64*)&targetfield[icx];
		for (int v = 0; v < 4; v+=1) {
			for (int c = 0; c < 3; c+=1) {
				COMPLEX_PRECISION val = bgq_spinorfield_get(sourcefield,isOdd,t,x,y,z,v,c);
				sp->v[v].c[c] = val;
			}
		}
	}
}


double bgq_spinorfield_compare(const bool isOdd, bgq_spinorfield const bgqfield, spinor * const reffield, bool silent) {
	assert(bgqfield);
	assert(reffield);

#if BGQ_FIELD_COORDCHECK
	bgq_spinorfield_resetcoord(bgqfield, isOdd, -1, -1, -1, -1);
#endif

	double norm1_sum = 0;
	double norm1_max = 0;
	double norm2_sum = 0;
	double norm2_max = 0;
	int count = 0;

//#pragma omp parallel
	{
		double worker_norm1_sum = 0;
		double worker_norm2_sum = 0;
		double worker_norm1_max = 0;
		double worker_norm2_max = 0;
		int worker_count = 0;

#pragma omp for schedule(static) nowait
		for (int txyz = 0; txyz < VOLUME; txyz += 1) {
			WORKLOAD_DECL(txyz, VOLUME);
			const int t = WORKLOAD_CHUNK(LOCAL_LT);
			const int x = WORKLOAD_CHUNK(LOCAL_LX);
			const int y = WORKLOAD_CHUNK(LOCAL_LY);
			const int z = WORKLOAD_CHUNK(LOCAL_LZ);
			WORKLOAD_CHECK

			if ( ((t+x+y+z)&1) != isOdd )
				continue;

			const int ix = g_ipt[t][x][y][z]; /* lexic coordinate */
			assert(ix == Index(t,x,y,z));

			int iy = g_lexic2eo[ix]; /* even/odd coordinate (even and odd sites in two different fields of size VOLUME/2, first even field followed by odd) */
			assert(0 <= iy && iy < (VOLUME+RAND));
			int icx = g_lexic2eosub[ix]; /*  even/odd coordinate relative to field base */
			assert(0 <= icx && icx < VOLUME/2);
			assert(icx == iy - (isOdd ? (VOLUME+RAND)/2 : 0));

			bool first = true;
			spinor_array64 *sp = (spinor_array64*) &reffield[icx];
			for (int v = 0; v < 4; v += 1) {
				for (int c = 0; c < 3; c += 1) {
					COMPLEX_PRECISION refvalue = sp->v[v].c[c];
					COMPLEX_PRECISION bgqvalue = bgq_spinorfield_get(bgqfield, isOdd, t, x, y, z, v, c);
					double norm1_val = cabs(refvalue - bgqvalue);
					double norm2_val = norm1_val*norm1_val;

					if (norm1_val > 0.01) {
						if (first) {
							//if (!silent)
							//	master_print("Coordinate (%d,%d,%d,%d)(%d,%d): ref=(%f + %fi) != bgb=(%f + %fi) off by %f\n", t,x,y,z,v,c,creal(refvalue), cimag(refvalue),creal(bgqvalue),cimag(bgqvalue),norm1_val);
							worker_count += 1;
							//fprintf(stderr, "Coordinate (%d,%d,%d,%d)(%d,%d): ref=(%f + %fi) != bgb=(%f + %fi) off by %f\n", t,x,y,z,v,c,creal(refvalue), cimag(refvalue),creal(bgqvalue),cimag(bgqvalue),norm1_val);
							//__asm__("int3");
						}
						first = false;
					}

					worker_norm1_sum += norm1_val;
					worker_norm2_sum += norm2_val;
					worker_norm1_max = max(worker_norm1_max, norm1_val);
					worker_norm2_max = max(worker_norm2_max, norm2_val);
				}
			}
		}

#if BGQ && defined(XLC)
#pragma omp tm_atomic
#else
#pragma omp critical
#endif
		{
			norm1_sum += worker_norm1_sum;
			norm2_sum += worker_norm2_sum;
			norm1_max = max(norm1_max, worker_norm1_max);
			norm2_max = max(norm2_max, worker_norm2_max);
			count += worker_count;
		}
	}

#if BGQ_FIELD_COORDCHECK
	bgq_spinorfield_resetcoord(bgqfield, isOdd, 1, 1, 0, 0);
#endif

	double send_sum[] = { norm1_sum, norm2_sum };
	double recv_sum[] = { -1, -1 };
	MPI_Allreduce(send_sum, recv_sum, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	norm1_sum = recv_sum[0];
	norm2_sum = recv_sum[1];

	double send_max[] = { norm1_max, norm2_max };
	double recv_max[] = { -1, -1 };
	MPI_Allreduce(send_max, recv_max, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
	norm1_max = recv_max[0];
	norm2_max = recv_max[1];

	//double avg1 = norm1_sum / (VOLUME * 4 * 3);
	//double avg2 = sqrt(norm2_sum) / (VOLUME * 4 * 3);
	//double norminf = norm1_max;

	if (count>0) {
		if (!silent)
			master_print("%d sites of %d wrong\n", count, VOLUME);
	}

	return norm1_max;
}

int bgq_spinorfield_find_index(bgq_spinorfield spinorfield) {
	const int fieldsize = sizeof(bgq_spinorsite) * VOLUME_SITES; // alignment???
	long long offset = (char*)spinorfield - (char*)g_spinorfields_data;
	assert(offset >= 0);
	assert(mod(offset, fieldsize) == 0);
	int index = offset / fieldsize;
	assert(index < g_num_spinorfields);
	return index;
}

bool bgq_spinorfield_isOdd(bgq_spinorfield spinorfield) {
	const size_t fieldsize = sizeof(bgq_spinorsite) * VOLUME_SITES; // alignment???
	assert(spinorfield >= g_spinorfields_data);
	const size_t offset = (char*)spinorfield - (char*)g_spinorfields_data;
	assert(offset % fieldsize == 0);
	size_t index = offset / fieldsize;
	assert(index < g_num_total_spinorfields);

	int isOdd = g_spinorfield_isOdd[index];
	assert(isOdd != -1);
	return isOdd;
}

void bgq_spinorfield_setOdd(bgq_spinorfield spinorfield, bool isOdd, bool overwrite) {
	assert(spinorfield);

	const size_t fieldsize = sizeof(bgq_spinorsite) * VOLUME_SITES; // alignment???
	assert(spinorfield >= g_spinorfields_data);
	const size_t offset = (char*)spinorfield - (char*)g_spinorfields_data;
	assert(offset % fieldsize == 0);
	size_t index = offset / fieldsize;
	assert(index < g_num_total_spinorfields);

	if (!overwrite) {
		int oldIsOdd = g_spinorfield_isOdd[index];
		if (oldIsOdd == -1) {
			assert(!"If setting oddness, it should already have been set");
		} else if (oldIsOdd == isOdd) {
			// Ok
		} else {
			assert(!"Oddness not matching\n");
			return;
		}
	}
	g_spinorfield_isOdd[index] = isOdd;
}


bool assert_spinorfield_coord(bgq_spinorfield spinorfield, bool isOdd, int t, int x, int y, int z, int tv, int k, int v, int c, bool isRead, bool isWrite) {
	if (!isRead && !isWrite && (z == LOCAL_LZ) ) {
		// Allow for prefetching
		return true;
	}

	assert(spinorfield);
	assert(false <= isOdd && isOdd <= true);
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= k && k < PHYSICAL_LK);

	// Get the index of the used spinorfield
	const size_t fieldsize = sizeof(bgq_spinorsite) * VOLUME_SITES; // alignment???
	const size_t offset = (char*)spinorfield - (char*)g_spinorfields_data;
	assert(offset >= 0);
	assert(mod(offset, fieldsize) == 0);
	int index = offset / fieldsize;
	assert(index < g_num_total_spinorfields);

	// Check that the coordinate is really an odd/even coordinate
	assert(((t+x+y+z)&1) == isOdd);

#if BGQ_FIELD_COORDCHECK
	// Check that field is used as an odd/even field
	if (g_spinorfield_isOdd[index] == -1) {
		// not yet defined, just ensure that at all following uses are the same
		g_spinorfield_isOdd[index] = isOdd;
	} else {
		assert(g_spinorfield_isOdd[index] == isOdd);
	}
#endif

	const int teo = t/PHYSICAL_LP;

	// Check that zv and k match the coordinate
	assert(teo/PHYSICAL_LK == tv);
	assert(mod(teo,PHYSICAL_LK) == k);

	// Get the memory address it points to
	const int idx = ((tv*PHYSICAL_LX + x)*PHYSICAL_LY + y)*PHYSICAL_LZ + z;
	assert(0 <= idx && idx < PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LY*PHYSICAL_LZ);

	// Validate the address
	bgq_spinorsite *site = &spinorfield[idx];
	assert(g_spinorfields[index] <= site && site < g_spinorfields[index] + VOLUME_SITES);
	assert(mod((size_t)site,PRECISION_VECTOR_ALIGNMENT)==0);

	COMPLEX_PRECISION *val = &site->s[v][c][k];
	assert( (&g_spinorfields[index]->s[0][0][0] <= val) && (val < &(g_spinorfields[index]+VOLUME_SITES)->s[0][0][0]) );
	assert(mod((size_t)site,16)==0);
	assert(mod((size_t)&site->s[v][c][0],PRECISION_VECTOR_ALIGNMENT)==0);

#if BGQ_FIELD_COORDCHECK
	// Get the equivalent address in debug address space, and check it
	bgq_spinorfield_checkcoord(spinorfield,isOdd,t,x,y,z,tv,k,v,c,val);
	bgq_spinorcoord *coord = bgq_spinorfield_coordref(val);

	// Mark uses of this value
	if (isRead)
		coord->coord.reads += 1;
	if (isWrite)
		coord->coord.writes += 1;
#endif

	return true; // All checks passed
}


bool assert_spinorcoord(bgq_spinorfield spinorfield, bool isOdd, int t, int x, int y, int z, int tv, int k, bool isRead, bool isWrite) {
	for (int v = 0; v < 4; v += 1) {
			for (int c = 0; c < 3; c += 1) {
				assert_spinorfield_coord(spinorfield, isOdd, t, x, y, z, tv, k, v, c, isRead, isWrite);
			}
	}

	return true;
}





////////////////////////////////////////////////////////////////////////////////
// Gaugefield

bgq_gaugefield g_gaugefield;
static bgq_gaugeeodir g_gaugefield_data;
#if BGQ_FIELD_COORDCHECK
static bgq_gaugeeodir g_gaugefield_data_debug;
#endif




static COMPLEX_PRECISION *bgq_gaugefield_ref(bgq_gaugefield gaugefield, bool isOdd, int t, int x, int y, int z, direction d, int i, int l, bool isRead, bool isWrite) {
	assert(gaugefield);
	assert(false <= isOdd && isOdd <= true);
	assert(isOdd == ((t+x+y+z)&1));
	assert(-1 <= t && t < LOCAL_LT);
	assert(-1 <= x && x < LOCAL_LX);
	assert(-1 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	assert(TUP <= d && d <= ZDOWN);
	assert(0 <= i && i < 3);
	assert(0 <= l && l < 3);
	assert((d&1)==0); // Must be up-direction

	const int teo = divdown(t, PHYSICAL_LP);
	const int tv = divdown(teo, PHYSICAL_LK);
	const int k = moddown(teo, PHYSICAL_LK);


	assert(assert_gaugeval(gaugefield,isOdd,t,x,y,z,tv,k,d,i,l,isRead,isWrite));
	bgq_gaugesite *site = BGQ_GAUGESITE_ACCESS(gaugefield, isOdd, tv,x,y,z,d);
	COMPLEX_PRECISION *val = &site->c[i][l][k];

	return val;
}


COMPLEX_PRECISION bgq_gaugefield_get(bgq_gaugefield gaugefield, bool isOdd, int t, int x, int y, int z, direction d, int i, int l) {
	bool adjoint = false;
	switch (d) {
	case TDOWN:
	case TDOWN_SHIFT:
		t -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = TUP;
		break;
	case XDOWN:
		x -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = XUP;
		break;
	case YDOWN:
		y -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = YUP;
		break;
	case ZDOWN:
		z -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = ZUP;
		break;
	case TUP_SHIFT:
		d = TUP;
		break;
	default:
		break;
	}

	if (adjoint) {
		const int tmp = i;
		i = l;
		l = tmp;
	}

	COMPLEX_PRECISION result = *bgq_gaugefield_ref(gaugefield, isOdd, t,x,y,z,d,i,l, true, false);
	if (adjoint) {
		result = conj(result);
	}
	return result;
}


void bgq_gaugefield_set(bgq_gaugefield gaugefield, bool isOdd, int t, int x, int y, int z, direction d, int i, int l, double _Complex value) {
	assert(gaugefield);
	assert(false <= isOdd && isOdd <= true);
	assert(isOdd == ((t+x+y+z)&1));
	assert(-1 <= t && t < LOCAL_LT);
	assert(-1 <= x && x < LOCAL_LX);
	assert(-1 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	assert(TUP <= d && d <= ZDOWN);
	assert(0 <= i && i < 3);
	assert(0 <= l && l < 3);

	bool adjoint = false;
	switch (d) {
	case TDOWN:
	case TDOWN_SHIFT:
		t -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = TUP;
		break;
	case XDOWN:
		x -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = XUP;
		break;
	case YDOWN:
		y -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = YUP;
		break;
	case ZDOWN:
		z -= 1;
		isOdd = !isOdd;
		adjoint = true;
		d = ZUP;
		break;
	case TUP_SHIFT:
		d = TUP;
		break;
	default:
		break;
	}

	const int teo = divdown(t, PHYSICAL_LP);
	const int tv = divdown(teo, PHYSICAL_LK);
	const int k = moddown(teo, PHYSICAL_LK);

	if (adjoint) {
		int tmp = i;
		i = l;
		l = tmp;
		crealf(value);
		value = conj(value);
	}

	{
	COMPLEX_PRECISION *valptr = BGQ_GAUGEVAL(gaugefield, isOdd, t, x, y, z, tv, k, d, i, l, false, true);
	*valptr = value;
	}

	if (d == ZUP) {
		// For z-direction...
		if (z==LOCAL_LZ-1) {
			// wraparound, also store as z=-1 coordinate
			COMPLEX_PRECISION *wrapptr = BGQ_GAUGEVAL(gaugefield, isOdd, t, x, y, -1, tv, k, d, i, l, false, true);
			*wrapptr = value;
		}
	}

	if (d == TUP) {
			// For t-direction...
			// also store as shifted

			// Move one to the right
		    const int teo_shift = moddown(teo+2, 1+LOCAL_LT/PHYSICAL_LP)-1;
			const int tv_shift = divdown(teo_shift, PHYSICAL_LK);
			const int k_shift = moddown(teo_shift, PHYSICAL_LK);

			if ( (t == LOCAL_LT-1) || (t == LOCAL_LT-2) ) {
				assert(tv_shift == -1);
				assert(k_shift == 1);
			} else {
				assert( (tv_shift==tv && k_shift==1) || (tv_shift==tv+1 && k_shift==0) );
			}

			COMPLEX_PRECISION *shiftptr = BGQ_GAUGEVAL(gaugefield, isOdd, t, x, y, z, tv_shift, k_shift, TUP_SHIFT, i, l, false, true);
			*shiftptr = value;
		}
}



#if BGQ_FIELD_COORDCHECK
typedef
	union {
		COMPLEX_PRECISION val;
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


static bgq_gaugecoord *bgq_gaugefield_coordref(COMPLEX_PRECISION *site) {
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (direction dir = TUP; dir <= TUP_SHIFT; dir += 2) {
			bgq_gaugesite *eodata = g_gaugefield_data.eodir[isOdd][dir/2];
			if ((void*)eodata <= (void*)site && (void*)site < (void*)(eodata + GAUGE_EOVOLUME)) {
				// found!, get equivalent debug site
				size_t offset = (char*)site - (char*)eodata;
				char *debugdata = (char*)g_gaugefield_data_debug.eodir[isOdd][dir/2];
				bgq_gaugecoord *result = (bgq_gaugecoord*)(debugdata + offset);
				return result;
			}
		}
	}

	assert(!"Pointer does not point to gaugefield");
	return NULL;
}


static void bgq_gaugefield_checkcoord(bgq_gaugefield gaugefield, bool isOdd, int t, int x, int y, int z, int tv, int k, direction dir, int i, int l, COMPLEX_PRECISION *value) {
	assert( sizeof(COMPLEX_PRECISION) == sizeof(bgq_gaugecoord) );

	bgq_gaugecoord *coord = bgq_gaugefield_coordref(value);
	if (coord->coord.init) {
		assert(coord->coord.isOdd == isOdd);
		assert(coord->coord.t == t+1);
		assert(coord->coord.x == x+1);
		assert(coord->coord.y == y+1);
		assert(coord->coord.z == moddown(z,LOCAL_LZ));
		assert(coord->coord.dir == dir);
		assert(coord->coord.i == i);
		assert(coord->coord.l == l);
	} else {
		assert(coord->val == 0);

		coord->coord.init = true;
		coord->coord.isOdd = isOdd;
		coord->coord.t = t+1;
		coord->coord.x = x+1;
		coord->coord.y = y+1;
		coord->coord.z = moddown(z,LOCAL_LZ);
		coord->coord.dir = dir;
		coord->coord.i = i;
		coord->coord.l = l;
		coord->coord.writes = 0;
		coord->coord.reads = 0;
	}
}

static void bgq_gaugefield_resetcoord_checkval(bgq_gaugefield gaugefield, bool isOdd, int t, int x, int y, int z, int tv, int k, direction dir, int i, int l,
		COMPLEX_PRECISION *val,
		int expected_reads_min, int expected_reads_max, int expected_writes_min, int expected_writes_max) {

	bgq_gaugecoord *coord = bgq_gaugefield_coordref(val);
	int reads = 0;
	int writes = 0;
	if (coord->coord.init) {
		reads = coord->coord.reads;
		writes = coord->coord.writes;
	}
	if (expected_reads_min >= 0)
		assert(reads >= expected_reads_min);
	if (expected_reads_max >= 0)
		assert(expected_reads_max >= reads);
	if (expected_writes_min >= 0)
		assert(writes >= expected_writes_min);
	if (expected_writes_max >= 0)
		assert(expected_writes_max >= writes);

	bgq_gaugefield_checkcoord(gaugefield,isOdd,t,x,y,z,tv,k,dir,i,l,val);
	coord->coord.writes = 0;
	coord->coord.reads = 0;
}

void bgq_gaugefield_resetcoord(bgq_gaugefield gaugefield, int expected_reads_min, int expected_reads_max, int expected_writes_min, int expected_writes_max) {
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
		const int teo = divdown(t, PHYSICAL_LP);
		const int tv = divdown(teo, PHYSICAL_LK);
		const int k = moddown(teo, PHYSICAL_LK);

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
						COMPLEX_PRECISION *value = BGQ_GAUGEVAL(gaugefield,isOdd,t,x,y,z,tv,k, dir,i,l,false,false);
						bgq_gaugefield_resetcoord_checkval(gaugefield, isOdd, t,x,y,z, tv,k, dir, i, l, value, expected_reads_min, expected_reads_max, expected_writes_min, expected_writes_max);
					}

					if (dir == ZUP && z == LOCAL_LZ-1) {
						COMPLEX_PRECISION *wrapvalue = BGQ_GAUGEVAL(gaugefield,isOdd,t,x,y,-1,tv,k, dir,i,l,false,false);
						bgq_gaugefield_resetcoord_checkval(gaugefield, isOdd, t,x,y,-1, tv,k, dir, i, l, wrapvalue, expected_reads_min, expected_reads_max, expected_writes_min, expected_writes_max);
					}

					if (dir == TUP) {
						const int teo_shift = moddown(teo+2, 1+LOCAL_LT/PHYSICAL_LP) - 1;
						const int tv_shift = divdown(teo_shift,PHYSICAL_LK);
						const int k_shift = moddown(teo_shift, PHYSICAL_LK);

						COMPLEX_PRECISION *shiftvalue = BGQ_GAUGEVAL(gaugefield, isOdd, t, x, y, z, tv_shift, k_shift, TUP_SHIFT, i, l, false, false);
						bgq_gaugefield_resetcoord_checkval(gaugefield, isOdd, t, x, y, z, tv_shift, k_shift, TUP_SHIFT, i, l, shiftvalue, expected_reads_min, expected_reads_max, expected_writes_min, expected_writes_max);
					}
				}
			}
		}
	}
}
#endif


void bgq_init_gaugefield() {
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (direction dir = TUP; dir <= TUP_SHIFT; dir += 2) {
			g_gaugefield_data.eodir[isOdd][dir/2] = malloc_aligned(GAUGE_EOVOLUME * sizeof(bgq_gaugesite), 128 /*L2 cache line size*/);
		}
	}

	g_gaugefield = &g_gaugefield_data;

#if BGQ_FIELD_COORDCHECK
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (direction dir = TUP; dir <= TUP_SHIFT; dir += 2) {
			g_gaugefield_data_debug.eodir[isOdd][dir/2] = malloc_aligned(GAUGE_EOVOLUME * sizeof(bgq_gaugesite), 128 /*L2 cache line size*/);
		}
	}

	bgq_gaugefield_resetcoord(g_gaugefield, -1,-1,-1,-1);
#endif
}




void bgq_free_gaugefield() {
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (int dir = TUP; dir <= TUP_SHIFT; dir += 2) {
			free(g_gaugefield_data.eodir[isOdd][dir/2]);
			g_gaugefield_data.eodir[isOdd][dir/2] = NULL;
		}
	}
	g_gaugefield = NULL;

#if BGQ_FIELD_COORDCHECK
	for (int isOdd = false; isOdd <= true; isOdd += 1) {
		for (int dir = TUP; dir <= TUP_SHIFT; dir += 2) {
			free(g_gaugefield_data_debug.eodir[isOdd][dir/2]);
			g_gaugefield_data_debug.eodir[isOdd][dir/2] = NULL;
		}
	}
#endif
}
typedef struct {
	double _Complex c[3][3];
} su3_array64;

void bgq_transfer_gaugefield(bgq_gaugefield const targetfield, su3 ** const sourcefield) {
	assert(targetfield);
	assert(sourcefield);

#if BGQ_FIELD_COORDCHECK
	bgq_gaugefield_resetcoord(targetfield, -1,-1,-1,-1);
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
					bgq_gaugefield_set(targetfield, isOdd, t, x, y, z, dir, i, l, m->c[i][l]);
				}
			}
		}
	}

#if BGQ_FIELD_COORDCHECK
	bgq_gaugefield_resetcoord(targetfield, 0,0,1,1);
#endif
}



bool assert_gaugesite(bgq_gaugefield gaugefield, bool isOdd, int t, int x, int y, int z, int tv, int k, direction dir, bool isRead, bool isWrite) {
	for (int i = 0; i < 3; i += 1) {
		for (int l = 0; l < 3; l += 1) {
			assert_gaugeval(gaugefield, isOdd, t, x, y, z, tv, k, dir, i, l, isRead, isWrite);
		}
	}

	return true;
}


bool assert_gaugeval(bgq_gaugefield gaugefield, bool isOdd, int t, int x, int y, int z, int tv, int k, direction dir, int i, int l, bool isRead, bool isWrite) {
	assert(gaugefield);
	assert(false <= isOdd && isOdd <= true);
	assert(-1 <= t && t < LOCAL_LT);
	assert(-1 <= x && x < LOCAL_LX);
	assert(-1 <= y && y < LOCAL_LY);
	assert(-1 <= z && z < PHYSICAL_LZ);
	assert(-1 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= k && k < PHYSICAL_LK);
	assert(TUP <= dir && dir <= TDOWN_SHIFT);
	assert(0 <= i && i < 3);
	assert(0 <= l && l < 3);

	// There is really just one gaugefield
	assert( gaugefield == g_gaugefield );

	// We really only store the up-fields
	assert((dir&1)==0);

	// Check that the coordinate is really an odd/even coordinate
	assert( ((t+x+y+z)&1) == isOdd );
	int teo = divdown(t,PHYSICAL_LP);
	if (dir == TUP_SHIFT) {
		teo = mod(teo+2, 1+LOCAL_LT/PHYSICAL_LP)-1;
	}

	// Check that zv and k match the coordinate
	assert(divdown(teo,PHYSICAL_LK) == tv);
	assert(moddown(teo,PHYSICAL_LK) == k);

	// Get the index
	const int idx = ((((tv)+1)*(PHYSICAL_LX+1) + ((x)+1))*(PHYSICAL_LY+1) + ((y)+1))*(PHYSICAL_LZ+1) + ((z)+1);
	assert((0 <= idx) && (idx < GAUGE_EOVOLUME));

	// Get the address
	bgq_gaugesite *eofield = gaugefield->eodir[isOdd][dir/2];
	bgq_gaugesite *site = &eofield[idx];
	//COMPLEX_PRECISION *address = &site->c[i][l][k];
	assert(mod((size_t)&site->c[i][l][0],PRECISION_VECTOR_ALIGNMENT)==0);


#if BGQ_FIELD_COORDCHECK
	// get the debug data
	bgq_gaugefield_checkcoord(gaugefield,isOdd,t,x,y,moddown(z,LOCAL_LZ),tv,k,dir,i,l,address);
	bgq_gaugecoord *debugdata = bgq_gaugefield_coordref(address);
	if (isRead)
		debugdata->coord.reads += 1;
	if (isWrite)
		debugdata->coord.writes += 1;
#endif


	// All checks passed
	return true;
}




////////////////////////////////////////////////////////////////////////////////
// Weylfields

typedef union {
	COMPLEX_PRECISION val;
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

bgq_weylfield weylxchange_recv[6];
bgq_weylfield weylxchange_send[6];
size_t weylxchange_size[3];
int weylexchange_destination[6];
MPI_Request weylexchange_request_recv[6];
MPI_Request weylexchange_request_send[6];

#if BGQ_FIELD_COORDCHECK
static bgq_weylcoord *weylxchange_recv_debug[PHYSICAL_LP][6];
static bgq_weylcoord *weylxchange_send_debug[PHYSICAL_LP][6];
#endif

//static bool l1p_first = true;



//double recvbuf;
MPI_Request recvrequest;
//double sendbuf;
MPI_Request sendrequest;


void bgq_hm_init() {
	weylxchange_size[TUP / 2] = PHYSICAL_LXV * PHYSICAL_LY * PHYSICAL_LZ * sizeof(bgq_weylsite);
	weylxchange_size[XUP / 2] = PHYSICAL_LTV * PHYSICAL_LY * PHYSICAL_LZ * sizeof(bgq_weylsite);
	weylxchange_size[YUP / 2] = PHYSICAL_LTV * PHYSICAL_LX * PHYSICAL_LZ * sizeof(bgq_weylsite);

	weylexchange_destination[TUP] = g_nb_t_up;
	weylexchange_destination[TDOWN] = g_nb_t_dn;
	weylexchange_destination[XUP] = g_nb_x_up;
	weylexchange_destination[XDOWN] = g_nb_x_dn;
	weylexchange_destination[YUP] = g_nb_y_up;
	weylexchange_destination[YDOWN] = g_nb_y_dn;

	//master_print("g_proc_id=%d g_cart_id=%d\n", g_proc_id, g_cart_id);
	for (direction d = TUP; d <= YDOWN; d += 1) {
		size_t size = weylxchange_size[d/2];
		weylxchange_recv[d] = (bgq_weylfield)malloc_aligned(size, 128);
		weylxchange_send[d] = (bgq_weylfield)malloc_aligned(size, 128);
		#if BGQ_FIELD_COORDCHECK
				for (int isOdd = false; isOdd <= true; isOdd += 1) {
					weylxchange_recv_debug[isOdd][d] = (bgq_weylcoord*) malloc_aligned(size, 128);
					weylxchange_send_debug[isOdd][d] = (bgq_weylcoord*) malloc_aligned(size, 128);
				}
		#endif

		MPI_CHECK(MPI_Recv_init(weylxchange_recv[d], size / PRECISION_BYTES, MPI_PRECISION, weylexchange_destination[d], d^1, g_cart_grid, &(weylexchange_request_recv[d])));
		MPI_CHECK(MPI_Send_init(weylxchange_send[d], size / PRECISION_BYTES, MPI_PRECISION, weylexchange_destination[d], d, g_cart_grid, &(weylexchange_request_send[d])));
	}

	//MPI_CHECK(MPI_Recv_init(weylxchange_recv[0], weylxchange_size[0] / PRECISION_BYTES, MPI_PRECISION, g_nb_t_dn, 0, g_cart_grid, &recvrequest));
	//MPI_CHECK(MPI_Send_init(weylxchange_send[0], weylxchange_size[0] / PRECISION_BYTES, MPI_PRECISION, g_nb_t_up, 0, g_cart_grid, &sendrequest));
}


void bgq_hm_free() {
	for (direction d = TUP; d <= YDOWN; d += 1) {
		free(weylxchange_recv[d]);
		weylxchange_recv[d] = NULL;
		free(weylxchange_send[d]);
		weylxchange_send[d] = NULL;
#if BGQ_FIELD_COORDCHECK
		for (int isOdd = false; isOdd <= true; isOdd += 1) {
			free(weylxchange_recv_debug[isOdd][d]);
			weylxchange_recv_debug[isOdd][d] = NULL;
			free(weylxchange_send_debug[isOdd][d]);
			weylxchange_send_debug[isOdd][d] = NULL;
		}
#endif

		MPI_Request_free(&weylexchange_request_recv[d]);
		MPI_Request_free(&weylexchange_request_send[d]);
	}

#if BGQ_PREFETCH_LIST
#pragma omp parallel
	{
		L1P_PatternUnconfigure();
	}
#endif
}


#if BGQ_FIELD_COORDCHECK
static bgq_weylcoord *bgq_weylfield_coordref(bool isOdd, COMPLEX_PRECISION *val) {
	assert( sizeof(bgq_weylcoord) == sizeof(COMPLEX_PRECISION) );
	// Logically, odd and even sites are different locations
	// But in HoppingMatrix, only one of them is used in a time, odd and even computation/communication do not interleave
	// Therefore we coalesce even and odd fields to the same physical memory location
	// But logically these sites are still distinct, so we need different debug infos for them

	for (direction d = TUP; d <= YDOWN; d += 1) {
		if ( ((char*)weylxchange_recv[d] <= (char*)val) && ((char*)val < (char*)weylxchange_recv[d] + weylxchange_size[d/2]) ) {
			size_t offset = (char*)val - (char*)(weylxchange_recv[d]);
			return (bgq_weylcoord*)((char*)weylxchange_recv_debug[isOdd][d] + offset);
		}

		if ( ((char*)weylxchange_send[d] <= (char*)val) && ((char*)val < (char *)weylxchange_send[d] + weylxchange_size[d/2]) ) {
			size_t offset = (char*)val - (char*)(weylxchange_send[d]);
			return (bgq_weylcoord*)((char*)weylxchange_send_debug[isOdd][d] + offset);
		}
	}

	master_error(1, "Unknown Weylfield\n");
	return NULL;
}

static void bgq_weylfield_checkcoord(bgq_weylfield weylfield, bool isOdd, int t, int x, int y, int z, int v, int c, COMPLEX_PRECISION *value) {
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
		assert(coord->val == 0);

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
#endif




void bgq_weylfield_t_resetcoord(bgq_weylfield weylfield, int t, bool isOdd, int expected_reads_min, int expected_reads_max, int expected_writes_min, int expected_writes_max) {
	assert(weylfield);

#pragma omp parallel for schedule(static)
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
				/*COMPLEX_PRECISION *value = */BGQ_WEYLVAL_T(weylfield,isOdd,t,x,y,z,xv,k,v,c,false,false);

#if BGQ_FIELD_COORDCHECK
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
#endif
			}
		}
	}
}



bool assert_weylval_t(bgq_weylfield weylfield, bool isOdd, int t, int x, int y, int z, int xv, int k, int v, int c, bool isRead, bool isWrite) {
	if (!isRead && !isWrite && (z == LOCAL_LZ)) {
		// Allow for prefetching
		return true;
	}

	assert(weylfield);
	assert(false <= isOdd && isOdd <= true); // bogus
	assert( (t == -1) || (t == 0) || (t == LOCAL_LT-1) || (t == LOCAL_LT) ); /* We are one out of the volume, either up or down */
	assert(0 <= x && x < LOCAL_LX);
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	assert(0 <= xv && xv < PHYSICAL_LXV);
	assert(0 <= k && k < PHYSICAL_LK);

	// Only 4 fields of this type actually exist
	assert( (weylfield == weylxchange_recv[TUP]) || (weylfield == weylxchange_recv[TDOWN])
			|| (weylfield == weylxchange_send[TUP]) || (weylfield == weylxchange_send[TDOWN]) );

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
	bgq_weylsite *site = &weylfield[idx];
	assert(&weylfield[0] <= site && site < &weylfield[PHYSICAL_LXV*PHYSICAL_LY*PHYSICAL_LZ]);
	assert(mod((size_t)&site->s[v][c][0],PRECISION_VECTOR_ALIGNMENT)==0);

	COMPLEX_PRECISION *val = &site->s[v][c][k];
	assert(&weylfield[0].s[0][0][0] <= val && val < &weylfield[PHYSICAL_LXV*PHYSICAL_LY*PHYSICAL_LZ].s[0][0][0]);

#if BGQ_FIELD_COORDCHECK
	bgq_weylcoord *coord = bgq_weylfield_coordref(isOdd, val);
	bgq_weylfield_checkcoord(weylfield, isOdd, t, x, y, z, v, c, val);

	assert(coord->coord.init);
	if (isRead)
		coord->coord.reads += 1;
	if (isWrite)
		coord->coord.writes += 1;
#endif

	return true;
}


bool assert_weylfield_t(bgq_weylfield weylfield, bool isOdd, int t, int x, int y, int z, int xv, int k, bool isRead, bool isWrite) {
	for (int v = 0; v < 2; v += 1) {
		for (int c = 0; c < 3; c += 1) {
			assert_weylval_t(weylfield,isOdd,t,x,y,z,xv,k,v,c,isRead,isWrite);
		}
	}
	return true;
}


void bgq_weylfield_x_resetcoord(bgq_weylfield weylfield, int x, bool isOdd, int expected_reads_min, int expected_reads_max, int expected_writes_min, int expected_writes_max) {
	assert(weylfield);

#pragma omp parallel for schedule(static)
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
				/*COMPLEX_PRECISION *value = */BGQ_WEYLVAL_X(weylfield,isOdd,t,x,y,z,tv,k,v,c,false,false);
#if BGQ_FIELD_COORDCHECK
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
#endif
			}
		}
	}
}



bool assert_weylval_x(bgq_weylfield weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k, int v, int c, bool isRead, bool isWrite) {
	assert(weylfield);
	assert(false <= isOdd && isOdd <= true);
	assert(0 <= t && t < LOCAL_LT);
	assert( (x == -1) || (x == 0) || (x == LOCAL_LX-1) || (x == LOCAL_LX) );
	assert(0 <= y && y < LOCAL_LY);
	assert(0 <= z && z < LOCAL_LZ);
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= k && k < PHYSICAL_LK);

	// Only 4 fields of this type actually exist
	assert((weylfield == weylxchange_recv[XUP]) || (weylfield == weylxchange_recv[XDOWN])
			|| (weylfield == weylxchange_send[XUP]) || (weylfield == weylxchange_send[XDOWN]) );

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
	bgq_weylsite *site = &weylfield[idx];
	assert(&weylfield[0] <= site && site < &weylfield[PHYSICAL_LTV*PHYSICAL_LY*PHYSICAL_LZ]);
	assert(mod((size_t)&site->s[v][c][0],PRECISION_VECTOR_ALIGNMENT)==0);

	COMPLEX_PRECISION *val = &site->s[v][c][k];
	assert(&weylfield[0].s[0][0][0] <= val && val < &weylfield[PHYSICAL_LTV*PHYSICAL_LY*PHYSICAL_LZ].s[0][0][0]);

#if BGQ_FIELD_COORDCHECK
	bgq_weylcoord *coord = bgq_weylfield_coordref(isOdd, val);
	bgq_weylfield_checkcoord(weylfield, isOdd, t, x, y, z, v, c, val);

	assert(coord->coord.init);
	if (isRead)
		coord->coord.reads += 1;
	if (isWrite)
		coord->coord.writes += 1;
#endif

	return true;
}


bool assert_weylfield_x(bgq_weylfield weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k, bool isRead, bool isWrite) {
	for (int v = 0; v < 2; v += 1) {
		for (int c = 0; c < 3; c += 1) {
			assert_weylval_x(weylfield,isOdd,t,x,y,z,tv,k,v,c,isRead,isWrite);
		}
	}
	return true;
}



void bgq_weylfield_y_resetcoord(bgq_weylfield weylfield, int y, bool isOdd, int expected_reads_min, int expected_reads_max, int expected_writes_min, int expected_writes_max) {
	assert(weylfield);

#pragma omp parallel for schedule(static)
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
				/*COMPLEX_PRECISION *value = */BGQ_WEYLVAL_Y(weylfield,isOdd,t,x,y,z,tv,k,v,c,false,false);
#if BGQ_FIELD_COORDCHECK
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
#endif
			}
		}
	}
}


bool assert_weylval_y(bgq_weylfield weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k, int v, int c, bool isRead, bool isWrite) {
	assert(weylfield);
	assert(false <= isOdd && isOdd <= true);
	assert(0 <= t && t < LOCAL_LT);
	assert(0 <= x && x < LOCAL_LX);
	assert( (-1 == y) || (y == 0) || (y == LOCAL_LY-1) || (y == LOCAL_LY) );
	assert(0 <= z && z < LOCAL_LZ);
	assert(0 <= tv && tv < PHYSICAL_LTV);
	assert(0 <= k && k < PHYSICAL_LK);

	// Only 4 fields of this type actually exist
	assert((weylfield == weylxchange_recv[YUP]) || (weylfield == weylxchange_recv[YDOWN])
			|| (weylfield == weylxchange_send[YUP]) || (weylfield == weylxchange_send[YDOWN]));

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
	bgq_weylsite *site = &weylfield[idx];
	assert(&weylfield[0] <= site && site < &weylfield[PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LZ]);
	assert(mod((size_t)&site->s[v][c][0],PRECISION_VECTOR_ALIGNMENT)==0);

	COMPLEX_PRECISION *val = &site->s[v][c][k];
	assert(&weylfield[0].s[0][0][0] <= val && val < &weylfield[PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LZ].s[0][0][0]);

#if BGQ_FIELD_COORDCHECK
	bgq_weylcoord *coord = bgq_weylfield_coordref(isOdd, val);
	bgq_weylfield_checkcoord(weylfield, isOdd, t, x, y, z, v, c, val);

	assert(coord->coord.init);
	if (isRead)
		coord->coord.reads += 1;
	if (isWrite)
		coord->coord.writes += 1;
#endif

	return true;
}


bool assert_weylfield_y(bgq_weylfield weylfield, bool isOdd, int t, int x, int y, int z, int tv, int k, bool isRead, bool isWrite) {
	for (int v = 0; v < 2; v += 1) {
		for (int c = 0; c < 3; c += 1) {
			assert_weylval_y(weylfield,isOdd,t,x,y,z,tv,k,v,c,isRead,isWrite);
		}
	}
	return true;
}


static void bgq_weylfield_t_foreach(bgq_weylfield weylfield, direction dir, bool isSend, bool isOdd, bgq_weylsite_callback callback, int tag) {
	assert(weylfield);
	assert( (dir==TUP) || (dir==TDOWN) );
	assert(callback);

	int t;
	int tsrc;
	switch (dir) {
	case TUP:
		t = isSend ? LOCAL_LT : LOCAL_LT-1;
		tsrc = isSend ? LOCAL_LT-1 : LOCAL_LT;
		break;
	case TDOWN:
		t = isSend ? -1 : 0;
		tsrc = isSend ? 0 : -1;
		break;
	default:
		assert(!"Wrong function for direction");
		return;
	}

	#pragma omp for schedule(static)
	for (int tyz = 0; tyz < PHYSICAL_LXV*PHYSICAL_LY*PHYSICAL_LZ; tyz+=1) {
		WORKLOAD_DECL(tyz,PHYSICAL_LXV*PHYSICAL_LY*PHYSICAL_LZ);
		const int xv = WORKLOAD_PARAM(PHYSICAL_LXV);
		const int y = WORKLOAD_PARAM(PHYSICAL_LY);
		const int z = WORKLOAD_PARAM(PHYSICAL_LZ);
		WORKLOAD_CHECK

		const int x1 = ((isOdd+t+y+z)&1)+xv*PHYSICAL_LP*PHYSICAL_LK;
		const int x2 = x1 + 2;

		bgq_weylsite *weylsite = BGQ_WEYLSITE_T(weylfield, !isOdd, tsrc, xv, y, z, x1, x2, false, false);

		for (int v = 0; v < 2; v+=1) {
			for (int c = 0; c < 3; c+=1) {
				(*callback)(weylfield, dir, isSend, isOdd, t,x1,y,z,v,c, &weylsite->s[v][c][0], tag);
				(*callback)(weylfield, dir, isSend, isOdd, t,x2,y,z,v,c, &weylsite->s[v][c][1], tag);
			}
		}
	}
}

static void bgq_weylfield_x_foreach(bgq_weylfield weylfield, direction dir, bool isSend, bool isOdd, bgq_weylsite_callback callback, int tag) {
	assert(weylfield);
	assert( (dir==XUP) || (dir==XDOWN) );
	assert(callback);

	int x;
	int xsrc;
	switch (dir) {
	case XUP:
		x = isSend ? LOCAL_LX : LOCAL_LX-1;
		xsrc = isSend ? LOCAL_LX-1 : LOCAL_LX;
		break;
	case XDOWN:
		x = isSend ? -1 : 0;
		xsrc = isSend ? 0 : -1;
		break;
	default:
		assert(!"Wrong function for direction");
		return;
	}

	#pragma omp for schedule(static)
	for (int tyz = 0; tyz < PHYSICAL_LTV*PHYSICAL_LY*PHYSICAL_LZ; tyz+=1) {
		WORKLOAD_DECL(tyz,PHYSICAL_LTV*PHYSICAL_LY*PHYSICAL_LZ);
		const int tv = WORKLOAD_PARAM(PHYSICAL_LTV);
		const int y = WORKLOAD_PARAM(PHYSICAL_LY);
		const int z = WORKLOAD_PARAM(PHYSICAL_LZ);
		WORKLOAD_CHECK

		const int t1 = ((isOdd+x+y+z)&1)+tv*PHYSICAL_LP*PHYSICAL_LK;
		const int t2 = t1 + 2;

		bgq_weylsite *weylsite = BGQ_WEYLSITE_X(weylfield, !isOdd, tv, xsrc, y, z, t1, t2, false, false);

		for (int v = 0; v < 2; v+=1) {
			for (int c = 0; c < 3; c+=1) {
				(*callback)(weylfield, dir, isSend, isOdd, t1, x,y,z,v,c, &weylsite->s[v][c][0], tag);
				(*callback)(weylfield, dir, isSend, isOdd, t2, x,y,z,v,c, &weylsite->s[v][c][1], tag);
			}
		}
	}
}

static void bgq_weylfield_y_foreach(bgq_weylfield weylfield, direction dir, bool isSend, bool isOdd, bgq_weylsite_callback callback, int tag) {
	assert(weylfield);
	assert( (dir==YUP) || (dir==YDOWN) );
	assert(callback);

	int y;
	int ysrc;
	switch (dir) {
	case YUP:
		y = isSend ? LOCAL_LY : LOCAL_LY-1;
		ysrc = isSend ? LOCAL_LY-1 : LOCAL_LY;
		break;
	case YDOWN:
		y = isSend ? -1 : 0;
		ysrc = isSend ? 0 : -1;
		break;
	default:
		assert(!"Wrong function for direction");
		return;
	}

	#pragma omp for schedule(static)
	for (int txz = 0; txz < PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LZ; txz+=1) {
		WORKLOAD_DECL(txz,PHYSICAL_LTV*PHYSICAL_LX*PHYSICAL_LZ);
		const int tv = WORKLOAD_PARAM(PHYSICAL_LTV);
		const int x = WORKLOAD_PARAM(PHYSICAL_LX);
		const int z = WORKLOAD_PARAM(PHYSICAL_LZ);
		WORKLOAD_CHECK

		const int t1 = ((isOdd+x+y+z)&1)+tv*PHYSICAL_LP*PHYSICAL_LK;
		const int t2 = t1 + 2;

		bgq_weylsite *weylsite = BGQ_WEYLSITE_Y(weylfield, !isOdd, tv, x, ysrc, z, t1, t2, false, false);

		for (int v = 0; v < 2; v+=1) {
			for (int c = 0; c < 3; c+=1) {
				(*callback)(weylfield, dir, isSend, isOdd, t1, x,y,z,v,c, &weylsite->s[v][c][0], tag);
				(*callback)(weylfield, dir, isSend, isOdd, t2, x,y,z,v,c, &weylsite->s[v][c][1], tag);
			}
		}
	}
}

void bgq_weylfield_foreach(bgq_weylfield weylfield, direction dir, bool isSend, bool isOdd, bgq_weylsite_callback callback, int tag) {
	switch (dir) {
	case TUP:
	case TDOWN:
		bgq_weylfield_t_foreach(weylfield,dir,isSend,isOdd,callback,tag);
		break;
	case XUP:
	case XDOWN:
		bgq_weylfield_x_foreach(weylfield,dir,isSend,isOdd,callback,tag);
		break;
	case YUP:
	case YDOWN:
		bgq_weylfield_y_foreach(weylfield,dir,isSend,isOdd,callback,tag);
		break;
	default:
		assert(!"Not yet implemented");
		break;
	}
}


////////////////////////////////////////////////////////////////////////////////


