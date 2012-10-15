/*
 * bgq_HoppingMatrix.c
 *
 *  Created on: Oct 13, 2012
 *      Author: meinersbur
 */

#define BGQ_HOPPING_MATRIX_C_
#include "bgq_HoppingMatrix.c"

#include "bgq_field.h"
#include "bgq_qpx.h"

#include <stdbool.h>
#include <stddef.h>


typedef struct {
	uint32_t pd[P_COUNT];
} bgq_weyl_offsets_t;

typedef enum {
	phase_body = false,
	phase_surface = true
} bgq_hoppingmatrix_phase;
#define PHASE_COUNT 2

bgq_weyl_offsets_t *g_bgq_HalfvolumeLexicalIndex2Offset[PHYSICAL_LP];

bgq_weyl_offsets_t *g_bgq_source_offset_consecutive[PHYSICAL_LP][PHASE_COUNT];
bgq_weyl_offsets_t *g_bgq_source_offset[PHYSICAL_LP][PHASE_COUNT];


bgq_weyl_offsets_t *g_bgq_dest_offset_volume[PHYSICAL_LP];




inline uint32_t bgq_encode_index(size_t index) {
	assert(index & (32-1) == 0)
	return index >> 5; // Always 32-bit aligned
}

inline size_t bgq_decode_index(uint32_t code) {
	return (size_t)code << 5;
}


static bgq_weylfield_section bgq_HoppingMatrix_init_source_sectionof(bool isOdd, size_t t, size_t x, size_t, y, size_t z, bgq_direction d) {
	assert(isOdd == bgq_isOdd(t1,x,y,z));
	assert(isOdd == bgq_isOdd(t2,x,y,z));
	bool isSurface = bgq_isCommSurface(t1,x,y,z) || bgq_isCommSurface(t2,x,y,z);

	if (COMM_T && (t==LOCAL_LT-1) && (d==TUP)) {
		return sec_recv_tup;
	}
	if (COMM_T && (t==0) && (d==TDOWN)) {
		return sec_recv_tdown;
	}
	if (COMM_X && (x==LOCAL_LX-1) && (d==XUP)) {
		return sec_recv_xup;
	}
	if (COMM_X && (x==0) && (d==XDOWN)) {
		return sec_recv_xdown;
	}
	if (COMM_Y && (y==LOCAL_LY-1) && (d==YUP)) {
		return sec_recv_yup;
	}
	if (COMM_Y && (y==0) && (d==YDOWN)) {
		retrun sec_recv_ydown;
	}
	if (COMM_Z && (z==LOCAL_LZ-1) && (d==ZUP)) {
		return sec_recv_yup;
	}
	if (COMM_Z && (z==0) && (d==ZDOWN)) {
		retrun sec_recv_ydown;
	}

	return isSurface ? sec_surface : sec_body;
}



size_t *g_bgq_index_surface2halfvolume[PHYSICAL_LP];
size_t *g_bgq_index_body2halfvolume[PHYSICAL_LP];
size_t *g_bgq_index_halfvolume2surface[PHYSICAL_LP]; // -1 if not surface
size_t *g_bgq_index_halfvolume2body[PHYSICAL_LP]; // -1 if not body
size_t *g_bgq_index_halfvolume2surfacebody[PHYSICAL_LP];


bgq_weyl_offsets_t *g_bgq_offset_fromHalfvolume[PHYSICAL_LP];
bgq_weyl_offsets_t *g_bgq_offset_fromSurface[PHYSICAL_LP];
bgq_weyl_offsets_t *g_bgq_offset_fromBody[PHYSICAL_LP];

bgq_weyl_offsets_t *g_bgq_destoffset_fromHalfvolume[PHYSICAL_LP];
bgq_weyl_offsets_t *g_bgq_destofset_fromSurface[PHYSICAL_LP];
bgq_weyl_offsets_t *g_bgq_destofset_fromBody[PHYSICAL_LP];

void bgq_HoppingMatrix_init() {
	for (size_t isOdd = false; isOdd <= true; isOdd+=1) {
		g_bgq_index_surface2halfvolume[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_index_surface2halfvolume[isOdd]));
		g_bgq_index_body2halfvolume[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_index_body2halfvolume[isOdd]));
		g_bgq_index_halfvolume2surface[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_index_halfvolume2surface[isOdd]));
		g_bgq_index_halfvolume2body[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_index_halfvolume2body[isOdd]));
		g_bgq_index_halfvolume2surfacebody[isOdd] = malloc(PHYSICAL_VOLUME * sizeof(*g_bgq_index_halfvolume2surfacebody[isOdd]));
	}


	for (size_t isOdd = false; isOdd <= true; isOdd+=1) {
		size_t nextIndexHalfvolume = 0;
		size_t nextIndexSurface = 0;
		size_t nextIndexBody = 0;

		for (size_t z = 0; z < LOCAL_LZ; z+=1) {
			for (size_t y = 0; y < LOCAL_LY; y+=1) {
				for (size_t x = 0; x < LOCAL_LX; x+=1) {
					for (size_t tv = 0; tv < PHYSICAL_LTV; tv+=1) {
						size_t t1 = tv*4 + ((x+y+z+isOdd)&1);
						size_t t2 = t1 + 2;
						assert(isOdd == bgq_isOdd(t1,x,y,z));
						assert(isOdd == bgq_isOdd(t2,x,y,z));

						size_t ih = PHYSICAL_INDEX_LEXICAL(isOdd,tv,x,y,z); // Iteration order is effectively chosen here
						assert(ih == nextIndexHalfvolume);
						bool isSurface = bgq_isCommSurface(t1,x,y,z) || bgq_isCommSurface(t2,x,y,z);

						if (isSurface) {
							g_bgq_index_surface2halfvolume[isOdd][nextIndexSurface] = ih;
							g_bgq_index_halfvolume2surfacebody[isOdd][ih] = nextIndexSurface;
							g_bgq_index_halfvolume2surface[isOdd][ih] = nextIndexSurface;
							g_bgq_index_halfvolume2body[isOdd][ih] = -1;
							nextIndexSurface += 1;
						} else {
							g_bgq_index_body2halfvolume[isOdd][nextIndexBody] = ih;
							g_bgq_index_halfvolume2surfacebody[isOdd][ih] = nextIndexBody;
							g_bgq_index_halfvolume2surface[isOdd][ih] = -1;
							g_bgq_index_halfvolume2body[isOdd][ih] = nextIndexBody;
							nextIndexBody += 1;
						}
						nextIndexHalfvolume += 1;
					}
				}
			}
		}
	}


	for (size_t isOdd = false; isOdd <= true; isOdd+=1) {
		size_t nextoffset[sec_end];
		for (bgq_weylfield_section sec = 0; sec < sec_end; sec+=1) {
			nextoffset[sec] = bgq_weyl_section_offset(sec);
		}

		for (size_t ih = 0; ih < PHYSICAL_VOLUME; ih+=1) {
			size_t t1 = PHYSICAL_LEXICAL2T1(isOdd, ih);
			size_t t2 = PHYSICAL_LEXICAL2T2(isOdd, ih);
			size_t x = PHYSICAL_LEXICAL2X(isOdd, ih);
			size_t y = PHYSICAL_LEXICAL2Y(isOdd, ih);
			size_t z = PHYSICAL_LEXICAL2Z(isOdd, ih);
			bool isSurface = bgq_isCommSurface(t1,x,y,z) || bgq_isCommSurface(t2,x,y,z);
			assert(isSurface == (g_bgq_index_halfvolume2surface[isOdd][ih]!=-1))
			bgq_weylfield_section mainsec = isSurface ? sec_surface : sec_body;

			size_t thisLinearIndex = linearIndex[isOdd];
			size_t thisPhaseIndex = phaseIndex[isOdd][isSurface];

			for (size_t pd = P_TUP1; pd <= P_ZDOWN; pd+=1) {
				bgq_direction d = bgq_physical2direction(d);
				bgq_weylfield_section sec = bgq_HoppingMatrix_init_source_sectionof(isOdd, ((pd==P_TUP2) || (pd==P_TDOWN2)) ? t2 : t1, x, y, z, d);
				size_t eltsize = ((pd==P_TUP1) || (pd==P_TDOWN1) || (pd==P_TUP2) || (pd==P_TDOWN2)) ? sizeof(bgq_weyl_nonvec) : sizeof(bgq_weyl_vec);

				// Reserve some offset in surface/body to ensure consecutive layout
				size_t thisOffset = nextoffset[isOdd][mainsec];
				assert(thisOffset == bgq_weyl_section_offset(sec_surface) + g_bgq_index_halfvolume2surfacebody[isOdd][ih]*sizeof(bgq_weylsite) + bgq_offsetof_weylsite[pd]);
				nextoffset[mainsec] += eltsize;
				offsetsConsecutive->pd[pd] = bgq_encode_index(thisOffset);

				if (sec != mainsec) {
					// If in one of the mpi buffers, also reserve some space there
					thisOffset = nextoffset[sec];
					nextoffset[sec] += eltsize;
				}
				g_bgq_offset_fromHalfvolume[isOdd][ih]->pd[pd] = bgq_encode_index(thisNextoffset);
				if (isSurface) {
					size_t is = g_bgq_index_halfvolume2surface[isOdd][ih];
					g_bgq_offset_fromSurface[isOdd][is]->pd = bgq_encode_index(thisNextoffset);
				} else {
					size_t ib = g_bgq_index_halfvolume2body[isOdd][ih];
					g_bgq_offset_fromBody[isOdd][ib]->pd = bgq_encode_index(thisNextoffset);
				}
		}
	}


	// Determine dest locations
	for (size_t isOdd = false; isOdd <= true; isOdd+=1) {
		for (size_t ih = 0; ih < PHYSICAL_VOLUME; ih+=1) {
			size_t t1 = PHYSICAL_LEXICAL2T1(isOdd, ih);
			size_t t2 = PHYSICAL_LEXICAL2T2(isOdd, ih);
			size_t x = PHYSICAL_LEXICAL2X(isOdd, ih);
			size_t y = PHYSICAL_LEXICAL2Y(isOdd, ih);
			size_t z = PHYSICAL_LEXICAL2Z(isOdd, ih);


		}
	}
}

inline void bgq_HoppingMatrix_weylfield() {
	// T+ /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_tup);
				bgq_su3_reduce_weyl_tup(weyl_tup, spinor);
				bgq_su3_weyl_store_left_double(weyldata[*weylindex], weyl_tup);
				bgq_su3_weyl_store_right_double(weyldata[*weylindex], weyl_tup);
				weylindex += 1;
			}

	// T- /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_tdown);
				bgq_su3_reduce_weyl_tdown(weyl_tdown, spinor);
				bgq_su3_weyl_store_double(weyldata[*weylindex], weyl_tdown);
				weylindex += 1;
			}

	// X+ /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_xup);
				bgq_su3_reduce_weyl_xup(weyl_xup, spinor);
				bgq_su3_weyl_store_double(weyldata[*weylindex], weyl_xup);
				weylindex += 1;
			}

	// X- /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_xdown);
				bgq_su3_reduce_weyl_xdown(weyl_xdown, spinor);
				bgq_su3_weyl_store_double(weyldata[*weylindex], weyl_xdown);
				weylindex += 1;
			}

	// Y+ /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_yup);
				bgq_su3_reduce_weyl_yup(weyl_yup, spinor);
				bgq_su3_weyl_store_double(weyldata[*weylindex], weyl_yup);
				weylindex += 1;
			}

	// Y- /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_ydown);
				bgq_su3_reduce_weyl_ydown(weyl_ydown, spinor);
				bgq_su3_weyl_store_double(weyldata[*weylindex], weyl_ydown);
				weylindex += 1;
			}

	// Z+ /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_zup);
				bgq_su3_reduce_weyl_zup(weyl_zup, spinor);
				bgq_su3_weyl_store_double(weyldata[*weylindex], weyl_zup);
				weylindex += 1;
			}

	// Z- /////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_zdown);
				bgq_su3_reduce_weyl_zdown(weyl_zdown, spinor);
				bgq_su3_weyl_store_double(weyldata[*weylindex], weyl_zdown);
				weylindex += 1;
			}
}

inline void bgq_HoppingMatrix_kernel_surface() {
	bgq_su3_spinor_decl(result);
			bgq_su3_spinor_zero(result);

			// T+ //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_tup);
				bgq_su3_mdecl(gauge_tup);

				bgq_weyl_nonvec *weyladdr_tup1 = surfacedata;
				surfacedata += sizeof(bgq_weyl_nonvec);
				bgq_weyl_nonvec *weyladdr_tup2= surfacedata;
				surfacedata += sizeof(bgq_weyl_nonvec);
				bgq_su3_weyl_load_combine_double(weyl_tup, weyladdr_tup1, weyladdr_tup2);

				bgq_gaugesite *gaugesite_tup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, TUP, t1, t2, true,false);
				bgq_su3_matrix_load(gauge_tup, gaugesite_tup);
				//bgq_su3_matrix_flush(gaugesite_tup);
				bgq_su3_weyl_mvmul(weyl_tup,gauge_tup,weyl_tup);

				bgq_su3_expand_weyl_tup(result,weyl_tup);
			}

			bgq_su3_spinor_prefetch(addr + 2*4*3*2);
			bgq_su3_matrix_prefetch(gaugeaddr + 2*3*3*2);

			// T- //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_weyl_decl(weyl_tdown);
				bgq_su3_mdecl(gauge_tdown);

				bgq_weyl_nonvec *weyladdr_tdown1 = surfacedata;
				surfacedata += sizeof(bgq_weyl_nonvec);
				bgq_weyl_nonvec *weyladdr_tdown2= surfacedata;
				surfacedata += sizeof(bgq_weyl_nonvec);
				bgq_su3_weyl_load_combine_double(weyl_tdown, weyladdr_tdown1, weyladdr_tdown2);

				bgq_gaugesite *gaugesite_tdown = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, TUP, t1-1, t2-1, true,false);
				bgq_su3_matrix_load(gauge_tdown, gaugesite_tdown);
				//bgq_su3_matrix_flush(gaugesite_tdown);
				bgq_su3_weyl_mvinvmul(weyl_tdown, gauge_tdown, weyl_tdown);

				bgq_su3_accum_weyl_tdown(result, weyl_tdown);
			}

			bgq_su3_spinor_prefetch(addr + 2*4*3*2);
			bgq_su3_matrix_prefetch(gaugeaddr + 2*3*3*2);

			// X+ //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_spinor_decl(spinor_xup);
				bgq_su3_weyl_decl(weyl_xup);
				bgq_su3_mdecl(gauge_xup);

				bgq_spinorsite *spinorsite_xup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x+1, y, z, t1, t2, true,false);
				bgq_su3_spinor_load(spinor_xup, spinorsite_xup);
				//bgq_su3_spinor_flush(spinorsite_xup);
				//bgq_su3_spinor_invalidate(spinorsite_xup);
				bgq_su3_reduce_weyl_xup(weyl_xup, spinor_xup);

				bgq_gaugesite *gaugesite_xup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, XUP, t1, t2, true,false);
				bgq_su3_matrix_load(gauge_xup, gaugesite_xup);
				//bgq_su3_matrix_flush(gaugesite_xup);
				//bgq_su3_matrix_invalidate(gaugesite_xup);
				bgq_su3_weyl_mvmul(weyl_xup, gauge_xup, weyl_xup);

				bgq_su3_accum_weyl_xup(result, weyl_xup);
			}

			bgq_su3_spinor_prefetch(addr + 2*4*3*2);
			bgq_su3_matrix_prefetch(gaugeaddr + 2*3*3*2);

			// X- //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_spinor_decl(spinor_xdown);
				bgq_su3_weyl_decl(weyl_xdown);
				bgq_su3_mdecl(gauge_xdown);

				bgq_spinorsite *spinorsite_xdown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x-1, y, z, t1, t2, true,false);
				bgq_su3_spinor_load(spinor_xdown, spinorsite_xdown);
				//bgq_su3_spinor_flush(spinorsite_xdown);
				//bgq_su3_spinor_invalidate(spinorsite_xdown);
				bgq_su3_reduce_weyl_xdown(weyl_xdown, spinor_xdown);

				bgq_gaugesite *gaugesite_xdown = BGQ_GAUGESITE(gaugefield, isOdd, tv, x-1, y, z, XUP, t1, t2, true,false);
				bgq_su3_matrix_load(gauge_xdown, gaugesite_xdown);
				//bgq_su3_matrix_flush(gaugesite_xdown);
				//bgq_su3_matrix_invalidate(gaugesite_xdown);
				bgq_su3_weyl_mvinvmul(weyl_xdown, gauge_xdown, weyl_xdown);

				bgq_su3_accum_weyl_xdown(result, weyl_xdown);
			}

			bgq_su3_spinor_prefetch(addr + 2*4*3*2);
			bgq_su3_matrix_prefetch(gaugeaddr + 2*3*3*2);

			// Y+ //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_spinor_decl(spinor_yup);
				bgq_su3_weyl_decl(weyl_yup);
				bgq_su3_mdecl(gauge_yup);

				bgq_spinorsite *spinorsite_yup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y+1, z, t1, t2, true,false);
				bgq_su3_spinor_load(spinor_yup, spinorsite_yup);
				//bgq_su3_spinor_flush(spinorsite_yup);
				//bgq_su3_spinor_invalidate(spinorsite_yup);
				bgq_su3_reduce_weyl_yup(weyl_yup, spinor_yup);

				bgq_gaugesite *gaugesite_yup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, YUP, t1, t2, true,false);
				bgq_su3_matrix_load(gauge_yup, gaugesite_yup);
				//bgq_su3_matrix_flush(gaugesite_yup);
				//bgq_su3_matrix_invalidate(gaugesite_yup);
				bgq_su3_weyl_mvmul(weyl_yup, gauge_yup, weyl_yup);

				bgq_su3_accum_weyl_yup(result, weyl_yup);
			}

			bgq_su3_spinor_prefetch(addr + 2*4*3*2);
			bgq_su3_matrix_prefetch(gaugeaddr + 2*3*3*2);

			// Y- //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_spinor_decl(spinor_ydown);
				bgq_su3_weyl_decl(weyl_ydown);
				bgq_su3_mdecl(gauge_ydown);

				bgq_spinorsite *spinorsite_ydown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y-1, z, t1, t2, true,false);
				bgq_su3_spinor_load(spinor_ydown, spinorsite_ydown);
				//bgq_su3_spinor_flush(spinorsite_ydown);
				//bgq_su3_spinor_invalidate(spinorsite_ydown);
				bgq_su3_reduce_weyl_ydown(weyl_ydown, spinor_ydown);

				bgq_gaugesite *gaugesite_ydown = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y-1, z, YUP, t1, t2, true,false);
				bgq_su3_matrix_load(gauge_ydown, gaugesite_ydown);
				//bgq_su3_matrix_flush(gaugesite_ydown);
				//bgq_su3_matrix_invalidate(gaugesite_ydown);
				bgq_su3_weyl_mvinvmul(weyl_ydown, gauge_ydown, weyl_ydown);

				bgq_su3_accum_weyl_ydown(result, weyl_ydown);
			}

			bgq_su3_spinor_prefetch(addr + 2*4*3*2);
			bgq_su3_matrix_prefetch(gaugeaddr + 2*3*3*2);

			// Z+ //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_spinor_decl(spinor_zup);
				bgq_su3_weyl_decl(weyl_zup);
				bgq_su3_mdecl(gauge_zup);

				bgq_spinorsite *spinorsite_zup = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z+1, t1, t2, true,false);
				bgq_su3_spinor_load(spinor_zup, spinorsite_zup);
				//bgq_su3_spinor_flush(spinorsite_zup);
				//bgq_su3_spinor_invalidate(spinorsite_zup);
				bgq_su3_reduce_weyl_zup(weyl_zup, spinor_zup);

				bgq_gaugesite *gaugesite_zup = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z, ZUP, t1, t2, true,false);
				bgq_su3_matrix_load(gauge_zup, gaugesite_zup);
				//bgq_su3_matrix_flush(gaugesite_zup);
				//bgq_su3_matrix_invalidate(gaugesite_zup);
				bgq_su3_weyl_mvmul(weyl_zup, gauge_zup, weyl_zup);

				bgq_su3_accum_weyl_zup(result, weyl_zup);
			}

			bgq_su3_spinor_prefetch(addr + 2*4*3*2);
			bgq_su3_matrix_prefetch(gaugeaddr + 2*3*3*2);

			// Z- //////////////////////////////////////////////////////////////////////////
			{
				bgq_su3_spinor_decl(spinor_zdown);
				bgq_su3_weyl_decl(weyl_zdown);
				bgq_su3_mdecl(gauge_zdown);

				bgq_spinorsite *spinorsite_zdown = BGQ_SPINORSITE(spinorfield, !isOdd, tv, x, y, z-1, t1, t2, true, false);
				bgq_su3_spinor_load(spinor_zdown, spinorsite_zdown);
				//bgq_su3_spinor_flush(spinorsite_zdown);
				//bgq_su3_spinor_invalidate(spinorsite_zdown);
				bgq_su3_reduce_weyl_zdown(weyl_zdown, spinor_zdown);

				bgq_gaugesite *gaugesite_zdown = BGQ_GAUGESITE(gaugefield, isOdd, tv, x, y, z-1, ZUP, t1, t2, true, false);
				bgq_su3_matrix_load(gauge_zdown, gaugesite_zdown);
				//bgq_su3_matrix_flush(gaugesite_zdown);
				//bgq_su3_matrix_invalidate(gaugesite_zdown);
				bgq_su3_weyl_mvinvmul(weyl_zdown, gauge_zdown, weyl_zdown);

				bgq_su3_accum_weyl_zdown(result, weyl_zdown);
			}

			// Store the result
			//for (int i = 0; i < 4; i+=1) {
			bgq_spinorsite *targetsite = BGQ_SPINORSITE(targetfield, isOdd, tv, x, y, z, t1,t2, false, true);
			bgq_su3_spinor_zeroload(targetsite); /* no impact on performance */
			bgq_su3_spinor_store(targetsite, result);
			//bgq_su3_spinor_flush(targetsite); /* reduces performance */
			//}
}


void bgq_HoppingMatrix(bool isOdd, int targetfield_index, int spinorfield_index, int tid, int threads) {
	bgq_weylfield weyldata;

// 1. Distribute
	{
	size_t workload = PHYSICAL_VOLUME;
	size_t threadload = (PHYSICAL_VOLUME+threads-1)/threads;
	size_t begin = tid*threadload;
	size_t end = min(workload, (tid+1)*threadload);
	uint32_t *weylindex = bgq_distr_indices[begin*8];

	for (int xyztv = begin; xyztv < end; xyztv+=1) {
		bgq_su3_spinor_decl(spinor);
		// Get the spinor from somewhere
		bgq_HoppingMatrix_weylfield();
	}
	}



// 2. Start communication

// 3. Compute the body
	{
		size_t workload = PHYSICAL_BODY;
		size_t threadload = (workload+threads-1)/threads;
		size_t begin = tid*threadload;
		size_t end = min(workload,begin+threadload);
		COMPLEX_PRECISION *surfacedata = (char*)weyldata + bgq_weyl_section_offset(sec_surface);
		for (int i = begin; i < end; i+=1) {
			bgq_HoppingMatrix_kernel_surface();
		}
	}

// 4. Wait for communication to finish

// 5. Compute the surface

}
