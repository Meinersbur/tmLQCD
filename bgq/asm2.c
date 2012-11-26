

#include "bgq_qpx.h"
#include "bgq_field.h"

typedef struct {
	bgq_vector4double_decl(ks); // sum
	bgq_vector4double_decl(kc); // compensation or carried error (i.e. ks+kc yields better approximation than just ks)
} bgq_vectorkahan_t;



static inline bgq_vectorkahan_t bgq_reduce_initkahan() {
	bgq_vectorkahan_t result;
	bgq_zero(result.ks);
	bgq_zero(result.kc);
	return result;
}


#define bgq_kahan_add(ks, kc, val) bgq_kahan_add_raw(bgq_vars(ks), bgq_vars(kc), bgq_vars(val))
static inline void bgq_kahan_add_raw(bgq_params(*accumulator_ks), bgq_params(*accumulator_kc), bgq_params(val)) {
	bgq_vector4double_decl(ks);
	bgq_vector4double_decl(kc);
	bgq_mov(ks, *accumulator_ks);
	bgq_mov(kc, *accumulator_kc);

	bgq_vector4double_decl(tr); // val+compensation
	bgq_vector4double_decl(ts); // next sum
	bgq_vector4double_decl(tt); // val+next compensation

	bgq_add(tr, val, kc); // val + carried_error
	bgq_add(ts, tr, ks); // val + carried_error + sum - error
	bgq_sub(tt, ts, ks); // val + carried_error - error
	bgq_mov(ks, ts);
	bgq_sub(kc, tr, tt); // error

    bgq_mov(*accumulator_ks, ks);
    bgq_mov(*accumulator_kc, kc);
}

void bgq_reduce_norm(bgq_vectorkahan_t *accumulator, bgq_su3_spinor_params(spinor), ucoord ic) {
	bgq_vector4double_decl(siteresult);
	bgq_mul(siteresult, spinor_v0_c0, spinor_v0_c0);
	bgq_madd(siteresult, spinor_v0_c1, spinor_v0_c1, siteresult);
	bgq_madd(siteresult, spinor_v0_c2, spinor_v0_c2, siteresult);
	bgq_madd(siteresult, spinor_v1_c0, spinor_v1_c0, siteresult);
	bgq_madd(siteresult, spinor_v1_c1, spinor_v1_c1, siteresult);
	bgq_madd(siteresult, spinor_v1_c2, spinor_v1_c2, siteresult);
	bgq_madd(siteresult, spinor_v2_c0, spinor_v2_c0, siteresult);
	bgq_madd(siteresult, spinor_v2_c1, spinor_v2_c1, siteresult);
	bgq_madd(siteresult, spinor_v2_c2, spinor_v2_c2, siteresult);
	bgq_madd(siteresult, spinor_v3_c0, spinor_v3_c0, siteresult);
	bgq_madd(siteresult, spinor_v3_c1, spinor_v3_c1, siteresult);
	bgq_madd(siteresult, spinor_v3_c2, spinor_v3_c2, siteresult);

	bgq_kahan_add(&accumulator->ks, &accumulator->kc, siteresult);
}
