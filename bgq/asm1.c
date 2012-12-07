
//#include "bgq_qpx.h"
//#include "bgq_field.h"

vector4double x;
vector4double y;
vector4double z;
vector4double w;

void test() {
	vector4double a = vec_mul(x, y);
	vector4double b = vec_add(a, z);
	w = b;
}
