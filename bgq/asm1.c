

#include <bgq_utils.h>

typedef void (*calledfunc)(int);

inline void callme1(int arg) {
	printf("Callme1 %d\n", arg);
}

inline void callme2(int arg) {
	printf("Callme2 %d\n", arg);
}




inline void caller(calledfunc func, int arg) {
	(*func)(arg);
}


void inst1() {
	caller(&callme1, 1);
}

void inst2() {
	caller(&callme2, 2);
}

int main(int argc, char **argv) {
	inst1();
	return 0;
}
