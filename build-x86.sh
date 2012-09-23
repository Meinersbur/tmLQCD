#! /bin/sh
export ac_cv_func_malloc_0_nonnull=yes

ROOTPATH=`pwd`


cd $ROOTPATH/lime

echo
echo Configure Lime..................
./configure CC="mpicc" CCFLAGS="-g -O0"
## --host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile CC=mpixlc_r CCFLAGS="-g -qfullpath -I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include" build_alias=ppc64-ibm-linux host_alias=ppc-ibm-bprts 

echo
echo Make Lime..........
make




cd $ROOTPATH

echo
echo Configure tmLQCD.....

LAPACK=""
LAPACK="${LAPACK} -llapack"


#CFLAGS="-g -O0 -ffast-math -fopenmp -Wall -Wundef -DBGQ=1 -DBGQ_QPX=0 -DBGQ_FIELD_COORDCHECK=0 -I/usr/lib/gcc/x86_64-linux-gnu/4.6/include -I/usr/include/openmpi -DBGQ_PREFETCH_EXPLICIT=0 -DBGQ_PREFETCH_STREAM=0 -DBGQ_PREFETCH_LIST=0 -DMPI=1 -DBGQ_HM_CARRY=1 -DBGQ_REPLACE=0" 
CFLAGS=""
CFLAGS="${CFLAGS} -I/usr/lib/gcc/x86_64-linux-gnu/4.6/include"
CFLAGS="${CFLAGS} -I/usr/include/openmpi"
CFLAGS="${CFLAGS} -I/usr/include/openmpi"
CFLAGS="${CFLAGS} -O0"
CFLAGS="${CFLAGS} -g"
CFLAGS="${CFLAGS} -ffast-math"
CFLAGS="${CFLAGS} -fopenmp"
CFLAGS="${CFLAGS} -Wall"
CFLAGS="${CFLAGS} -Wundef"


CFLAGS="${CFLAGS} -DBGQ=1"
#CFLAGS="${CFLAGS} -DXLC=1"
#CFLAGS="${CFLAGS} -DNDEBUG=1"
#CFLAGS="${CFLAGS} -DBGQ_HM_NOKAMUL=1"
CFLAGS="${CFLAGS} -DBGQ_QPX=0"
CFLAGS="${CFLAGS} -DBGQ_PREFETCH_EXPLICIT=0"
CFLAGS="${CFLAGS} -DBGQ_PREFETCH_STREAM=0"
CFLAGS="${CFLAGS} -DBGQ_PREFETCH_LIST=0"
CFLAGS="${CFLAGS} -DBGQ_FIELD_COORDCHECK=0"
CFLAGS="${CFLAGS} -DMPI=1"
#CFLAGS="${CFLAGS} -DXLC=1"
CFLAGS="${CFLAGS} -DBGQ_HM_CARRY=1"
CFLAGS="${CFLAGS} -DBGQ_REPLACE=1"


#LDFLAGS="-lgomp -lblas"
LDFLAGS=""
LDFLAGS="${LDFLAGS} -lgomp"
LDFLAGS="${LDFLAGS} -lblas"

#-./configure --with-lapack="-llapack" --with-limedir="${ROOTPATH}/lime" --enable-mpi --with-mpidimension=3 --disable-halfspinor CC="mpicc" CFLAGS="-g -O0 -ffast-math -fopenmp -Wall -Wundef -DBGQ=1 -DBGQ_QPX=0 -DBGQ_FIELD_COORDCHECK=0 -I/usr/lib/gcc/x86_64-linux-gnu/4.6/include -I/usr/include/openmpi -DBGQ_PREFETCH_EXPLICIT=0 -DBGQ_PREFETCH_STREAM=0 -DBGQ_PREFETCH_LIST=0 -DMPI=1 -DBGQ_HM_CARRY=1 -DBGQ_REPLACE=0" OPTARGS="-O0 -g eart " SOPTARGS="-O0 -g badfg" --enable-optimize=no LDFLAGS="-lgomp -lblas"
CONFIGURE=""
#CONFIGURE="${CONFIGURE} --with-alignment=32"
CONFIGURE="${CONFIGURE} --without-bgldram"
CONFIGURE="${CONFIGURE} --with-limedir=${ROOTPATH}/lime"
CONFIGURE="${CONFIGURE} --enable-mpi"
#CONFIGURE="${CONFIGURE} --enable-qpx"
CONFIGURE="${CONFIGURE} --with-mpidimension=XYT"
#CONFIGURE="${CONFIGURE} --enable-omp"
CONFIGURE="${CONFIGURE} --enable-gaugecopy"
#CONFIGURE="${CONFIGURE} --enable-halfspinor"
CONFIGURE="${CONFIGURE} --disable-halfspinor"
CONFIGURE="${CONFIGURE} --enable-largefile"
CONFIGURE="${CONFIGURE} --with-lapack="\"'${LAPACK}'\"
CONFIGURE="${CONFIGURE} CC=mpicc"
CONFIGURE="${CONFIGURE} CFLAGS="\"'${CFLAGS}'\"
#CONFIGURE="${CONFIGURE} F77=bgf77"
CONFIGURE="${CONFIGURE} LDFLAGS="\"'${LDFLAGS}'\"
#CONFIGURE="${CONFIGURE} FC=bgxlf_r"

CONFIGURE="${CONFIGURE} --enable-optimize=no"
#CONFIGURE="${CONFIGURE} --with-lemondir=${ROOTPATH}/lemon"
CONFIGURE="${CONFIGURE} --disable-sse2"
CONFIGURE="${CONFIGURE} --disable-sse3"



echo ${CONFIGURE}
echo
eval ./decho ${CONFIGURE}
echo
eval ./configure ${CONFIGURE}


echo
echo Making tmLQCD
make -j4 bgqbench benchmark invert hmc_tm

