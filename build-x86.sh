#! /bin/sh
ROOTPATH=`pwd`
cd $ROOTPATH

echo
echo Configure tmLQCD.....


CPPFLAGS=""
CPPFLAGS="${CPPFLAGS} -I/usr/include/openmpi"
#CPPFLAGS="${CPPFLAGS} -DNDEBUG=1"
CPPFLAGS="${CPPFLAGS} -DBGQ_QPX=0"
#CPPFLAGS="${CPPFLAGS} -DBGQ=1"
#CPPFLAGS="${CPPFLAGS} -DXLC=1"
#CPPFLAGS="${CPPFLAGS} -DBGQ_HM_NOKAMUL=1"
CPPFLAGS="${CPPFLAGS} -DBGQ_PREFETCH_EXPLICIT=0"
CPPFLAGS="${CPPFLAGS} -DBGQ_PREFETCH_STREAM=0"
CPPFLAGS="${CPPFLAGS} -DBGQ_PREFETCH_LIST=0"
#CPPFLAGS="${CPPFLAGS} -DBGQ_COORDCHECK=1"
CPPFLAGS="${CPPFLAGS} -DMPI=1"
#CPPFLAGS="${CPPFLAGS} -DXLC=1"
#CPPFLAGS="${CPPFLAGS} -DBGQ_HM_CARRY=1"
CPPFLAGS="${CPPFLAGS} -DBGQ_REPLACE=0"
CPPFLAGS="${CPPFLAGS} -DPAPI=0"
CPPFLAGS="${CPPFLAGS} -DBGQ_UNVECTORIZE=1"


CFLAGS=""
CFLAGS="${CFLAGS} -g"
CFLAGS="${CFLAGS} -O0"
CFLAGS="${CFLAGS} -ffast-math"
CFLAGS="${CFLAGS} -fopenmp"
CFLAGS="${CFLAGS} -Wall"
CFLAGS="${CFLAGS} -Wundef"


LIBS=""
LIBS="${LIBS} -lgomp"
LIBS="${LIBS} -lblas"
#LIBS="${LIBS} -lefence"


LDFLAGS=""


#module load lapack
LAPACK=""
LAPACK="${LAPACK} -llapack"


CONFIGURE=""
#CONFIGURE="${CONFIGURE} --with-alignment=32"
#CONFIGURE="${CONFIGURE} --with-fixedvolume"
CONFIGURE="${CONFIGURE} --without-bgldram"
CONFIGURE="${CONFIGURE} --with-limedir=${ROOTPATH}/lime"
#CONFIGURE="${CONFIGURE} --with-lemondir=${ROOTPATH}/lemon"
CONFIGURE="${CONFIGURE} --enable-mpi"
CONFIGURE="${CONFIGURE} --with-mpidimension=XYZT"
#CONFIGURE="${CONFIGURE} --enable-omp"
CONFIGURE="${CONFIGURE} --enable-gaugecopy"
#CONFIGURE="${CONFIGURE} --enable-halfspinor"
CONFIGURE="${CONFIGURE} --disable-halfspinor"
CONFIGURE="${CONFIGURE} --enable-largefile"
CONFIGURE="${CONFIGURE} --with-lapack="\"'${LAPACK}'\"
CONFIGURE="${CONFIGURE} CC=mpicc"
#CONFIGURE="${CONFIGURE} F77=bgf77"
#CONFIGURE="${CONFIGURE} FC=bgxlf_r
CONFIGURE="${CONFIGURE} CPPFLAGS="\"'${CPPFLAGS}'\"
CONFIGURE="${CONFIGURE} CFLAGS="\"'${CFLAGS}'\"
CONFIGURE="${CONFIGURE} LIBS="\"'${LIBS}'\"
CONFIGURE="${CONFIGURE} LDFLAGS="\"'${LDFLAGS}'\"
CONFIGURE="${CONFIGURE} --enable-optimize=no"
CONFIGURE="${CONFIGURE} --disable-qpx"
CONFIGURE="${CONFIGURE} --disable-spi"
CONFIGURE="${CONFIGURE} --disable-sse2"
CONFIGURE="${CONFIGURE} --disable-sse3"



eval ./decho ${CONFIGURE}
echo
eval ./configure ${CONFIGURE}


echo
echo Making tmLQCD
make -j4 bgqbench benchmark invert hmc_tm

