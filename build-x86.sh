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

./configure --with-lapack="-llapack" --with-limedir="${ROOTPATH}/lime" --enable-mpi --with-mpidimension=3 --disable-halfspinor CC="mpicc" CFLAGS="-g -O0 -ffast-math -fopenmp -Wall -Wundef -DBGQ=1 -DBGQ_QPX=0 -DBGQ_FIELD_COORDCHECK=0 -I/usr/lib/gcc/x86_64-linux-gnu/4.6/include -I/usr/include/openmpi -DBGQ_PREFETCH_EXPLICIT=0 -DBGQ_PREFETCH_STREAM=0 -DBGQ_PREFETCH_LIST=0 -DMPI=1 -DBGQ_HM_NOKAMUL=0" OPTARGS="-O0 -g eart " SOPTARGS="-O0 -g badfg" --enable-optimize=no LDFLAGS="-lgomp -lblas"
make -j4 bgqbench benchmark hmc_tm invert

# mpirun -np 16 hmc_tm -v -f input.inp

