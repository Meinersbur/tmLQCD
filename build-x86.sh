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

./configure --with-lapack="-llapack" --with-limedir="${ROOTPATH}/lime" --enable-mpi --with-mpidimension=XYZ --enable-halfspinor CC="mpicc" CCFLAGS="-g -O0"
make benchmark hmc_tm invert

# mpirun -np 16 hmc_tm -v -f input.inp

