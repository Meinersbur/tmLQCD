#! /bin/sh
export ac_cv_func_malloc_0_nonnull=yes

./configure --with-lapack="-llapack" --enable-mpi --with-mpidimension=4 CFLAGS="-g -O3 -I/usr/include/mpich2" INCLUDES="-I/usr/include/mpich2"  LIBS="-L/usr/lib64/ -lopa -lmpl -lrt -lpthread -lblas -lmpich"
# CFLAGS=" -I/usr/include/mpich2"  
# LIBS="-L/usr/lib -lopa -lmpl  -lrt -lcr -lpthread" INCLUDES="-I/usr/include/mpich2"  
make benchmark hmc_tm

# mpirun -np 16 hmc_tm -v -f input.inp

