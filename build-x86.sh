#! /bin/sh
export ac_cv_func_malloc_0_nonnull=yes

./configure --with-lapack="-llapack" --enable-mpi --with-mpidimension=2 CC="mpicc" CCFLAGS="-g -O3"
make benchmark hmc_tm

# mpirun -np 16 hmc_tm -v -f input.inp

