#! /bin/sh

./configure --with-lapack="-llapack" --enable-mpi --with-mpidimension=4 CC="mpicc" CCFLAGS="-g -O3"
make benchmark hmc_tm

# mpirun -np 16 hmc_tm -v -f input.inp

