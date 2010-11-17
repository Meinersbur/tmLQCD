#! /bin/sh


echo Making Lime..................
mkdir ../build-lime-bg
cd ../build-lime-bg

# ./configure --prefix=/homea/hch02/hch02d/tmLQCD/lime --host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile CC=mpixlc_r CCFLAGS="-I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include build_alias=ppc64-ibm-linux host_alias=ppc-ibm        -bprts --no-create --no-recursion
../tmLQCD/lime/configure --prefix=/homea/hch02/hch023/tmLQCD-5.1.1/lime --host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile CC=mpixlc_r CCFLAGS="-I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include build_alias=ppc64-ibm-linux host_alias=ppc-ibm-bprts --no-create --no-recursion"
make



#echo Making lemon..................
#mkdir ../build-lemon-bg
#cd ../build-lemon-bg

#../tmLQCD/lemon/configure --srcdir=/homea/hch02/hch02d/tmLQCD/lemon --prefix=/homea/hch02/hch02d/lemon-build-bg --enable-largefile CC=mpixlc_r CCFLAGS="-I/homea/hch02/hch02d/tmLQCD/lemon/include -I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include"
#make 



#echo Making hmc_tm
#mkdir ../build-bg
#cd ../build-bg

#../tmLQCD/configure --enable-mpi --with-mpidimension=4 --enable-gaugecopy --enable-halfspinor --without-gprof --without-bgldram  -with-limedir=/homea/hch02/hch02d/build-lime-bg --with-lemondir=/homea/hch02/hch02d/tmLQCD-5.1.5/lemon-build --host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile --with-lapack="-L/bgsys/local/lapack/lib -L/opt/ibmmath/essl/4.4/lib -lesslbg -llapack -lesslbg -lxlf90_r" CC=mpixlc_r CCFLAGS="-I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include F77=mpixlf77_r" | tee configure.log

#make hmc_tm | tee make.log


