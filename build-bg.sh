#! /bin/sh

ROOTPATH=`pwd`

#mkdir ../build-lemon-bg
#cd ../build-lemon-bg

#echo
#echo Configure lemon..................
# ../lemon/configure --srcdir=/homea/hch02/hch023/tmLQCD-5.1.5/lemon --prefix=/homea/hch02/hch023/tmLQCD-5.1.5/lemon-build --enable-la        rgefile CC=mpixlc_r

cd $ROOTPATH/lemon

echo 
echo Configure Lemon....
./configure --enable-largefile CC=mpixlc_r

echo
echo Make Lemon......
make 

 


#mkdir ../build-lime-bg
#cd ../build-lime-bg

cd $ROOTPATH/lime

echo
echo Configure Lime..................
# ./configure --prefix=/homea/hch02/hch023/tmLQCD-5.1.1/lime --host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile CC=mpixlc_r CCFLAGS="-I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include build_alias=ppc64-ibm-linux host_alias=ppc-ibm        -bprts --no-create --no-recursion
./configure --host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile CC=mpixlc_r CCFLAGS="-I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include" build_alias=ppc64-ibm-linux host_alias=ppc-ibm-bprts 

echo
echo Make Lime..........
make



cd $ROOTPATH

echo
echo Configure tmLQCD.....

./configure --enable-mpi --with-mpidimension=4 --enable-gaugecopy --enable-halfspinor --without-gprof --without-bgldram --with-limedir=/homea/hch02/hch023/tmLQCD-5.1.1/lime --with-lemondir=/homea/hch02/hch023/tmLQCD-5.1.5/lemon-build --host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile --with-lapack="-L/bgsys/local/lapack/lib -L/opt/ibmmath/essl/4.4/lib -lesslbg -llapack -lesslbg -lxlf90_r" CC=mpixlc_r CCFLAGS="-I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include" F77=mpixlf77_r --with-lapackdir=/usr/local/bg_soft/lapack

#./configure --enable-mpi --with-mpidimension=4 --enable-gaugecopy --enable-halfspinor --without-gprof --without-bgldram -with-limedir=/homea/hch02/hch02d/tmLQCD-5.1.1/lime -with-lemondir=/homea/hch02/hch02d/tmLQCD-5.1.5/lemon-build host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile --with-lapack="-L/bgsys/local/lapack/lib -L/opt/ibmmath/essl/4.4/lib -lesslbg -llapack -lesslbg -lxlf90_r" CC=mpixlc_r CCFLAGS="-I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include F77=mpixlf77_r" | tee configure.log

echo 
echo Making tmLQCD
make hmc_tm | tee make.log


