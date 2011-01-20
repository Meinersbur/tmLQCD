#! /bin/sh

ROOTPATH=`pwd`

#echo
#echo Configure lemon..................

cd $ROOTPATH/lemon

echo 
echo Configure Lemon....
#./configure --enable-largefile CC=mpixlc_r CCFLAGS="-g -qfullpath"

echo
echo Make Lemon......
make 

 



cd $ROOTPATH/lime

echo
echo Configure Lime..................
#./configure --host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile CC=mpixlc_r CCFLAGS="-g -qfullpath -I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include" build_alias=ppc64-ibm-linux host_alias=ppc-ibm-bprts 

echo
echo Make Lime..........
make



cd $ROOTPATH

echo
echo Configure tmLQCD.....

# Heavy optimization
# -O5 -qnoautoconfig -qprefetch -qarch=450d -qtune=450 -qcache=level=1:type=i:size=32:line=32:assoc=64:cost=8 -qcache=level=1:type=d:size=32:line=32:assoc=64:cost=8 -qcache=level=2:type=c:size=4096:line=128:assoc=8:cost=40

# Normal
./configure --enable-mpi --with-mpidimension=4 --enable-gaugecopy --enable-halfspinor --without-gprof --without-bgldram --with-limedir=/homea/hch02/hch023/tmLQCD-5.1.1/lime --with-lemondir=/homea/hch02/hch023/tmLQCD-5.1.5/lemon-build --host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile --with-lapack="-L/bgsys/local/lapack/lib -L/opt/ibmmath/essl/4.4/lib -lesslbg -llapack -lesslbg -lxlf90_r" CC=mpixlc_r CCFLAGS="-g -qfullpath -I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include" F77=mpixlf77_r --with-lapackdir=/usr/local/bg_soft/lapack | tee configure.log

# With scalasca
#module load scalasca
#./configure --enable-mpi --with-mpidimension=4 --enable-gaugecopy --enable-halfspinor --without-gprof --without-bgldram --with-limedir=/homea/hch02/hch023/tmLQCD-5.1.1/lime --with-lemondir=/homea/hch02/hch023/tmLQCD-5.1.5/lemon-build --host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile --with-lapack="-L/bgsys/local/lapack/lib -L/opt/ibmmath/essl/4.4/lib -lesslbg -llapack -lesslbg -lxlf90_r" CC="skin mpixlc_r" CCFLAGS="-g -qfullpath -I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include" F77=mpixlf77_r --with-lapackdir=/usr/local/bg_soft/lapack | tee configure.log

# With mpitrace
#./configure --enable-mpi --with-mpidimension=4 --enable-gaugecopy --enable-halfspinor --without-gprof --without-bgldram --with-limedir=/homea/hch02/hch023/tmLQCD-5.1.1/lime --with-lemondir=/homea/hch02/hch023/tmLQCD-5.1.5/lemon-build --host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile --with-lapack="-L/bgsys/local/lapack/lib -L/opt/ibmmath/essl/4.4/lib -lesslbg -llapack -lesslbg -lxlf90_r" CC=mpixlc_r CCFLAGS="-g -qfullpath -lmpitrace -llicense -I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include" F77=mpixlf77_r --with-lapackdir=/usr/local/bg_soft/lapack | tee configure.log
# LIBS = -lmpitrace -llicense -lhmc -lsolver -llinalg -lhmc -lio -L/bgsys/local/lapack/lib -L/opt/ibmmath/essl/4.4/lib -lesslbg -llapack -lesslbg -lxlf90_r -llemon -llime   -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/comm/default/lib -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/comm/sys/lib -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/runtime/SPI -L/opt/ibmcmp/xlsmp/bg/1.7/bglib -L/opt/ibmcmp/xlmass/bg/4.4/bglib -L/opt/ibmcmp/xlf/bg/11.1/bglib -R/opt/ibmcmp/lib/bg/bglib -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/gnu-linux/lib/gcc/powerpc-bgp-linux/4.1.2 -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/gnu-linux/lib/gcc/powerpc-bgp-linux/4.1.2/../../../../powerpc-bgp-linux/lib -lmpich.cnk -lopa -ldcmf.cnk -ldcmfcoll.cnk -lpthread -lSPI.cna -lrt -lxlf90_r -lxlopt -lxlomp_ser -lxl -lxlfmath -ldl -lm -lgcc_eh -lm -L/usr/local/bg_soft/ihpct-2.2.2/lib -lgfortran -L/usr/local/bg_soft/toolchain/gnu-linux/powerpc-bgp-linux/lib


echo 
echo Making tmLQCD
make hmc_tm | tee make.log


