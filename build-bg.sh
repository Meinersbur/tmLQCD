#! /bin/sh

ROOTPATH=`pwd`

#echo
#echo Configure lemon..................

cd $ROOTPATH/lemon

echo 
echo Configure Lemon....
#./configure --enable-largefile CC=mpixlc_r CCFLAGS="-g -qfullpath" --prefix=`pwd`

echo
echo Make Lemon......
#make 
#mkdir lib
#cp src/liblemon.a lib

 



cd $ROOTPATH/lime

echo
echo Configure Lime..................
#./configure --host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile CC=mpixlc_r CCFLAGS="-g -qfullpath -I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include" build_alias=ppc64-ibm-linux host_alias=ppc-ibm-bprts 

echo
echo Make Lime..........
#make



cd $ROOTPATH

echo
echo Configure tmLQCD.....

module load lapack/3.3.0
module load scalasca

BGP_FLOOR="/bgsys/drivers/ppcfloor"
BGP_IDIRS="-I${BGP_FLOOR}/arch/include -I${BGP_FLOOR}/comm/include"
BGP_LIBS="-L${BGP_FLOOR}/comm/lib -L${BGP_FLOOR}/runtime/SPI -lmpich.cnk -ldcmfcoll.cnk -ldcmf.cnk -lrt -lSPI.cna -lpthread"

CC="mpixlc_r"
CCLD="${CC}"
F77="bgf77"
CFLAGS="-I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include -DPREFETCH"
INCLUDES="-I/homea/hch02/hch02d/papi/include/"

LDFLAGS="-qsmp=omp"
LIBS="-L/homea/hch02/hch02d/papi/lib -lpapi ${BGP_LIBS}"

# Enable OpenMP
#LDFLAGS="${LDFLAGS} -qsmp"
#CFLAGS="${CFLAGS} -DOMP"
#CFLAGS="${CFLAGS} -qsmp=noauto"
#CFLAGS="${CFLAGS} -qsmp=noauto -DOMP"
#SOPTARGS="${SOPTARGS} -qsmp=auto"

# How many dimensions?
DIM=4

# Reports on optimizations
CFLAGS="${CFLAGS} -qreport -qxflag=diagnostic -qphsinfo"
#CFLAGS="${CFLAGS} -qreport -qxflag=diagnostic -qphsinfo -v -V"

# More detailed error messages
CFLAGS="${CFLAGS} -qsrcmsg"

# Fast compilation for non-critical code
OPTARGS="${OPTARGS} -qnosmp -O2"

# Extreme optimization for field operations (Not all switches supported by all versions of xlc)
SOPTARGS="${SOPTARGS} -O5 -qnoautoconfig -O5 -qmaxmem=-1 -qprefetch -qarch=450d -qtune=450 -qcache=level=1:type=i:size=32:line=32:assoc=64:cost=8 -qcache=level=1:type=d:size=32:line=32:assoc=64:cost=8 -qcache=level=2:type=c:size=4096:line=128:assoc=8:cost=40 -qalign=natural -qhot=simd -qhot=vector"

# Half- or full spinor exchange?
SPINOR=--enable-halfspinor
#SPINOR=--disable-halfspinor

#Scalasca
module load scalasca 
CC="skin -v ${CC}"
CCLD="skin -v ${CCLD}"
CFLAGS="${CFLAGS} -DNOPAPI"

# Lime or lemon?
LIME="--with-limedir=`pwd`/lime"
#LEMON="--with-lemondir=`pwd`/lemon"

#./configure --enable-mpi --with-mpidimension=${DIM} --enable-gaugecopy ${SPINOR} --without-gprof --without-bgldram --host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile --with-lapack="-L$LAPACK_LIB -L/opt/ibmmath/essl/4.4/lib -lesslbg -llapack -lesslbg -lxlf90_r -lxlfmath" ${LIME} ${LEMON} CFLAGS="${CFLAGS}" INCLUDES="${INCLUDES}" LIBS="${LIBS}" LDFLAGS="${LDFLAGS}" F77="${F77}" CC="${CC}" CCLD="${CCLD}"




# LEMON
#./configure --enable-mpi --with-mpidimension=4 --enable-gaugecopy --enable-halfspinor --without-gprof --without-bgldram --host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile --with-lapack="-L$LAPACK_LIB -L/opt/ibmmath/essl/4.4/lib -lesslbg -llapack -lesslbg -lxlf90_r" --with-limedir=`pwd`/lime --with-lemondir=`pwd`/lemon CC="mpixlc_r" CFLAGS="-I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include -DPREFETCH" INCLUDES="-I/homea/hch02/hch02d/papi/include/" LIBS="-L/homea/hch02/hch02d/papi/lib -lpapi" F77="bgf77"

# LIME
#./configure --enable-mpi --with-mpidimension=4 --enable-gaugecopy --enable-halfspinor --without-gprof --without-bgldram --host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile --with-lapack="-L$LAPACK_LIB -L/opt/ibmmath/essl/4.4/lib -lesslbg -llapack -lesslbg -lxlf90_r" --with-limedir=`pwd`/lime CC="mpixlc_r" CCFLAGS="-I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include -DPREFETCH" INCLUDES="-I/homea/hch02/hch02d/papi/include/" LIBS="-L/homea/hch02/hch02d/papi/lib -lpapi" F77="bgf77"

# PAPI 4.1.2
#./configure --enable-mpi --with-mpidimension=4 --enable-gaugecopy --enable-halfspinor --without-gprof --without-bgldram --host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile --with-lapack="-L$LAPACK_LIB -L/opt/ibmmath/essl/4.4/lib -lesslbg -llapack -lesslbg -lxlf90_r" --with-limedir=`pwd`/lime CC="mpixlc_r" CCFLAGS="-I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include" INCLUDES="-I/homea/hch02/hch02d/papi/include/" LIBS="-L/homea/hch02/hch02d/papi/lib -lpapi" F77="bgf77"

# Works
#./configure --enable-mpi --with-mpidimension=4 --enable-gaugecopy --enable-halfspinor --without-gprof --without-bgldram --host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile --with-lapack="-L$LAPACK_LIB -L/opt/ibmmath/essl/4.4/lib -lesslbg -llapack -lesslbg -lxlf90_r" --with-limedir=`pwd`/lime CC="mpixlc_r" CCFLAGS="-I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include -I/bgsys/local/papi/papi-c-3.9.0/include" INCLUDES="-I/bgsys/local/papi/papi-c-3.9.0/include" LIBS="-L/bgsys/local/papi/papi-c-3.9.0/lib -lpapi" F77="bgf77"




# Heavy optimization
# -O5 -qnoautoconfig -qprefetch -qarch=450d -qtune=450 -qcache=level=1:type=i:size=32:line=32:assoc=64:cost=8 -qcache=level=1:type=d:size=32:line=32:assoc=64:cost=8 -qcache=level=2:type=c:size=4096:line=128:assoc=8:cost=40

# Optimization report
# -qreport -qxflag=diagnostic

# OpenMP
# CFLAGS="-qsmp=omp"
# LDFLAGS="-qsmp=omp"

# Scalasca
# CC="skin mpixlc_r"


# With scalasca
#module load scalasca
#./configure --enable-mpi --with-mpidimension=4 --enable-gaugecopy --enable-halfspinor --without-gprof --without-bgldram --with-limedir=/homea/hch02/hch023/tmLQCD-5.1.1/lime --with-lemondir=/homea/hch02/hch023/tmLQCD-5.1.5/lemon-build --host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile --with-lapack="-L/bgsys/local/lapack/lib -L/opt/ibmmath/essl/4.4/lib -lesslbg -llapack -lesslbg -lxlf90_r" CC="skin -pomp mpixlc_r" CCFLAGS="-g -qfullpath -I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include" F77=mpixlf77_r --with-lapackdir=/usr/local/bg_soft/lapack | tee configure.log

# With mpitrace
#./configure --enable-mpi --with-mpidimension=4 --enable-gaugecopy --enable-halfspinor --without-gprof --without-bgldram --with-limedir=/homea/hch02/hch023/tmLQCD-5.1.1/lime --with-lemondir=/homea/hch02/hch023/tmLQCD-5.1.5/lemon-build --host=ppc-ibm-bprts --build=ppc64-ibm-linux --enable-largefile --with-lapack="-L/bgsys/local/lapack/lib -L/opt/ibmmath/essl/4.4/lib -lesslbg -llapack -lesslbg -lxlf90_r" CC=mpixlc_r CCFLAGS="-g -qfullpath -lmpitrace -llicense -I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/include" F77=mpixlf77_r --with-lapackdir=/usr/local/bg_soft/lapack | tee configure.log
# LIBS = -lmpitrace -llicense -lhmc -lsolver -llinalg -lhmc -lio -L/bgsys/local/lapack/lib -L/opt/ibmmath/essl/4.4/lib -lesslbg -llapack -lesslbg -lxlf90_r -llemon -llime   -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/comm/default/lib -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/comm/sys/lib -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/runtime/SPI -L/opt/ibmcmp/xlsmp/bg/1.7/bglib -L/opt/ibmcmp/xlmass/bg/4.4/bglib -L/opt/ibmcmp/xlf/bg/11.1/bglib -R/opt/ibmcmp/lib/bg/bglib -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/gnu-linux/lib/gcc/powerpc-bgp-linux/4.1.2 -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/gnu-linux/lib/gcc/powerpc-bgp-linux/4.1.2/../../../../powerpc-bgp-linux/lib -lmpich.cnk -lopa -ldcmf.cnk -ldcmfcoll.cnk -lpthread -lSPI.cna -lrt -lxlf90_r -lxlopt -lxlomp_ser -lxl -lxlfmath -ldl -lm -lgcc_eh -lm -L/usr/local/bg_soft/ihpct-2.2.2/lib -lgfortran -L/usr/local/bg_soft/toolchain/gnu-linux/powerpc-bgp-linux/lib


echo 
echo Making tmLQCD
echo \#! /bin/sh > make.sh
echo make benchmark hmc_tm OPTARGS=\"${OPTARGS}\" SOPTARGS=\"${SOPTARGS}\" >> make.sh
make benchmark hmc_tm OPTARGS="${OPTARGS}" SOPTARGS="${SOPTARGS}" | tee make.log


