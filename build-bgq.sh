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
make 
mkdir lib
cp src/liblemon.a lib

 



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



module load lapack

#--with-lapack="-L/bgsys/local/lib/ -L/usr/local/bg_soft/lapack/3.3.0/lib -lesslbg -llapack -lesslbg -lxlf90_r -L/opt/ibmcmp/xlf/bg/14.1/lib64 -lxl -lxlopt -lxlf90_r -lxlfmath 
# -L/opt/ibmcmp/xlsmp/bg/3.1/bglib64 -lxlsmp -lpthread"
LAPACK=""
LAPACK="${LAPACK} -L/bgsys/local/lib/"
LAPACK="${LAPACK} -L/usr/local/bg_soft/lapack/3.3.0/lib"
LAPACK="${LAPACK} -lesslbg"
LAPACK="${LAPACK} -llapack"
LAPACK="${LAPACK} -lesslbg"
LAPACK="${LAPACK} -lxlf90_r"
LAPACK="${LAPACK} -L/opt/ibmcmp/xlf/bg/14.1/lib64"
LAPACK="${LAPACK} -lxl"
LAPACK="${LAPACK} -lxlopt"
LAPACK="${LAPACK} -lxlf90_r"
LAPACK="${LAPACK} -lxlfmath"
LAPACK="${LAPACK} -L/opt/ibmcmp/xlsmp/bg/3.1/bglib64"
LAPACK="${LAPACK} -lxlsmp"
LAPACK="${LAPACK} -lpthread"

# CFLAGS="-I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/xl/include -O5 -qprefetch=aggressive -qarch=qp -qtune=qp -qmaxmem=-1 -qsimd=noauto -qsmp=noauto -qstrict=all -DBGQ" 
CFLAGS=""
CFLAGS="${CFLAGS} -I/bgsys/drivers/ppcfloor/arch/include"
CFLAGS="${CFLAGS} -I/bgsys/drivers/ppcfloor/comm/xl/include"
CFLAGS="${CFLAGS} -O3"
CFLAGS="${CFLAGS} -qprefetch=aggressive"
CFLAGS="${CFLAGS} -qarch=qp"
CFLAGS="${CFLAGS} -qtune=qp"
CFLAGS="${CFLAGS} -qmaxmem=-1"
#CFLAGS="${CFLAGS} -qsimd=noauto"
CFLAGS="${CFLAGS} -qsmp=noauto"
#CFLAGS="${CFLAGS} -qstrict=all"
CFLAGS="${CFLAGS} -DBGQ=1"
CFLAGS="${CFLAGS} -DXLC=1"
#CFLAGS="${CFLAGS} -qipa=level=2"

CFLAGS="${CFLAGS} -g"
CFLAGS="${CFLAGS} -DNDEBUG=1"
#CFLAGS="${CFLAGS} -DBGQ_HM_NOKAMUL=1"
#CFLAGS="${CFLAGS} -qsimd=auto"
CFLAGS="${CFLAGS} -qstrict=none"
CFLAGS="${CFLAGS} -DBGQ_QPX=1"
CFLAGS="${CFLAGS} -DBGQ_PREFETCH_EXPLICIT=0"
CFLAGS="${CFLAGS} -DBGQ_PREFETCH_STREAM=0"
CFLAGS="${CFLAGS} -DBGQ_PREFETCH_LIST=0"
CFLAGS="${CFLAGS} -DBGQ_FIELD_COORDCHECK=0"
CFLAGS="${CFLAGS} -DMPI=1"
CFLAGS="${CFLAGS} -DBGQ_HM_CARRY=0"
CFLAGS="${CFLAGS} -DBGQ_REPLACE=0"
#CFLAGS="${CFLAGS} -qsmp=noauto"

# LDFLAGS="-L/opt/ibmcmp/xlf/bg/14.1/lib64 -L/usr/local/bg_soft/lapack/3.3.0 -lxl -lxlopt -lxlf90_r -L/bgsys/drivers/ppcfloor/bgpm/lib/ -lxlfmath -L/opt/ibmcmp/xlsmp/bg/3.1/bglib64 
# -lxlsmp -lpthread -L/bgsys/ibm_essl/prod/opt/ibmmath/lib64"
LDFLAGS=""
LDFLAGS="${LDFLAGS} -L/opt/ibmcmp/xlf/bg/14.1/lib64"
LDFLAGS="${LDFLAGS} -L/usr/local/bg_soft/lapack/3.3.0"
LDFLAGS="${LDFLAGS} -lxl"
LDFLAGS="${LDFLAGS} -lxlopt"
LDFLAGS="${LDFLAGS} -lxlf90_r"
LDFLAGS="${LDFLAGS} -L/bgsys/drivers/ppcfloor/bgpm/lib/"
LDFLAGS="${LDFLAGS} -lxlfmath"
LDFLAGS="${LDFLAGS} -L/opt/ibmcmp/xlsmp/bg/3.1/bglib64"
LDFLAGS="${LDFLAGS} -lxlsmp"
LDFLAGS="${LDFLAGS} -lpthread"
LDFLAGS="${LDFLAGS} -L/bgsys/ibm_essl/prod/opt/ibmmath/lib64"

LDFLAGS="${LDFLAGS} -lSPI_l1p"
LDFLAGS="${LDFLAGS} -L/bgsys/drivers/ppcfloor/spi/lib"
#LDFLAGS="${LDFLAGS} -qipa=level=2"
#LDFLAGS="${LDFLAGS} -L/opt/ibmcmp/xlsmp/bg/3.1/bglib64"
#LDFLAGS="${LDFLAGS} -L/bgsys/ibm_essl/prod/opt/ibmmath/lib64"
#LDFLAGS="${LDFLAGS} -L/usr/local/bg_soft/lapack/3.3.0/lib"
#LDFLAGS="${LDFLAGS} -L/bgsys/local/lib/"

#--with-alignment=32 --without-bgldram --with-limedir=/work/pra067/pra06700/juqueen/programs/lime_c --enable-mpi --enable-qpx --with-mpidimension=4 --enable-omp --enable-gaugecopy
# --disable-halfspinor --enable-largefile 
#--with-lapack="-L/bgsys/local/lib/ -L/usr/local/bg_soft/lapack/3.3.0/lib -lesslbg -llapack -lesslbg -lxlf90_r -L/opt/ibmcmp/xlf/bg/14.1/lib64 -lxl -lxlopt -lxlf90_r -lxlfmath -L/opt/ibmcmp/xlsmp/bg/3.1/bglib64 -lxlsmp -lpthread" 
# CC=/bgsys/drivers/ppcfloor/comm/xl/bin/mpixlc_r 
# CFLAGS="-I/bgsys/drivers/ppcfloor/arch/include/ -I/bgsys/drivers/ppcfloor/comm/xl/include -O5 -qprefetch=aggressive -qarch=qp -qtune=qp -qmaxmem=-1 -qsimd=noauto -qsmp=noauto -qstrict=all -DBGQ" 
#n F77=bgf77 
# LDFLAGS="-L/opt/ibmcmp/xlf/bg/14.1/lib64 -L/usr/local/bg_soft/lapack/3.3.0 -lxl -lxlopt -lxlf90_r -L/bgsys/drivers/ppcfloor/bgpm/lib/ -lxlfmath -L/opt/ibmcmp/xlsmp/bg/3.1/bglib64 -lxlsmp -lpthread -L/bgsys/ibm_essl/prod/opt/ibmmath/lib64"
# FC=bgxlf_r
CONFIGURE=""
#CONFIGURE="${CONFIGURE} --with-alignment=32"
CONFIGURE="${CONFIGURE} --without-bgldram"
CONFIGURE="${CONFIGURE} --with-limedir=${ROOTPATH}/lime"
CONFIGURE="${CONFIGURE} --enable-mpi"
#CONFIGURE="${CONFIGURE} --enable-qpx"
CONFIGURE="${CONFIGURE} --with-mpidimension=XYT"
#CONFIGURE="${CONFIGURE} --enable-omp"
CONFIGURE="${CONFIGURE} --enable-gaugecopy"
#CONFIGURE="${CONFIGURE} --enable-halfspinor"
CONFIGURE="${CONFIGURE} --disable-halfspinor"
CONFIGURE="${CONFIGURE} --enable-largefile"
CONFIGURE="${CONFIGURE} --with-lapack="\"'${LAPACK}'\"
CONFIGURE="${CONFIGURE} CC=/bgsys/drivers/ppcfloor/comm/xl/bin/mpixlc_r"
CONFIGURE="${CONFIGURE} CFLAGS="\"'${CFLAGS}'\"
CONFIGURE="${CONFIGURE} F77=bgf77"
CONFIGURE="${CONFIGURE} LDFLAGS="\"'${LDFLAGS}'\"
CONFIGURE="${CONFIGURE} FC=bgxlf_r"

CONFIGURE="${CONFIGURE} --enable-optimize=no"
#CONFIGURE="${CONFIGURE} --with-lemondir=${ROOTPATH}/lemon"



echo ${CONFIGURE}
echo
eval ./decho ${CONFIGURE}
echo
eval ./configure ${CONFIGURE}


echo 
echo Making tmLQCD
make -j32 bgqbench benchmark invert hmc_tm


