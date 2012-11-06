#! /bin/sh
ROOTPATH=`pwd`
cd $ROOTPATH

echo
echo Configure tmLQCD.....


CPPFLAGS=""
CPPFLAGS="${CPPFLAGS} -I/bgsys/drivers/ppcfloor/arch/include"
CPPFLAGS="${CPPFLAGS} -I/bgsys/drivers/ppcfloor/comm/xl/include"
CPPFLAGS="${CPPFLAGS} -I${HOME}/usr/include"
CPPFLAGS="${CPPFLAGS} -DNDEBUG=1"
CPPFLAGS="${CPPFLAGS} -DBGQ_QPX=1"
CPPFLAGS="${CPPFLAGS} -DXLC=1"
#CPPFLAGS="${CPPFLAGS} -DBGQ=1"
CPPFLAGS="${CPPFLAGS} -DPAPI=1"


CFLAGS=""
CFLAGS="${CFLAGS} -g"
CFLAGS="${CFLAGS} -O5"
CFLAGS="${CFLAGS} -qlibmpi"
#CFLAGS="${CFLAGS} -qprefetch=aggressive"
CFLAGS="${CFLAGS} -qarch=qp"
CFLAGS="${CFLAGS} -qtune=qp"
CFLAGS="${CFLAGS} -qmaxmem=-1"
CFLAGS="${CFLAGS} -qtm" # Enable transactional memory
#CFLAGS="${CFLAGS} -qsimd=noauto"
CFLAGS="${CFLAGS} -qsmp=noauto"
#CFLAGS="${CFLAGS} -qstrict=all"
#CFLAGS="${CFLAGS} -qstrict=none"
CFLAGS="${CFLAGS} -qstrict=order"
#CFLAGS="${CFLAGS} -qipa=level=2"
#CFLAGS="${CFLAGS} -qsimd=auto"
#CFLAGS="${CFLAGS} -qsmp=noauto"


LIBS=""
LIBS="${LIBS} -lxl"
LIBS="${LIBS} -lxlopt"
LIBS="${LIBS} -lxlf90_r"
LIBS="${LIBS} -lxlfmath"
LIBS="${LIBS} -lxlsmp"
LIBS="${LIBS} -lpthread"
LIBS="${LIBS} -lbgq"
LIBS="${LIBS} -lSPI_l1p"
LIBS="${LIBS} -lSPI"
LIBS="${LIBS} -lbgpm"

LDFLAGS=""
LDFLAGS="${LDFLAGS} -qipa=level=2"
LDFLAGS="${LDFLAGS} -L/opt/ibmcmp/xlf/bg/14.1/lib64"
LDFLAGS="${LDFLAGS} -L/opt/ibmcmp/xlsmp/bg/3.1/bglib64"
LDFLAGS="${LDFLAGS} -L/usr/local/bg_soft/lapack/3.3.0"
LDFLAGS="${LDFLAGS} -L/bgsys/ibm_essl/prod/opt/ibmmath/lib64"
LDFLAGS="${LDFLAGS} -L/bgsys/drivers/ppcfloor/spi/lib"
LDFLAGS="${LDFLAGS} -L/bgsys/drivers/ppcfloor/bgpm/lib/"
#LDFLAGS="${LDFLAGS} -qipa=level=2"
#LDFLAGS="${LDFLAGS} -L/opt/ibmcmp/xlsmp/bg/3.1/bglib64"
#LDFLAGS="${LDFLAGS} -L/bgsys/ibm_essl/prod/opt/ibmmath/lib64"
#LDFLAGS="${LDFLAGS} -L/usr/local/bg_soft/lapack/3.3.0/lib"
#LDFLAGS="${LDFLAGS} -L/bgsys/local/lib/"
#LDFLAGS="${LDFLAGS} -L$HOME/usr/lib -lpapi"


module load lapack
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


CONFIGURE=""
#CONFIGURE="${CONFIGURE} --with-alignment=32"
#CONFIGURE="${CONFIGURE} --with-fixedvolume"
CONFIGURE="${CONFIGURE} --without-bgldram"
CONFIGURE="${CONFIGURE} --with-limedir=${HOME}/lime"
#CONFIGURE="${CONFIGURE} --with-lemondir=${ROOTPATH}/lemon"
CONFIGURE="${CONFIGURE} --enable-mpi"
CONFIGURE="${CONFIGURE} --with-mpidimension=XYZT"
#CONFIGURE="${CONFIGURE} --enable-omp"
CONFIGURE="${CONFIGURE} --enable-gaugecopy"
CONFIGURE="${CONFIGURE} --enable-halfspinor"
#CONFIGURE="${CONFIGURE} --disable-halfspinor"
CONFIGURE="${CONFIGURE} --enable-largefile"
CONFIGURE="${CONFIGURE} --with-lapack="\"'${LAPACK}'\"
CONFIGURE="${CONFIGURE} CC=/bgsys/drivers/ppcfloor/comm/xl/bin/mpixlc_r"
#CONFIGURE="${CONFIGURE} CC=/bgsys/drivers/ppcfloor/comm/xl.ndebug/bin/mpixlc_r"
CONFIGURE="${CONFIGURE} F77=bgf77"
CONFIGURE="${CONFIGURE} FC=bgxlf_r"
CONFIGURE="${CONFIGURE} CPPFLAGS="\"'${CPPFLAGS}'\"
CONFIGURE="${CONFIGURE} CFLAGS="\"'${CFLAGS}'\"
CONFIGURE="${CONFIGURE} LIBS="\"'${LIBS}'\"
CONFIGURE="${CONFIGURE} LDFLAGS="\"'${LDFLAGS}'\"
CONFIGURE="${CONFIGURE} --enable-optimize=no"
CONFIGURE="${CONFIGURE} --enable-qpx"
CONFIGURE="${CONFIGURE} --enable-spi"
CONFIGURE="${CONFIGURE} --disable-sse2"
CONFIGURE="${CONFIGURE} --disable-sse3"


eval ./decho ${CONFIGURE}
echo
eval ./configure ${CONFIGURE}


echo 
echo Making tmLQCD
make -j32 bgqbench benchmark invert hmc_tm


