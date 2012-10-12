#! /bin/sh

ROOTPATH=`pwd`





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
CFLAGS="${CFLAGS} -I${HOME}/usr/include"
CFLAGS="${CFLAGS} -O5"
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
CFLAGS="${CFLAGS} -g"
#CFLAGS="${CFLAGS} -qsimd=auto"
#CFLAGS="${CFLAGS} -qsmp=noauto"

#CFLAGS="${CFLAGS} -DXLC=1"
#CFLAGS="${CFLAGS} -DNDEBUG=1"
#CFLAGS="${CFLAGS} -DBGQ=1"
#CFLAGS="${CFLAGS} -DBGQ_QPX=1"
#CFLAGS="${CFLAGS} -DPAPI=1"


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
LDFLAGS="${LDFLAGS} -L$HOME/usr/lib -lpapi"

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
#CONFIGURE="${CONFIGURE} --with-fixedvolume"
CONFIGURE="${CONFIGURE} --without-bgldram"
CONFIGURE="${CONFIGURE} --with-limedir=${HOME}/lime"
CONFIGURE="${CONFIGURE} --enable-mpi"
#CONFIGURE="${CONFIGURE} --enable-qpx"
CONFIGURE="${CONFIGURE} --with-mpidimension=XYZT"
#CONFIGURE="${CONFIGURE} --enable-omp"
CONFIGURE="${CONFIGURE} --enable-gaugecopy"
CONFIGURE="${CONFIGURE} --enable-halfspinor"
#CONFIGURE="${CONFIGURE} --disable-halfspinor"
CONFIGURE="${CONFIGURE} --enable-largefile"
CONFIGURE="${CONFIGURE} --with-lapack="\"'${LAPACK}'\"
CONFIGURE="${CONFIGURE} CC=/bgsys/drivers/ppcfloor/comm/xl/bin/mpixlc_r"
#CONFIGURE="${CONFIGURE} CC=/bgsys/drivers/ppcfloor/comm/xl.ndebug/bin/mpixlc_r"
#CONFIGURE="${CONFIGURE} CCDEP=mpixlc_r"
CONFIGURE="${CONFIGURE} CFLAGS="\"'${CFLAGS}'\"
CONFIGURE="${CONFIGURE} F77=bgf77"
CONFIGURE="${CONFIGURE} LDFLAGS="\"'${LDFLAGS}'\"
CONFIGURE="${CONFIGURE} FC=bgxlf_r"
#CONFIGURE="${CONFIGURE} --with-lemondir=${ROOTPATH}/lemon"

CONFIGURE="${CONFIGURE} --enable-optimize=no"
CONFIGURE="${CONFIGURE} --enable-qpx"
CONFIGURE="${CONFIGURE} --enable-spi"



eval ./decho ${CONFIGURE}
echo
eval ./configure ${CONFIGURE}


echo 
echo Making tmLQCD
make -j benchmark invert hmc_tm


