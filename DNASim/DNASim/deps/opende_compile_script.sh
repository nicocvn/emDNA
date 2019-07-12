#!/bin/sh

# compilation script for ODE static library

# script argument
LIB_MODE=${1}

# current directory
CURRENT_DIR=`pwd`

# opende directories
ODE_SRC_DIR="${CURRENT_DIR}/ode-0.13/"
ODE_INSTALL_DIR="${CURRENT_DIR}/ode_${LIB_MODE}"

# check for opende directory
if [ ! -d "${ODE_SRC_DIR}" ]; then
	echo "ODE sources directory not found, cannot build ODE library";
	echo "DNASim compilation aborted";
	exit 1;
fi

# check if new build is needed
if [ -f "${ODE_INSTALL_DIR}/lib/libode.a" ]; then
    exit 0;
fi

# opende install directory
rm -Rf ${ODE_INSTALL_DIR}
mkdir ${ODE_INSTALL_DIR}

# compiler flags
export ODE_CXXFLAGS="-O2 -ffast-math -funroll-loops -msse3 -mssse3 \
						-fstrict-aliasing \
						-m64 -w"
export ODE_CFLAGS="-O2 -ffast-math -funroll-loops -msse3 -mssse3 \
						-fstrict-aliasing \
						-m64 -w"
if [ "${LIB_MODE}" = "shared" ];
then
    export ODE_CXXFLAGS="${ODE_CXXFLAGS} -fPIC"
    export ODE_CFLAGS="${ODE_CFLAGS} -fPIC"
fi
if [ `uname` = "Darwin" ];
then
    export ODE_CXX=/usr/bin/c++
    export ODE_CXXFLAGS="${ODE_CXXFLAGS} -std=gnu++11 -stdlib=libc++"
fi

# switch to src directory
cd ${ODE_SRC_DIR}

# configure stage
export CXXFLAGS="${ODE_CXXFLAGS}"
export CFLAGS="${ODE_CFLAGS}"
export CXX="${ODE_CXX}"
sh configure --prefix="${ODE_INSTALL_DIR}" \
    --disable-demos --without-x \
	--with-trimesh=none \
    --with-drawstuff=none \
    --enable-double-precision \
	--enable-libccd \
    --disable-threading-intf \
    --enable-builtin-threading-imp \
    2>&1 > ../opende_configure_log

# make directives
/usr/bin/make all 2>&1 > ../opende_build_log
/usr/bin/make install >> ../opende_build_log 2>&1
cd ..

# unset
unset CXXFLAGS
unset CFLAGS
unset CXX
