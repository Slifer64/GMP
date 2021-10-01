#! /bin/sh

# colors for output
COLOR_RED="\033[1;31m"
COLOR_GREEN="\033[1;32m"
COLOR_RESET="\033[0m"

success=0

MAIN_PATH=${PWD}

# path to the osqp_lib ros package
OSQP_LIB_PATH=${MAIN_PATH}"/src/lib/osqp_lib"

cd ${MAIN_PATH}"/install/" && \

# download and build osqp
git clone --recursive https://github.com/osqp/osqp && \

cd osqp/ && \

mkdir osqp_lib && \

mkdir build && \
cd build && \

cmake .. -DCMAKE_INSTALL_PREFIX:PATH=../osqp_lib && \
make install && \

cd ../osqp_lib/ && \

# move folders and libs to osqp_lib ros package
mv include/osqp include/osqp_lib && \
rsync -a include/ ${OSQP_LIB_PATH}"/include/" && \
mv lib/libosqp.so ${OSQP_LIB_PATH}"/lib/libosqp.so" && \
mv lib/libqdldl.so ${OSQP_LIB_PATH}"/lib/libqdldl.so" && \

# clean up
cd ${MAIN_PATH}/install/ && \
rm -rf osqp/ && \
success=1

# output message
if [ $success -eq 1 ]; then
  echo $COLOR_GREEN"OSQP installed successfully!"$COLOR_RESET	
else
  echo $COLOR_RED"Failed to complete OSQP installation..."$COLOR_RESET
fi

