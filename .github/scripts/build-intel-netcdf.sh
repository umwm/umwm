#!/usr/bin/env bash
set -euo pipefail

prefix=${1:?usage: build-intel-netcdf.sh PREFIX}
jobs=$(nproc)
workdir=$(mktemp -d)
trap 'rm -rf "$workdir"' EXIT

zlib_version=1.3.2
hdf5_version=1.14.6
netcdf_c_version=4.10.1
netcdf_fortran_version=4.6.4

download_extract() {
  local url=$1
  local archive=$2
  curl --fail --location --retry 3 "$url" --output "$archive"
  tar -xf "$archive"
}

mkdir -p "$prefix"
cd "$workdir"

download_extract \
  "https://zlib.net/fossils/zlib-${zlib_version}.tar.gz" \
  zlib.tar.gz
cmake -S "zlib-${zlib_version}" -B build-zlib \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=icx \
  -DCMAKE_INSTALL_PREFIX="$prefix"
cmake --build build-zlib --parallel "$jobs"
cmake --install build-zlib

download_extract \
  "https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.6/hdf5-1.14.6.tar.gz" \
  hdf5.tar.gz
cmake -S "hdf5-${hdf5_version}" -B build-hdf5 \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=icx \
  -DCMAKE_INSTALL_PREFIX="$prefix" \
  -DBUILD_TESTING=OFF \
  -DHDF5_BUILD_HL_LIB=ON \
  -DHDF5_BUILD_TOOLS=ON \
  -DHDF5_ENABLE_Z_LIB_SUPPORT=ON \
  -DZLIB_ROOT="$prefix"
cmake --build build-hdf5 --parallel "$jobs"
cmake --install build-hdf5

download_extract \
  "https://github.com/Unidata/netcdf-c/archive/refs/tags/v${netcdf_c_version}.tar.gz" \
  netcdf-c.tar.gz
cmake -S "netcdf-c-${netcdf_c_version}" -B build-netcdf-c \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=icx \
  -DCMAKE_INSTALL_PREFIX="$prefix" \
  -DCMAKE_PREFIX_PATH="$prefix" \
  -DBUILD_TESTING=OFF \
  -DNETCDF_BUILD_UTILITIES=ON \
  -DNETCDF_ENABLE_DAP=OFF \
  -DNETCDF_ENABLE_HDF5=ON \
  -DNETCDF_ENABLE_NCZARR=OFF
cmake --build build-netcdf-c --parallel "$jobs"
cmake --install build-netcdf-c

download_extract \
  "https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v${netcdf_fortran_version}.tar.gz" \
  netcdf-fortran.tar.gz
cmake -S "netcdf-fortran-${netcdf_fortran_version}" -B build-netcdf-fortran \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=icx \
  -DCMAKE_Fortran_COMPILER=ifx \
  -DCMAKE_INSTALL_PREFIX="$prefix" \
  -DCMAKE_PREFIX_PATH="$prefix" \
  -DBUILD_TESTING=OFF
cmake --build build-netcdf-fortran --parallel "$jobs"
cmake --install build-netcdf-fortran

test "$("$prefix/bin/nc-config" --has-nc4)" = yes
test "$("$prefix/bin/nf-config" --has-nc4)" = yes
