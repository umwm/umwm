#!/usr/bin/env bash
set -euo pipefail

# netCDF-Fortran built with nvfortran, for the netCDF-C library of the system.

prefix=${1:?usage: build-nvhpc-netcdf.sh PREFIX}
jobs=$(nproc)
workdir=$(mktemp -d)
trap 'rm -rf "$workdir"' EXIT

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
  "https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v${netcdf_fortran_version}.tar.gz" \
  netcdf-fortran.tar.gz
cmake -S "netcdf-fortran-${netcdf_fortran_version}" -B build-netcdf-fortran \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_Fortran_COMPILER=nvfortran \
  -DCMAKE_INSTALL_PREFIX="$prefix" \
  -DBUILD_TESTING=OFF
cmake --build build-netcdf-fortran --parallel "$jobs"
cmake --install build-netcdf-fortran

test "$("$prefix/bin/nf-config" --has-nc4)" = yes
