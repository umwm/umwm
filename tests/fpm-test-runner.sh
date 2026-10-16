#!/bin/sh
set -eu

test_executable=$(realpath "$1")
shift
cd "$(dirname "$0")"
python_command=${PYTHON:-python3}

if [ ! -f input/.rankine_forcing.stamp ]; then
  mkdir -p input
  "$python_command" forcing/idealized_rankine_forcing.py \
    --output-dir input \
    --start-time 2026-01-01T00:00:00 \
    --last-hour 1 \
    --nx 21 \
    --ny 11 \
    --dx-m 10000 \
    --dy-m 10000 \
    --vmax 30 \
    --rmw-m 50000 \
    --initial-center-x-m 0 \
    --initial-center-y-m 0 \
    --translation-x-m-per-hour 10000 \
    --translation-y-m-per-hour 0 \
    --uc 0.25 \
    --vc -0.5 \
    --rhoa 1.3 \
    --rhow 1025
  touch input/.rankine_forcing.stamp
fi

exec "$test_executable" "$@"
