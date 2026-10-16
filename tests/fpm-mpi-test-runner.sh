#!/bin/sh
set -eu

script_dir=$(CDPATH= cd -- "$(dirname "$0")" && pwd)
exec mpiexec -n "${FPM_TEST_PROCESSES:-1}" "$script_dir/fpm-test-runner.sh" "$@"
