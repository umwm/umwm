#!/bin/bash
#
# Removes the output and restart files that the model writes to an output
# directory, the current directory by default, and nothing else: the copies of
# the namelists, umwmout_main.nml and the others, are kept.
#
#   ./clean.sh [OUTPUT]

dir=${1:-.}
rm -vf "$dir"/umwmout_*.nc "$dir"/umwmout.grid "$dir"/umwmspc_*.nc "$dir"/umwmrst_*.nc
