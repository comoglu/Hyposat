#!/bin/bash
# Hyposat ExternalLocator wrapper — model: sp6_A
export HYPOSAT_MODEL=sp6_A
exec python3 /home/ubuntu/Projects/Hyposat/seiscomp/hyposat_wrapper.py \
  --hypomod-log=/home/ubuntu/Projects/Hyposat/logs/hypomod_sp6.log \
  "$@"
