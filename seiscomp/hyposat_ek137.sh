#!/bin/bash
# Hyposat ExternalLocator wrapper — model: ek137_A
export HYPOSAT_MODEL=ek137_A
exec python3 /home/ubuntu/Projects/Hyposat/seiscomp/hyposat_wrapper.py \
  --hypomod-log=/home/ubuntu/Projects/Hyposat/logs/hypomod_ek137.log \
  "$@"
