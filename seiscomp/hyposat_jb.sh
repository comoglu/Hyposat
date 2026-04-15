#!/bin/bash
# Hyposat ExternalLocator wrapper — model: jb_A
export HYPOSAT_MODEL=jb_A
exec python3 /home/ubuntu/Projects/Hyposat/seiscomp/hyposat_wrapper.py \
  "$@"
