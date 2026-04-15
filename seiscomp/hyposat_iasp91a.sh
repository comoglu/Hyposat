#!/bin/bash
# Hyposat ExternalLocator wrapper — model: iasp91a_A
export HYPOSAT_MODEL=iasp91a_A
exec python3 /home/ubuntu/Projects/Hyposat/seiscomp/hyposat_wrapper.py \
  "$@"
