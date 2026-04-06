#!/bin/bash
# Hyposat ExternalLocator wrapper — model: iasp91_A
export HYPOSAT_MODEL=iasp91_A
exec python3 /home/ubuntu/Projects/Hyposat/seiscomp/hyposat_wrapper.py \
  --hypomod-log=/home/ubuntu/Projects/Hyposat/logs/hypomod_iasp91.log \
  "$@"
