#!/bin/bash
# Hyposat ExternalLocator wrapper — model: ak135_A
export HYPOSAT_MODEL=ak135_A
exec python3 /home/ubuntu/Projects/Hyposat/seiscomp/hyposat_wrapper.py "$@"
