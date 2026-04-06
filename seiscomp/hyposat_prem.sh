#!/bin/bash
# Hyposat ExternalLocator wrapper — model: prem_A
export HYPOSAT_MODEL=prem_A
exec python3 /home/ubuntu/Projects/Hyposat/seiscomp/hyposat_wrapper.py "$@"
