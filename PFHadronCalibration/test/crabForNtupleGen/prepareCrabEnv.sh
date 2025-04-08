#!/bin/bash

echo "--------------------------------------------------"
echo "Step 1: Setting up the CRAB3 environment"
echo "--------------------------------------------------"

if source /cvmfs/cms.cern.ch/crab3/crab.sh; then
    echo "[INFO] CRAB3 environment variables successfully loaded."
else
    echo "[ERROR] Failed to source CRAB3 environment. Exiting."
    exit 1
fi

echo
echo "--------------------------------------------------"
echo "Step 2: Initializing VOMS proxy (valid for 72 hours)"
echo "--------------------------------------------------"

if voms-proxy-init --voms cms --valid 72:00; then
    echo "[INFO] VOMS proxy successfully initialized."
else
    echo "[ERROR] VOMS proxy initialization failed. Exiting."
    exit 1
fi

echo
echo "--------------------------------------------------"
echo "CRAB environment setup complete. Ready to submit jobs."
echo "--------------------------------------------------"
