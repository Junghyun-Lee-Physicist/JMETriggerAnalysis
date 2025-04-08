#!/bin/bash

echo
echo "--------------------------------------------------"
echo "Initializing VOMS proxy (valid for 72 hours)"
echo "--------------------------------------------------"

condorPath="${CMSSW_BASE}/src/JMETriggerAnalysis/JESCorrections/test/condorForJEC/"

if voms-proxy-init --voms cms --valid 72:00 --out "${condorPath}/proxy.cert"; then
    # Check if the proxy file was actually created
    if [ -f "${condorPath}/proxy.cert" ]; then
        echo "[INFO] VOMS proxy successfully initialized to ${condorPath}."
    else
        echo "[ERROR] VOMS proxy initialization reported success, but ${condorPath}/proxy.cert does not exist."
        exit 1
    fi
else
    echo "[ERROR] VOMS proxy initialization failed. Exiting."
    exit 1
fi

echo
echo "--------------------------------------------------"
echo "CRAB environment setup complete. Ready to submit jobs."
echo "--------------------------------------------------"
