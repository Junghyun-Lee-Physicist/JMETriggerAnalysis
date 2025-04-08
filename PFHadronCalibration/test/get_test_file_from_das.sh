
#!/bin/bash

# Configuration
CLIENT="dasgoclient"
DATASET="/Pi_Par-E-0p2to200_PGun/Run3Winter25Digi-NoPU_142X_mcRun3_2025_realistic_v7-v1/GEN-SIM-RAW"
# GT : 142X_mcRun3_2025_realistic_v7

echo "------------------------------------------------------------"
echo "Step 1: DAS Query to extract a sample input file"
echo "------------------------------------------------------------"
echo "[INFO] Target dataset:"
echo "       $DATASET"
echo
echo "[INFO] Querying DAS for available files..."

FILE=$($CLIENT -query="file dataset=$DATASET" | head -n 1)

if [[ -z "$FILE" ]]; then
    echo "[ERROR] No file retrieved. Please check the dataset name or DAS availability."
    exit 1
fi

echo
echo "------------------------------------------------------------"
echo "Step 2: Retrieved file"
echo "------------------------------------------------------------"
echo "$FILE"
echo
echo "[INFO] You can now use this file as input for local cmsRun testing."

echo "------------------------------------------------------------"
echo "Done."
echo "------------------------------------------------------------"
