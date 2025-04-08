
#!/bin/bash

CFG_FILE="pfHadCalibNTuple_cfg.py"
TMP_DUMP="testCfg.py"

echo "Step 1: Checking Global Tag in configuration"
echo "-----------------------------------------------"
echo "Target config file: $CFG_FILE"
echo "Dumping the expanded configuration using python3 $CFG_FILE dumpPython=$TMP_DUMP ..."
echo

python3 "$CFG_FILE" dumpPython="$TMP_DUMP"

if [[ ! -f "$TMP_DUMP" ]]; then
    echo "Failed to generate $TMP_DUMP. Please check if $CFG_FILE runs correctly."
    exit 1
fi

echo
echo "Step 2: Searching for Global Tag in dumped config"
echo "----------------------------------------------------"
GTTAG=$(grep -i globaltag "$TMP_DUMP")

if [[ -z "$GTTAG" ]]; then
    echo "No Global Tag found in $TMP_DUMP. Check if it is set programmatically."
else
    echo -e "Found Global Tag:\n"
    echo "$GTTAG"
fi

echo
echo "Step 3: Cleaning up temporary file"
echo "-------------------------------------"
rm -f "$TMP_DUMP"
echo "Removed $TMP_DUMP"

echo
echo "Done!"
