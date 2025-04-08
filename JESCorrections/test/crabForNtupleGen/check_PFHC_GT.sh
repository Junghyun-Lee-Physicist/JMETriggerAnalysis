#!/bin/bash

DB_FILE="PFCalibration.db"
TAG_KEYWORD="PFCalibration"

cd ..

echo "------------------------------------------------------------"
echo "Step 1: Checking if database file exists"
echo "------------------------------------------------------------"

if [[ ! -f "$DB_FILE" ]]; then
    echo "[ERROR] Database file '$DB_FILE' not found in current directory."
    echo "[INFO] Please ensure the file exists before running this script."
    exit 1
fi

echo "[INFO] Found local database file: $DB_FILE"
echo

echo "------------------------------------------------------------"
echo "Step 2: Searching for tags related to '$TAG_KEYWORD'"
echo "------------------------------------------------------------"

conddb --db "$DB_FILE" search "$TAG_KEYWORD"

STATUS=$?

echo
if [[ $STATUS -eq 0 ]]; then
    echo "[INFO] Search completed successfully."
else
    echo "[ERROR] conddb command failed with exit code $STATUS."
    exit $STATUS
fi

echo
echo "------------------------------------------------------------"
echo "Done."
echo "------------------------------------------------------------"

cd -
