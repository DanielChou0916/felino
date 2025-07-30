#!/bin/bash

# --------------------------------------------
# Felino installation test script
# Usage: ./felino_test.sh
# --------------------------------------------

# Find executable
EXE=$(find . -maxdepth 1 -type f -name "felino-opt" -perm -111 | head -n 1)

if [ -z "$EXE" ]; then
    echo "[Felino] ❌ Error: felino-opt executable not found in the project root."
    echo "Please build the project first: make -j$(nproc)"
    exit 1
fi

# Use the fixed test input path
TEST_INPUT="check/test/installation_check_disp0.i"

if [ ! -f "$TEST_INPUT" ]; then
    echo "[Felino] ❌ Error: test input file not found at $TEST_INPUT"
    exit 1
fi

echo "[Felino] (='-'=) Running installation test..."
"$EXE" -i "$TEST_INPUT" --n-threads=1 > run.log 2>&1

if [ $? -eq 0 ]; then
    echo "[Felino] (='v'=) Installation test PASSED."
    rm -f run.log   # remove log file if test passed
    exit 0
else
    echo "[Felino] (='^'=) Installation test FAILED.  Check run.log for details"
    exit 1
fi
