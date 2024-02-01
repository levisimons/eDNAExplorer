#!/bin/bash

# Check if two arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: ./report.sh SCRIPT_NAME REPORT_ID"
    exit 1
fi

# Extract arguments
SCRIPT_NAME=$1
REPORT_ID=$2

# Execute Rscript with the arguments
conda run -n reports Rscript --vanilla $SCRIPT_NAME $REPORT_ID