#!/bin/bash
source /usr/local/miniconda/etc/profile.d/conda.sh
conda activate reports
echo "Executing: $@"
exec "$@"
