#!/bin/bash
source /usr/local/miniconda/etc/profile.d/conda.sh
conda activate GBIF_env
echo "Executing: $@"
exec "$@"
