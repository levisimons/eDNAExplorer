#!/bin/bash
source /opt/conda/etc/profile.d/conda.sh
conda activate GBIF_env
echo "Executing: $@"
exec "$@"