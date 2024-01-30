#!/bin/bash

if [ -z "$1" ]; then
  echo "No command provided. Exiting."
  exit 1
fi

COMMAND="$@"
SCRIPT_DIR=$(dirname "$0")
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)

# DOCKER_CMD="docker run --platform linux/amd64 -v $PROJECT_ROOT:/project jimjeffers/edna-explorer:v3 $COMMAND"
DOCKER_CMD="docker run -v $PROJECT_ROOT:/project reports-dev:tag $COMMAND"

echo "Executing: $DOCKER_CMD"
eval $DOCKER_CMD
