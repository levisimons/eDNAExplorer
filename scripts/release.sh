#!/bin/bash

if [ -z "$1" ]; then
  echo "No tag provided. Exiting."
  exit 1
fi

TAG=$1

docker tag edna-explorer:latest "jimjeffers/edna-explorer:${TAG}"
docker push "jimjeffers/edna-explorer:${TAG}"
