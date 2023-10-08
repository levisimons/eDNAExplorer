#!/bin/bash

# Load .env file
if [ -f ".env" ]; then
  source .env
else
  echo ".env file not found"
  exit 1
fi

# Check if projectID is provided
if [ -z "$1" ]; then
  echo "Usage: ./stage.sh <projectID>"
  exit 1
fi

# Set projectID from the first argument
projectID=$1

# Set custom region and endpoint
BUCKET_REGION="RegionOne"
BUCKET_ENDPOINT="https://js2.jetstream-cloud.org:8001"

# List contents of the directory before copying
# echo "Key: ${AWS_ACCESS_KEY_ID}"
# echo "Copying contents of s3://ednaexplorer/tronko_output/${projectID}/"
aws s3 ls "s3://ednaexplorer/" --endpoint-url ${BUCKET_ENDPOINT}

# Copy recursively using AWS CLI
aws s3 cp "s3://ednaexplorer/tronko_output/${projectID}/" "s3://ednaexplorer_staging/tronko_output/${projectID}/" --recursive \
  --endpoint-url ${BUCKET_ENDPOINT}
aws s3 cp "s3://ednaexplorer/projects/${projectID}/" "s3://ednaexplorer_staging/projects/${projectID}/" --recursive \
  --endpoint-url ${BUCKET_ENDPOINT}