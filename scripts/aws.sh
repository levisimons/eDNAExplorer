aws ecr get-login-password --region us-west-1 --profile edna-explorer | docker login --username AWS --password-stdin 576661951794.dkr.ecr.us-west-1.amazonaws.com
docker pull 576661951794.dkr.ecr.us-west-1.amazonaws.com/reports:latest