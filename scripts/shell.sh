#bin/bash
docker run --platform linux/amd64 --entrypoint /bin/bash -v $PWD:/project -it jimjeffers/edna-explorer:v2