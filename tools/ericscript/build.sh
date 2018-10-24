#!/bin/bash

tool=$(basename $(pwd))
version=$(grep $(head -1 environment.yml | cut -d' ' -f2)'=' environment.yml | cut -d'=' -f2)
container_name=nfcore/rnafusion:${tool}_v${version}
docker build . -t $tool
docker tag $tool $container_name
docker push $container_name