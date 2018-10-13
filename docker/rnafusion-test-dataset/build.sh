#!/bin/bash

docker build . -t "rnafusion-test-dataset"
docker tag rnafusion-test-dataset nfcore/rnafusion-test-dataset
