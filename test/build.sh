#!/bin/bash

cp ../environment.yml .
cp ../bin/get_genome.py .
docker build . -t "rnafusion-test"
