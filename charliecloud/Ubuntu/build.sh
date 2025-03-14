#!/bin/bash

export CH_IMAGE_STORAGE=`pwd`/images
ch-image build -t ubuntu_base -f DockerfileBase .
ch-image build --rebuild  -t ubuntu_sodasignals -f Dockerfile .
