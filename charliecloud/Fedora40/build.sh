#!/bin/bash

export CH_IMAGE_STORAGE=`pwd`/images
ch-image build -t fedora40_sodasignals -f Dockerfile .
