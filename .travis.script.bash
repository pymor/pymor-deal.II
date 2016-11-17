#!/bin/bash

set -e
set -x

cd /home/dealii/src
mkdir build && cd build
cmake .
make

