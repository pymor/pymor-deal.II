#!/bin/bash

set -e
set -x

cd /home/dealii/src
cd lib
cmake .
make

