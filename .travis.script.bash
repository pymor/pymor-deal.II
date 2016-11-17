#!/bin/bash

set -e
set -x

BASE_DIR=/home/dealii/src
export PYTHONPATH=${PYTHONPATH}:${BASE_DIR}/lib:${BASE_DIR}/src

cd ${BASE_DIR}
mkdir build && cd build
cmake ..
make

cd ${BASE_DIR}/src/pydealii/pymor
python demo.py
