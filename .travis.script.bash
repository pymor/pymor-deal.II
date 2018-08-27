#!/bin/bash

set -e
set -x

BASE_DIR=/home/dealii/src
BUILD_DIR=/tmp/build
export PYTHONPATH=${PYTHONPATH}:${BUILD_DIR}/lib:${BASE_DIR}/src

mkdir ${BUILD_DIR} && cd ${BUILD_DIR}
cmake ${BASE_DIR} 
make

cd ${BASE_DIR}/src/pydealii/pymor
python demo.py
