#!/bin/bash

set -e
set -x

BASE_DIR=${HOME}/src
BUILD_DIR=${BASE_DIR}/build
export PYTHONPATH=${PYTHONPATH}:${BUILD_DIR}/lib:${BASE_DIR}/src

mkdir ${BUILD_DIR} && cd ${BUILD_DIR}
cmake ${BASE_DIR}
make

cd ${BASE_DIR}/src/pydealii/pymor
python demo.py
