#!/bin/bash

set -e
set -x

BASE_DIR=/src
BUILD_DIR=/tmp/build
export PYTHONPATH=${BUILD_DIR}/lib
# Using
#   export PYTHONPATH=${BUILD_DIR}/lib:${PYTHONPATH}
# will fail with empty PYTHONPATH as the empty string
# after ':' will be treated as the current directory
# in which the operator module is imported instead of
# the operator module from stdlib.

mkdir ${BUILD_DIR} && cd ${BUILD_DIR}
cmake ${BASE_DIR}
make

cd ${BASE_DIR}/pymor_src
pip install -e.

cd ${BASE_DIR}/src/pydealii/pymor
python demo.py
