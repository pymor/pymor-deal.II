#!/bin/bash

set -e
set -x

BASE_DIR=/home/pymor/src
export PYTHONPATH=${BASE_DIR}/lib:${BASE_DIR}/src
# Using
#   export PYTHONPATH=${BUILD_DIR}/lib:${BASE_DIR}/src:${PYTHONPATH}
# will fail with empty PYTHONPATH as the empty string
# after ':' will be treated as the current directory
# in which the operator module is imported instead of
# the operator module from stdlib.
cd ${BASE_DIR}
python setup.py build_ext -i
cd ${BASE_DIR}/pymor_src
pip install -e .
pip install -r requirements-travis.txt
pip install -r requirements-optional.txt

cd ${BASE_DIR}/src/
xvfb-run -a pytest test/demo.py
