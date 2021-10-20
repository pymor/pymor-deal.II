#!/bin/bash

set -e
set -x

BASE_DIR="$(cd "$(dirname ${BASH_SOURCE[0]})" ; cd ../ ; pwd -P )"

# Using
#   export PYTHONPATH=${BUILD_DIR}/lib:${BASE_DIR}/src:${PYTHONPATH}
# will fail with empty PYTHONPATH as the empty string
# after ':' will be treated as the current directory
# in which the operator module is imported instead of
# the operator module from stdlib.
cd ${BASE_DIR}
export CCACHE_DIR=${BASE_DIR}/cache
git submodule update --init

python -m venv ~/venv
source ~/venv/bin/activate
python -m pip install -U pip pytest wheel

python -m pip install git+https://github.com/pymor/pymor.git#egg=pymor
