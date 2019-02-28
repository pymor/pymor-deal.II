#!/bin/bash

set -e
set -x

BASE_DIR=/home/pymor/src

# Using
#   export PYTHONPATH=${BUILD_DIR}/lib:${BASE_DIR}/src:${PYTHONPATH}
# will fail with empty PYTHONPATH as the empty string
# after ':' will be treated as the current directory
# in which the operator module is imported instead of
# the operator module from stdlib.
cd ${BASE_DIR}

pip install pytest

pip install git+https://github.com/pymor/pymor.git[full]#egg=pymor
# if we're in a versioned branch pip will downgrade pymor here from pypi
pip install .

cd ${BASE_DIR}/src/
xvfb-run -a pytest -s -r sxX  test/demo.py
