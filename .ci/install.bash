#!/bin/bash

set -e
set -x

BASE_DIR="$(cd "$(dirname ${BASH_SOURCE[0]})" ; cd ../ ; pwd -P )"

source ${BASE_DIR}/.ci/common_setup.bash

python -m pip install .

python -m pip uninstall -y pymor_dealii

python -m pip install -e .

