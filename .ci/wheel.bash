#!/bin/bash

set -e
set -x

BASE_DIR="$(cd "$(dirname ${BASH_SOURCE[0]})" ; cd ../ ; pwd -P )"

source ${BASE_DIR}/.ci/common_setup.bash

python -m build
