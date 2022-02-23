#!/bin/bash

set -e
set -x

BASE_DIR="$(cd "$(dirname ${BASH_SOURCE[0]})" ; cd ../ ; pwd -P )"

source ${BASE_DIR}/.ci/common_setup.bash

# if we're in a versioned branch pip will downgrade pymor here from pypi
python -m pip install .

xvfb-run -a pytest -s -r sxX  test/*py
