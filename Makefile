#!/usr/bin/env make

# customization points via makefile key-value arguments
#
# interpreter in images: 3.{6,7,8} currently available
# DOCKER_BASE_PYTHON=3.7
#

THIS_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
DOCKER_RUN=docker run -v $(THIS_DIR):/pymor --env-file  $(THIS_DIR)/.env
DOCKER_COMPOSE=CI_COMMIT_SHA=$(shell git log -1 --pretty=format:"%H") \
	NB_USER=$(NB_USER) $(COMPOSE_SUDO) docker-compose -f .binder/docker-compose.yml -p pymor
NB_DIR=notebooks
NB_USER:=${USER}
ifeq ($(PYMOR_SUDO), 1)
	COMPOSE_SUDO:=sudo -E
else
	COMPOSE_SUDO:=
endif

# this loads $(ENV_FILE) as both makefile variables and into shell env
ENV_FILE?=.env
include $(ENV_FILE)
export $(shell sed 's/=.*//' $(ENV_FILE))

.PHONY: docker_file

docker_file:
	 sed -e "s;CI_IMAGE_TAG;$(CI_IMAGE_TAG);g" -e "s;DOCKER_BASE_PYTHON;$(DOCKER_BASE_PYTHON);g" \
		 .binder/Dockerfile.in > .binder/Dockerfile

docker_image: docker_file
	$(DOCKER_COMPOSE) build

docker_docs: docker_image
	NB_DIR=notebooks $(DOCKER_COMPOSE) run docs ./.ci/gitlab/test_docs.bash

docker_run: docker_image
	$(DOCKER_COMPOSE) run --service-ports jupyter bash

docker_exec: docker_image
	$(DOCKER_COMPOSE) run --service-ports jupyter bash -l -c "${DOCKER_CMD}"

docker_tutorials: NB_DIR=docs/_build/html
docker_tutorials: docker_docs docker_jupyter

docker_test: docker_image
	$(DOCKER_COMPOSE) run --service-ports jupyter /src/.ci/travis.script.bash

docker_jupyter: docker_image
	NB_DIR=$(NB_DIR) $(DOCKER_COMPOSE) up jupyter
docker_wheel_check: docker_image
	PYMOR_TEST_OS=$(PYMOR_TEST_OS) $(DOCKER_COMPOSE) run --service-ports wheel_check bash

docs:
	PYTHONPATH=${PWD}/src/:${PYTHONPATH} make -C docs html
