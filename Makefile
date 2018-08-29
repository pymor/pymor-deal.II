.PHONY: dockertest dockerrun

DEAL=9.0.0
PY=3.5
IMAGE="pymor/dealii:v$(DEAL)_py$(PY)"

dockertest:
	docker run -v $(shell pwd):/src $(IMAGE) /src/.ci/travis.script.bash

dockerrun:
	docker run -it -v $(shell pwd):/src $(IMAGE) bash
