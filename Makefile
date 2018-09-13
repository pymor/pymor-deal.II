.PHONY: dockertest dockerrun

DEAL=9.0.0
PY=3.5
IMAGE="pymor/dealii:v$(DEAL)_py$(PY)"

dockertest:
	docker run -v $(shell pwd):/home/pymor/src $(IMAGE) /home/pymor/src/.ci/travis.script.bash

dockerrun:
	docker run -it -v $(shell pwd):/home/pymor/src $(IMAGE) bash

alldockertest:
	for deal in "9.0.0" "8.5.1" ; do \
	for py in 5 6 7 ; do \
	    docker run -v $(shell pwd):/home/pymor/src pymor/dealii:v$${deal}_py3.$${py} \
		/home/pymor/src/.ci/travis.script.bash ; \
	done \
	done
