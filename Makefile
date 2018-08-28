.PHONY: dockertest

DEAL=9.0.0
IMAGE="pymor/dealii-testing_pymor_0.4:dealii_$(DEAL)"

dockertest:
	make -C docker $(DEAL)
	docker run -v $(shell pwd):/home/pymor/src $(IMAGE) /home/pymor/src/.ci/travis.script.bash

dockerrun:
	make -C docker $(DEAL)
	docker run -it -v $(shell pwd):/home/pymor/src $(IMAGE) bash
