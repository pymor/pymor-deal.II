DEALII_VERSIONS = 9.0.0 8.5.1

.PHONY: dealiis $(DEALII_VERSIONS) push

dealiis: $(DEALII_VERSIONS)


$(DEALII_VERSIONS):
	cd "demo" && \
	docker build --build-arg DEALII_VERSION=$@ -t "pymor/dealii-testing_pymor_master:dealii_$@" .
push:
	docker push pymor/dealii-testing_pymor_master

all: dealiis
