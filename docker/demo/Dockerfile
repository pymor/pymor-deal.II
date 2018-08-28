ARG PYVER=3.5
FROM pymor/python:$PYVER

MAINTAINER rene.milk@wwu.de

ARG DEALII_VERSION

USER root
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update
#for add-apt-repo
RUN apt-get install -y software-properties-common \
    gfortran cmake \
    git \
    gsl-bin \
    libblas-dev \
    libbz2-dev \
    libgsl-dev \
    liblapack-dev \
    libnetcdf-c++4-dev \
    libnetcdf-cxx-legacy-dev \
    libnetcdf-dev \
    ninja-build \
    numdiff \
    unzip \
    wget \
    zlib1g-dev \
    libopenmpi-dev \
    openmpi-bin

RUN pip install -U pip

RUN cd /tmp && git clone https://github.com/dealii/dealii.git dealii-${DEALII_VERSION}-src && \
    cd dealii-${DEALII_VERSION}-src && \
    git checkout v${DEALII_VERSION} && \
    mkdir build && cd build && \
    cmake -DDEAL_II_WITH_MPI=ON \
          -DDEAL_II_COMPONENT_EXAMPLES=OFF \
          -DCMAKE_INSTALL_PREFIX=/usr/local \
          -DCMAKE_BUILD_TYPE=Release \
          -GNinja \
          ../ && \
    ninja && ninja install && \
    cd / && rm -rf /tmp/deal*

# THE END
ENV DEBIAN_FRONTEND teletype

RUN apt-get install -y gosu

ADD entrypoint.sh /usr/local/bin/
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
