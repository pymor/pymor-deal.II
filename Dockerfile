FROM dealii/base:gcc-mpi
MAINTAINER rene.milk@wwu.de

ARG DEALII_VERSION

USER root
ENV DEBIAN_FRONTEND noninteractive
RUN echo "deb mirror://mirrors.ubuntu.com/mirrors.txt xenial main restricted universe multiverse" > /etc/apt/sources.list && \
    echo "deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-updates main restricted universe multiverse" >> /etc/apt/sources.list && \
    echo "deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-security main restricted universe multiverse" >> /etc/apt/sources.list && \
    echo "deb-src mirror://mirrors.ubuntu.com/mirrors.txt xenial main restricted universe multiverse" >> /etc/apt/sources.list && \
    echo "deb-src mirror://mirrors.ubuntu.com/mirrors.txt xenial-updates main restricted universe multiverse" >> /etc/apt/sources.list && \
    echo "deb-src mirror://mirrors.ubuntu.com/mirrors.txt xenial-security main restricted universe multiverse" >> /etc/apt/sources.list
RUN apt-get update
#for add-apt-repo
RUN apt-get install -y software-properties-common python-pip
RUN add-apt-repository -y ppa:pymor/stable
RUN apt-get update
#RUN apt-get build-dep -y --no-install-recommends python-pymor-demos
RUN pip install -U pip
RUN pip install pymor==0.4.1
USER ${USER}

RUN git clone https://github.com/dealii/dealii.git dealii-${DEALII_VERSION}-src && \
    cd dealii-${DEALII_VERSION}-src && \
    git checkout ${DEALII_VERSION} && \
    mkdir build && cd build && \
    cmake -DDEAL_II_WITH_MPI=ON \
          -DDEAL_II_COMPONENT_EXAMPLES=OFF \
          -DCMAKE_INSTALL_PREFIX=~/dealii-${DEALII_VERSION} \
          -DCMAKE_BUILD_TYPE=Release \
          -GNinja \
          ../

# THE END
ENV DEBIAN_FRONTEND teletype
