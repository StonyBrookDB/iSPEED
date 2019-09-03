FROM ubuntu:latest

# install dependencies
RUN apt-get update && \
    apt-get install -y gcc g++ cmake libboost-all-dev libgeos++-dev freeglut3-dev libcgal-dev software-properties-common wget git && \
    add-apt-repository ppa:ubuntugis/ppa && \
    apt-get update && \
    apt-get install -y gdal-bin

# install eigen
RUN wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.bz2 && \
    tar xvf 3.3.7.tar.bz2 && \
    mv eigen-eigen-323c052e1731/Eigen /usr/local/include && \
    rm -rf eigen-eigen-323c052e1731 && \
    rm 3.3.7.tar.bz2

# install MPFR
RUN wget https://www.mpfr.org/mpfr-current/mpfr-4.0.2.tar.xz && \
    tar xvf mpfr-4.0.2.tar.xz && \
    cd mpfr-4.0.2/ && \
    ./configure --prefix=/usr/local/ && \
    make -j20 && \
    make install && \
    cd ../ && \
    rm -rf mpfr-4.0.2*

# install spatialindex
RUN wget http://download.osgeo.org/libspatialindex/spatialindex-src-1.8.5.tar.gz && \
    tar zxvf spatialindex-src-1.8.5.tar.gz && \
    cd spatialindex-src-1.8.5 && \
    ./configure --prefix=/usr/local && \
    make -j20 && \
    make install && \
    cd ../ && \
    rm -rf spatialindex*

# setting up the environment
ENV C_INCLUDE_PATH=/iSPEED/src/:/usr/local/include \
    CPLUS_INCLUDE_PATH=/iSPEED/src/:/usr/local/include \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

# install iSPEED
COPY src /iSPEED/src
RUN cd /iSPEED && \
    mkdir -p build/bin && \
    cd src && \
    make
