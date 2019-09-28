FROM ubuntu

MAINTAINER Ivan Ermilov <ivan.s.ermilov@gmail.com>
MAINTAINER Giannis Mouchakis <gmouchakis@iit.demokritos.gr>

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      openjdk-8-jdk \
      net-tools \
      curl \
      netcat \
      gnupg \
      apt-utils \
      gcc g++ cmake libboost-all-dev libgeos++-dev freeglut3-dev libcgal-dev software-properties-common wget git \
    && rm -rf /var/lib/apt/lists/*
      
ENV JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/

RUN curl -O https://dist.apache.org/repos/dist/release/hadoop/common/KEYS

RUN gpg --import KEYS

ENV HADOOP_VERSION 3.1.2
ENV HADOOP_URL http://mirror.cc.columbia.edu/pub/software/apache/hadoop/common/hadoop-$HADOOP_VERSION/hadoop-$HADOOP_VERSION.tar.gz

RUN set -x \
    && curl -fSL "$HADOOP_URL" -o /tmp/hadoop.tar.gz \
    && tar -xvf /tmp/hadoop.tar.gz -C /opt/ \
    && rm /tmp/hadoop.tar.gz*

RUN ln -s /opt/hadoop-$HADOOP_VERSION/etc/hadoop /etc/hadoop

RUN mkdir /opt/hadoop-$HADOOP_VERSION/logs

RUN mkdir /hadoop-data

ENV HADOOP_HOME=/opt/hadoop-$HADOOP_VERSION
ENV HADOOP_CONF_DIR=/etc/hadoop
ENV MULTIHOMED_NETWORK=1
ENV USER=root
ENV PATH $HADOOP_HOME/bin/:$PATH

ADD entrypoint.sh /entrypoint.sh

RUN chmod a+x /entrypoint.sh


# to install the running environment of ispeed

# install dependencies
RUN add-apt-repository ppa:ubuntugis/ppa && \
    apt-get update && \
    apt-get install -y gdal-bin \
    && rm -rf /var/lib/apt/lists/*

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

ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
# setting up the environment

ENTRYPOINT ["/entrypoint.sh"]
