FROM hadoop-base:latest

# setting up the environment
ENV C_INCLUDE_PATH=/iSPEED/src/:/usr/local/include \
    CPLUS_INCLUDE_PATH=/iSPEED/src/:/usr/local/include \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib \
    PATH=$PATH:/iSPEED/build/bin \
    DEBUG=TRUE

# copy in and compile iSPEED
COPY src /iSPEED/src
COPY script /iSPEED/script

RUN cd /iSPEED/src && \
    make install && \
    make

WORKDIR /iSPEED/build/bin


