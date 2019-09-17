prefix=/usr/local

libbase=${HADOOPGIS_LIB_PATH}
incbase=${HADOOPGIS_INC_PATH}

builddir=../build
bindir=${builddir}/bin
libdir=${builddir}/lib



CPPFLAGS	= -std=c++0x -g -O2 -Wno-unused-result
INCFLAGS	= -I ./ 
LIBS		= -lgeos -lspatialindex -lboost_program_options -lstdc++ -lgmp -lmpfr
CGALFLAGS	= -lCGAL -Wl,-rpath

ifdef DEBUG
    CPPFLAGS += -DDEBUG
else
    CPPFLAGS += -DNDEBUG
endif
 
   
CC = gcc
CXX = g++


#compile all the cpp files
SRCS := $(wildcard */*.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

PPMC_SRCS := $(wildcard PPMC/*.cpp)
PPMC_OBJS := $(patsubst %.cpp,%.o,$(PPMC_SRCS))

%.o: %.cpp
	$(CXX) $(INCFLAGS) $(CPPFLAGS) -c $? -o $@


#compile additional range coder .c files
RC_SRCS := $(wildcard PPMC/rangeCoder/*.c)
RC_OBJS := $(patsubst %.c,%.o,$(RC_SRCS))

%.o: %.c
	$(CC) $(INCFLAGS) -c $? -o $@


all: queryprocessor_3d  stats_extract_space_dims_3d  manipulate_3d  ppmc  \
	compress_combine  compress_load  fg_3d ot_3d resque_3d sampler

queryprocessor_3d: framework/queryprocessor_3d.o mapreducecpp/mapreduce_cpp.o progparams/queryprocessor_params_3d.o
	$(CXX) $? $(INCFLAGS) $(LIBS) $(CPPFLAGS) -o ${bindir}/$@ -Wno-unused-result

manipulate_3d: transform/manipulate_3d.o
	$(CXX) -frounding-math -DCGAL_EIGEN3_ENABLED $? $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ${bindir}/$@

ppmc: compression/ppmc.o $(PPMC_OBJS) $(RC_OBJS)
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $? $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ${bindir}/$@

resque_3d: resque/resque_3d.o resque/spjoin.o progparams/resque_params_3d.o $(PPMC_OBJS) $(RC_OBJS)
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math -DCGAL_EIGEN3_ENABLED $? $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ${bindir}/$@

#some other simple programs

stats_extract_space_dims_3d: utilities/stats_extract_space_dims_3d.o
	$(CXX) $? $(INCFLAGS) $(LIBS) $(CPPFLAGS) -o ${bindir}/$@

compress_combine: compression/compress_combine.o
	$(CXX) $? $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ${bindir}/$@

compress_load: compression/compress_load.o
	$(CXX) $? $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ${bindir}/$@

skeleton_3d: indices/skeleton_3d.o
	$(CXX) -DCGAL_EIGEN3_ENABLED $? $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ${bindir}/$@

voronoi_3d: indices/voronoi_3d.o
	$(CXX) $? $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ${bindir}/$@

sampler: transform/sampler.o
	$(CXX) $? $(INCFLAGS) $(CPPFLAGS) -o ${bindir}/$@

fg_3d: partitionalgo/fg/fg_3d.o progparams/partition_params_3d.o
	$(CXX) $? $(INCFLAGS) $(LIBS) $(CPPFLAGS) -o ${bindir}/$@

ot_3d: partitionalgo/ot/ot_3d.o progparams/partition_params_3d.o
	$(CXX) $? $(INCFLAGS) $(LIBS) $(CPPFLAGS) -o ${bindir}/$@


.PHONY: all clean

install:
	rm -rf ${builddir} && \
	mkdir ${builddir} && \
	mkdir ${bindir} && \
	mkdir ${libdir} && \
	cp mapreducejava/myCustomLibs.jar ${libdir} && \
	cp ../script/* ${bindir}/

clean:
	rm -rf ${builddir} *.o */*.o */*/*.o

