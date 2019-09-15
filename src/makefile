prefix=/usr/local

libbase=${HADOOPGIS_LIB_PATH}
incbase=${HADOOPGIS_INC_PATH}

builddir=../build/bin

ifdef VALGRINDDEBUG
    CPPFLAGS += -O0 -g -DDEBUG#for memory leak checking
else
    CPPFLAGS += -O2 
endif

CPPFLAGS	= -std=c++0x -DCOMPRESSED
INCFLAGS	= -I ./ 
LIBS		= -lgeos -lspatialindex -lboost_program_options -lstdc++ -lgmp -lmpfr
CGALFLAGS	= -lCGAL -Wl,-rpath

ifdef DEBUG
    CPPFLAGS += -DDEBUGTIME -DEBUG
else
    CPPFLAGS += -DNDEBUGTIME -DNDEBUG
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
	compress_combine  compress_load  fg_3d ot_3d resque_3d

queryprocessor_3d: framework/queryprocessor_3d.o mapreducecpp/mapreduce_cpp.o
	$(CXX) $? $(INCFLAGS) $(LIBS) $(CPPFLAGS) -o ${builddir}/$@ -Wno-unused-result

manipulate_3d: transform/manipulate_3d.o
	$(CXX) -frounding-math -DCGAL_EIGEN3_ENABLED $? $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ${builddir}/$@

ppmc: compression/ppmc.o $(PPMC_OBJS) $(RC_OBJS)
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $? $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ${builddir}/$@

resque_3d: resque/resque_3d.o
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math -DCGAL_EIGEN3_ENABLED $? $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ${builddir}/$@

#some other simple programs


stats_extract_space_dims_3d: utilities/stats_extract_space_dims_3d.o
	$(CXX) $? $(INCFLAGS) $(LIBS) $(CPPFLAGS) -o ${builddir}/$@

compress_combine: compression/compress_combine.o
	$(CXX) $? $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ${builddir}/$@

compress_load: compression/compress_load.o
	$(CXX) $? $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ${builddir}/$@

skeleton_3d: indices/skeleton_3d.o
	$(CXX) -DCGAL_EIGEN3_ENABLED $? $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ${builddir}/$@

voronoi_3d: indices/voronoi_3d.o
	$(CXX) $? $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ${builddir}/$@

sampler: transform/sampler.o
	$(CXX) $? $(INCFLAGS) $(CPPFLAGS) -o ${builddir}/$@

fg_3d: partitionalgo/fg/fg_3d.o
	$(CXX) $? $(INCFLAGS) $(LIBS) $(CPPFLAGS) -o ${builddir}/$@

ot_3d: partitionalgo/ot/ot_3d.o
	$(CXX) $? $(INCFLAGS) $(LIBS) $(CPPFLAGS) -o ${builddir}/$@


.PHONY: all clean

install:
	mkdir -p ${builddir}

clean:
	@rm -f ${builddir}/* *.o */*.o */*/*.o

