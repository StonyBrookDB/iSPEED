prefix=/usr/local

libbase=${HADOOPGIS_LIB_PATH}
incbase=${HADOOPGIS_INC_PATH}

builddir=../build/bin


ifdef VALGRINDCHECK
    OPTFLAGS	+= -O0 -g -W #for memory leak checking
else
    OPTFLAGS += -O2
endif

INCFLAGS	= -I $(incbase) -I ./ 
LIBS		= -lgeos -lspatialindex -lboost_program_options -L $(libbase)
CGALFLAGS	= -lCGAL -Wl,-rpath,${libbase}

ifdef DEBUG
    OPTFLAGS += -DNDEBUG
endif


ifdef NOTCOMPRESSED
    OPTFLAGS += -DNCOMPRESSED
else
    OPTFLAGS += -DCOMPRESSED 
endif


ifdef DEBUGTIME
    OPTFLAGS += -DDEBUGTIME
else
    OPTFLAGS += -DNDEBUGTIME
endif
 
   
CC = gcc
CXX = g++

all: queryprocessor_3d  stats_extract_space_dims_3d  manipulate_3d  ppmc  \
	compress_combine  compress_load  fg_3d ot_3d resque_3d


queryprocessor_3d: framework/queryprocessor_3d.cpp
	$(CXX) -std=c++0x $< $(INCFLAGS) $(LIBS) $(OPTFLAGS) -o ${builddir}/$@

stats_extract_space_dims_3d: utilities/stats_extract_space_dims_3d.cpp
	$(CXX) -std=c++0x $< $(INCFLAGS) $(LIBS) $(OPTFLAGS) -o ${builddir}/$@

manipulate_3d: transform/manipulate_3d.cpp
	$(CXX)  -frounding-math -std=c++0x  -DCGAL_EIGEN3_ENABLED $< $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(OPTFLAGS) -o ${builddir}/$@

ppmc: PPMC/ppmc.cpp
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math -O3 $< $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(OPTFLAGS) -o ${builddir}/$@

compress_combine: PPMC/compress_combine.cpp
	$(CXX) -std=c++0x $< $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(OPTFLAGS) -o ${builddir}/$@

compress_load: PPMC/compress_load.cpp
	$(CXX) -std=c++0x $< $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(OPTFLAGS) -o ${builddir}/$@

skeleton_3d: indices/skeleton_3d.cpp
	$(CXX) -std=c++0x  -DCGAL_EIGEN3_ENABLED $< $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(OPTFLAGS) -o ${builddir}/$@

voronoi_3d: indices/voronoi_3d.cpp
	$(CXX) -std=c++0x $< $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(OPTFLAGS) -o ${builddir}/$@



fg_3d: partitionalgo/fg/fg_3d.cpp
	$(CXX) -std=c++0x $< $(INCFLAGS) $(LIBS) $(OPTFLAGS) -o ${builddir}/$@

ot_3d: partitionalgo/ot/ot_3d.cpp
	$(CXX) -std=c++0x $< $(INCFLAGS) $(LIBS) $(OPTFLAGS) -o ${builddir}/$@

resque_3d: resque/resque_3d.cpp
	$(CXX)  -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math  -DCGAL_EIGEN3_ENABLED $< $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(OPTFLAGS) -o ${builddir}/$@


install:
	mkdir -p ${builddir}

clean:
	@rm -f ${builddir}/*

