
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <map>
#include <cstdlib>
#include <getopt.h>
#include <time.h>
#include <float.h>
#include <sys/time.h>



#include <spatialindex/SpatialIndex.h>

#include <utilities/tokenizer.h>
#include "../progparams/global_define.h"
//#include <progparams/resque_params_3d.hpp>

// Constants used for building the R-tree
#define FillFactor 0.9
#define IndexCapacity 10 
#define LeafCapacity 50
#define COMPRESS true
#include <indices/rtree_builder_3d.hpp>

using namespace std;
using namespace SpatialIndex;


/* Function protoypes */
//bool build_index_tiles(SpatialIndex::IStorageManager * &storage, SpatialIndex::ISpatialIndex * &spidx);
//bool process_input(const int join_idx, const int geom_idx, SpatialIndex::IStorageManager * &storage, SpatialIndex::ISpatialIndex * &spidx);

bool build_index_tiles(char *cachefile, IStorageManager* &storage, ISpatialIndex * &spidx,
	std::map<SpatialIndex::id_type, std::string> *id_tiles);
void process_input(IStorageManager * &storage, ISpatialIndex * &spidx,
		std::map<id_type, string> *id_tiles);
