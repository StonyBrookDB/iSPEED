
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

#include <progparams/string_constants.h>
#include <utilities/tokenizer.h>

// Constants
#include <progparams/resque_constants_3d.h>

// Program parameters
#include <progparams/resque_params_3d.hpp>

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

bool build_index_tiles(struct query_op &stop, struct query_temp &sttemp,
	IStorageManager* &storage, ISpatialIndex * &spidx,
	std::map<SpatialIndex::id_type, std::string> *id_tiles);
void process_input(struct query_op &stop, struct query_temp &sttemp,
		const int join_idx, const int geom_idx, 
		IStorageManager * &storage, ISpatialIndex * &spidx,
		std::map<id_type, string> *id_tiles);
void init(struct query_op &stop, struct query_temp &sttemp);
