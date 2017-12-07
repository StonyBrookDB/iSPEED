#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <map>
#include <vector>
#include <cstdlib> 

#include <progparams/string_constants.h>
#include <utilities/tokenizer.h>

#include <progparams/partition_params_3d.hpp>

#include <boost/program_options.hpp>

static int GLOBAL_MAX_LEVEL = 10000000;

using namespace std;

class SpatialObject {
	public:
		double low[3];
		double high[3];
		SpatialObject(double min_x, double min_y, double min_z, double max_x, double max_y, double max_z);
};

class OctreeNode {
	public:
		double low[3];
		double high[3];
		int level;
		bool isLeaf;
		bool canBeSplit;
		int size;
		OctreeNode* children[8];
		vector<SpatialObject*> objectList;

		OctreeNode(double min_x, double min_y, double min_z, double max_x, 
				double max_y, double max_z, int level);

		~OctreeNode();

		bool addObject(SpatialObject *object);
		bool intersects(SpatialObject *object);
		bool addObjectIgnore(SpatialObject *object);    
};

#include <partitionalgo/ot/OctreeNode.hpp>
