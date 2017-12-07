#ifndef RESQUE_DATASTRUCTS_3D
#define RESQUE_DATASTRUCTS_3D

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <sstream>
//#include <unordered_set>
//#include <unordered_map>

#include <CGAL/Random.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/intersections.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/convex_decomposition_3.h> 
#include <CGAL/Tetrahedron_3.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <fstream>


#include <boost/foreach.hpp>

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Delaunay_triangulation_3.h>

//#include <progparams/string_constants.h>
//#include <utilities/tokenizer.h>

//typedef CGAL::Exact_integer  NT;
//instead of
//typedef CGAL::Extended_homogeneous<NT>  Kernel;
// workaround for VC++
//struct Kernel : public CGAL::Extended_homogeneous<NT> {};


#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <spatialindex/SpatialIndex.h>

const int OSM_SRID = 4326;
const int ST_INTERSECTS = 1;
const int ST_TOUCHES = 2;
const int ST_CROSSES = 3;
const int ST_CONTAINS = 4;
const int ST_ADJACENT = 5;
const int ST_DISJOINT = 6;
const int ST_EQUALS = 7;
const int ST_DWITHIN = 8;
const int ST_WITHIN = 9;
const int ST_OVERLAPS = 10;
const int ST_NEAREST = 11;
const int ST_NEAREST_2 = 12;
const int ST_NN_VORONOI = 13;
const int ST_NN_RTREE = 14;

const std::string PARAM_PREDICATE_INTERSECTS = "st_intersects";
const std::string PARAM_PREDICATE_NEAREST_NEIGHBOR_VORONOI = "st_nn_voronoi"; // nn for 3d
const std::string PARAM_PREDICATE_NEAREST_NEIGHBOR_RTREE = "st_nn_rtree"; // nn for 3d
const std::string PARAM_PREDICATE_TOUCHES = "st_touches";
const std::string PARAM_PREDICATE_CROSSES = "st_crosses";
const std::string PARAM_PREDICATE_CONTAINS = "st_contains";
const std::string PARAM_PREDICATE_ADJACENT = "st_adjacent";
const std::string PARAM_PREDICATE_DISJOINT = "st_disjoint";
const std::string PARAM_PREDICATE_EQUALS = "st_equals";
const std::string PARAM_PREDICATE_DWITHIN = "st_dwithin";
const std::string PARAM_PREDICATE_WITHIN = "st_within";
const std::string PARAM_PREDICATE_OVERLAPS = "st_overlaps";
const std::string PARAM_PREDICATE_NEAREST = "st_nearest";
const std::string PARAM_PREDICATE_NEAREST_NO_BOUND = "st_nearest2";

const int SID_1 = 1;
const int SID_2 = 2;
const int SID_NEUTRAL = 0;

const int STATS_AREA_1 = -1;
const int STATS_AREA_2 = -2;
const int STATS_UNION_AREA = -3;
const int STATS_INTERSECT_AREA = -4;
const int STATS_JACCARD_COEF = -5;
const int STATS_DICE_COEF = -6;
const int STATS_TILE_ID = -7;
const int STATS_MIN_DIST = -8;
//for 3d
const int STATS_VOLUME_1 = -9;
const int STATS_VOLUME_2 = -10;
const int STATS_INTERSECT_VOLUME = -11;
const int STATS_NN_DISTANCE = -12;

const std::string PARAM_STATS_AREA_1 = "area1";
const std::string PARAM_STATS_AREA_2 = "area2";
const std::string PARAM_STATS_UNION_AREA = "union";
const std::string PARAM_STATS_INTERSECT_AREA = "intersect";
const std::string PARAM_STATS_JACCARD_COEF = "jaccard";
const std::string PARAM_STATS_DICE_COEF = "dice";
const std::string PARAM_STATS_TILE_ID= "tileid";
const std::string PARAM_STATS_MIN_DIST = "mindist";

//for 3d
const std::string PARAM_STATS_VOLUME_1 = "volume1";
const std::string PARAM_STATS_VOLUME_2 = "volume2";
const std::string PARAM_STATS_INTERSECT_VOLUME = "intersect_volume";
const std::string PARAM_STATS_NN_DISTANCE = "nn_distance";

const std::string STR_SET_DELIM = ":";
const std::string STR_OUTPUT_DELIM = ",";

// for sp join yes or no intersection
bool intersection_flag = false;

char * resque_decomp_buffer = NULL;
char * shm_ptr = NULL;
//long maxoffset = 0;
std::string dummyoutputname = "/tmp/nonsense";

// for 3d mbb
struct mbb_3d {
	double low[3];
	double high[3];
};

// for 3d skeleton vertex
/*struct id_point{
	int id;
	Point p;
};*/

/* Placeholder for nearest neighbor rankings */
struct query_nn_dist {
	int object_id;
	double distance;
};


/* Struct to hold temporary values during spatial operations */
struct query_temp {
	double area1;
	double area2;
	double union_area;
	double intersect_area;
	double dice;
	double jaccard;
	double distance;

	// for 3d
	double volume1;
	double volume2;
	double intersect_volume;
	//Delaunay* voronoi;
	double nn_distance;

	std::stringstream stream;
	std::string tile_id;

	/* Bucket information */
	double bucket_min_x;
	double bucket_min_y;
	double bucket_max_x;
	double bucket_max_y;

	/* Data current to the tile/bucket */
	//std::map<int, std::vector< std::vector<std::string> > > rawdata;
	//std::map<int, std::vector<geos::geom::Geometry*> > polydata;
	
	/* Data current to the tile/bucket */
	std::map<int, std::vector< std::vector<std::string> > > rawdata;
	std::map<int, std::vector<struct mbb_3d *> > mbbdata;
	//std::map<std::vector<Point*>, int > idpoint;
//	std::vector<double> mbbdata_vector;
	//std::map<int, std::vector<Polyhedron *> > polydata;
	std::map<int, std::vector<CGAL::Polyhedron_3<CGAL::Exact_predicates_exact_constructions_kernel> *> > polydata;
	
	// for comprssed data
	std::map<int, std::vector<long> > offsetdata;
	std::map<int, std::vector<long> > lengthdata;
	std::istringstream poly_str[2]; // 0 for data set 1, 1 for data set 2
	//std::stringstream poly_str[2]; // 0 for data set 1, 1 for data set 2

	/* Nearest neighbor temporary placeholders */
	std::list<struct query_nn_dist*> nearest_distances;
};


/* Query operator */
struct query_op {
	bool use_cache_file;
	bool reading_mbb; // input is mbb
	char* cachefilename;
	

	int join_predicate; /* Join predicate - see resquecommon.h for the full list*/
	int shape_idx_1; /* geometry field number of the 1st set */
	int shape_idx_2; /* geometry field number of the 2nd set */
	int join_cardinality; /* Join cardinality */
	double expansion_distance; /* Distance used in dwithin and k-nearest query */
	int k_neighbors; /* k - the number of neighbors used in k-nearest neighbor */
	/* the number of additional meta data fields
	 i.e. those first fields not counting towards object data
	*/

	int sid_second_set; // set id for the "second" dataset.
	// Value == 2 for joins between 2 data sets
	// Value == 1 for selfjoin

	/* Mapping-specific */
	bool extract_mbb;
	bool collect_mbb_stat;
	bool use_sampling;
	double sample_rate;
	bool drop_join_idx; // Does not append join index when emitting object to tile

	//map<int, geos::geom::Geometry*> geom_tiles; // mapping of actual geomery for tiles
	std::map<int, long> count_tiles; // mapping of actual geomery index for tiles
	char* prefix_1; // directory prefix for the files from dataset 1
	char* prefix_2; // directory prefix for the files from dataset 2

	/* Reducing-specific - RESQUE */
	int offset; // offset/adjustment for meta data field in RESQUE
	std::vector<int> proj1; /* Output fields for 1st set  */
	std::vector<int> proj2; /* Output fields for 2nd set */

	/* Output fields -  parallel arrays/lists */
	std::vector<int> output_fields; // fields to the output
	// e.g. 1 is field #1, 2 is fiend #2, the dataset they belong to 
	// is stored in corresponding position 
	//in the the output_fields_set_id list 
	std::vector<int> output_fields_set_id; // meta information to indicate fields to output

	/* Indicate whether symmetric result pair should be included or not */	
	bool result_pair_duplicate;

	/* Indicate whether to use the geographical distance for points on earth */
	bool use_earth_distance;

	/* Special bits to indicate whether statistics/fields to compute */
	bool needs_area_1;
	bool needs_area_2;
	bool needs_union;
	bool needs_intersect;
	bool needs_dice;
	bool needs_jaccard;
	bool needs_min_distance;

	// for 3d
	bool needs_volume_1;
	bool needs_volume_2;
	bool needs_intersect_volume;
	bool needs_nn_distance;

	// for 3d compression
	size_t shm_max_size; // size of the shared mem segment
	unsigned decomp_lod; // decompression level/level of details 0 to 100

};

// Function prototypes defined in resque_params
void usage();
void set_projection_param(char * arg);
int get_join_predicate(char * predicate_str);
bool extract_params(int argc, char** argv, struct query_op & stop, struct query_temp &sttemp);
#endif
