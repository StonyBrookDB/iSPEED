#ifndef RESQUE_DATASTRUCTS_3D
#define RESQUE_DATASTRUCTS_3D

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include <unistd.h>
#include <getopt.h>

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

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <spatialindex/SpatialIndex.h>
#include <progparams/string_constants.h>

const int OSM_SRID = 4326;
enum Jointype{
	ST_ERROR = 0,
	ST_INTERSECTS = 1,
	ST_TOUCHES = 2,
	ST_CROSSES = 3,
	ST_CONTAINS = 4,
	ST_ADJACENT = 5,
	ST_DISJOINT = 6,
	ST_EQUALS = 7,
	ST_DWITHIN = 8,
	ST_WITHIN = 9,
	ST_OVERLAPS = 10,
	ST_NEAREST = 11,
	ST_NEAREST_2 = 12,
	ST_NN_VORONOI = 13,
	ST_NN_RTREE = 14
};

const std::string join_type_str[15] = {
		"st_error", "st_intersects", "st_touches", "st_crosses", "st_contains"
		"st_adjacent", "st_disjoint", "st_equals", "st_dwithin", "st_within",
		"st_overlaps", "st_nearest", "st_nearest2", "st_nn_voronoi", "st_nn_rtree"};

const int SID_1 = 1;
const int SID_2 = 2;
const int SID_NEUTRAL = 0;


/*
 * Struct to hold temporary values during spatial operations
 * for one tile, it contains the mbb, offset, length of objects which
 * are parsed from lines in stdin
 * */
struct query_temp {

	std::string tile_id;
	std::map<int, std::vector<struct mbb_3d *> > mbbdata;
	// for compressed data
	std::map<int, std::vector<long> > offsetdata;
	std::map<int, std::vector<long> > lengthdata;

};


/* Query operator */
struct query_op {

	Jointype join_predicate = ST_INTERSECTS; /* Join predicate*/
	int join_cardinality = 2; /* Join cardinality, can only be 1 or 2*/
	double expansion_distance = 0; /* Distance used in dwithin and k-nearest query */
	int k_neighbors = 1; /* k - the number of neighbors used in k-nearest neighbor */

	// for 3d
	bool needs_volume_1 = true;
	bool needs_volume_2 = true;
	bool needs_intersect_volume = true;
	bool needs_nn_distance = true;

	// for 3d compression
	size_t shm_max_size = 0; // size of the shared mem segment
	unsigned decomp_lod = 100; // decompression level/level of details 0 to 100

};

// Function prototypes defined in resque_params
void usage();
bool extract_params(int argc, char** argv, struct query_op & stop, struct query_temp &sttemp);
#endif
