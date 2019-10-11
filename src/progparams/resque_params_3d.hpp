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

#include "global_define.h"

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
