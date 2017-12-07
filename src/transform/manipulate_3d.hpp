
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

typedef CGAL::Bbox_3                                     Bbox;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef CGAL::Exact_predicates_exact_constructions_kernel  Kernel;
typedef Kernel::Triangle_3                               Triangle;
typedef Kernel::Segment_3                                Segment;
typedef CGAL::Polyhedron_3<Kernel>                       Polyhedron;
typedef Polyhedron::Halfedge_const_handle                Halfedge_const_handle;
typedef Polyhedron::Facet_const_iterator                 Facet_const_iterator;
typedef Polyhedron::Facet_const_handle                   Facet_const_handle;
//typedef CGAL::Box_intersection_d::Box_with_handle_d<double, 3, Facet_const_handle>               Box;
//typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,Iterator> Box;
typedef Kernel::FT FT;

typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Kernel::Vector_3  Vector_3;
typedef CGAL::Triangulation_3<Kernel> Triangulation;
typedef Triangulation::Point        CGAL_Point;
typedef Triangulation::Tetrahedron 	Tetrahedron;
typedef Nef_polyhedron::Volume_const_iterator Volume_const_iterator;


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
bool process_input(struct query_op &stop, struct query_temp &sttemp,
		const int join_idx, const int geom_idx, 
		IStorageManager * &storage, ISpatialIndex * &spidx,
		std::map<id_type, string> *id_tiles);
void init(struct query_op &stop, struct query_temp &sttemp);
