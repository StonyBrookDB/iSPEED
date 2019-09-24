/* This contains headers (dependencies) of RESQUE spatial processing engine */
#ifndef RESQUE_3D_HPP
#define RESQUE_3D_HPP


#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <map>
#include <cstdlib>
#include <getopt.h>
#include <time.h>
#include <sys/shm.h>
#include <unordered_map>

#include <boost/algorithm/string/replace.hpp>

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

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <fstream>

#include <boost/foreach.hpp>
#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/OFF_to_nef_3.h>



//compressed data

#include <spatialindex/SpatialIndex.h>

#include <progparams/string_constants.h>
#include <utilities/tokenizer.h>

//#include <extensions/rtree3d/rtree_traversal_3d.h>

// Constants
#include <progparams/resque_constants_3d.h>
// Program parameters
#include <PPMC/configuration.h>
#include <PPMC/mymesh.h>
#include "../progparams/resque_params_3d.hpp"

typedef CGAL::Bbox_3                                     Bbox;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef CGAL::Exact_predicates_exact_constructions_kernel  Kernel;
typedef Kernel::Triangle_3                               Triangle;
typedef std::vector<Triangle>                               Triangles;
typedef Triangles::iterator                                   Iterator;
typedef Kernel::Segment_3                                Segment;
typedef CGAL::Polyhedron_3<Kernel>                       Polyhedron;
typedef Polyhedron::Halfedge_const_handle                Halfedge_const_handle;
typedef Polyhedron::Facet_const_iterator                 Facet_const_iterator;
typedef Polyhedron::Facet_const_handle                   Facet_const_handle;

//typedef CGAL::Box_intersection_d::Box_with_handle_d<double, 3, Facet_const_handle>               Box;
typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,Iterator> Box;

typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Kernel::Vector_3  Vector_3;
typedef CGAL::Triangulation_3<Kernel> Triangulation;
typedef Triangulation::Point        CGAL_Point;
typedef Triangulation::Tetrahedron 	Tetrahedron;
typedef Nef_polyhedron::Volume_const_iterator Volume_const_iterator;

typedef Kernel::Point_3                                       CGAL_Point3;
typedef CGAL::Delaunay_triangulation_3<Kernel, CGAL::Fast_location> Delaunay;


// Skeleton related
typedef CGAL::Simple_cartesian<double>		Sk_Kernel;
typedef Sk_Kernel::Point_3					Sk_Point;
typedef CGAL::Polyhedron_3<Sk_Kernel>		Sk_Polyhedron;

typedef boost::graph_traits<Sk_Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef CGAL::Mean_curvature_flow_skeletonization<Sk_Polyhedron> Skeletonization;
typedef Skeletonization::Skeleton                             	 Skeleton;
typedef Skeleton::vertex_descriptor                              Skeleton_vertex;
typedef CGAL::Delaunay_triangulation_3<Sk_Kernel, CGAL::Fast_location> Sk_Delaunay;


// for AABB tree distance calculation
typedef CGAL::Triangulation_3<Sk_Kernel> Sk_Triangulation;
typedef Sk_Triangulation::Point        Sk_CGAL_Point;

typedef Sk_Kernel::FT Sk_FT;
typedef CGAL::AABB_face_graph_triangle_primitive<Sk_Polyhedron> SK_Primitive;
typedef CGAL::AABB_traits<Sk_Kernel, SK_Primitive> Sk_Traits;
typedef CGAL::AABB_tree<Sk_Traits> Sk_Tree;
typedef Sk_Tree::Point_and_primitive_id Point_and_primitive_id;

// Skeleton related
//typedef CGAL::Simple_cartesian<double>                        Sc_Kernel;
//typedef Sc_Kernel::Point_3                                       Sc_Point;
//typedef CGAL::Polyhedron_3<Sc_Kernel>                            Sc_Polyhedron;
//typedef Sc_Kernel::Triangle_3                               Sc_Triangle;
//typedef std::vector<Sc_Triangle>                               Sc_Triangles;

//typedef Kernel::FT FT;
// for AABB tree distance calculation

/* // Haders for Yanhui's change 
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
*/

//
// // Headers for Hoang's change
//typedef std::list<Triangle>::iterator Cgal_Iterator;
//typedef CGAL::AABB_triangle_primitive<Kernel,Cgal_Iterator> Primitive;
//typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
////typedef Tree::Point_and_primitive_id Point_and_primitive_id;
////typedef CGAL::AABB_traits<Sc_Kernel, Sc_Primitive> AABB_triangle_traits;
////typedef CGAL::AABB_tree<AABB_triangle_traits> Sc_Tree;
//
//// typedef CGAL::AABB_face_graph_triangle_primitive<Sc_Polyhedron> Sc_Primitive;
////typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
//typedef CGAL::AABB_tree<Traits> Tree;
////typedef std::list<Triangle>::iterator Iterator;
//typedef Tree::Point_and_primitive_id Point_and_primitive_id;

// Constants used for building the R-tree
#define FillFactor 0.9
#define IndexCapacity 10 
#define LeafCapacity 50
#define COMPRESS true
#include <indices/rtree_builder_3d.hpp>

// Default key
#define COMPRESSION_KEY 5678
#define SHMSZ 10000000000
#define NUMBER_DIMENSIONS 3


//clock variables
extern clock_t start_reading_data;
extern clock_t start_query_exec;

extern clock_t total_reading;
extern clock_t total_query_exec;
//using namespace SpatialIndex;

/* Function prototypes */

int read_cache_file(struct query_op &stop, struct query_temp &sttemp);
void release_mem(struct query_op &stop, struct query_temp &sttemp, int maxCard);
void obtain_field(struct query_op &stop, struct query_temp &sttemp, 
	int position, int pos1, int pos2);
void obtain_field(struct query_op &stop, struct query_temp &sttemp, 
	int position, std::vector<std::string> &set1fields, int pos2);
void update_nn(struct query_op &stop, struct query_temp &sttemp, 
	int object_id, double distance);
//void update_bucket_dimension(const geos::geom::Envelope * env);
//double get_distance(const geos::geom::Point* p1, const geos::geom::Point* p2);
//double get_distance_earth(const geos::geom::Point* p1,const geos::geom::Point* p2);
//bool build_index_geoms();
//using namespace SpatialIndex;
//using namespace geos;
//using namespace geos::geom;

//bool build_index_geoms(map<int,Geometry*> & geom_polygons, ISpatialIndex* & spidx, IStorageManager* & storage);
bool build_index_geoms(std::vector<struct mbb_3d *> & geom_mbbs, SpatialIndex::ISpatialIndex* & spidx, SpatialIndex::IStorageManager* & storage);

bool join_with_predicate(struct query_op &stop, struct query_temp &sttemp,
		Polyhedron &geom1 , Polyhedron &geom2,
		const struct mbb_3d * env1, const struct mbb_3d * env2, const int jp); // for 3d spatial join
Polyhedron extract_geometry(long offset, long length, unsigned i_decompPercentage,
	struct query_op &stop, struct query_temp &sttemp, int dataset_id); // to extract geometry from compressed data
Sk_Polyhedron sk_extract_geometry(long offset, long length, unsigned i_decompPercentage,
	struct query_op &stop, struct query_temp &sttemp, int dataset_id);

void report_result(struct query_op &stop, struct query_temp &sttemp, int i, int j);
void report_result(struct query_op &stop, struct query_temp &sttemp, 
	std::vector<std::string> &set1fields, int j, bool skip_window_data);
int join_bucket_spjoin(struct query_op &stop, struct query_temp &sttemp);
int join_bucket_nn_voronoi(struct query_op &stop, struct query_temp &sttemp);
int join_bucket_nn_rtree(struct query_op &stop, struct query_temp &sttemp);
int join_bucket_spjoin(struct query_op &stop, struct query_temp &sttemp);

int join_bucket(struct query_op &stop, struct query_temp &sttemp);

int execute_query(struct query_op &stop, struct query_temp &sttemp);

//for 3d spatial join
double get_volume(Nef_polyhedron &inputpoly);
void get_triangle(Polyhedron P, std::vector<Triangle>& triangles,  std::vector<Box>& boxes, std::vector<Box*>& ptr);
bool intersects(Polyhedron &P1, Polyhedron &P2, const struct mbb_3d * env1, const struct mbb_3d * env2);


struct Report {
  Triangles* triangles;
  Triangles* cell_triangles;

  Report(Triangles& triangles, Triangles& cell_triangles)
    : triangles(&triangles), cell_triangles(&cell_triangles)
  {}

  // callback functor that reports all truly intersecting triangles
  void operator()(const Box* a, const Box* b) const
  {
    if (intersection_flag) {
    	return;
    }
    if ( ! a->handle()->is_degenerate() && ! b->handle()->is_degenerate()
         && CGAL::do_intersect( *(a->handle()), *(b->handle()))) {
      intersection_flag = true;
     // std::cerr << "Intersection? " << intersection_flag << std::endl;
    }
  }
};

#endif

