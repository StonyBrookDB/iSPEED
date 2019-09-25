/*
 * test.cpp
 *
 *  Created on: Sep 25, 2019
 *      Author: teng
 */
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/mesh_segmentation.h>
#include <fstream>
#include <PPMC/configuration.h>
#include <PPMC/mymesh.h>
#include <sys/shm.h>
#include <iostream>

typedef CGAL::Simple_cartesian<double>                               Kernel;
typedef Kernel::Point_3                                              CGAL_Point;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor           vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::halfedge_descriptor         halfedge_descriptor;
typedef boost::graph_traits<Polyhedron>::face_descriptor             face_descriptor;
typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron>        Skeletonization;
typedef Skeletonization::Skeleton                                    Skeleton;
typedef Skeleton::vertex_descriptor                                  Skeleton_vertex;
// Property map associating a facet with an integer as id to an
// element in a vector stored internally
using namespace std;

template<class ValueType>
struct Facet_with_id_pmap
    : public boost::put_get_helper<ValueType&,
             Facet_with_id_pmap<ValueType> >
{
    typedef face_descriptor key_type;
    typedef ValueType value_type;
    typedef value_type& reference;
    typedef boost::lvalue_property_map_tag category;
    Facet_with_id_pmap(
      std::vector<ValueType>& internal_vector
    ) : internal_vector(internal_vector) { }
    reference operator[](key_type key) const
    { return internal_vector[key->id()]; }
private:
    std::vector<ValueType>& internal_vector;
};
int main()
{
	int shmid;
	char * shm_ptr = NULL;
	std::string dummyoutputname = "/tmp/nonsense";
	char *decomp_buffer =  new char[BUFFER_SIZE];

	//use the same key to locate the segment.
	// Getting access to shared memory with all compressed objects stored in there
	if ((shmid = shmget(5678, 0, 0666)) < 0) {
		perror("shmget");
		exit(1);
	}
	// Now we attach the segment to our data space.
	if ((shm_ptr = (char *) shmat(shmid, (const void *)NULL, 0)) == (char *) -1) {
		perror("shmat");
		exit(1);
	}
	int i_mode = DECOMPRESSION_MODE_ID; // compression mode


	// Codec features status.
	bool b_useAdaptiveQuantization = false;
	bool b_useLiftingScheme = true;
	bool b_useCurvaturePrediction = true;
	bool b_useConnectivityPredictionFaces = true;
	bool b_useConnectivityPredictionEdges = true;
	bool b_allowConcaveFaces = true;
	bool b_useTriangleMeshConnectivityPredictionFaces = true;
	unsigned i_quantBit = 12;
	//unsigned i_decompPercentage = 100;

	// Init the random number generator.
	srand(4212);
	MyMesh *currentMesh = new MyMesh(NULL,// dummyoutputname,
				100,
					 i_mode, i_quantBit, b_useAdaptiveQuantization,
					 b_useLiftingScheme, b_useCurvaturePrediction,
					 b_useConnectivityPredictionFaces, b_useConnectivityPredictionEdges,
					 b_allowConcaveFaces, b_useTriangleMeshConnectivityPredictionFaces,
				dummyoutputname,
				(char*)(shm_ptr), 16759, decomp_buffer);


	assert(currentMesh!=NULL);
	currentMesh->completeOperation();
	Polyhedron geom;
	std::stringstream os;
	os << *currentMesh;
	os >>geom;
	cout << os.str();
	if (!CGAL::is_triangle_mesh(geom))
	{
		std::cout << "Input geometry is not triangulated." << std::endl;
		return EXIT_FAILURE;
	}
	std::cout<<"start extraction"<<std::endl;
	// extract the skeleton
	Skeleton skeleton;
	try{
		CGAL::extract_mean_curvature_flow_skeleton(geom, skeleton);
	}catch(std::exception &exc){
		std::cerr<<exc.what()<<std::endl;
	}
	std::cout<<"extracted"<<std::endl;

	// init the polyhedron simplex indices
	CGAL::set_halfedgeds_items_id(geom);
	//for each input vertex compute its distance to the skeleton
	std::vector<double> distances(num_vertices(geom));
	BOOST_FOREACH(Skeleton_vertex v, vertices(skeleton) )
	{
		const CGAL_Point& skel_pt = skeleton[v].point;
		BOOST_FOREACH(vertex_descriptor mesh_v, skeleton[v].vertices)
		{
				const CGAL_Point& mesh_pt = mesh_v->point();
					distances[mesh_v->id()] = std::sqrt(CGAL::squared_distance(skel_pt, mesh_pt));
		}
	}
	// create a property-map for sdf values
	std::vector<double> sdf_values( num_faces(geom) );
	Facet_with_id_pmap<double> sdf_property_map(sdf_values);
	// compute sdf values with skeleton
	BOOST_FOREACH(face_descriptor f, faces(geom))
	{
		double dist = 0;
		BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(f, geom), geom))
			dist+=distances[target(hd, geom)->id()];
		sdf_property_map[f] = dist / 3.;
	}
	// post-process the sdf values
	CGAL::sdf_values_postprocessing(geom, sdf_property_map);
	// create a property-map for segment-ids (it is an adaptor for this case)
	std::vector<std::size_t> segment_ids( num_faces(geom) );
	Facet_with_id_pmap<std::size_t> segment_property_map(segment_ids);
	// segment the mesh using default parameters
	std::cout << "Number of segments: "
		<< CGAL::segmentation_from_sdf_values(geom, sdf_property_map, segment_property_map) <<"\n";
  return EXIT_SUCCESS;
}
