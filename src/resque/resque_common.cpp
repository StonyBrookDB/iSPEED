/*
 * extract.cpp
 * with functions and variables which may be shared by different modules
 *
 *  Created on: Sep 24, 2019
 *      Author: teng
 */



#include "resque_3d.hpp"
using namespace std;
using namespace CGAL;



/* Performance metrics */
clock_t start_reading_data;
clock_t start_query_exec;

clock_t total_reading;
clock_t total_query_exec;



char * shm_ptr = NULL;
//long maxoffset = 0;

void attach_shm(struct query_op &stop){

	std::cerr<<"attaching to shared memory with id "<<COMPRESSION_KEY<<std::endl;
	// each instance attach to the shared memory for
	// the compressed objects
	int shmid;
	//use the same key to locate the segment.
	size_t maxoffset2 = stop.shm_max_size;
	// Getting access to shared memory with all compressed objects stored in there
	if ((shmid = shmget(COMPRESSION_KEY, maxoffset2, 0666)) < 0) {
		perror("shmget");
		exit(1);
	}
	// Now we attach the segment to our data space.
	if ((shm_ptr = (char *) shmat(shmid, (const void *)NULL, 0)) == (char *) -1) {
		perror("shmat");
		exit(1);
	}
	std::cerr<<"shared memory attached"<<std::endl;
}

/* Release objects in memory (for the current tile/bucket) */
void release_mem(struct query_op &stop, struct query_temp &sttemp, int maxCard) {
	if (stop.join_cardinality <= 0) {
		return ;
	}
	for (int j = 0; j < stop.join_cardinality && j < maxCard; j++ ) {
		int delete_index = j + 1; // index are adjusted to start from 1
		int len = sttemp.mbbdata[delete_index].size();
		for (int i = 0; i < len ; i++) {
				//delete sttemp.polydata[delete_index][i];
			delete sttemp.mbbdata[delete_index][i]; // release mbb
			//delete sstemp.offset[delete_index][i];
			//sttemp.offsetdata[delete_index][i].clear();
			//sttemp.lengthdata[delete_index][i].clear();
		}
    		//sttemp.polydata[delete_index].clear();
		sttemp.offsetdata[delete_index].clear();
    	sttemp.lengthdata[delete_index].clear();
		sttemp.mbbdata[delete_index].clear();
  	}
}

/* Create an R-tree index on a given set of polygons */
bool build_index_geoms(std::vector<struct mbb_3d *> & geom_mbbs, SpatialIndex::ISpatialIndex* & spidx,
		SpatialIndex::IStorageManager* & storage) {
	// build spatial index on tile boundaries
	SpatialIndex::id_type  indexIdentifier;
	CustomDataStream stream(&geom_mbbs);
	storage = SpatialIndex::StorageManager::createNewMemoryStorageManager();
	spidx   = SpatialIndex::RTree::createAndBulkLoadNewRTree(SpatialIndex::RTree::BLM_STR, stream, *storage,
			FillFactor,
			IndexCapacity,
			LeafCapacity,
			3,
			SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
	// Error checking
	return spidx->isIndexValid();
}



MyMesh *extract_mesh(long offset, long length, unsigned i_decompPercentage){

#ifdef DEBUG
	std::cerr << "extracting mesh with offset and length " << offset << TAB << length << std::endl;
#endif

	// Initialize parameters
	// Codec features status.
	int i_mode = DECOMPRESSION_MODE_ID; // compression mode
	bool b_useAdaptiveQuantization = false;
	bool b_useLiftingScheme = true;
	bool b_useCurvaturePrediction = true;
	bool b_useConnectivityPredictionFaces = true;
	bool b_useConnectivityPredictionEdges = true;
	bool b_allowConcaveFaces = true;
	bool b_useTriangleMeshConnectivityPredictionFaces = true;
	unsigned i_quantBit = 12;

	// Init the random number generator.
	srand(4212);
	// decompress the polyhedron from binary
	MyMesh *currentMesh = new MyMesh(
					 i_decompPercentage,
					 i_mode, i_quantBit, b_useAdaptiveQuantization,
					 b_useLiftingScheme, b_useCurvaturePrediction,
					 b_useConnectivityPredictionFaces, b_useConnectivityPredictionEdges,
					 b_allowConcaveFaces, b_useTriangleMeshConnectivityPredictionFaces,
					 (char*)(shm_ptr + offset), length);


	assert(currentMesh!=NULL);
	currentMesh->completeOperation();
#ifdef DEBUG
	std::cerr << "mesh decompressed" << std::endl;
#endif
	assert(currentMesh->is_closed()&&"polyhedron should be closed");
	//assert(currentMesh->is_pure_triangle()&&"polyhedron mesh should be pure triangle");
//	if(!currentMesh->is_pure_triangle()){
//		std::cout<<*currentMesh;
//		exit(1);
//	}

	return currentMesh;

}

/*
  Extract one polyhedron geometry from compressed data with given offset and length
*/
Polyhedron *extract_geometry(long offset, long length, unsigned i_decompPercentage) {
	Polyhedron *geom = new Polyhedron();
	MyMesh *currentMesh = extract_mesh(offset, length, i_decompPercentage);
	std::stringstream os;
	os << *currentMesh;
	os >> *geom;
	delete currentMesh;
	os.clear();
#ifdef DEBUG
	cerr << "geometry is extracted successfully" << endl;
#endif
	return geom;
}


/*
 * the kernel of the polyhedron extracted is Simple_Cartisian
 * */
Sc_Polyhedron *sc_extract_geometry(long offset, long length, unsigned i_decompPercentage) {
	Sc_Polyhedron *geom = new Sc_Polyhedron();
	MyMesh *currentMesh = extract_mesh(offset, length, i_decompPercentage);
	std::stringstream os;
	os << *currentMesh;
	os >> *geom;
	delete currentMesh;
	os.clear();
#ifdef DEBUG
	cerr << "geometry is extracted successfully" << endl;
#endif
	return geom;
}

Sc_Polyhedron *sc_extract_geometry_from_file(const char *path){
	assert(path!=NULL && "path cannot be NULL");
#ifdef DEBUG
	cerr<<"reading polyhedron from "<<path<<endl;
#endif
	std::ifstream input(path);
	Sc_Polyhedron *poly = new Sc_Polyhedron();
	input >> *poly;
	input.close();
	assert(CGAL::is_triangle_mesh(*poly) && "polyhedron should be in triangle mesh");
#ifdef DEBUG
	cerr << "geometry is extracted successfully" << endl;
#endif
	return poly;
}

void extract_skeleton(long offset, long length, unsigned i_decompPercentage, std::vector<Sc_Point> &P){
#ifdef DEBUG
	cerr << "extracting the Skeleton!" << endl;
#endif
	Sc_Skeleton skeleton;
	// todo this is still something wrong here while loading large polyhedron
	// thus we do not use this method, use the extract_skeleton_advance instead
	// de-compressed by the pmcc, fix it in the future
	Sc_Polyhedron *geom = sc_extract_geometry(offset, length, i_decompPercentage);
	//Sc_Polyhedron geom = sc_extract_geometry_from_file("/home/teng/gisdata/processed/good.off");

	CGAL::extract_mean_curvature_flow_skeleton(*geom, skeleton);

	BOOST_FOREACH(Sc_Skeleton_vertex v, vertices(skeleton)){
		Sc_Point p = skeleton[v].point;
		P.push_back(p);
	}
	delete geom;
#ifdef DEBUG
	cerr << "extracted one Skeleton!" << endl;
#endif
}


void extract_skeleton_advance(long offset, long length, unsigned i_decompPercentage, std::vector<Sc_Point> &P){
	// some local definition
	typedef CGAL::Surface_mesh<Point>                             Triangle_mesh;
	typedef boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;
	typedef CGAL::Mean_curvature_flow_skeletonization<Triangle_mesh> Skeletonization;
	typedef Skeletonization::Skeleton                             Skeleton;
	typedef Skeleton::vertex_descriptor                           Skeleton_vertex;
	typedef Skeleton::edge_descriptor                             Skeleton_edge;

#ifdef DEBUG
	cerr << "extracting the Skeleton (advanced)!" << endl;
#endif
	MyMesh *currentMesh = extract_mesh(offset, length, i_decompPercentage);
	std::stringstream os;
	os << *currentMesh;
	Triangle_mesh tmesh;
	os >> tmesh;
	if (!CGAL::is_triangle_mesh(tmesh)){
		std::cerr << "Input geometry is not triangulated." << std::endl;
		exit(-1);
	}
	try{
		Skeleton skeleton;
		Skeletonization mcs(tmesh);
		// 1. Contract the mesh by mean curvature flow.
		mcs.contract_geometry();
		// 2. Collapse short edges and split bad triangles.
		mcs.collapse_edges();
		mcs.split_faces();
		// 3. Fix degenerate vertices.
		mcs.detect_degeneracies();
		// Perform the above three steps in one iteration.
		mcs.contract();
		// Iteratively apply step 1 to 3 until convergence.
		mcs.contract_until_convergence();
		// Convert the contracted mesh into a curve skeleton and
		// get the correspondent surface points
		mcs.convert_to_skeleton(skeleton);

		BOOST_FOREACH(Skeleton_vertex v, boost::vertices(skeleton)){
			auto p = skeleton[v].point;
			P.push_back(Sc_Point(p.x(),p.y(),p.z()));
		}
	}catch(std::exception &exc){
		std::cerr<<exc.what()<<std::endl;
		exit(-1);
	}
	delete currentMesh;
#ifdef DEBUG
	cerr << "extracted one Skeleton!" << endl;
#endif
}


void get_triangle(Polyhedron *P, std::vector<Triangle>& triangles,  std::vector<Box>& boxes, std::vector<Box*>& ptr){
	for ( Facet_const_iterator i = P->facets_begin(); i != P->facets_end(); ++i){
		Triangle t(i->halfedge()->vertex()->point(),
				   i->halfedge()->next()->vertex()->point(),
				   i->halfedge()->next()->next()->vertex()->point());
		triangles.push_back(t);
	}
	// Create the corresponding std::vector of bounding boxes
	for ( Iterator i = triangles.begin(); i != triangles.end(); ++i){
		boxes.push_back(Box(i->bbox(), i));
	}

	for ( std::vector<Box>::iterator i = boxes.begin(); i != boxes.end(); ++i){
		ptr.push_back( &*i);
	}
}


bool intersects_mbb(const struct mbb_3d * m1, const struct mbb_3d *m2) {
	return !(m1->low[0] > m2->high[0] || m1->high[0] < m2->low[0]
	      || m1->low[1] > m2->high[1] || m1->high[1] < m2->low[1]
	      || m1->low[2] > m2->high[2] || m1->high[2] < m2->low[2] );
}

// testing intersecting between polyhedrons
bool intersection_flag = false;
//struct Report {
//  Triangles* triangles;
//  Triangles* cell_triangles;
//
//  Report(Triangles& triangles, Triangles& cell_triangles)
//    : triangles(&triangles), cell_triangles(&cell_triangles)
//  {}
//
//  // callback functor that reports all truly intersecting triangles
//  void operator()(const Box* a, const Box* b) const
//  {
//    if (intersection_flag) {
//    	return;
//    }
//    if ( ! a->handle()->is_degenerate() && ! b->handle()->is_degenerate()
//         && CGAL::do_intersect( *(a->handle()), *(b->handle()))) {
//      intersection_flag = true;
//    }
//  }
//};
struct Report {

  Report(){}

  // callback functor that reports all truly intersecting triangles
  void operator()(const Box* a, const Box* b) const
  {
    if (intersection_flag) {
    	return;
    }
    if ( ! a->handle()->is_degenerate() && ! b->handle()->is_degenerate()
         && CGAL::do_intersect( *(a->handle()), *(b->handle()))) {
      intersection_flag = true;
    }
  }
};
bool intersects(Polyhedron *P1, Polyhedron *P2) {

	// Use Nef_polyhedron for intersection detection
	Triangles triangles1, triangles2;
	std::vector<Box> boxes1, boxes2;
	std::vector<Box*> boxes1_ptr, boxes2_ptr;

	get_triangle(P1, triangles1, boxes1, boxes1_ptr);
	get_triangle(P2, triangles2, boxes2, boxes2_ptr);

	intersection_flag = false;
	CGAL::box_intersection_d( boxes1_ptr.begin(), boxes1_ptr.end(),
							  boxes2_ptr.begin(), boxes2_ptr.end(),
							  Report());
	return intersection_flag;
}

double get_volume(Nef_polyhedron &inputpoly) {
	// to check if the intersected object can be converted to polyhedron or not
	std::vector<Polyhedron> PList;
	if(inputpoly.is_simple()) {
		Polyhedron P;
		inputpoly.convert_to_polyhedron(P);
		PList.push_back(P);
	} else {
		// decompose non-convex volume to convex parts
		convex_decomposition_3(inputpoly);
		Volume_const_iterator ci = ++inputpoly.volumes_begin();
		for( ; ci != inputpoly.volumes_end(); ++ci) {
			if(ci->mark()) {
				Polyhedron P;
				inputpoly.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);
				PList.push_back(P);
			}
		}
	}
	// std::cerr<< "# of Polyhedrons: " << PList.size() <<std::endl;
	// triangulate the polyhedrons to generate mesh and use tetrahedron to calculate volumes
	Polyhedron poly;
	double total_volume = 0, hull_volume = 0;
	for(int i = 0; i < PList.size(); i++) {
		poly = PList[i];
		std::vector<CGAL_Point> L;
		for (Polyhedron::Vertex_const_iterator  it = poly.vertices_begin(); it != poly.vertices_end(); it++) {
			L.push_back(CGAL_Point(it->point().x(), it->point().y(), it->point().z()));
		}
		Triangulation T(L.begin(), L.end());
		hull_volume = 0;
		for(Triangulation::Finite_cells_iterator it = T.finite_cells_begin(); it != T.finite_cells_end(); it++) {
			Tetrahedron tetr = T.tetrahedron(it);
			hull_volume += to_double(tetr.volume());
		}
		total_volume += hull_volume;
	}
	return total_volume;
}



