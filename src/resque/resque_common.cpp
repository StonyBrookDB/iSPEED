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
		sttemp.rawdata[delete_index].clear();
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
	// Initialize parameters
	int i_mode = DECOMPRESSION_MODE_ID; // compression mode

#ifdef DEBUG
	std::cerr << "attempting to extract " << offset << TAB << length << std::endl;
#endif

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
				i_decompPercentage,
					 i_mode, i_quantBit, b_useAdaptiveQuantization,
					 b_useLiftingScheme, b_useCurvaturePrediction,
					 b_useConnectivityPredictionFaces, b_useConnectivityPredictionEdges,
					 b_allowConcaveFaces, b_useTriangleMeshConnectivityPredictionFaces,
				dummyoutputname,
				(char*)(shm_ptr + offset), length, resque_decomp_buffer);


	assert(currentMesh!=NULL);
	currentMesh->completeOperation();
#ifdef DEBUG
	std::cerr << "mesh extracted" << std::endl;
#endif
	return currentMesh;

}

/*
  Extract one polyhedron geometry from compressed data with given offset and length
*/
Polyhedron extract_geometry(long offset, long length, unsigned i_decompPercentage,
	struct query_op &stop, struct query_temp &sttemp, int dataset_id) {
	Polyhedron geom;
	MyMesh *currentMesh = extract_mesh(offset, length, i_decompPercentage);
	std::stringstream os;
	os << *currentMesh;
	#ifdef DEBUG
	std::cerr << "done decomp" << std::endl;
	#endif
	os >> geom;
	// only when volume is needed
	if (stop.needs_intersect_volume) {
		sttemp.poly_str[dataset_id].str(os.str());
	}

	delete currentMesh;
	return geom;
}


/*
 * the kernel of the polyhedron extracted is Simple_Cartisian
 * */
Sc_Polyhedron sc_extract_geometry(long offset, long length, unsigned i_decompPercentage,
	struct query_op &stop, struct query_temp &sttemp, int dataset_id) {
	Sc_Polyhedron geom;
	MyMesh *currentMesh = extract_mesh(offset, length, i_decompPercentage);
	std::stringstream os;
	os << *currentMesh;
	os >> geom;
	// only when volume is needed
	if (stop.needs_intersect_volume) {
		sttemp.poly_str[dataset_id].str(os.str());
	}
	delete currentMesh;
#ifdef DEBUG
	cerr << "geometry is extracted successfully" << endl;
#endif
	os.clear();
	return geom;
}

Sc_Polyhedron sc_extract_geometry_from_file(const char *path){
	assert(path!=NULL && "path cannot be NULL");
	std::ifstream input(path);
	Sc_Polyhedron poly;
	input >> poly;
	input.close();
	return poly;
}

Sc_Skeleton extract_skeleton(long offset, long length, unsigned i_decompPercentage){
	Sc_Skeleton skeleton;

//	MyMesh *currentMesh = extract_mesh(offset, length, i_decompPercentage);
//	std::stringstream os;
//	os << *currentMesh;
//	Sc_Triangle_mesh tmesh;
//	os >> tmesh;
//	if (!CGAL::is_triangle_mesh(tmesh)){
//		std::cerr << "Input geometry is not triangulated." << std::endl;
//		exit(-1);
//	}
//	try{
//		Sc_Skeletonization mcs(tmesh);
//		// 1. Contract the mesh by mean curvature flow.
//		mcs.contract_geometry();
//		// 2. Collapse short edges and split bad triangles.
//		mcs.collapse_edges();
//		mcs.split_faces();
//		// 3. Fix degenerate vertices.
//		mcs.detect_degeneracies();
//		// Perform the above three steps in one iteration.
//		mcs.contract();
//		// Iteratively apply step 1 to 3 until convergence.
//		mcs.contract_until_convergence();
//		// Convert the contracted mesh into a curve skeleton and
//		// get the correspondent surface points
//		mcs.convert_to_skeleton(skeleton);
//	}catch(std::exception &exc){
//		std::cerr<<exc.what()<<std::endl;
//		exit(-1);
//	}
//	delete currentMesh;

	return skeleton;
}


void get_triangle(Polyhedron P, std::vector<Triangle>& triangles,  std::vector<Box>& boxes, std::vector<Box*>& ptr){

	 // std::vector<Box> boxes;
	 for ( Facet_const_iterator i = P.facets_begin(); i != P.facets_end(); ++i){
		triangles.push_back(
		    Triangle( i->halfedge()->vertex()->point(),
			i->halfedge()->next()->vertex()->point(),
			i->halfedge()->next()->next()->vertex()->point()));

	  }
	 // Create the corresponding std::vector of bounding boxes
	  for ( Iterator i = triangles.begin(); i != triangles.end(); ++i)
	    boxes.push_back(Box( i->bbox(), i));

	  for ( std::vector<Box>::iterator i = boxes.begin(); i != boxes.end(); ++i)
	    ptr.push_back( &*i);
}


bool intersects_mbb(const struct mbb_3d * m1, const struct mbb_3d *m2) {
	return !(m1->low[0] > m2->high[0] || m1->high[0] < m2->low[0]
	      || m1->low[1] > m2->high[1] || m1->high[1] < m2->low[1]
	      || m1->low[2] > m2->high[2] || m1->high[2] < m2->low[2] );
}


bool intersects(Polyhedron &P1, Polyhedron &P2, const struct mbb_3d * env1, const struct mbb_3d * env2) {
	//return true;
	// Use Nef_polyhedron for intersection detection
	if(!intersects_mbb(env1, env2))
		return false;
	else{
		Triangles triangles1, triangles2;
		std::vector<Box> boxes1, boxes2;
		std::vector<Box*> boxes1_ptr, boxes2_ptr;

		get_triangle(P1, triangles1, boxes1, boxes1_ptr);
		get_triangle(P2, triangles2, boxes2, boxes2_ptr);

		intersection_flag = false;

		CGAL::box_intersection_d( boxes1_ptr.begin(), boxes1_ptr.end(), boxes2_ptr.begin(), boxes2_ptr.end(), Report(triangles1, triangles2));

		/*std::cout << "yes or no: " << intersection_flag << std::endl; // to compare * operator and yes or no question

		// use * operator for intersection detection
		if(P1->is_closed() && P2->is_closed()) {
			//std::cerr << "got here 1" << std::endl;
			Nef_polyhedron N1(*P1);
			Nef_polyhedron N2(*P2);
			//std::cerr << "got here 2" << std::endl;
		//	return true;
			Nef_polyhedron inputpoly = N1 * N2;

			bool star_flag = false;
			if(inputpoly.number_of_vertices() > 0) { star_flag = true; }
			else { star_flag = false; }

			std::cout << "*: " << star_flag << std::endl;
		}
		else
			std::cerr << "ERROR: Polyhedron is not closed!" << std::endl;*/


		/*if(inputpoly.number_of_vertices() > 0) { return true; }
		else { return false; }
		}
		else
			std::cerr << "ERROR: Polyhedron is not closed!" << std::endl;*/

		return intersection_flag;
	}

	return false;
}

double get_volume(Nef_polyhedron &inputpoly) {
	// to check if the intersected object can be converted to polyhedron or not
	std::vector<Polyhedron> PList;
	if(inputpoly.is_simple()) {
		Polyhedron P;
		inputpoly.convert_to_polyhedron(P);
		PList.push_back(P);
	}
	else {
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
	//std::cout<< "# of Polyhedrons: " << PList.size() <<std::endl;
	// triangulate the polyhedrons to generate mesh and use terahedron to calculate volumes
	Polyhedron poly;
	double total_volume = 0, hull_volume = 0;
	for(int i = 0; i < PList.size(); i++) {
		poly = PList[i];
		std::vector<CGAL_Point> L;
		for (Polyhedron::Vertex_const_iterator  it = poly.vertices_begin(); it != poly.vertices_end(); it++) {
			CGAL_Point p(it->point().x(), it->point().y(), it->point().z());
			L.push_back(p);
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



