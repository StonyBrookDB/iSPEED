/*
 * spjoin.cpp
 *
 *  Created on: Sep 15, 2019
 *      Author: teng
 */

#include <resque/resque_3d.hpp>

using namespace std;

// performs spatial join on the current tile (bucket)
int join_bucket_spjoin(struct query_op &stop, struct query_temp &sttemp) {
	SpatialIndex::IStorageManager *storage = NULL;
	SpatialIndex::ISpatialIndex *spidx = NULL;
	/* Indicates where original data is mapped to */
	int idx1 = SID_1;
	int idx2 = SID_2;

	int pairs = 0; // number of satisfied results

	double low[3], high[3];  // Temporary value placeholders for MBB

	try {
		#ifdef DEBUG
		std::cerr<<"start building r-tree index"<<std::endl;
		time_t rtree_st, rtree_et;
		double rtree_tt;
		time(&rtree_st);
		#endif
		/* Handling for special nearest neighbor query */
		// build the actual spatial index for input polygons from idx2
		if (! build_index_geoms(sttemp.mbbdata[idx2], spidx, storage)) {
			#ifdef DEBUG
			std::cerr << "Building index on geometries from set 2 has failed" << std::endl;
			#endif
			return -1;
		}

		#ifdef DEBUG
		time(&rtree_et);
		rtree_tt = difftime(rtree_et,rtree_st);
		std::cerr << "********************************************" << std::endl;
		std::cerr << "R-tree construction total execution time: "
			<< rtree_tt
			<< " seconds." << std::endl;
		std::cerr << "********************************************" << std::endl;
		#endif

		/*#ifdef DEBUG
		time_t mbb_st, mbb_et;
		double mbb_tt;
		time(&mbb_st);
		#endif*/

		std::vector<struct mbb_3d *> geom_mbb1 = sttemp.mbbdata[idx1];
		for (int i = 0; i < geom_mbb1.size(); i++) {
			/* Extract minimum bounding box */
			struct mbb_3d * env1 = geom_mbb1[i];
			low[0] = env1->low[0];
			low[1] = env1->low[1];
			low[2] = env1->low[2];
			high[0] = env1->high[0];
			high[1] = env1->high[1];
			high[2] = env1->high[2];

			if (stop.join_predicate == ST_DWITHIN) {
				low[0] -= stop.expansion_distance;
				low[1] -= stop.expansion_distance;
				low[2] -= stop.expansion_distance;
				high[0] += stop.expansion_distance;
				high[1] += stop.expansion_distance;
				high[2] += stop.expansion_distance;
			}
			//std::cerr << " done Got MBB1!" << std::endl;
			/* Regular handling */
			SpatialIndex::Region r(low, high, 3);
			MyVisitor vis;
			vis.matches.clear();
			/* R-tree intersection check */
			spidx->intersectsWithQuery(r, vis);
			if(vis.matches.size()==0){
				continue;
			}

			// checking true intersection
			Polyhedron *geom1 = extract_geometry(sttemp.offsetdata[idx1][i],
					sttemp.lengthdata[idx1][i], stop.decomp_lod);
			for (uint32_t j = 0; j < vis.matches.size(); j++){
				Polyhedron *geom2 = extract_geometry(sttemp.offsetdata[idx2][vis.matches[j]],
					sttemp.lengthdata[idx2][vis.matches[j]], stop.decomp_lod);
				struct mbb_3d * env2 = sttemp.mbbdata[idx2][vis.matches[j]];
				pairs += join_with_predicate(stop, geom1, geom2);
				delete geom2;
			}
			delete geom1;
		}
	} catch (Tools::Exception& e) {
		std::cerr << "******ERROR******" << std::endl;
		std::cerr << e.what() << std::endl;
		return -1;
	} // end of catch

	delete spidx;
	delete storage;
	return pairs ;
}


/* Perform (Refine) spatial computation between 2 geometry objects */
bool join_with_predicate(struct query_op &stop,Polyhedron *geom1 , Polyhedron *geom2){
	bool flag = false;
	// we currently support only intersects
	// which can be easily extended to other operations
	// like touch, cross, contains, within etc.
	switch (stop.join_predicate){
		case ST_INTERSECTS:
			flag = intersects(geom1, geom2);
			break;
		default:
			std::cerr << "ERROR: unknown spatial predicate " << std::endl;
			break;
	}
	/* Spatial computation is only performed once for a result pair */
	if (flag) {
		Nef_polyhedron *N1 = NULL;
		Nef_polyhedron *N2 = NULL;

		if (stop.needs_volume_1) {
			N1 = new Nef_polyhedron(*geom1);
			cerr<<get_volume(*N1)<<endl;
		}

		if (stop.needs_volume_2) {
			N2 = new Nef_polyhedron(*geom2);
			cerr<<get_volume(*N2)<<endl;
		}

		// Slow because of this
		if (stop.needs_intersect_volume) {
			if(N1==NULL){
				N1 = new Nef_polyhedron(*geom1);
			}
			if(N2==NULL){
				N2 = new Nef_polyhedron(*geom2);
			}
			Nef_polyhedron Inter = (*N1) * (*N2);
			assert(Inter.number_of_vertices()>0
					&&"Those two polygedrons should be intersected");
			cerr<<get_volume(Inter)<<endl;
			Inter.clear();
			delete N1;
			delete N2;
		}
	}
	return flag;
}



