/*
 * spjoin.cpp
 *
 *  Created on: Sep 15, 2019
 *      Author: teng
 */

#include <resque/resque_3d.hpp>

using namespace std;

/*
  Extract one polyhedron geometry from compressed data with given offset and length
*/
Polyhedron* extract_geometry(long offset, long length, unsigned i_decompPercentage,
	struct query_op &stop, struct query_temp &sttemp, int dataset_id) {
	Polyhedron* geom;

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
				// fbuffer, length, resque_decomp_buffer);
		            // b_allowConcaveFaces, b_useTriangleMeshConnectivityPredictionFaces, NULL);
	assert(currentMesh!=NULL);
	currentMesh->completeOperation();

	//std::cerr << "current mesh: " << *currentMesh << std::endl;
	std::stringstream os;
	os << *currentMesh;
	//os.clear();
	#ifdef DEBUG
	std::cerr << "done decomp" << std::endl;
	#endif
	geom = new Polyhedron();
	os >> *geom;
	//std::cerr << "os: " << os.str() << std::endl;

	// only when volume is needed
	//if (stop.needs_intersect_volume) {
		sttemp.poly_str[dataset_id].str(os.str());
	//}

	//delete[] fbuffer;
	//std::cerr << "constructing poly" << std::endl;
	//std::cerr << "geom: " << *geom << std::endl;
	delete currentMesh;
	return geom;
}




// performs spatial join on the current tile (bucket)
int join_bucket_spjoin(struct query_op &stop, struct query_temp &sttemp) {
	SpatialIndex::IStorageManager *storage = NULL;
	SpatialIndex::ISpatialIndex *spidx = NULL;
	bool selfjoin = stop.join_cardinality == 1  ? true : false;
	/* Indicates where original data is mapped to */
	int idx1 = SID_1;
	int idx2 = selfjoin ? SID_1 : SID_2;

	int pairs = 0; // number of satisfied results

	double low[3], high[3];  // Temporary value placeholders for MBB


	try {
		//std::cerr << "shm add: " << (long) shm << std::endl;

		//std::cerr << "idx1: " << idx1 << std::endl;
		//std::cerr << "idx2: " << idx2 << std::endl;
		/* Build index on the "second data set */
		//std::vector<struct mbb_3d *> geom_mbb2 = sttemp.mbbdata[idx2];
		//Polyhedron* geom3 = extract_geometry(sttemp.offsetdata[SID_1][0], sttemp.lengthdata[SID_1][0]);
		//delete geom3;

		int len1 = sttemp.mbbdata[idx1].size();
		int len2 = sttemp.mbbdata[idx2].size();

	 	if (len1 <= 0 || len2 <= 0) {
			 return 0;
		}

		#ifdef DEBUG
		std::cerr<<"\nstart one round of joining"<<std::endl;
		std::cerr << "Length of data"<<idx1<<": " << len1 << std::endl;
		std::cerr << "Length of data"<<idx2<<": " << len2 << std::endl;
		#endif

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

			// Extract geometry from compressed data in the shared memeory segment
			//std::cerr << "obj: " << i << std::endl;
			//std::cerr << "offset: " << sttemp.offsetdata[idx1][i] << std::endl;
			//std::cerr << "length: " << sttemp.lengthdata[idx1][i] << std::endl;
			//Polyhedron* geom1 = extract_geometry(sttemp.offsetdata[idx1][i], sttemp.lengthdata[idx1][i],
			//		stop.decomp_lod, stop, sttemp, 0);

			/* Extract minimum bounding box */
			//Polyhedron* geom1 = poly_set_one[i];
			//std::cerr << "before Got MBB1!" << std::endl;
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
			#ifdef DEBUG
			std::cerr<<vis.matches.size()<<" objects in dataset 2 is matched"<<std::endl;
			#endif

			/*#ifdef DEBUG
			time(&mbb_et);
			rtree_tt = difftime(mbb_et,mbb_st);
			std::cerr << "********************************************" << std::endl;
			std::cerr << "One MBB filtering total execution time: "
				<< mbb_tt
				<< " seconds." << std::endl;
			std::cerr << "********************************************" << std::endl;
			#endif*/

			// This is where iSPEED difference from Hadoopgis starts:

			Polyhedron* geom1 = extract_geometry(sttemp.offsetdata[idx1][i],
					sttemp.lengthdata[idx1][i], stop.decomp_lod, stop, sttemp, 0);

			for (uint32_t j = 0; j < vis.matches.size(); j++){

				/* Skip results that have been seen before (self-join case) */
				/*
				if (selfjoin && ((vis.matches[j] == i) ||  // same objects in self-join
				    (!stop.result_pair_duplicate && vis.matches[j] <= i))) { // duplicate pairs
					#ifdef DEBUG
					std::cerr << "skipping (selfjoin): " << j << " " << vis.matches[j] << std::endl;
					#endif
					continue;
				}
				*/

				Polyhedron* geom2 = extract_geometry(sttemp.offsetdata[idx2][vis.matches[j]],
					sttemp.lengthdata[idx2][vis.matches[j]], stop.decomp_lod, stop, sttemp, 1);
				struct mbb_3d * env2 = sttemp.mbbdata[idx2][vis.matches[j]];

				// for now only decomp time
				//std::cout << i << TAB << vis.matches[j] << std::endl;
				#ifdef DEBUG
				std::cerr << "Checking actual intersection between " << i << TAB << vis.matches[j] << std::endl;
				time_t geometry_st, geometry_et;
				double geometry_tt;
				time(&geometry_st);
				#endif

				if (join_with_predicate(stop, sttemp, geom1, geom2, env1, env2,
							stop.join_predicate))  {
					//	report_result(stop, sttemp, i, vis.matches[j]);
					#ifdef DEBUG
					std::cerr << "Actual intersected with each other " << i << TAB << vis.matches[j] << std::endl;
					#endif
					pairs++;
				}else{
					#ifdef DEBUG
					std::cerr << "not intersected with each other " << i << TAB << vis.matches[j] << std::endl;
					#endif
				}

				#ifdef DEBUG
				time(&geometry_et);
				geometry_tt = difftime(geometry_et,geometry_st);
				std::cerr << "********************************************" << std::endl;
				std::cerr << "One geometry spatial refinement total execution time: "
					<< geometry_tt
					<< " seconds." << std::endl;
				std::cerr << "********************************************" << std::endl;
				#endif
				delete geom2;

			}
			delete geom1;
		}
	} catch (Tools::Exception& e) {
		std::cerr << "******ERROR******" << std::endl;
		#ifdef DEBUG
		std::cerr << e.what() << std::endl;
		#endif
		return -1;
	} // end of catch

	delete spidx;
	delete storage;
	return pairs ;
}


/* Perform (Refine) spatial computation between 2 geometry objects */
bool join_with_predicate(
		struct query_op &stop, struct query_temp &sttemp,
		Polyhedron * geom1 , Polyhedron * geom2,
		const struct mbb_3d * env1, const struct mbb_3d * env2,
		const int jp){
	bool flag = false; // flag == true means the predicate is satisfied

	//BufferOp * buffer_op1;
	//BufferOp * buffer_op2;
	Polyhedron* geom_buffer1;
	Polyhedron* geom_buffer2;
	Polyhedron* geomUni;
	Polyhedron* geomIntersect;


	#ifdef DEBUG

	#endif

	switch (jp){
		case ST_INTERSECTS:
			flag = intersects(geom1, geom2, env1, env2);
			break;

		/*#ifdef SKIP2D // skip other spatial operations for now
		case ST_TOUCHES:
			flag = geom1->touches(geom2);
			break;

		case ST_CROSSES:
			flag = geom1->crosses(geom2);
			break;

		case ST_CONTAINS:
			flag = env1->contains(env2) && geom1->contains(geom2);
			break;

		case ST_ADJACENT:
			flag = !geom1->disjoint(geom2);
			break;

		case ST_DISJOINT:
			flag = geom1->disjoint(geom2);
			break;

		case ST_EQUALS:
			flag = env1->equals(env2) && geom1->equals(geom2);
			break;
		*/
		//case ST_DWITHIN:
			/* Special spatial handling for the point-point case */
			//if (geom1->getGeometryTypeId() == geos::geom::GEOS_POINT
			//	&& geom2->getGeometryTypeId() == geos::geom::GEOS_POINT) 				{
				/* Replace with spherical distance computation if points are on eath */
			/*	if (stop.use_earth_distance) {
					flag = get_distance_earth(
						dynamic_cast<const geos::geom::Point*>(geom1),
						dynamic_cast<const geos::geom::Point*>(geom2))
						<= stop.expansion_distance;
				} else {
					flag = DistanceOp::distance(geom1, geom2)
						<= stop.expansion_distance;
				}*/

				/* flag = distance(
					dynamic_cast<const geos::geom::Point*>(geom1),
					dynamic_cast<const geos::geom::Point*>(geom2) )
					 <= stop.expansion_distance; */
			//}
			/*else {
				// Regular handling for other object types
				buffer_op1 = new BufferOp(geom1);
				// buffer_op2 = new BufferOp(geom2);
				if (NULL == buffer_op1)
					std::cerr << "NULL: buffer_op1" <<std::endl;
				geom_buffer1 = buffer_op1->getResultGeometry(stop.expansion_distance);
				env1 = geom_buffer1->getEnvelopeInternal();
				// geom_buffer2 = buffer_op2->getResultGeometry(expansion_distance);
				//Envelope * env_temp = geom_buffer1->getEnvelopeInternal();
				if (NULL == geom_buffer1) {
					std::cerr << "NULL: geom_buffer1" << std::endl;
				}
				flag = join_with_predicate(stop, sttemp, geom_buffer1, geom2,
					env1, env2, ST_INTERSECTS);
				delete geom_buffer1;
				delete buffer_op1;
			}
			break;

		case ST_WITHIN:
			flag = geom1->within(geom2);
			break;

		case ST_OVERLAPS:
			flag = geom1->overlaps(geom2);
			break;*/
		/*
		case ST_NEAREST:
		case ST_NEAREST_2:
			// Execution only reaches here if this is already the nearest neighbor
			flag = true;
			break;
		*/
		//#endif

		default:
			std::cerr << "ERROR: unknown spatial predicate " << std::endl;
			break;
	}
	/* Spatial computation is only performed once for a result pair */
	if (flag) {

		if (stop.needs_volume_1) {
			Nef_polyhedron N1(*geom1);
			sttemp.volume1 = get_volume(N1);
		}

		if (stop.needs_volume_2) {
			Nef_polyhedron N2(*geom2);
			sttemp.volume2 = get_volume(N2);
		}

		// Slow because of this
		if (stop.needs_intersect_volume) {
			if((*geom1).is_closed() && (*geom2).is_closed()) {
		//	if(p1.is_closed() && p2.is_closed()){

				/*
				istringstream poly1str;
				poly1str << *geom1;

     				Nef_polyhedron NP;
     				CGAL::OFF_to_nef_3(poly1str, N1);

				istringstream poly2str;
				poly1str << *geom2;

     				Nef_polyhedron NP;
     				CGAL::OFF_to_nef_3(poly2str, N2);
				*/
				Nef_polyhedron N1;
				Nef_polyhedron N2;
				CGAL::OFF_to_nef_3(sttemp.poly_str[0], N1);
				CGAL::OFF_to_nef_3(sttemp.poly_str[1], N2);
				//Nef_polyhedron N1(*geom1);
				//Nef_polyhedron N2(*geom2);

				//Nef_polyhedron N1(p1);
				//Nef_polyhedron N2(p2);
				Nef_polyhedron Inter = N1 * N2;
				//Nef_polyhedron Inter = N1 * N2;

				if(Inter.number_of_vertices() > 0) {
				   	sttemp.intersect_volume = get_volume(Inter);
				}
				else {
					std::cerr << "ERROR: Polyhedrons are not intersected!" << std::endl;
 				}
				Inter.clear();
				N1.clear();
				N2.clear();

			}
			else {
				std::cerr << "ERROR: Polyhedron is not closed!" << std::endl;
			}
		}

		/*#ifdef SKIP2D // skip 2d cases
		if (stop.needs_area_1) {
			sttemp.area1 = geom1->getArea();
		}
		if (stop.needs_area_2) {
			sttemp.area2 = geom2->getArea();
		}
		if (stop.needs_union) {
			Geometry * geomUni = geom1->Union(geom2);
			sttemp.union_area = geomUni->getArea();
			delete geomUni;
		}
		if (stop.needs_intersect) {
			Geometry * geomIntersect = geom1->intersection(geom2);
			sttemp.intersect_area = geomIntersect->getArea();
			delete geomIntersect;
		}*/
		/* Statistics dependent on previously computed statistics */
		/*if (stop.needs_jaccard) {
			sttemp.jaccard = compute_jaccard(sttemp.union_area, sttemp.intersect_area);
		}

		if (stop.needs_dice) {
			sttemp.dice = compute_dice(sttemp.area1, sttemp.area2, sttemp.intersect_area);
		}

		if (stop.needs_min_distance) {
			if (stop.use_earth_distance
				&& geom1->getGeometryTypeId() == geos::geom::GEOS_POINT
				&& geom2->getGeometryTypeId() == geos::geom::GEOS_POINT) 				{
				sttemp.distance = get_distance_earth(
						dynamic_cast<const geos::geom::Point*>(geom1),
						dynamic_cast<const geos::geom::Point*>(geom2));
			}
			else {
				sttemp.distance = DistanceOp::distance(geom1, geom2);
			}
		}
		#endif*/
	}
	return flag;
}



