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



bool intersects(Polyhedron *P1, Polyhedron *P2, const struct mbb_3d * env1, const struct mbb_3d * env2) {
	//return true;
	// Use Nef_polyhedron for intersection detection
	if(!intersects(env1, env2))
		return false;
	else{
		Triangles triangles1, triangles2;
		std::vector<Box> boxes1, boxes2;
		std::vector<Box*> boxes1_ptr, boxes2_ptr;
		
		get_triangle(*P1, triangles1, boxes1, boxes1_ptr);
		get_triangle(*P2, triangles2, boxes2, boxes2_ptr);
	
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

bool intersects(const struct mbb_3d * m1, const struct mbb_3d *m2) {
	return !(m1->low[0] > m2->high[0] || m1->high[0] < m2->low[0] 
	      || m1->low[1] > m2->high[1] || m1->high[1] < m2->low[1]
	      || m1->low[2] > m2->high[2] || m1->high[2] < m2->low[2] );
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


	//#ifdef COMPRESSED
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

		//#ifdef DEBUG
		std::cerr << "Length of data1: " << len1 << std::endl;
		std::cerr << "Length of data2: " << len2<< std::endl;
		//#endif

		#ifdef DEBUGTIME
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

		#ifdef DEBUGTIME
		time(&rtree_et);
		rtree_tt = difftime(rtree_et,rtree_st);
		std::cerr << "********************************************" << std::endl;
		std::cerr << "R-tree construction total execution time: " 
			<< rtree_tt
			<< " seconds." << std::endl;
		std::cerr << "********************************************" << std::endl;
		#endif
		
//if (0){
		/*#ifdef DEBUGTIME
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


			/*#ifdef DEBUGTIME
			time(&mbb_et);
			rtree_tt = difftime(mbb_et,mbb_st);
			std::cerr << "********************************************" << std::endl;
			std::cerr << "One MBB filtering total execution time: " 
				<< mbb_tt
				<< " seconds." << std::endl;
			std::cerr << "********************************************" << std::endl;
			#endif*/

			
			
	//if(0){
			Polyhedron* geom1 = extract_geometry(sttemp.offsetdata[idx1][i], sttemp.lengthdata[idx1][i],
					stop.decomp_lod, stop, sttemp, 0);

			for (uint32_t j = 0; j < vis.matches.size(); j++) 
			{
				
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
					sttemp.lengthdata[idx2][vis.matches[j]], stop.decomp_lod,
					stop, sttemp, 1);
				struct mbb_3d * env2 = sttemp.mbbdata[idx2][vis.matches[j]];
				//std::cout << "already got geom1 and geom2: " << std::endl; // to compare * operator and yes or no question				

				// for now only decomp time
				//std::cout << i << TAB << vis.matches[j] << std::endl;
				#ifdef DEBUG
				std::cerr << "Checking actual intersection between " << i << TAB << vis.matches[j] << std::endl;
				#endif
				
				#ifdef DEBUGTIME
				time_t geometry_st, geometry_et;
				double geometry_tt;
				time(&geometry_st);
				#endif
				
				if (join_with_predicate(stop, sttemp, geom1, geom2, env1, env2,
							stop.join_predicate))  {
					//std::cerr << i << TAB << vis.matches[j] << std::endl;
				//	report_result(stop, sttemp, i, vis.matches[j]);
				//	std::cout << i << TAB << vis.matches[j] << << sttemp. <<std::endl;				

				#ifdef DEBUG
				std::cerr << "Actual intersected with each other " << i << TAB << vis.matches[j] << std::endl;
				#endif
					pairs++;
				}
				
				#ifdef DEBUGTIME
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
//}
		}
		//shmdt(shm);  // detach shared memory segment

		

	} // end of try

	catch (Tools::Exception& e) {
	//catch (...) {
		std::cerr << "******ERROR******" << std::endl;
		#ifdef DEBUG
		std::cerr << e.what() << std::endl;
		#endif
		return -1;
	} // end of catch
	
	

	std::cout << pairs << std::endl;

	delete spidx;
	delete storage;
	return pairs ; 
//}
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

