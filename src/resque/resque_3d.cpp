
#include <resque/resque_3d.hpp>

using namespace std;

/* 
 * RESQUE processing engine v4.0
 *   It supports spatial join and nearest neighbor query with different predicates
 *   1) parseParameters
 *   2) readCacheFile - metadata such as partition schemata
 *   3) for every input line in the current tile
 *         an input line represents an object
 *         save geometry and original data in memory
 *         execute join operation when finish reading a tile
 *   4) Join operation between 2 sets or a single set
 *         build Rtree index on the second data set
 *         for every object in the first data set
 *            using Rtree index of the second data set
 *              check for MBR/envelope intersection
 *              output the pair result or save pair statistics
 *   5) Output final statistics (opt)
 *   Requirement (input files): see the Wiki
 * */

/* Performance metrics */
clock_t start_reading_data;
clock_t start_query_exec;

clock_t total_reading;
clock_t total_query_exec;

Sc_Polyhedron* sc_extract_geometry(long offset, long length, unsigned i_decompPercentage,
	struct query_op &stop, struct query_temp &sttemp, int dataset_id) {
	Sc_Polyhedron* geom;
	/*
	// moved to global vars
	int shmid;
    	key_t key;
    	char *shm, *s;
	*/
	/*
	int shmid;
	//use the same key to locate the segment.
	if ((shmid = shmget(COMPRESSION_KEY, SHMSZ, 0666)) < 0) {
		perror("shmget");
		exit(1);
	}
	    
	// Now we attach the segment to our data space.
	if ((shm = (char *) shmat(shmid, NULL, 0)) == (char *) -1) {
		perror("shmat");
		exit(1);
	}
	*/
	// Initialize parameters
	int i_mode = DECOMPRESSION_MODE_ID; // compression mode

	#ifdef DEBUG
	std::cout << "sk_geometry: attempting to extract " << offset << TAB << length << std::endl;
	#endif
	std::cerr << "sk_geometry: attempting to extract " << offset << TAB << length << std::endl;
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
	//std::cerr << "got to this part 1" << std::endl;
	//read the mesh:

	// Never pass NULL to a std::string type parameter that uses a lazy initialization constructor :(

	/*
	ifstream file("tmpData/failcompressed.pp3d", std::ios::binary | std::ios::ate);
	file.seekg(0, std::ios::beg);

	char * fbuffer = new char[2000000];
        file.read(fbuffer,length);
    	for (size_t i = 0; i < length; ++i) {
		printf("%02X ", fbuffer[i]);
	}
	printf("\n");
	*/
	/*
	printf("resque shared buffer: ");
    	for (size_t i = 0; i < length; ++i) {
		printf("%02X ", resque_decomp_buffer[i]);
	}
	printf("\n");
	printf("shm current at offset  ");
    	for (size_t i = 0; i < length; ++i) {
		printf("%02X ", ((char*)(shm_ptr + offset))[i]);
	}
	printf("\n");
	*/
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
	
	currentMesh->completeOperation();
	
	// debug

	//std::cerr << "current mesh: " << *currentMesh << std::endl;
	std::stringstream os;
	os << *currentMesh;
	//os.clear();
	
	std::cerr << "done decomp" << std::endl;
	geom = new Sc_Polyhedron();
	os >> *geom;
	//std::cerr << "os: " << os.str() << std::endl;
	
	// only when volume is needed
	//if (stop.needs_intersect_volume) {
		sttemp.poly_str[dataset_id].str(os.str());
	//}

	//delete[] fbuffer;
	std::cerr << "constructing sk poly" << std::endl;
	//std::cerr << "geom: " << *geom << std::endl;
	delete currentMesh;
	return geom;
}

// Performs spatial query on data stored in query_temp using operator query_op
int execute_query(struct query_op &stop, struct query_temp &sttemp)
{
	// Processing variables
	std::string input_line; // Temporary line
	std::vector<std::string> fields; // Temporary fields
	int sid = 0; // Join index ID for the current object
	std::string tile_id = ""; // The current tile_id
	std::string previd = ""; // the tile_id of the previously read object
	int tile_counter = 0; // number of processed tiles

	/* CGAL variables for spatial computation */
	Polyhedron *poly = NULL;
	
	long offset = 0, length = 0;	
	struct mbb_3d *mbb_ptr = NULL;

	/* Define the resource when using cache-file  */
	//int maxCardRelease = std::min(stop.join_cardinality, 2);
	int maxCardRelease = 2;

	#ifdef DEBUG
	std::cerr << "Bucket info:[ID] |A|x|B|=|R|" <<std::endl;
	#endif  

	#ifdef DEBUGTIME
	start_reading_data = clock();
	#endif

	#ifdef DEBUGTIME
	time_t data_st, data_et;
	double data_tt;
	time(&data_st);
	#endif

	std::stringstream ss;
	// with compressed data
	int shmid;		
	//use the same key to locate the segment.
	size_t maxoffset2 = stop.shm_max_size;
	//std::cerr << "max offset" << maxoffset2 << std::endl;
	// Getting access to shared memory with all compressed objects stored in there
	if ((shmid = shmget(COMPRESSION_KEY, maxoffset2, 0666)) < 0) {
		perror("shmget");
		exit(1);
	}
	//std::cerr << "maxoffset: " << maxoffset << std::endl;
	// Now we attach the segment to our data space.
	if ((shm_ptr = (char *) shmat(shmid, (const void *)NULL, 0)) == (char *) -1) {
		perror("shmat");
		exit(1);
	}

	/*printf("shm current at offset  ");
        for (size_t i = 0; i < 1000; ++i) {
                 printf("%02X ",  shm_ptr[i]);
         }
         printf("\n");*/

	char *decomp_buffer =  new char[BUFFER_SIZE];
	resque_decomp_buffer = decomp_buffer;
	//std::cerr << "decomp_buffer" << (long) decomp_buffer << TAB << BUFFER_SIZE << std::endl;
	//std::cerr << "decomp_buffer" << (long) resque_decomp_buffer << TAB << BUFFER_SIZE << std::endl;
	
	// Read line by line inputs
	int tt = 0;
	while (std::cin && getline(std::cin, input_line) && !std::cin.eof()) {
#ifdef DEBUG
		std::cerr<<tt++<<" line content:"<<input_line<<std::endl;
#endif
		tokenize(input_line, fields, TAB, true);
		if(fields.size()==0){
			break;
		}
	
		tile_id = fields[0];
		sid = atoi(fields[1].c_str());

		/* Parsing fields from input */
		try { 
			// Parsing MBB
			mbb_ptr = new struct mbb_3d();
			for (int k = 0; k < NUMBER_DIMENSIONS; k++) {
				//mbb_ptr->low[k] = stod(fields[3 + k]);
				mbb_ptr->low[k] = std::atof(fields[3 + k].c_str());
			}
			for (int k = 0; k < NUMBER_DIMENSIONS; k++) {
				//mbb_ptr->high[k] = stod(fields[6 + k]);
				mbb_ptr->high[k] = std::atof(fields[6 + k].c_str());
			}
			#ifdef DEBUG
			//compare with the input line to make sure the parsing process is correct
			std::cerr << "MBB: ";
			for (int k = 0; k < NUMBER_DIMENSIONS; k++) {
				std::cerr << TAB << mbb_ptr->low[k];
			}
			for (int k = 0; k < NUMBER_DIMENSIONS; k++) {
				std::cerr << TAB << mbb_ptr->high[k];
			}
			std::cerr << std::endl;
			#endif
		}
		catch (...) {
			std::cerr << "******MBB Parsing Error******" << std::endl;
			return -1;
		}
		
		try { 
			// Parsing Geometry from shared memory segment
			//poly = new Polyhedron();
			//boost::replace_all(fields[9], BAR, "\n");
			//ss.str(fields[9]);	
			//std::cout << input_line << std::endl;
			//ss >> *poly;
//			std::cout << ss.str() << std::endl;
			// extracting MBB information
			// mbb_ptr = get_mbb(poly);
			offset = atol(fields[9].c_str());
			length = atol(fields[10].c_str());
			//maxoffset = std::max(offset + length, maxoffset); // track the maximum offset that has been seen so far
		}
		catch (...) {
			std::cerr << "******Offset and Length Parsing Error******" << std::endl;
			return -1;
		}

		/* Process the current tile (bucket) when finishing reading all objects belonging
		   to the current tile */
		if (previd.compare(tile_id) != 0 && previd.size() > 0 ) {
			#ifdef DEBUGTIME
			total_reading += clock() - start_reading_data;
			start_query_exec = clock();
			#endif
			sttemp.tile_id = previd;
			// Process the current tile in memory
			int pairs = join_bucket(stop, sttemp); // number of satisfied predicates
			#ifdef DEBUG
			std::cerr <<"Special T[" << previd << "] |" << sttemp.mbbdata[SID_1].size()
				<< "|x|" << sttemp.mbbdata[stop.sid_second_set].size()
				<< "|=|" << pairs << "|" << std::endl;
			#endif


			#ifdef DEBUGTIME
			total_query_exec += clock() - start_query_exec;
			start_reading_data = clock();
			#endif


			tile_counter++; 
			release_mem(stop, sttemp, maxCardRelease);
		}

		// populate the bucket for join 
		//sttemp.polydata[sid].push_back(poly);
		sttemp.offsetdata[sid].push_back(offset);
		sttemp.lengthdata[sid].push_back(length);
		sttemp.mbbdata[sid].push_back(mbb_ptr);
		fields.pop_back(); // Remove the last geometry field to save space
		fields.pop_back();		
		
		sttemp.rawdata[sid].push_back(fields);

		/* Update the field */
		previd = tile_id; 
		fields.clear();
	}

	#ifdef DEBUGTIME
	time(&data_et);
	data_tt = difftime(data_et,data_st);
	std::cerr << "********************************************" << std::endl;
	std::cerr << "Data loading and parsing total execution time: " 
		<< data_tt
		<< " seconds." << std::endl;
	std::cerr << "********************************************" << std::endl;
	#endif

	#ifdef DEBUGTIME
	total_reading += clock() - start_reading_data;
	start_query_exec = clock();
	#endif

	// Process the last tile (what remains in memory)
	sttemp.tile_id = tile_id;
	int pairs = join_bucket(stop, sttemp); // number of satisfied predicates

	#ifdef DEBUGTIME
	total_query_exec += clock() - start_query_exec;
	start_reading_data = clock();
	#endif


	#ifdef DEBUG
	std::cerr <<"Special 2 T[" << previd << "] |" << sttemp.mbbdata[SID_1].size() << "|x|" 
		<< sttemp.mbbdata[stop.sid_second_set].size() 
		<< "|=|" << pairs << "|" << std::endl;

	shmdt(shm_ptr);
	#endif

	tile_counter++;

	release_mem(stop, sttemp, stop.join_cardinality);
	
	// clean up newed objects
	delete[] decomp_buffer;

	return tile_counter;
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
		sttemp.rawdata[delete_index].clear();
  	}
}

/* Create an R-tree index on a given set of polygons */
bool build_index_geoms(std::vector<struct mbb_3d *> & geom_mbbs, SpatialIndex::ISpatialIndex* & spidx, SpatialIndex::IStorageManager* & storage) {
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


/* 
 *  Perform spatial computation on a given tile with data 
 *   located in polydata and rawdata
 *   
 */
int join_bucket(struct query_op &stop, struct query_temp &sttemp)
{	
	// return join_bucket_voronoi
	// return join_bucket_knn
	return join_bucket_spjoin(stop, sttemp);
}

/* main body of the engine */
int main(int argc, char** argv)
{
	int c = 0; // Number of results satisfying the predicate

	struct query_op stop;
	struct query_temp sttemp;

	if (!extract_params(argc, argv, stop, sttemp)) { // Function is located in params header file
		#ifdef DEBUG 
		std::cerr <<"ERROR: query parameter extraction error." << std::endl 
		     << "Please see documentations, or contact author." << std::endl;
		#endif
		usage();
		return 1;
	}
	// Query execution	
	// Spatial join and nearest neighbors from joint datasets (stdin)
	switch (stop.join_cardinality) {
		case 1:
		case 2:
			// adjusting set id
			stop.sid_second_set = stop.join_cardinality == 1 ? SID_1 : SID_2;
			#ifdef DEBUG
			std::cerr <<"sid_second_set" << stop.sid_second_set << std::endl ;
			#endif
			c = execute_query(stop, sttemp);
			break;
		default:
			#ifdef DEBUG
			std::cerr <<"ERROR: join cardinality does not match engine capacity." << std::endl ;
			#endif
			return 1;
	}

	if (c >= 0 ) {
		#ifdef DEBUG 
		std::cerr << "Query Load: [" << c << "]" << std::endl;
		#endif
	} else {
		#ifdef DEBUG 
		std::cerr <<"Error: ill formatted data. Terminating ....... " << std::endl;
		#endif
		return 1;
	}

	#ifdef DEBUGTIME
	std::cerr << "Total reading time: " 
		<< (double) total_reading / CLOCKS_PER_SEC 
		<< " seconds." << std::endl;
	std::cerr << "Total query exec time: " 
		<< (double) total_query_exec / CLOCKS_PER_SEC 
		<< " seconds." << std::endl;
	#endif

	std::cout.flush();
	std::cerr.flush();
	return 0;
}



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


bool intersects_mbb(const struct mbb_3d * m1, const struct mbb_3d *m2) {
	return !(m1->low[0] > m2->high[0] || m1->high[0] < m2->low[0]
	      || m1->low[1] > m2->high[1] || m1->high[1] < m2->low[1]
	      || m1->low[2] > m2->high[2] || m1->high[2] < m2->low[2] );
}


bool intersects(Polyhedron *P1, Polyhedron *P2, const struct mbb_3d * env1, const struct mbb_3d * env2) {
	//return true;
	// Use Nef_polyhedron for intersection detection
	if(!intersects_mbb(env1, env2))
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


