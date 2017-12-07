
#include <resque/resque_3d.hpp>

// Spatial join operations
#include <resque/spjoin_3d.hpp>


/* 
 * RESQUE processing engine v3.0
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

/* Initialize default values in query structs (operator and temporary placeholders) */
/* To be potentially removed to adjust for initialization already 
 * 	been done in param extraction method */
void init(struct query_op &stop, struct query_temp &sttemp){
	stop.offset = 2; // default format or value for offset
}

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
/*
  Extract one polyhedron geometry from compressed data with given offset and length
*/
Polyhedron* extract_geometry(long offset, long length, unsigned i_decompPercentage,
	struct query_op &stop, struct query_temp &sttemp, int dataset_id) {
	Polyhedron* geom;
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
	std::cout << "attempting to extract " << offset << TAB << length << std::endl;
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
	
	//std::cerr << "done decomp" << std::endl;
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


/* report result separated by separator */
void report_result(struct query_op &stop, struct query_temp &sttemp, int i, int j)
{
	sttemp.stream.str("");
	sttemp.stream.clear();
	/* ID used to access rawdata for the "second" data set */
	#ifdef DEBUG
	std::cerr << "size of raw data: " << sttemp.rawdata[SID_1].size() << TAB << sttemp.rawdata[stop.sid_second_set].size() << std::endl;

	#endif
	if (stop.output_fields.size() == 0) {
		/* No output fields have been set. Print all fields read */
		for (int k = 0; k < sttemp.rawdata[SID_1][i].size(); k++) {
			sttemp.stream << sttemp.rawdata[SID_1][i][k] << SEP;
		}
		for (int k = 0; k < sttemp.rawdata[stop.sid_second_set][j].size(); k++) {
			sttemp.stream << SEP << sttemp.rawdata[stop.sid_second_set][j][k];
		}
	}
	else {
		/* Output fields are listed */
		int k = 0;
		for (; k < stop.output_fields.size() - 1; k++) {
	//		std::cerr << "outputting fields " << stop.output_fields[k];
			obtain_field(stop, sttemp, k, i, j);
			sttemp.stream << SEP;
		}
		obtain_field(stop, sttemp, k, i, j);

	//		std::cerr << "outputting fields " << stop.output_fields[k];
	}

	sttemp.stream << std::endl;
	std::cout << sttemp.stream.str();
}

/* Reporting result for the case when processing 1 by 1 object from data set 1
 *  skip_window_data == true when there is simply a single window query (data set 2)
 *      only fields from data set 1  will be output
 *  skip_window_data == false when there are more than one objects in data set 2
 * */
void report_result(struct query_op &stop, struct query_temp &sttemp, 
	std::vector<std::string> &set1fields, int j, bool skip_window_data)
{
	sttemp.stream.str("");
	sttemp.stream.clear();
	/* ID used to access rawdata for the "second" data set */

	if (stop.output_fields.size() == 0) {
		/* No output fields have been set. Print all fields read */
		for (int k = 0; k < set1fields.size(); k++) {
			sttemp.stream << set1fields[k] << SEP;
		}

		if (!skip_window_data) {
			for (int k = 0; k < sttemp.rawdata[SID_2][j].size(); k++) {
				sttemp.stream << SEP << sttemp.rawdata[SID_2][j][k];
			}
		}
	}
	else {
		/* Output fields are listed */
		int k = 0;
		for (; k < stop.output_fields.size() - 1; k++) {
	//		std::cerr << "outputting fields " << stop.output_fields[k];
			obtain_field(stop, sttemp, k, set1fields, j);
			sttemp.stream << SEP;
		}
		obtain_field(stop, sttemp, k, set1fields, j);

	//		std::cerr << "outputting fields " << stop.output_fields[k];
	}

	sttemp.stream << std::endl;
	std::cout << sttemp.stream.str();
}


/* Performs a spatial query processing where set 2 is obtained from the cache file */
int execute_query_cache_file(struct query_op &stop, struct query_temp &sttemp) {
	int num_obj_file;
	int count = 0; // Returns the number

	// Processing variables
	std::string input_line; // Temporary line
	std::vector<std::string> fields; // Temporary fields
	int sid = 0; // Join index ID for the current object
	int index = -1;  // Geometry field position for the current object
	std::string tile_id = ""; // The current tile_id
	std::string previd = ""; // the tile_id of the previously read object
	int tile_counter = 0; // number of processed tiles

	/* GEOS variables for spatial computation */
	SpatialIndex::IStorageManager *storage = NULL;
	SpatialIndex::ISpatialIndex *spidx = NULL;

	Polyhedron *poly = NULL;
	struct mbb_3d *mbb_ptr = NULL;
	Polyhedron *windowpoly = NULL;
	struct mbb_3d *windowmbb_ptr = NULL;

	std::ifstream input(stop.cachefilename);

	sid = SID_2;
	index = stop.shape_idx_2 ; 
	num_obj_file = 0;

	std::stringstream ss;
	// Reading from the cache file
	while(!input.eof() && getline(input, input_line)) {
		tokenize(input_line, fields, TAB, true);

		/* Handling of objects with missing geometry */
		if (fields[index].size() <= 0) 
			continue ; //skip empty spatial object 
		
		#ifdef DEBUG
		std::cerr << "geometry: " << fields[stop.shape_idx_2]<< std::endl;
		#endif  
		
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

			/* Below code achieves the same thing
			mbb_ptr->low[0] = atoi(fields[3]);
			mbb_ptr->low[1] = atoi(fields[4]);
			mbb_ptr->low[2] = atoi(fields[5]);
			mbb_ptr->high[0] = atoi(fields[6]);
			mbb_ptr->high[1] = atoi(fields[7]);
			mbb_ptr->high[2] = atoi(fields[8]);
			*/
			#ifdef DEBUG
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
		
		/* Parsing polyhedron input */
		try { 
			// Parsing Geometry
			poly = new Polyhedron();
			boost::replace_all(fields[9], BAR, "\n");
			ss.str(fields[9]);	
			//std::cout << input_line << std::endl;
			ss >> *poly;
		}
		catch (...) {
			std::cerr << "******Geometry Parsing Error******" << std::endl;
			return -1;
		}
		sttemp.polydata[sid].push_back(poly);
		sttemp.mbbdata[sid].push_back(mbb_ptr);
		fields.pop_back(); // Remove the last geometry field to save space
		sttemp.rawdata[sid].push_back(fields);

		num_obj_file++;

		fields.clear();
	}
	#ifdef DEBUG
	std::cerr << "Read " << num_obj_file << " from the cache file." << std::endl;
	#endif
	if (num_obj_file <= 0) {
		#ifdef DEBUG
		std::cerr << "No object in cache file." << std::endl;
		#endif
		return -1; 
	}
	
	if (num_obj_file == 1) {
		// Single window range query
		windowpoly = poly;
		windowmbb_ptr = mbb_ptr;
	} else {
		// Build R*-tree index

		/* Build index on the "second data set */
		std::vector<struct mbb_3d *> geom_mbb2;
		geom_mbb2.clear();

		int len2 = sttemp.mbbdata[SID_2].size();
		// Make a copy of the std::vector to map to build index (API restriction)
		for (int j = 0; j < len2; j++) {
			geom_mbb2.push_back(sttemp.mbbdata[SID_2][j]);
		}

		/* Handling for special nearest neighbor query */	
		// build the actual spatial index for input polygons from idx2
		if (!build_index_geoms(geom_mbb2, spidx, storage)) {
			#ifdef DEBUG
			std::cerr << "Building index on geometries from set 2 has failed" << std::endl;
			#endif
			return -1;
		}

		// must clear memory of storage and spidx at the end
	}

	index = stop.shape_idx_1 ; 
	// Process standard input (dataset 1)
	while (std::cin && getline(std::cin, input_line) && !std::cin.eof()) {
		tokenize(input_line, fields, TAB, true);
		/* Handling of objects with missing geometry */
		if (fields[index].size() <= 0) 
			continue ; //skip empty spatial object 
		
		#ifdef DEBUG
		std::cerr << "geometry: " << fields[stop.shape_idx_1]<< std::endl;
		#endif  

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
			// Parsing Geometry
			poly = new Polyhedron();
			boost::replace_all(fields[9], BAR, "\n");
			ss.str(fields[9]);	
			//std::cout << input_line << std::endl;
			ss >> *poly;
			
		}
		catch (...) {
			std::cerr << "******Polyhedron Parsing Error******" << std::endl;
			return -1;
		}

		if (num_obj_file == 1) {
			// Uses a function from spjoin file here
			if (intersects(poly, windowpoly, mbb_ptr, windowmbb_ptr) && join_with_predicate(stop, sttemp, poly, 
				windowpoly, mbb_ptr, windowmbb_ptr, stop.join_predicate)) {
				//report_result(stop, sttemp, fields, 0, true); // the index when there is only 1 object is 0
			}
			
		} else {
		}
		
		delete poly;
		delete windowpoly;
		delete mbb_ptr;
		delete windowmbb_ptr;

		fields.clear();
	}

	// clean up newed objects
	if (num_obj_file > 1) {
		delete spidx;
		delete storage;
	}
	
	return count;
}

// Performs spatial query on data stored in query_temp using operator query_op
int execute_query(struct query_op &stop, struct query_temp &sttemp)
{
	// Processing variables
	std::string input_line; // Temporary line
	std::vector<std::string> fields; // Temporary fields
	int sid = 0; // Join index ID for the current object
	int index = -1;  // Geometry field position for the current object
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
	//#ifdef COMPRESSED  
	// with compressed data
	int shmid;		
	//use the same key to locate the segment.
	size_t maxoffset2 = stop.shm_max_size;
	//std::cerr << "max offset" << maxoffset2 << std::endl;
	if ((shmid = shmget(COMPRESSION_KEY, maxoffset2, 0666)) < 0) {
		perror("shmget");
		exit(1);
	}
	//std::cerr << "maxoffset: " << maxoffset << std::endl;
	// Now we attach the segment to our data space.
	if ((shm_ptr = (char *) shmat(shmid, NULL, 0)) == (char *) -1) {
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
	
	while (std::cin && getline(std::cin, input_line) && !std::cin.eof()) {
		tokenize(input_line, fields, TAB, true);
	
		tile_id = fields[0];
		
		sid = atoi(fields[1].c_str());
		
		#ifdef DEBUG
		std::cerr << input_line << std::endl;
		#endif 
		switch (sid) {
			case SID_1:
				index = stop.shape_idx_1 ; 
				break;
			case SID_2:
				index = stop.shape_idx_2 ; 
				break;
			default:
				std::cerr << "wrong sid : " << sid << std::endl;
				return false;
		}
		

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
			int pairs = join_bucket(stop, sttemp); // number of satisfied predicates

			#ifdef DEBUGTIME
			total_query_exec += clock() - start_query_exec;
			start_reading_data = clock();
			#endif


			#ifdef DEBUG
			std::cerr <<"Special T[" << previd << "] |" << sttemp.mbbdata[SID_1].size() 
				<< "|x|" << sttemp.mbbdata[stop.sid_second_set].size() 
				<< "|=|" << pairs << "|" << std::endl;
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

	//#ifdef COMPRESSED

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
	#ifdef COMPRESSED
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
	#else
  	for (int j = 0; j < stop.join_cardinality && j < maxCard; j++ ) {
    		int delete_index = j + 1; // index are adjusted to start from 1
    		int len = sttemp.polydata[delete_index].size();
    		for (int i = 0; i < len ; i++) {
      			delete sttemp.polydata[delete_index][i];
			delete sttemp.mbbdata[delete_index][i]; // release mbb
			sttemp.rawdata[delete_index][i].clear();
		}
    		sttemp.polydata[delete_index].clear();
    		sttemp.rawdata[delete_index].clear();
		sttemp.mbbdata[delete_index].clear();
  	}
	#endif
}


/* Compute distance between two points using Euclidian distance */
/*double get_distance(const geos::geom::Point * p1, const geos::geom::Point * p2) 
{	return sqrt(pow(p1->getX() - p2->getX(), 2) 
			+ pow(p1->getY() - p2->getY(), 2));
}*/

/* Compute geographical distance between two points on earth */
/*double get_distance_earth(const geos::geom::Point * p1, const geos::geom::Point * p2) 
{
	return earth_distance(p1->getX(), p1->getY(), p2->getX(), p2->getY());
}
*/

/* Output the field at given position  */
void obtain_field(struct query_op &stop, struct query_temp &sttemp, 
	int position, int pos1, int pos2)
{
	//std::cerr << "Set id" << stop.output_fields_set_id[position] << std::endl;
	if (stop.output_fields_set_id[position] == SID_1) {
		sttemp.stream << sttemp.rawdata[SID_1][pos1][stop.output_fields[position]];
			
	}
	else if (stop.output_fields_set_id[position] == SID_2) {
		sttemp.stream << sttemp.rawdata[stop.sid_second_set][pos2][stop.output_fields[position]];	
	}
	else if (stop.output_fields_set_id[position] == SID_NEUTRAL) {
		switch (stop.output_fields[position]) {
			case STATS_AREA_1:
				sttemp.stream << sttemp.area1;	
				break;
			case STATS_AREA_2:
				sttemp.stream << sttemp.area2;	
				break;
			case STATS_UNION_AREA:
				sttemp.stream << sttemp.union_area;
				break;
			case STATS_INTERSECT_AREA:
				sttemp.stream << sttemp.intersect_area;
				break;
			case STATS_JACCARD_COEF:
				sttemp.stream << sttemp.jaccard;
				break;
			case STATS_DICE_COEF:
				sttemp.stream << sttemp.dice;
				break;
			case STATS_TILE_ID:
				sttemp.stream << sttemp.tile_id;
				break;
			case STATS_MIN_DIST:
				sttemp.stream << sttemp.distance;
				break;
			case STATS_VOLUME_1:  // for 3d 
				sttemp.stream << sttemp.volume1;
				break;
			case STATS_VOLUME_2:
				sttemp.stream << sttemp.volume2;
				break;
			case STATS_INTERSECT_VOLUME:
				sttemp.stream << sttemp.intersect_volume;
				break;
			case STATS_NN_DISTANCE:
				sttemp.stream << sttemp.nn_distance;
				break;
			default:
				return;
		}					
	}
}

void obtain_field(struct query_op &stop, struct query_temp &sttemp, 
	int position, std::vector<std::string> &set1fields, int pos2)
{
	//std::cerr << "Set id" << stop.output_fields_set_id[position] << std::endl;
	if (stop.output_fields_set_id[position] == SID_1) {
		sttemp.stream << set1fields[stop.output_fields[position]];	
	}
	else if (stop.output_fields_set_id[position] == SID_2) {
		sttemp.stream << sttemp.rawdata[SID_2][pos2][stop.output_fields[position]];	
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
	return join_bucket_spjoin(stop, sttemp);
}

/* main body of the engine */
int main(int argc, char** argv)
{
	int c = 0; // Number of results satisfying the predicate

	struct query_op stop;
	struct query_temp sttemp;

	init(stop, sttemp); // setting the query operator and temporary variables to default

	if (!extract_params(argc, argv, stop, sttemp)) { // Function is located in params header file
		#ifdef DEBUG 
		std::cerr <<"ERROR: query parameter extraction error." << std::endl 
		     << "Please see documentations, or contact author." << std::endl;
		#endif
		usage();
		return 1;
	}

	#ifdef DEBUG 
	std::cerr <<"use cache file: " << stop.use_cache_file << std::endl ;
	#endif
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

