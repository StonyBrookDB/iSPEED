/* 
 * This program has many functionalities depending on the parameters
 *   set in progparams/resque_params_2d.hpp
 *
 *   Flags used in processing are located in stop (spatial operator struct)
 *
 *   It can map objects to their respective tiles given 
 *   a partition schema which is read from a disk
 * or it can be used to extract MBBs from objects
 *
 *  It can also use sampling rate. See the usage() in progparams/resque_params_2d.hpp
 *        and progparams/resque_datastructs.h for the list functionalities
 * */

#include <transform/manipulate_3d.hpp>

using namespace std;

/* Performance metrics */
clock_t start_reading_data;
clock_t start_query_exec;

clock_t total_reading;
clock_t total_query_exec;

const string STR_3D_HEADER = "OFF";

// Initialize default values in spatial operator structure
void init(struct query_op &stop, struct query_temp &sttemp) {
	stop.extract_mbb = false;
	stop.collect_mbb_stat = false;
	stop.use_sampling = false;
	stop.sample_rate = 1.0;
	stop.offset = 0;

	stop.prefix_1 = NULL;
	stop.prefix_2 = NULL;
	stop.shape_idx_1 = 0;
	stop.shape_idx_2 = 0;
}


/* Build indexing on tile boundaries from cache file */
bool build_index_tiles(struct query_op &stop, struct query_temp &sttemp, 
	IStorageManager* &storage, ISpatialIndex * &spidx, 
	std::map<SpatialIndex::id_type, std::string> *id_tiles) {
	// build spatial index on tile boundaries 
	id_type  indexIdentifier;
	GEOSDataStreamFileTile stream(stop.cachefilename, id_tiles); // input from cache file
	storage = StorageManager::createNewMemoryStorageManager();
	spidx   = RTree::createAndBulkLoadNewRTree(RTree::BLM_STR, stream, *storage, 
			FillFactor,
			IndexCapacity,
			LeafCapacity,
			3, 
			RTree::RV_RSTAR, indexIdentifier);

	// Error checking 
	return spidx->isIndexValid();
}

/* Process standard input and emit objects to their respective partitions */
bool process_input(struct query_op &stop, struct query_temp &sttemp,
		const int join_idx, const int geom_idx, 
		IStorageManager * &storage, ISpatialIndex * &spidx,
		std::map<id_type, string> *id_tiles) {
	
	bool firstLineRead = false;

	MyVisitor vis;

	/* Space info */
	double space_low[3];
	double space_high[3];

	long count_objects = 0;
	#ifdef DEBUG
	long count_emitted = 0;
	long count_bad = 0;
	#endif


	#ifdef DEBUGTIME
	start_reading_data = clock();
	#endif

	// Variables to track
	int random_id = 0;
	if (stop.extract_mbb) {
		struct timeval time; 
		gettimeofday(&time,NULL);
		srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
		random_id = rand() % 10000; // Can use the number of reducers here for better load balancing/key hashing
	}
	
	string input_line;
	vector<string> fields;
	vector<string> inner_fields;
	double tmp_x, tmp_y, tmp_z;
	double low[3];
	double high[3];
	double xx = 20, yy = 20, zz = 20; // buffer size: the maximum MBB of blood vessels

	char* stdin_file_name = NULL; // name of the input file
	stdin_file_name = getenv("mapreduce_map_input_file");
	string prefix(stdin_file_name);
	// srand (static_cast <unsigned>(time(NULL)));

	/* Handling standard input */
	while(cin && getline(cin, input_line) && !cin.eof()){
		/* Sampling (optional parameter) */
		if (stop.use_sampling &&
			(static_cast<float>(rand()) / static_cast<float>(RAND_MAX)) > stop.sample_rate) {
			continue; // skip the record
		}

		// Removal of \r symbol on Windows
		if (input_line.at(input_line.size() - 1) == '\r') {
			input_line = input_line.substr(0, input_line.size() - 1);
		}

		#ifdef COMPRESSED
		tokenize(input_line, fields, TAB, true);
		#else
		tokenize(input_line, fields, BAR, true);
		#endif
		count_objects++;
		try {

			#ifdef COMPRESSED
			// for compressed data
			low[0] = stod(fields[3]);
			low[1] = stod(fields[4]);
			low[2] = stod(fields[5]);
			high[0] = stod(fields[6]);
			high[1] = stod(fields[7]);
			high[2] = stod(fields[8]);
			if(stop.sid_second_set == stod(fields[2])){
				low[0] -= xx;
				low[1] -= yy;
				low[2] -= zz;
				high[0] += xx;
				high[1] += yy;
				high[2] += zz;			
			}
			#endif

			Region r(low, high, 3);
			/* Find objects matching with intersecting tiles*/	
			spidx->intersectsWithQuery(r, vis);

			#ifdef DEBUG
			cerr << "intersecting with: " << vis.matches.size() << endl;
			#endif
			/* Emit objects to intersecting tiles */
			for (uint32_t i = 0 ; i < vis.matches.size(); i++ ) {
				if (stop.reading_mbb) {
					cout << (*id_tiles)[vis.matches[i]]
						<< TAB << low[0] << TAB << low[1] << TAB << low[2]
						<< TAB << high[0] << TAB << high[1] << TAB << high[2] << endl;

				} else if (stop.drop_join_idx) {
					// Append tile id only
					cout <<  (*id_tiles)[vis.matches[i]]
						<< TAB << input_line <<  endl ;
				} else {
					// Append tile id and join index

					#ifdef COMPRESSED
					cout <<  (*id_tiles)[vis.matches[i]]
						<< TAB << fields[2]
						//							<< object_id
						<< TAB << fields[1] 
						// MBB Info	
						<< TAB << low[0] << TAB << low[1] << TAB << low[2]
						<< TAB << high[0] << TAB << high[1] << TAB << high[2]
						<< TAB << fields[9] << TAB << fields[10] << endl; // offset is the beginning point of this object
					#else
					cout <<  (*id_tiles)[vis.matches[i]]
						<< TAB << join_idx
						//							<< object_prefix 

						<< TAB << prefix << BAR << count_objects // Object "ID"	
						// MBB Info	
						<< TAB << low[0] << TAB << low[1] << TAB << low[2]
						<< TAB << high[0] << TAB << high[1] << TAB << high[2]
						<< TAB << input_line <<  endl ; // geometry
					#endif

				}
				#ifdef DEBUG
				count_emitted++;
				#endif
			}
			vis.matches.clear();

		}
		catch (...)
		{
#ifdef DEBUG
			count_bad++;
			cerr << "WARNING: Record is not well formatted " << input_line << endl;
#endif
			continue; // ignore bad formatting and continue
		}

		fields.clear();

#ifdef DEBUGTIME
		total_query_exec += clock() - start_query_exec;
		start_reading_data = clock();
#endif
	}
	/* Output dimensions of space */
	if (stop.extract_mbb) {
		cout << "-1" << TAB << space_low[0] << TAB << space_low[1] << TAB << space_low[2] 
			<< TAB << space_high[0]	<< TAB << space_high[1]	<< TAB << space_high[2] << TAB << count_objects << endl;
	}

#ifdef DEBUG
	/* Output useful statistics */
	cerr << "Number of processed objects: " << count_objects << endl;
	cerr << "Number of times objects were emitted: " << count_emitted << endl;
	cerr << "Number of not well formatted objects: " << count_bad << endl;
#endif
}


int main(int argc, char **argv) {	
	struct query_op stop;
	struct query_temp sttemp;
	std::map<id_type, string> id_tiles;

	init(stop, sttemp);


	if (!extract_params(argc, argv, stop, sttemp)) {
#ifdef DEBUG 
		cerr <<"ERROR: query parameter extraction error." << endl 
			<< "Please see documentations, or contact author." << endl;
#endif
		usage();
		return -1;
	}

	char* stdin_file_name = NULL; // name of the input file
	int join_idx = -1; // index of the current file (0 or 1) matching to dataset 1 or 2
	int geom_idx = -1; // geometry field index

	IStorageManager * storage = NULL;
	ISpatialIndex * spidx = NULL;
	if (!stop.extract_mbb) {
		if( !build_index_tiles(stop, sttemp, storage, spidx, &id_tiles)) {
#ifdef DEBUG
			cerr << "ERROR: Index building on tile structure has failed ." << std::endl;
#endif
			return 1 ;
		}
		else { 
#ifdef DEBUG  
			cerr << "GRIDIndex Generated successfully." << endl;
#endif
		}
	}

	// Process input line by line
	process_input(stop, sttemp, join_idx, geom_idx, storage, spidx, &id_tiles);

	/* Clean up indices */
	delete spidx;
	delete storage;

	id_tiles.clear();

	cout.flush();
	cerr.flush();

	return 0;

}


