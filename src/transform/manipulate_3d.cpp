/* 
 *
 *   It can map objects to their respective tiles given 
 *   a partition schema which is read from a disk
 * or it can be used to extract MBBs from objects
 *
 * */

#include <transform/manipulate_3d.hpp>

using namespace std;

/* Performance metrics */
clock_t start_reading_data_manipulate;
clock_t start_query_exec_manipulate;

clock_t total_reading_manipulate;
clock_t total_query_exec_manipulate;

const string STR_3D_HEADER = "OFF";

/* Build indexing on tile boundaries from cache file */
bool build_index_tiles(char *cachefilename,
			IStorageManager* &storage, ISpatialIndex * &spidx,
			std::map<SpatialIndex::id_type, std::string> *id_tiles) {
	// build spatial index on tile boundaries
	id_type  indexIdentifier;
	GEOSDataStreamFileTile stream(cachefilename, id_tiles); // input from cache file
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
void process_input(IStorageManager * &storage, ISpatialIndex * &spidx,
		std::map<id_type, string> *id_tiles) {

	bool firstLineRead = false;

	MyVisitor vis;

	long count_objects = 0;
	#ifdef DEBUG
	long count_emitted = 0;
	long count_bad = 0;
	start_reading_data_manipulate = clock();
	#endif

	string input_line;
	vector<string> fields;
	vector<string> inner_fields;
	double tmp_x, tmp_y, tmp_z;
	double low[3];
	double high[3];
	double xx = 20, yy = 20, zz = 20; // buffer size: the maximum MBB of blood vessels

	/* Handling standard input */
#ifdef DEBUG
	std::cerr<<"mapper_id obj_id dataset_id 6*mbbs global_offset length"<<std::endl;
#endif
	while(cin && getline(cin, input_line) && !cin.eof()){

		// Removal of \r symbol on Windows
		if (input_line.at(input_line.size() - 1) == '\r') {
			input_line = input_line.substr(0, input_line.size() - 1);
		}

		tokenize(input_line, fields, TAB, true);
		count_objects++;
		try {
			// for compressed data
			low[0] = stod(fields[3]);
			low[1] = stod(fields[4]);
			low[2] = stod(fields[5]);
			high[0] = stod(fields[6]);
			high[1] = stod(fields[7]);
			high[2] = stod(fields[8]);
			//the second dataset will be buffered to assign it to multiple
			//tiles if needed
			if(stod(fields[2])==2){
				low[0] -= xx;
				low[1] -= yy;
				low[2] -= zz;
				high[0] += xx;
				high[1] += yy;
				high[2] += zz;
			}

			Region r(low, high, 3);
			/* Find objects matching with intersecting tiles*/
			spidx->intersectsWithQuery(r, vis);
			/* Emit objects to intersecting tiles */
			for (uint32_t i = 0 ; i < vis.matches.size(); i++ ) {
				cout<< (*id_tiles)[vis.matches[i]]	//tile id
					<< TAB << fields[1]				//object id
					<< TAB << fields[2]				//dataset id
					// MBB Info
					<< TAB << low[0] << TAB << low[1] << TAB << low[2]
					<< TAB << high[0] << TAB << high[1] << TAB << high[2]
					<< TAB << fields[9] << TAB << fields[10] << endl; // offset is the beginning point of this object
				#ifdef DEBUG
				count_emitted++;
				#endif
			}
			vis.matches.clear();
		}catch (...) {
			#ifdef DEBUG
			count_bad++;
			cerr << "WARNING: Record is not well formatted " << input_line << endl;
			#endif
			continue; // ignore bad formatting and continue
		}

		fields.clear();

		#ifdef DEBUG
		total_query_exec_manipulate += clock() - start_query_exec_manipulate;
		start_reading_data_manipulate = clock();
		#endif
	}

	#ifdef DEBUG
	/* Output useful statistics */
	cerr << "Number of objects processed:\t" << count_objects << endl;
	cerr << "Number of objects emitted:\t" << count_emitted << endl;
	cerr << "Number of objects malformed:\t" << count_bad << endl;
	#endif
}


int main(int argc, char **argv) {
	std::map<id_type, string> id_tiles;
	char *cachefilename = NULL;
	if(argc < 2){
		std::cerr<<"usage: manipulate /path/to/partition/file"<<std::endl;
		return 0;
	}
	cachefilename = argv[1];
	IStorageManager * storage = NULL;
	ISpatialIndex * spidx = NULL;
	if( !build_index_tiles(cachefilename, storage, spidx, &id_tiles)) {
		cerr << "ERROR: Index building on tile structure has failed ." << std::endl;
		exit(0) ;
	} else {
		#ifdef DEBUG
		cerr << "GRIDIndex Generated successfully." << endl;
		#endif
	}

	// Process input line by line
	process_input(storage, spidx, &id_tiles);

	/* Clean up indices */
	delete spidx;
	delete storage;

	id_tiles.clear();

	cout.flush();
	cerr.flush();

	return 0;

}


