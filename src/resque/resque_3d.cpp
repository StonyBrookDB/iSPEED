
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


/* main body of the engine */
int main(int argc, char** argv)
{
	int c = 0; // Number of results satisfying the predicate

	struct query_op stop;
	struct query_temp sttemp;
	if (!extract_params(argc, argv, stop, sttemp)) { // Function is located in params header file
		std::cerr <<"ERROR: query parameter extraction error." << std::endl 
		     << "Please see documentations, or contact author." << std::endl;
		usage();
		return 1;
	}
	// Query execution
	// Spatial join and nearest neighbors from joint datasets (stdin)

	assert(stop.join_cardinality==1||stop.join_cardinality==2);
	stop.sid_second_set = stop.join_cardinality == 1 ? SID_1 : SID_2;

	c = execute_query(stop, sttemp);
	if (c >= 0 ) {
		std::cerr << "Query Load: [" << c << "]" << std::endl;
	} else {
		std::cerr <<"Error: ill formatted data. Terminating ....... " << std::endl;
		return 1;
	}

	std::cerr << "Total reading time: " 
		<< (double) total_reading / CLOCKS_PER_SEC 
		<< " seconds." << std::endl;
	std::cerr << "Total query exec time: " 
		<< (double) total_query_exec / CLOCKS_PER_SEC 
		<< " seconds." << std::endl;

	std::cout.flush();
	std::cerr.flush();
	return 0;
}

/*dispatch joins to target module*/
int join_bucket(struct query_op &stop, struct query_temp &sttemp){

	assert(stop.join_cardinality == 2 && "cannot conduct nn on same data set");

	//return if either one dataset is empty
 	if (sttemp.mbbdata[SID_1].size() <= 0 || sttemp.mbbdata[SID_2].size() <= 0) {
		return 0;
	}

	/* Define the resource when using cache-file  */
	//int maxCardRelease = std::min(stop.join_cardinality, 2);
	int maxCardRelease = 2;
	total_reading += clock() - start_reading_data;
	start_query_exec = clock();
	// Process the current tile in memory
	int pairs = 0;
	switch(stop.join_predicate){
	case ST_NN_VORONOI:
		pairs = join_bucket_nn_voronoi(stop, sttemp);
		break;
	case ST_NN_RTREE:
		pairs = join_bucket_nn_rtree(stop, sttemp);
		break;
	default:
		pairs = join_bucket_spjoin(stop, sttemp);
	}
	std::cerr <<"Special T[" << sttemp.tile_id << "] |" << sttemp.mbbdata[SID_1].size()
			  << "|x|" << sttemp.mbbdata[stop.sid_second_set].size()
		      << "|=|" << pairs << "|" << std::endl;
	total_query_exec += clock() - start_query_exec;
	start_reading_data = clock();
	release_mem(stop, sttemp, maxCardRelease);
	return pairs;
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


	start_reading_data = clock();

	attach_shm(stop);
#ifdef DEBUG
	std::cerr<<"line content: partition_id object_id dataset_id mbbs*6 offset length"<<std::endl;
#endif
	// Read line by line inputs
	while (std::cin && getline(std::cin, input_line) && !std::cin.eof()) {

		tokenize(input_line, fields, TAB, true);
		if(fields.size()!=11){//skip the corrupted lines
			std::cerr<<"malformed: "<<fields.size()<<endl;
						std::cerr<<input_line;
			continue;
		}

#ifdef DEBUG
		// the input is in format (11 fields):
		// partition_id object_id dataset_id mbbs*6 offset length
		//std::cerr<<input_line<<std::endl;
#endif
		/* Parsing fields from input */
		tile_id = fields[0];
		// dataset id
		sid = atoi(fields[2].c_str());
		try {
			// Parsing MBB
			mbb_ptr = new struct mbb_3d();
			for (int k = 0; k < NUMBER_DIMENSIONS; k++) {
				mbb_ptr->low[k] = std::atof(fields[3 + k].c_str());
			}
			for (int k = 0; k < NUMBER_DIMENSIONS; k++) {
				mbb_ptr->high[k] = std::atof(fields[6 + k].c_str());
			}
		} catch (...) {
			std::cerr << "******MBB Parsing Error******" << std::endl;
			return -1;
		}

		try {
			offset = atol(fields[9].c_str());
			length = atol(fields[10].c_str());
		}catch (...) {
			std::cerr << "******Offset and Length Parsing Error******" << std::endl;
			return -1;
		}

		/* Process the current tile (bucket) when finishing reading all objects belonging
		   to the current tile */
		if (previd.compare(tile_id) != 0 && previd.size() > 0 ) {
			sttemp.tile_id = previd;
			join_bucket(stop, sttemp);
			tile_counter++;
		}

		// populate the bucket for join
		sttemp.offsetdata[sid].push_back(offset);
		sttemp.lengthdata[sid].push_back(length);
		sttemp.mbbdata[sid].push_back(mbb_ptr);

		/* Update the field */
		previd = tile_id;
		fields.clear();
	}
	// Process the last tile (what remains in memory)
	sttemp.tile_id = tile_id;
	join_bucket(stop, sttemp);
	tile_counter++;

	shmdt(shm_ptr);

	return tile_counter;
}




