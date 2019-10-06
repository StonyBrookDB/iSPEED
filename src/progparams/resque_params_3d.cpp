#include <progparams/resque_params_3d.hpp>
#include <boost/program_options.hpp>
#include <utilities/tokenizer.h>


// for sp join yes or no intersection
bool intersection_flag = false;

char * shm_ptr = NULL;
//long maxoffset = 0;

/* Containing methods to extract parameters and store them in query operator */

/* Display help message to users */
void usage(){
	std::cerr  << std::endl << "Usage: program_name [OPTIONS]" << std::endl << "OPTIONS:" << std::endl;
	std::cerr << TAB << "-p,  --predicate" << TAB <<  "The spatial join predicate for query processing. \
Acceptable values are [st_intersects, st_nn_voronoi, st_nn_rtree, st_disjoint, st_overlaps, st_within, st_equals,\
st_dwithin, st_crosses, st_touches, st_contains, st_nearest, st_nearest2]." << std::endl;
	std::cerr << TAB << "-d, --distance" << TAB << "Used together with st_dwithin predicate to \
indicate the join distance or used together with st_nearest to indicate the max distance \
to search for nearest neighbor. \
This field has no effect on other join predicates." << std::endl;	
	std::cerr << TAB << "-k, --knn" << TAB << "The number of nearest neighbor\
Only used in conjuction with the st_nearest or st_nearest2 predicate" << std::endl;
	std::cerr << TAB << "-c, --cachefile" << TAB << "The name of cache file. \
Dependent on operation/task" << std::endl;
	std::cerr << TAB << "-r, --replicate" << TAB <<  "Optional [true | false]. \
Indicates whether result pair are output twice for result pair that logically exists\
when the pair is swapped (position). e.g. intersection between object 1 and 2\
also indicates a logical intersection between object 2 and 1. \
Default value is false" << std::endl;
	std::cerr << TAB << "-e, --earth" << TAB << "Optional [true | false]\
Indicate wheather to compute spherical distance for point-point on earth." << std::endl;
}

/* This function converts expensive char comparison into fast int comparison 
 * during spatial processing to determine the spatial predicate being used*/
int get_join_predicate(char * predicate_str)
{
	if (strcmp(predicate_str, PARAM_PREDICATE_INTERSECTS.c_str()) == 0) {
		return ST_INTERSECTS ; 
	} 
	else if (strcmp(predicate_str, PARAM_PREDICATE_TOUCHES.c_str()) == 0) {
		return ST_TOUCHES;
	} 
	else if (strcmp(predicate_str, PARAM_PREDICATE_CROSSES.c_str()) == 0) {
		return ST_CROSSES;
	} 
	else if (strcmp(predicate_str, PARAM_PREDICATE_CONTAINS.c_str()) == 0) {
		return ST_CONTAINS;
	} 
	else if (strcmp(predicate_str, PARAM_PREDICATE_ADJACENT.c_str()) == 0) {
		return ST_ADJACENT;
	} 
	else if (strcmp(predicate_str, PARAM_PREDICATE_DISJOINT.c_str()) == 0) {
		return ST_DISJOINT;
	}
	else if (strcmp(predicate_str, PARAM_PREDICATE_EQUALS.c_str()) == 0) {
		return ST_EQUALS;
	}
	else if (strcmp(predicate_str, PARAM_PREDICATE_DWITHIN.c_str()) == 0) {
		return ST_DWITHIN;
	}
	else if (strcmp(predicate_str, PARAM_PREDICATE_WITHIN.c_str()) == 0) {
		return ST_WITHIN;
	}
	else if (strcmp(predicate_str, PARAM_PREDICATE_OVERLAPS.c_str()) == 0) {
		return ST_OVERLAPS;
	}
	else if (strcmp(predicate_str, PARAM_PREDICATE_NEAREST.c_str()) == 0) {
		return ST_NEAREST;
	}
	else if (strcmp(predicate_str, PARAM_PREDICATE_NEAREST_NO_BOUND.c_str()) == 0) {
		return ST_NEAREST_2;
	}
	else if (strcmp(predicate_str, PARAM_PREDICATE_NEAREST_NEIGHBOR_VORONOI.c_str()) == 0) {
		return ST_NN_VORONOI;
	}
	else if (strcmp(predicate_str, PARAM_PREDICATE_NEAREST_NEIGHBOR_RTREE.c_str()) == 0) {
		return ST_NN_RTREE;
	}
	else {
		#ifdef DEBUG
		std::cerr << "unrecognized join predicate " << std::endl;
		#endif
		return -1;
	}
}


// This function  extracts command line arguments
bool extract_params(int argc, char** argv, struct query_op &stop, struct query_temp &sttemp){ 
	// initlize query operator 
	stop.expansion_distance = 0.0;
	stop.k_neighbors = 0;
	stop.join_predicate = 0;


	stop.needs_area_1 = false;
	stop.needs_area_2 = false;
	stop.needs_union = false;
	stop.needs_intersect = false;
	stop.needs_dice = false;
	stop.needs_jaccard = false;

	
	//for 3d
	stop.needs_volume_1 = false;
	stop.needs_volume_2 = false;
	stop.needs_intersect_volume = false;
	stop.needs_nn_distance = false;

	stop.result_pair_duplicate = true;
	stop.use_earth_distance = false;

	
	// for compression
	stop.shm_max_size = 0;
	stop.decomp_lod = 100; // 100% decompression

	sttemp.nearest_distances.clear();	
	sttemp.area1 = -1;
	sttemp.area2 = -1;
	sttemp.union_area = -1;
	sttemp.intersect_area = -1;
	sttemp.dice = -1;
	sttemp.jaccard = -1;
	sttemp.distance = -1;

	//for 3d
	sttemp.volume1 = -1;
	sttemp.volume2 = -1;
	sttemp.intersect_volume = -1;
	sttemp.nn_distance = -1;


	int option_index = 0;
	/* getopt_long uses opterr to report error*/
	opterr = 0 ;
	struct option long_options[] =
	{
		{"distance",   required_argument, 0, 'd'},
		{"predicate",  required_argument, 0, 'p'},
		{"knn",  required_argument, 0, 'k'},
		{"earth",     required_argument, 0, 'e'},
		{"replicate",     required_argument, 0, 'r'},
		{"cachefile",     required_argument, 0, 'c'},
		{"lod",       required_argument, 0, 'y'},
		{"shmmaxsize", required_argument, 0, 'z'},
		{"join_cardinality",required_argument, 0, 'j'},

		// Specific to controller only
		{0, 0, 0, 0}
	};

	int c;
	while ((c = getopt_long (argc, argv, "c:j:p:d:k:r:e:y:z:", long_options, &option_index)) != -1){
		switch (c)
		{
			case 0:
				/* If this option set a flag, do nothing else now. */
				if (long_options[option_index].flag != 0)
					break;
				std::cout << "option " << long_options[option_index].name ;
				if (optarg)
					std::cout << "a with arg " << optarg ;
				std::cout << std::endl;
				break;
			case 'p':
				stop.join_predicate = get_join_predicate(optarg);
				#ifdef DEBUG
					std::cerr << "predicate: " << stop.join_predicate << std::endl;
       	                        #endif
				break;
			case 'd':
				stop.expansion_distance = atof(optarg);
				#ifdef DEBUG
					std::cerr << "Search distance/Within distance parameter " 
						<< stop.expansion_distance << std::endl;
				#endif
				break;

			case 'k':
				stop.k_neighbors = strtol(optarg, NULL, 10);
				#ifdef DEBUG
					std::cerr << "Number of neighbors: " 
						<< stop.k_neighbors << std::endl;
				#endif
				break;
			case 'r':
				stop.result_pair_duplicate = strcmp(optarg, "true") == 0;
				#ifdef DEBUG
					std::cerr << "Allows symmetric result pairs: "
						<< stop.result_pair_duplicate << std::endl;
				#endif
				break;

			case 'e':
				stop.use_earth_distance = strcmp(optarg, "true") == 0;
				#ifdef DEBUG
					std::cerr << "Using earth distance: "
						<< stop.use_earth_distance << std::endl;
				#endif
				break;
			/*
			case 'q':
				stop.sample_rate = atof(optarg);
				stop.use_sampling = true;
				#ifdef DEBUG
					cerr << "Sample rate: " << stop.sample_rate << std::endl;
				#endif
				break;
			*/
			case 'z':
				stop.shm_max_size = strtol(optarg, NULL, 10);
				#ifdef DEBUG
					std::cerr << "Shared memory max size: " << stop.shm_max_size << std::endl;
				#endif
				break;

			case 'y':
				stop.decomp_lod = strtol(optarg, NULL, 10);
				#ifdef DEBUG
					std::cerr << "Level of details for decompression: " << stop.decomp_lod << std::endl;
				#endif
				break;
			case 'j':

				stop.join_cardinality = strtol(optarg, NULL, 10);
				#ifdef DEBUG
				std::cerr << "join cardinality: " << stop.join_cardinality << std::endl;
				#endif
				if (stop.join_cardinality <= 0) {
					#ifdef DEBUG
					std::cerr << "Join cardinality are NOT set properly. Please refer to the documentation." << std::endl ;
					#endif
					return false;
				}
				break;
			case 'c':
				stop.cachefilename = optarg;
				#ifdef DEBUG
				std::cerr << "cache file name: " << stop.cachefilename << std::endl;
				#endif
				break;
			case '?':
				return false;
				/* getopt_long already printed an error message. */
				break;

			default:
				return false;
		}
	}

	// Adjusting the actual geometry field (shift) to account
	//   for tile_id and join_index
	#ifdef DEBUG
    std::cerr << "join cardinality: " << stop.join_cardinality << std::endl;
	#endif

	// query operator validation 
	/*
	if (!stop.drop_join_idx && stop.join_predicate <= 0 ) {
		#ifdef DEBUG 
		cerr << "Query predicate is NOT set properly. Please refer to the documentation." << std::endl ; 
		#endif
		return false;
	}*/
	// check if distance is set for dwithin predicate
	if (stop.join_predicate == ST_DWITHIN && stop.expansion_distance == 0.0) { 
		#ifdef DEBUG 
		std::cerr << "Distance parameter is NOT set properly. Please refer to the documentation." << std::endl ;
		#endif
		return false;
	}
	if ((stop.join_predicate == ST_NEAREST || stop.join_predicate == ST_NEAREST_2)
		 && stop.k_neighbors <= 0) {
		#ifdef DEBUG 
		std::cerr << "K-the number of nearest neighbors is NOT set properly. Please refer to the documentation." << std::endl ;
		#endif
		return false; 
	}

	return true;
}
