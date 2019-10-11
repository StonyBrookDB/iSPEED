#include <progparams/resque_params_3d.hpp>
#include <boost/program_options.hpp>
#include <utilities/tokenizer.h>

/* Containing methods to extract parameters and store them in query operator */

/* Display help message to users */
void usage(){
	std::cerr <<"Usage: program_name [OPTIONS]" << std::endl << "OPTIONS:" << std::endl;
	std::cerr << TAB << "-p,  --predicate" << TAB <<  "The spatial join predicate for query processing. "
			"Acceptable values are [st_intersects, st_nn_voronoi, st_nn_rtree, st_disjoint, st_overlaps, st_within, st_equals, "
			"st_dwithin, st_crosses, st_touches, st_contains, st_nearest, st_nearest2]." << std::endl;
	std::cerr << TAB << "-d, --distance" << TAB << "Used together with st_dwithin predicate to "
			"indicate the join distance or used together with st_nearest to indicate the max distance "
			"to search for nearest neighbor. This field has no effect on other join predicates." << std::endl;
	std::cerr << TAB << "-k, --knn" << TAB << "The number of nearest neighbor. "
			"Only used in conjuction with the st_nearest or st_nearest2 predicate" << std::endl;
	std::cerr << TAB << "-l, --lod" << TAB << "level of details, 100 by default" << std::endl;
	std::cerr << TAB << "-s, --shmmaxsize" << TAB << "size of the shared memory" << std::endl;
	std::cerr << TAB << "-j, --join_cardinality" << TAB << "1 for one dateset and 2 for two datasets" << std::endl;
	std::cerr << TAB << "-h, --help" << TAB << "print the usage" << std::endl;
}

/* This function converts expensive char comparison into fast int comparison 
 * during spatial processing to determine the spatial predicate being used*/
inline Jointype get_join_predicate(char * predicate_str)
{
	for(int i=1;i<15;i++){
		if (strcmp(predicate_str, join_type_str[i].c_str()) == 0) {
			return (Jointype)i ;
		}
	}
	std::cerr << "unrecognized join predicate " <<predicate_str<< std::endl;
	return ST_ERROR;
}


// This function  extracts command line arguments
bool extract_params(int argc, char** argv, struct query_op &stop, struct query_temp &sttemp){ 

	int option_index = 0;
	/* getopt_long uses opterr to report error*/
	opterr = 0 ;
	struct option long_options[] =
	{
		{"distance",   required_argument, 0, 'd'},
		{"predicate",  required_argument, 0, 'p'},
		{"knn",  required_argument, 0, 'k'},
		{"lod",       required_argument, 0, 'l'},
		{"shmmaxsize", required_argument, 0, 's'},
		{"join_cardinality",required_argument, 0, 'j'},
		{"help",no_argument, 0, 'h'},
		// Specific to controller only
		{0, 0, 0, 0}
	};

	int c;
	while ((c = getopt_long (argc, argv, "d:p:k:l:s:j:h", long_options, &option_index)) != -1){
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
				if(stop.join_predicate == ST_ERROR){
					std::cerr << "invalid predicate: " << optarg << std::endl;
					exit(-1);
				}
#ifdef DEBUG
				else{
					std::cerr << "predicate: " << join_type_str[stop.join_predicate] << std::endl;
				}
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
			case 's':
				stop.shm_max_size = strtol(optarg, NULL, 10);
				#ifdef DEBUG
					std::cerr << "Shared memory max size: " << stop.shm_max_size << std::endl;
				#endif
				break;

			case 'l':
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
				if (stop.join_cardinality != 1&&stop.join_cardinality != 2) {
					#ifdef DEBUG
					std::cerr << "Join cardinality are NOT set properly. Please refer to the documentation." << std::endl ;
					#endif
					return false;
				}
				break;
			case 'h':
				usage();
				exit(1);
			case '?':
				return false;
				/* getopt_long already printed an error message. */
				break;

			default:
				return false;
		}
	}

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
