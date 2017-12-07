
/* Containing methods to extract parameters and store them in query operator */
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <map>
#include <iterator>
#include <cstdlib> 
#include <vector>
#include <boost/program_options.hpp>
#include <progparams/partition_datastructs_3d.hpp>

#define NUMBER_DIMENSIONS 3

namespace po = boost::program_options;

void init_params_partitioning(struct partition_op & partop) {
	partop.to_be_normalized = false;
	partop.to_be_denormalized = false;
	partop.offset = 1; // Default offset
}

// Read the cache file containing MBB of regions
bool read_partition_file(struct partition_op &partop) {
	std::string inputline;
	std::string tile_id;
	std::ifstream infile(partop.file_name);

	while (getline(infile, inputline)) {
		std::istringstream ss(inputline);
		struct mbb_info *tmp = new struct mbb_info();
		if (!(ss >> tile_id >> tmp->low[0] >> tmp->low[1] >> tmp->low[2] >> tmp->high[0] >> tmp->high[1] >> tmp->high[2])) { 
			return false;
		}

		partop.region_mbbs[tile_id] = tmp;
	}
	// Update the global information
	if (partop.region_mbbs.size() == 1) {
		struct mbb_info *tmp = partop.region_mbbs.begin()->second;
		partop.low[0] = tmp->low[0];
		partop.low[1] = tmp->low[1];
		partop.low[2] = tmp->low[2];
		partop.high[0] = tmp->high[0];
		partop.high[1] = tmp->high[1];
		partop.high[2] = tmp->high[2];
	}

	return true;
}


void update_partop_info(struct partition_op & partop, std::string uppertileid, std::string newprefix) {
	partop.prefix_tile_id = newprefix;
	if (partop.region_mbbs.size() > 1) {
		struct mbb_info *tmp = partop.region_mbbs[uppertileid]; 
		partop.low[0] = tmp->low[0];
		partop.low[1] = tmp->low[1];
		partop.low[2] = tmp->low[2];
		partop.high[0] = tmp->high[0];
		partop.high[1] = tmp->high[1];
		partop.high[2] = tmp->high[2];
	}
}

/* Extract parameters for partitioning */
bool extract_params_partitioning(int ac, char** av, struct partition_op & partop){ 
	init_params_partitioning(partop);
	try {
		po::options_description desc("Options");
		desc.add_options()
			("help", "this help message")
			("norm", "Normalize the data")
			("denorm", "Denormalize the data")
			("bucket,b", po::value<long>(&partop.bucket_size), "Expected bucket size")
			("cachefilename,c", po::value<std::string>(&(partop.file_name)), "Name of cache file")
			("offset,f", po::value<int>(&partop.offset), "(Optional) Offset from where MBB starts")      
			("min_x,k", po::value<double>(&(partop.low[0])), "(Optional) Spatial min x")
			("min_y,l", po::value<double>(&(partop.low[1])), "(Optional) Spatial min y")
			("min_z,p", po::value<double>(&(partop.low[2])), "(Optional) Spatial min z")
			("max_x,m", po::value<double>(&(partop.high[0])), "(Optional) Spatial max x")
			("max_y,n", po::value<double>(&(partop.high[1])), "(Optional) Spatial max y")
			("max_z,q", po::value<double>(&(partop.high[2])), "(Optional) Spatial max z");
		po::variables_map vm;        
		po::store(po::parse_command_line(ac, av, desc), vm);
		po::notify(vm);   

		if (vm.count("help")) {
			std::cerr<< desc << std::endl;
			return 0;
		}

		if (vm.count("norm")) {
			partop.to_be_normalized = true;
		}
		if (vm.count("denorm")) {
			partop.to_be_denormalized = true;
		}

		if (vm.count("cachefilename")) {
			read_partition_file(partop);
		}

	} catch(std::exception& e) {
		std::cerr<< "error: " << e.what() << "\n";
		return false;
	}
	catch(...) {
		std::cerr<< "Exception of unknown type!\n";
		return false;
	}

	if (partop.to_be_denormalized && partop.to_be_normalized) {
		/* Mutually exclusive options */
		#ifdef DEBUG
		std::cerr<< "Mutually exclusive functions -n and -o" << std::endl;
		#endif
		return false;
	}

	return true;
}

void cleanup(struct partition_op & partop){
	std::cerr<< "cleaning up" << std::endl;
	for (std::map<std::string,struct mbb_info *>::iterator it= partop.region_mbbs.begin();
			it!= partop.region_mbbs.end(); ++it) {
		std::cerr<< "cleaning: " << it->first << "---" << it->second->low[0] << TAB << it->second->low[1] << TAB << it->second->high[0] << it->second->high[1] << std::endl;
		delete it->second;
	}
	partop.region_mbbs.clear();
}
