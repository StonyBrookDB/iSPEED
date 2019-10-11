#ifndef QUERY_PROCESSOR_DATASTRUCTS_3D_HPP
#define QUERY_PROCESSOR_DATASTRUCTS_3D_HPP

#include <unistd.h>
#include <getopt.h>
#include <boost/program_options.hpp>
#include <stdio.h>
#include <iostream>

#include "global_define.h"
#include "mapreducecpp/mapreduce_cpp.hpp"

/* Parameter placeholders */
// Holds spatial dimensions of the space
struct space_info {
	long num_objects;
	double space_low[3];
	double space_high[3];
	double total_size;
} ;

// Manages piped I/O
struct flow_info;
struct framework_vars;

// Combined framework variables
struct framework_vars {
	
	// MapReduce-related parameter
	std::string binary_prefix;
	std::string streaming_path; // Absolute path to the hadoop streaming jar
	std::string hadoopcmdpath; // Absolute path to the command line hadoop executable
	std::string hdfscmdpath; // Absolute path to the command line hadoop executable

	bool overwritepath = true; // Overwrite HDFS directory if it already exists
	int numreducers = 1; // number of reducer
	std::string compressed_data_path;

	struct space_info spinfo;

	// Input data variables	
	std::string input_path_1;
	std::string input_path_2;
	
	// Query variables
	ISPEEDOperation query_type;

	//for resque
	std::string predicate;
	double distance;
	int decomp_lod = 100;


	// Partitioning variables
	std::string partition_method = PARTITION_FG_3D; // Default selection
	double sampling_rate = 1; // to be changed for different sampling method
	long bucket_size = -1;

	// Temporary paths
	std::string output_path;

	std::string mbb_output;
	std::string resque_input;
	std::string partitionpath;
	std::string partitionpathout;
	std::string joinoutputpath;
};

bool extract_params(int argc, char **argv, struct framework_vars &fr_vars);

#endif
