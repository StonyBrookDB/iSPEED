#ifndef QUERY_PROCESSOR_DATASTRUCTS_3D_HPP
#define QUERY_PROCESSOR_DATASTRUCTS_3D_HPP

#include <unistd.h>
#include <getopt.h>
#include <boost/program_options.hpp>
#include <stdio.h>
#include <iostream>

#include <progparams/string_constants.h>
#include "mapreducecpp/mapreduce_cpp.hpp"

extern char nametemplate[];
extern char nametemplate2[];
extern char nametemplate3[];
extern char nametemplate4[];
extern char nametemplate5[];
extern char nametemplate6[];

// Supported query types
const std::string QUERYPROC_PARTITION = "partition";
const std::string QUERYPROC_CONTAINMENT = "containment";
const std::string QUERYPROC_JOIN = "spjoin";
const std::string QUERYPROC_NN_VORONOI = "spnn_voronoi";
const std::string QUERYPROC_NN_RTREE = "spnn_rtree";
// Names of binaries
//const std::string MANIPULATE = "manipulate_2d";
//const std::string SPACE_EXTRACTOR = "get_space_dimension";
const std::string SKELETON = "skeleton_3d";
const std::string VORONOI = "voronoi_3d";

// for data compression
const std::string COMPRESSION = "ppmc";

const std::string MANIPULATE = "manipulate_3d";
const std::string SPACE_EXTRACTOR = "stats_extract_space_dims_3d";
const std::string DUPLICATE_REMOVER = "duplicate_remover";
const std::string STAT_COLLECT_MAPPER = "collect_tile_stat";
const std::string STAT_COLLECT_REDUCER = "combine_stat";
const std::string MBB_SAMPLER = "sampler";

const std::string RESQUE = "resque_3d";

// for 3d
const std::string PARTITION_FG_3D = "fg_3d";
const std::string PARTITION_OT_3D = "ot_3d";

/* Parameter placeholders */
// Holds spatial dimensions of the space
struct space_info {
	long num_objects;
	double space_low[3];
	double space_high[3];
	double total_size;
} ;

enum Operation{
	COMPRESS = 0,
	PARTITION = 1,
	JOIN = 2,
	DUPLICATE_REMOVAL = 3,
	ERROR = 4
};

const string operation_str[5] = {"compress", "partition", "join", "duplicate_removal", "error"};


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
	long size_of_compressed_data = -1;

	struct space_info spinfo;

	// Input data variables	
	std::string input_path_1;
	std::string input_path_2;
	long size_1 = -1;
	long size_2 = -1;
	int join_cardinality = 1;
	
	// Query variables
	Operation query_type;

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
