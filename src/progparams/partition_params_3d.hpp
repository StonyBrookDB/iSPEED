#ifndef PARTITION_PARAMS_HPP
#define PARTITION_PARAMS_HPP

#include <map>
#include <unistd.h>
#include <getopt.h>

#include "global_define.h"

#define NUMBER_DIMENSIONS 3


/* Struct to hold MBB */
struct mbb_info {
	double low[3];
	double high[3];
};

/* Query operator */
struct partition_op {
	long bucket_size;
	int offset = 1; // Field offset
	std::string prefix_tile_id;
	int object_count;

	/* Space span  */
	double low[3];
	double high[3];

	/* only at most one can be true  */
	bool to_be_normalized = false;
	bool to_be_denormalized = false;

	bool parallel_partitioning;
	std::string file_name;
	std::map<std::string, struct mbb_info *> region_mbbs;
};

void cleanup(struct partition_op & partop);
bool extract_params_partitioning(int ac, char** av, struct partition_op & partop);
#endif
