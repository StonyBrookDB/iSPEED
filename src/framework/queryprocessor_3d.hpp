#ifndef QUERY_PROCESSOR_3D_HPP
#define QUERY_PROCESSOR_3D_HPP
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <map>
#include <cstdlib> 
#include <getopt.h>
#include <time.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/fcntl.h>
#include <fcntl.h>

#include <boost/program_options.hpp>
#ifdef DEBUGSTAT
	#include <boost/accumulators/accumulators.hpp>
	#include <boost/accumulators/statistics/stats.hpp>
	#include <boost/accumulators/statistics/mean.hpp>
	#include <boost/accumulators/statistics/variance.hpp>
	using namespace boost::accumulators;
#endif

#include <progparams/string_constants.h>
#include <utilities/tokenizer.h>



// MapReduce-C++ interaction
#include <mapreducecpp/mapreduce_cpp.hpp>

#include <progparams/queryprocessor_params_3d.hpp>

bool extract_params(int argc, char **argv, struct framework_vars &fr_vars);

void execute_spjoin(struct framework_vars &fr_vars);
void execute_partition(struct framework_vars &fr_vars);
void execute_compress(struct framework_vars &fr_vars);
void execute_duplicate_removal(struct framework_vars &fr_vars);

bool compress_data(struct framework_vars &fr_vars);
bool partition_data(struct framework_vars &fr_vars);
bool join_data(struct framework_vars &fr_vars, char *cachefilefullpath);


#ifdef DEBUGSTAT
void post_process_stat(char *tmpFile, stringstream &ss);
#endif

#endif
