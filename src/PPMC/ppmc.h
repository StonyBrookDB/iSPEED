
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <map>
#include <cstdlib>
#include <getopt.h>
#include <time.h>
#include <float.h>
#include <sys/time.h>

#include <stdio.h>
#include <assert.h>
#include <string.h>


#include <vector>
#include <istream>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/replace.hpp>


/// RESQUE related
#include <progparams/string_constants.h>
#include <utilities/tokenizer.h>

// Constants
#include <progparams/resque_constants_3d.h>

// Program parameters
#include <progparams/resque_params_3d.hpp>


#include <PPMC/compressionsupport.hpp>

void init(struct query_op &stop, struct query_temp &sttemp);

bool compress_data(char* stdin_file_name, std::string output_path, char* mapper_id, long join_id);
