
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
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#include <stdio.h>
#include <assert.h>
#include <string.h>


#include <vector>
#include <istream>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <GL/glut.h>
#include <CGAL/Timer.h>

/// RESQUE related
#include <progparams/string_constants.h>
#include <utilities/tokenizer.h>

// Constants
#include <progparams/resque_constants_3d.h>
#include "../progparams/resque_params_3d.cpp"

// Program parameters
#include "PPMC/mymesh.h"
#include "PPMC/configuration.h"
#define SHMSZ     10000000000


bool compress_data(std::string output_path, char* mapper_id, long join_id);
