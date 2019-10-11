
// separators for parsing
#ifndef RESQUE_STR_CONSTANTS
#define RESQUE_STR_CONSTANTS

#include <string>
using namespace std;

// name of binary tools
const std::string COMPRESSION = "ppmc";
const std::string MANIPULATE = "manipulate_3d";
const std::string DUPLICATE_REMOVER = "duplicate_remover";
const std::string MBB_SAMPLER = "sampler";
const std::string RESQUE = "resque_3d";

// for partitioning
const std::string PARTITION_FG_3D = "fg_3d";
const std::string PARTITION_OT_3D = "ot_3d";

// other constants
const std::string TAB = "\t";
const std::string SEP = TAB;
const std::string BAR = "|";
const std::string DASH = "-";
const std::string COMMA = ",";
const std::string SPACE = " ";
const std::string SLASH = "/";


enum ISPEEDOperation{
	COMPRESS = 0,
	PARTITION = 1,
	JOIN = 2,
	DUPLICATE_REMOVAL = 3,
	ERROR = 4
};

const string operation_str[5] = {"compress", "partition", "join", "duplicate_removal", "error"};


enum Jointype{
	ST_ERROR = 0,
	ST_INTERSECTS = 1,
	ST_TOUCHES = 2,
	ST_CROSSES = 3,
	ST_CONTAINS = 4,
	ST_ADJACENT = 5,
	ST_DISJOINT = 6,
	ST_EQUALS = 7,
	ST_DWITHIN = 8,
	ST_WITHIN = 9,
	ST_OVERLAPS = 10,
	ST_NEAREST = 11,
	ST_NEAREST_2 = 12,
	ST_NN_VORONOI = 13,
	ST_NN_RTREE = 14
};

const std::string join_type_str[15] = {
		"st_error", "st_intersects", "st_touches", "st_crosses", "st_contains",
		"st_adjacent", "st_disjoint", "st_equals", "st_dwithin", "st_within",
		"st_overlaps", "st_nearest", "st_nearest2", "st_nn_voronoi", "st_nn_rtree"};

const int SID_1 = 1;
const int SID_2 = 2;
const int SID_NEUTRAL = 0;
const int OSM_SRID = 4326;

const char nametemplate[] = "/tmp/hadoopgisXXXXXX";


inline ISPEEDOperation get_operation(string type){
	for(int i=0;i<4;i++){
		if(strcmp(type.c_str(),operation_str[i].c_str())==0){
			return (ISPEEDOperation)i;
		}
	}
	return ERROR;
}

inline Jointype get_join_predicate(const char * predicate_str){
	for(int i=1;i<15;i++){
		if (strcmp(predicate_str, join_type_str[i].c_str()) == 0) {
			return (Jointype)i ;
		}
	}
	std::cerr << "unrecognized join predicate " <<predicate_str<< std::endl;
	return ST_ERROR;
}

#endif
