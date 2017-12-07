
//  This program combines the compressed data stored in a location given by user input

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
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

#include <sys/types.h>
#include <boost/program_options.hpp>

#include "progparams/string_constants.h"
#include "utilities/tokenizer.h"

#define NUMBER_DIMENSIONS 3

// for the pipeline
#define MAPPER_ID_FIELD 0
#define OBJ_ID_FIELD 1
#define OFFSET_FIELD 9
#define LENGTH_FIELD 10
#define JOIN_IDX_FIELD 2
#define MBB_START_FIELD 3

#define COMP_FILE_PREFIX "/compressionoutput"

namespace po = boost::program_options;
using namespace std;

struct compress_struct {
	std::string inputindex;
	std::string inputbindir;
	std::string outputbin;
};

/* Extract parameters for query operator */
bool extract_params(int ac, char** av, struct compress_struct & compstruct){ 
	try {
		po::options_description desc("Options");
		desc.add_options()
			("help", "this help message")
			("inputindex,a", po::value<string>(&(compstruct.inputindex)), "Directory path of index file (Optional. By default it uses standard input)")
			("inputbin,b", po::value<string>(&(compstruct.inputbindir)), "Directory path of binary files (compressed data).")      
			("outputbin,p", po::value<string>(&(compstruct.outputbin)), "Output path for the final binary file");
		po::variables_map vm;        
		po::store(po::parse_command_line(ac, av, desc), vm);
		po::notify(vm);   

		if (vm.count("help")) {
			cerr << desc <<	endl;
			return 0;
		}


	} catch(exception& e) {
		cerr << "error: " << e.what() << "\n";
		return false;
	}
	catch(...) {
		cerr << "Exception of unknown type!\n";
		return false;
	}

	return true;
}


int main(int argc, char** argv) {
	struct compress_struct compstruct;
	if (!extract_params(argc, argv, compstruct)) {
		#ifdef DEBUG 
		cerr <<"ERROR: program parameter extraction error." << endl 
			<< "Please see documentations (runs -h), or contact author." << endl;
		#endif
		return -1;
	}

	string input_line; // Temporary line
	vector<string> fields; // Temporary fields

	/*
	// MBB and space info
	double tmp_x, tmp_y, tmp_z;
	double low[3];
	double high[3];
	bool firstLineRead = false;
	double space_low[3];
	double space_high[3];
	*/

	string mapper_id;
	
	string prev_id = "";
	size_t offset = 0;	
	int counter = 0;

	string file_prefix = COMP_FILE_PREFIX;

	long count_objects = 0;
	// output binary file for the compressed data
	ofstream finalFile((compstruct.outputbin).c_str(), ios::binary);
	// ios::out | ios::binary);
	stringstream ss;
	// Read standard input
	while (cin && getline(cin, input_line) && !cin.eof()) {
		
		// Removal of \r symbol on Windows
		// Here data is from HDFS, not from Windows
		/*if (input_line.at(input_line.size() - 1) == '\r') {
			input_line = input_line.substr(0, input_line.size() - 1);
		}*/
		tokenize(input_line, fields, TAB, true);
		mapper_id = fields[MAPPER_ID_FIELD];
	//	cerr << "offset: " << offset << endl;
		// output the fields
		cout << mapper_id 
			<< TAB << fields[OBJ_ID_FIELD]// << TAB << obj_id 
			<< TAB << fields[JOIN_IDX_FIELD];  // << TAB << join_idx
		// Output MBBs
		for (counter = 0; counter < 2 * NUMBER_DIMENSIONS; counter++) {
			cout << TAB << fields[MBB_START_FIELD + counter];
		}
		cout << TAB << offset
		     << TAB << fields[LENGTH_FIELD] << endl;

		// Adding to the running offset total. This will be the offset of the next object
		offset += stol(fields[LENGTH_FIELD]);

		// compute offset for next object

		#ifdef DEBUG
		cerr << "Processing " << count_objects << endl;
		#endif

		fields.clear();

		// Reading and writing the binary file if necessary
		if (prev_id.compare(mapper_id) != 0 && prev_id.size() > 0 ) {
			ss.str("");
			ss << compstruct.inputbindir << file_prefix << prev_id;
			string fname = ss.str();
			const char *fnamechar = fname.c_str();
			//string inputbinfilename = compstruct.inputbindir + file_prefix + prev_id;
			ifstream inFile(fnamechar, ios::in | ios::binary);
			finalFile << inFile.rdbuf();
		}
		prev_id = mapper_id;
	} 

	// Process last batch
	// Reading and writing the binary file if necessary
	if (prev_id.size() > 0 ) {
		#ifdef DEBUG
		cerr << "Got here" << endl;
		#endif
			ss.str("");
			ss << compstruct.inputbindir << file_prefix << prev_id;
			string fname = ss.str();
			const char *fnamechar = fname.c_str();
			//string inputbinfilename = compstruct.inputbindir + file_prefix + prev_id;
			ifstream inFile(fnamechar, ios::in | ios::binary);
			finalFile << inFile.rdbuf();
/*
		string inputbinfilename = compstruct.inputbindir + file_prefix + mapper_id;
		ifstream inFile(inputbinfilename.c_str(), ios::in | ios::binary);
		finalFile << inFile.rdbuf();
*/
			

	}

	finalFile.flush();
	finalFile.close();
	#ifdef DEBUG
	cerr << offset << endl; // the total size
	#endif
	return true;
}

