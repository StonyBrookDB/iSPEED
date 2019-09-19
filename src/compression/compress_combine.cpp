
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
	#ifdef DEBUG
	cerr<<"merging into file "<<compstruct.outputbin<<endl;
	#endif

	//for merging the space information
	//generated by different mappers
	double space_low[3] = {DBL_MAX,DBL_MAX,DBL_MAX};
	double space_high[3] = {0,0,0};
	double low_buf[3];
	double high_buf[3];

	string input_line; // Temporary line
	vector<string> fields; // Temporary fields

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
		cerr << input_line<<endl;
		//input format
		//mapper_id obj_id dataset_id mbb*6 length
		
		tokenize(input_line, fields, TAB, true);
		mapper_id = fields[MAPPER_ID_FIELD];
		//output the fields
		//mapper_id obj_id dataset_id mbb*6 offset length
		cout << fields[MAPPER_ID_FIELD]
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
		fields.clear();

		// Reading and writing the binary file if necessary
		if (prev_id.compare(mapper_id) != 0 && prev_id.size() > 0 ) {
			ss.str("");
			ss << compstruct.inputbindir << file_prefix << prev_id;
			string fname = ss.str();
			const char *fnamechar = fname.c_str();
			#ifdef DEBUG
			cerr<<"merging from file "<<ss.str()<<endl;
			#endif

			ifstream inFile(fnamechar, ios::in | ios::binary);
			inFile.seekg(0, ios::end);
			int size = inFile.tellg();
			inFile.seekg(0, ios::beg);
			char *rd_buffer = new char[size];
			inFile.read(rd_buffer, size-6*sizeof(double));
			//read and update the space low and high information
			inFile.read((char *)low_buf, 3*sizeof(double));
			inFile.read((char *)high_buf, 3*sizeof(double));
			for(int j=0;j<3;j++){
				if(low_buf[j]<space_low[j]){
					space_low[j] = low_buf[j];
				}
				if(high_buf[j]>space_high[j]){
					space_high[j] = high_buf[j];
				}
			}
			finalFile.write(rd_buffer, size-6*sizeof(double));
			delete(rd_buffer);
		}
		prev_id = mapper_id;
		count_objects++;
	} 

	// Process last batch
	// Reading and writing the binary file if necessary
	if (prev_id.size() > 0 ) {
		ss.str("");
		ss << compstruct.inputbindir << file_prefix << prev_id;
		string fname = ss.str();
		const char *fnamechar = fname.c_str();
		#ifdef DEBUG
		cerr<<"merging from file "<<ss.str()<<endl;
		#endif
		ifstream inFile(fnamechar, ios::in | ios::binary);
		inFile.seekg(0, ios::end);
		int size = inFile.tellg();
		inFile.seekg(0, ios::beg);
		char *rd_buffer = new char[size];
		inFile.read(rd_buffer, size-6*sizeof(double));
		//read and update the space low and high information
		inFile.read((char *)low_buf, 3*sizeof(double));
		inFile.read((char *)high_buf, 3*sizeof(double));
		for(int j=0;j<3;j++){
			if(low_buf[j]<space_low[j]){
				space_low[j] = low_buf[j];
			}
			if(high_buf[j]>space_high[j]){
				space_high[j] = high_buf[j];
			}
		}
		finalFile.write(rd_buffer, size-6*sizeof(double));
		delete(rd_buffer);
	}

	//attach the global space mbb to the file
	finalFile.write((char *)space_low, 3*sizeof(double));
	finalFile.write((char *)space_high, 3*sizeof(double));

	//cleaning
	finalFile.flush();
	finalFile.close();

	#ifdef DEBUG
	cerr << "processed:	"<<count_objects<<endl;
	cerr << "total size is:	"<<offset << endl; // the total size
	#endif
	return true;
}

