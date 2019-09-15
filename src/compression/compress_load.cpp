
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
#include <sys/types.h>
#include <sys/ipc.h> 
#include <sys/shm.h> 

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include <boost/program_options.hpp>

#include "progparams/string_constants.h"
#include "utilities/tokenizer.h"

#define NUMBER_DIMENSIONS 3


// Default key
#define COMPRESSION_KEY 5678

namespace po = boost::program_options;
using namespace std;

struct compress_struct {
	std::string binarypath;
	key_t key;
	bool cleanmode;
	long totalsize;
};

/* Extract parameters for query operator */
bool extract_params(int ac, char** av, struct compress_struct & compstruct){ 
	try {
		po::options_description desc("Options");
		desc.add_options()
			("help", "this help message")
			("inputbin,n", po::value<string>(&(compstruct.binarypath)), "In path for the final binary file")
			("key,k", po::value<key_t>(&(compstruct.key)), "(Optional) Key for the shared memory segment(compressed data).")      
			("size,s", po::value<long>(&(compstruct.totalsize)), "(Optional) Size of shared memory segment")
			("remove,r", "Remove all shared segment (cleaning mode)");
		po::variables_map vm;        
		po::store(po::parse_command_line(ac, av, desc), vm);
		po::notify(vm);   

		if (vm.count("help")) {
			cerr << desc <<	endl;
			return 0;
		}

		if (!vm.count("key")) {
			compstruct.key = COMPRESSION_KEY;
		}

		if (vm.count("remove")) {
			compstruct.cleanmode = true;
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
	#ifdef DEBUG
	cerr << "size_t max: " << SIZE_MAX << endl;
	
	#endif


	struct compress_struct compstruct;
	compstruct.cleanmode = false;
	streamsize size; // size of the memory segment

	if (!extract_params(argc, argv, compstruct)) {
		#ifdef DEBUG 
		cerr <<"ERROR: program parameter extraction error." << endl 
			<< "Please see documentations (runs -h), or contact author." << endl;
		#endif
		return -1;
	}
	
	// use shmget to create a memory segment
	int shmid;
	char *shm, *s;

	if (compstruct.cleanmode) {
		/*
		 * Load the segment
		 */
		if ((shmid = shmget(compstruct.key, compstruct.totalsize, 0666)) < 0) {
			perror("shmget error");
			exit(1);

			if ((shmctl(shmid,IPC_RMID, 0)) == -1) {
				cerr <<" ERROR(C++) with shmctl(IPC_RMID): " << strerror(errno) << endl;
				return -1;
			} else {
				return 0;
			}
		}
	}

	// Write
	ifstream file(compstruct.binarypath.c_str(), std::ios::binary | std::ios::ate);
	size = file.tellg(); // total size
	file.seekg(0, std::ios::beg);
	#ifdef DEBUG
	cerr << "size: " << size << " key: " << compstruct.key << endl;
	#endif

	// Create the segment.
	if ((shmid = shmget(compstruct.key, size, IPC_CREAT | 0666)) < 0) {
		perror("shmget");
		exit(1);
	}

	// Now we attach the segment to our data space.
	if ((shm = (char *) shmat(shmid, NULL, 0)) == (char *) -1) {
		perror("shmat");
		exit(1);
	}    

	cout << "Content of first byte: " << (* ((int *) shm)) << endl;
	// read in chunks
	

	// old read
	// file.read(shm, size);
	cerr << "Amount read = " << file.gcount() << endl;
	/*
	if (file.read(shm, size))
	{
		#ifdef DEBUG
		cerr << "Done loading objects into memory" << endl;
		#endif
	} else {
		cerr << "Cannot load" << endl;
		return 1;
	}*/

	cout << "Content of first byte: " << (* ((int *) shm)) << endl;
	shmdt(shm);

	#ifdef DEBUG
	cerr << "segment ID: " << shmid << endl; // the total size
	cerr << "total size: " << size  << endl;
	#endif
	return 0;
}

